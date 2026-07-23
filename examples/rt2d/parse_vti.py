#!/usr/bin/env python3
# 解析 FreeLB VTI 输出，检查磁场 H 是否在界面处产生畸变
import base64, struct, glob, re, os, sys
import numpy as np

def parse_dataarray(content, name):
    # 找到 <DataArray ... Name="name" ...>...</DataArray>
    pat = r'<DataArray[^>]*Name="%s"[^>]*>\s*(.*?)\s*</DataArray>' % name
    m = re.search(pat, content, re.DOTALL)
    if not m: return None
    raw_b64 = m.group(1).strip()
    # VTK base64 binary 格式: header块 + data块 (两个独立的base64编码块)
    # Header: uint32 count (4字节) = 8个base64字符 (含'='填充)
    # Data: 实际数据 (独立的base64块)
    # 必须分别解码, 不能作为一个base64字符串解码(否则'='后的data被忽略)
    # 找到第一个'='的位置, 跳过所有连续'=', 之后为data块
    first_pad = raw_b64.find('=')
    if first_pad == -1:
        # 无填充, 尝试整体解码(FreeLB可能用自定义格式)
        try:
            allbytes = base64.b64decode(raw_b64)
        except Exception:
            return None
        if len(allbytes) < 4: return None
        n = struct.unpack('<I', allbytes[:4])[0]
        data = allbytes[4:4+n]
        arr = np.frombuffer(data, dtype=np.float64)
        return arr
    # 跳过所有连续'='
    hdr_end = first_pad
    while hdr_end < len(raw_b64) and raw_b64[hdr_end] == '=':
        hdr_end += 1
    header_b64 = raw_b64[:first_pad + (4 - (first_pad % 4)) % 4]
    # 修正: base64块大小必须是4的倍数, 包含填充
    header_b64 = raw_b64[:hdr_end]
    data_b64 = raw_b64[hdr_end:].strip()
    try:
        header_bytes = base64.b64decode(header_b64)
        data_bytes = base64.b64decode(data_b64)
    except Exception:
        return None
    if len(header_bytes) < 4: return None
    n = struct.unpack('<I', header_bytes[:4])[0]  # 数据字节数
    arr = np.frombuffer(data_bytes[:n], dtype=np.float64)
    return arr

def parse_vti(path):
    with open(path,'r') as f:
        content = f.read()
    # extent
    m = re.search(r'<Piece Extent="(\d+) (\d+) (\d+) (\d+) (\d+) (\d+)"', content)
    if not m: return None
    x0,x1,y0,y1 = int(m.group(1)),int(m.group(2)),int(m.group(3)),int(m.group(4))
    nx, ny = x1-x0+1, y1-y0+1
    origin = re.search(r'Origin="([^"]+)"', content)
    fields = {}
    for name in ['PHI','HX','HY','HMAG','PSI','Velocity','Density','Force']:
        arr = parse_dataarray(content, name)
        if arr is not None:
            fields[name] = arr.reshape(ny, nx) if arr.size==nx*ny else arr
    # 如果没有HMAG, 从HX/HY计算
    if 'HMAG' not in fields and 'HX' in fields and 'HY' in fields:
        hx = fields['HX']; hy = fields['HY']
        if hx.shape == hy.shape:
            fields['HMAG'] = np.sqrt(hx*hx + hy*hy)
    fields['_extent'] = (x0,x1,y0,y1)
    fields['_origin'] = origin.group(1) if origin else ''
    return fields

step = sys.argv[1] if len(sys.argv)>1 else '100'
files = sorted(glob.glob('vtkoutput/vtidata/rt2d_T%s_B*.vti'%step))
print('Found %d VTI files at T%s'%(len(files), step))

# 收集所有块的 H 场统计，找跨越界面的块
interface_blocks = []
for fp in files:
    d = parse_vti(fp)
    if d is None: continue
    phi = d.get('PHI')
    hx = d.get('HX'); hy = d.get('HY'); hmag = d.get('HMAG')
    if phi is None or hy is None: continue
    x0,x1,y0,y1 = d['_extent']
    # 检查是否跨越界面 (phi 既有接近1也有接近0)
    pmin, pmax = phi.min(), phi.max()
    has_interface = (pmin < 0.3) and (pmax > 0.7)
    # H 统计
    hy_mean = hy.mean()
    info = (os.path.basename(fp), d['_extent'], pmin, pmax, has_interface,
            hy.min(), hy.max(), hy.mean())
    if has_interface:
        interface_blocks.append((fp, d))

# 全局 H 统计
all_hy = []; all_hx=[]; all_phi=[]; all_hmag=[]
for fp in files:
    d = parse_vti(fp)
    if d is None: continue
    if d.get('HY') is not None:
        all_hy.append(d['HY'].ravel())
        all_hx.append(d['HX'].ravel())
        all_hmag.append(d['HMAG'].ravel())
        all_phi.append(d['PHI'].ravel())
all_hy = np.concatenate(all_hy) if all_hy else np.array([])
all_hx = np.concatenate(all_hx) if all_hx else np.array([])
all_hmag = np.concatenate(all_hmag) if all_hmag else np.array([])
all_phi = np.concatenate(all_phi) if all_phi else np.array([])

H0_expected = 0.005672  # 计算的 H0
print('\n=== 全局 H 场统计 ===')
print('H0 (expected) = %.6f' % H0_expected)
print('HY: min=%.6e max=%.6e mean=%.6e' % (all_hy.min(), all_hy.max(), all_hy.mean()))
print('HX: min=%.6e max=%.6e mean=%.6e' % (all_hx.min(), all_hx.max(), all_hx.mean()))
print('HMAG: min=%.6e max=%.6e mean=%.6e' % (all_hmag.min(), all_hmag.max(), all_hmag.mean()))

# 关键诊断: 重相(phi>0.9) vs 轻相(phi<0.1) 的 H 值
heavy = all_phi > 0.9
light = all_phi < 0.1
print('\n=== 界面两侧 H 场对比 (关键诊断) ===')
print('重相(phi>0.9, 铁磁流体, 顶部):')
print('  HY mean=%.6e  HMAG mean=%.6e  (理论 H_h=0.667*H0=%.6e)' % (all_hy[heavy].mean(), all_hmag[heavy].mean(), 0.667*H0_expected))
print('轻相(phi<0.1, 非磁介质, 底部):')
print('  HY mean=%.6e  HMAG mean=%.6e  (理论 H_l=1.333*H0=%.6e)' % (all_hy[light].mean(), all_hmag[light].mean(), 1.333*H0_expected))
print('\n  比值 HMAG_heavy/HMAG_light = %.4f  (理论 0.667/1.333=0.500)' % (all_hmag[heavy].mean()/all_hmag[light].mean() if all_hmag[light].mean()!=0 else 0))

# 界面附近梯度
interface_region = (all_phi > 0.3) & (all_phi < 0.7)
print('\n界面区(0.3<phi<0.7): HMAG std=%.6e (梯度指示, 大=有畸变)' % all_hmag[interface_region].std())

if all_hmag[heavy].mean()/max(all_hmag[light].mean(),1e-20) > 0.85:
    print('\n>>> 诊断: 重相/轻相 H 值接近相同，磁场未产生界面畸变！MF求解器有问题 <<<')
else:
    print('\n>>> 诊断: 界面两侧 H 值不同，MF求解器产生了场畸变（正常）<<<')
