#!/usr/bin/env python3
# 诊断: 检查Force字段是否包含磁力贡献
# 对比: 重相(铁磁)磁力应向下(负y), 轻相磁力应为0
import base64, struct, glob, re, os, sys
import numpy as np

def parse_dataarray(content, name):
    pat = r'<DataArray[^>]*Name="%s"[^>]*>\s*(.*?)\s*</DataArray>' % name
    m = re.search(pat, content, re.DOTALL)
    if not m: return None
    raw_b64 = m.group(1).strip()
    first_pad = raw_b64.find('=')
    if first_pad == -1:
        try:
            allbytes = base64.b64decode(raw_b64)
        except Exception:
            return None
        if len(allbytes) < 4: return None
        n = struct.unpack('<I', allbytes[:4])[0]
        data = allbytes[4:4+n]
        arr = np.frombuffer(data, dtype=np.float64)
        return arr
    hdr_end = first_pad
    while hdr_end < len(raw_b64) and raw_b64[hdr_end] == '=':
        hdr_end += 1
    header_b64 = raw_b64[:hdr_end]
    data_b64 = raw_b64[hdr_end:].strip()
    try:
        header_bytes = base64.b64decode(header_b64)
        data_bytes = base64.b64decode(data_b64)
    except Exception:
        return None
    if len(header_bytes) < 4: return None
    n = struct.unpack('<I', header_bytes[:4])[0]
    arr = np.frombuffer(data_bytes[:n], dtype=np.float64)
    return arr

def parse_vti(path):
    with open(path,'r') as f:
        content = f.read()
    m = re.search(r'<Piece Extent="(\d+) (\d+) (\d+) (\d+) (\d+) (\d+)"', content)
    if not m: return None
    x0,x1,y0,y1 = int(m.group(1)),int(m.group(2)),int(m.group(3)),int(m.group(4))
    nx, ny = x1-x0+1, y1-y0+1
    fields = {}
    for name in ['PHI','HX','HY','HMAG','PSI','Velocity','Density','Force']:
        arr = parse_dataarray(content, name)
        if arr is not None:
            if name == 'Force' or name == 'Velocity':
                # Vector field: 2 components, stored as interleaved or contiguous
                if arr.size == nx*ny*2:
                    fields[name] = arr.reshape(ny, nx, 2)
                elif arr.size == nx*ny:
                    fields[name] = arr.reshape(ny, nx)
            else:
                fields[name] = arr.reshape(ny, nx) if arr.size==nx*ny else arr
    if 'HMAG' not in fields and 'HX' in fields and 'HY' in fields:
        hx = fields['HX']; hy = fields['HY']
        if hx.shape == hy.shape:
            fields['HMAG'] = np.sqrt(hx*hx + hy*hy)
    fields['_extent'] = (x0,x1,y0,y1)
    return fields

step = sys.argv[1] if len(sys.argv)>1 else '100'
files = sorted(glob.glob('vtkoutput/vtidata/rt2d_T%s_B*.vti'%step))
print('Found %d VTI files at T%s'%(len(files), step))

# 收集所有数据
all_phi=[]; all_hy=[]; all_hmag=[]; all_fx=[]; all_fy=[]
for fp in files:
    d = parse_vti(fp)
    if d is None: continue
    phi = d.get('PHI')
    hy = d.get('HY')
    hmag = d.get('HMAG')
    force = d.get('Force')
    if phi is None: continue
    all_phi.append(phi.ravel())
    if hy is not None: all_hy.append(hy.ravel())
    if hmag is not None: all_hmag.append(hmag.ravel())
    if force is not None:
        if force.ndim == 3:
            all_fx.append(force[:,:,0].ravel())
            all_fy.append(force[:,:,1].ravel())

all_phi = np.concatenate(all_phi) if all_phi else np.array([])
all_hy = np.concatenate(all_hy) if all_hy else np.array([])
all_hmag = np.concatenate(all_hmag) if all_hmag else np.array([])
all_fx = np.concatenate(all_fx) if all_fx else np.array([])
all_fy = np.concatenate(all_fy) if all_fy else np.array([])

H0 = 0.005672
g = 0.01**2 / 256  # gravity
rho_h = 3.0; rho_l = 1.0

print('\n=== Force 字段统计 ===')
if all_fx.size > 0:
    print('FX: min=%.6e max=%.6e mean=%.6e' % (all_fx.min(), all_fx.max(), all_fx.mean()))
    print('FY: min=%.6e max=%.6e mean=%.6e' % (all_fy.min(), all_fy.max(), all_fy.mean()))
else:
    print('No Force data found!')
    sys.exit(0)

# 分相统计
heavy = all_phi > 0.9
light = all_phi < 0.1
interface = (all_phi > 0.3) & (all_phi < 0.7)

print('\n=== 分相 Force 统计 ===')
print('重相(phi>0.9):')
if heavy.any() and all_fy.size == all_phi.size:
    print('  FX mean=%.6e  FY mean=%.6e' % (all_fx[heavy].mean(), all_fy[heavy].mean()))
    print('  (重力 FY_grav = rho_h*g = %.6e)' % (rho_h*g))
else:
    print('  (no data)')
print('轻相(phi<0.1):')
if light.any() and all_fy.size == all_phi.size:
    print('  FX mean=%.6e  FY mean=%.6e' % (all_fx[light].mean(), all_fy[light].mean()))
    print('  (重力 FY_grav = rho_l*g = %.6e)' % (rho_l*g))
else:
    print('  (no data)')
print('界面区(0.3<phi<0.7):')
if interface.any() and all_fy.size == all_phi.size:
    print('  FX mean=%.6e  FY mean=%.6e' % (all_fx[interface].mean(), all_fy[interface].mean()))
    print('  FY min=%.6e  FY max=%.6e' % (all_fy[interface].min(), all_fy[interface].max()))

# 关键诊断: 磁力方向
# 重相(铁磁)磁力应向下(FY<0, 同重力方向), 因为|H|从轻相H0到重相H0/2递减
# 轻相磁力应为0(chi=0)
print('\n=== 磁力诊断 ===')
if all_fy.size == all_phi.size and heavy.any():
    fy_heavy = all_fy[heavy].mean()
    fy_light = all_fy[light].mean() if light.any() else 0
    grav_heavy = rho_h * g
    grav_light = rho_l * g
    # 重相: FY_total = FY_grav + FY_mag → FY_mag = FY_total - FY_grav
    fy_mag_heavy = fy_heavy - grav_heavy
    fy_mag_light = fy_light - grav_light
    print('重相: FY_total=%.6e  FY_grav=%.6e  FY_mag(估算)=%.6e' % (fy_heavy, grav_heavy, fy_mag_heavy))
    print('轻相: FY_total=%.6e  FY_grav=%.6e  FY_mag(估算)=%.6e' % (fy_light, grav_light, fy_mag_light))
    if abs(fy_mag_heavy) < 1e-15:
        print('\n>>> 警告: 重相磁力接近0! 磁力未生效! <<<')
    elif fy_mag_heavy > 0:
        print('\n>>> 警告: 重相磁力为正(向上)! 方向错误, 应为负(向下)! <<<')
    else:
        print('\n>>> 磁力方向正确(向下, 促进RT), 大小=%.6e <<<' % abs(fy_mag_heavy))
        print('>>> 磁力/重力 = %.4f <<<' % (abs(fy_mag_heavy)/grav_heavy))
