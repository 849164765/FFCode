#!/usr/bin/env python3
"""对比 Bom=41.18 和 Bom=0 的RT不稳定性界面演化"""
import base64, struct, glob, re, sys, os
import numpy as np

def parse_dataarray(content, name):
    pat = r'<DataArray[^>]*Name="%s"[^>]*>\s*(.*?)\s*</DataArray>' % name
    m = re.search(pat, content, re.DOTALL)
    if not m: return None
    raw_b64 = m.group(1).strip()
    first_pad = raw_b64.find('=')
    if first_pad == -1:
        try: allbytes = base64.b64decode(raw_b64)
        except: return None
        if len(allbytes)<4: return None
        n = struct.unpack('<I', allbytes[:4])[0]
        return np.frombuffer(allbytes[4:4+n], dtype=np.float64)
    hdr_end = first_pad
    while hdr_end < len(raw_b64) and raw_b64[hdr_end]=='=': hdr_end+=1
    header_b64 = raw_b64[:hdr_end]; data_b64 = raw_b64[hdr_end:].strip()
    try:
        hb = base64.b64decode(header_b64); db = base64.b64decode(data_b64)
    except: return None
    if len(hb)<4: return None
    n = struct.unpack('<I', hb[:4])[0]
    return np.frombuffer(db[:n], dtype=np.float64)

def parse_vti(path):
    with open(path,'r') as f: content=f.read()
    m = re.search(r'<Piece Extent="(\d+) (\d+) (\d+) (\d+) (\d+) (\d+)"', content)
    if not m: return None
    x0,x1,y0,y1=int(m.group(1)),int(m.group(2)),int(m.group(3)),int(m.group(4))
    nx,ny=x1-x0+1,y1-y0+1
    fields={}
    for name in ['PHI','HMAG','Velocity']:
        arr=parse_dataarray(content,name)
        if arr is not None: fields[name]=arr
    fields['_nx']=nx; fields['_ny']=ny; fields['_x0']=x0; fields['_y0']=y0
    return fields

def collect_field(step, basedir, field_name='PHI'):
    """收集所有块的指定字段, 返回全局2D数组
    块布局: 8x16 (blockXNum=8, blockYNum=16), 每块内部32x64 + 1 ghost cell/边 = 34x66
    块ID → (bid%8, bid//8) → 全局偏移 (bx*32, by*64)
    内部数据在VTI中位于 [1:65, 1:33] (跳过ghost cell)
    """
    import re as _re
    files=sorted(glob.glob(f'{basedir}/vtkoutput/vtidata/rt2d_T{step}_B*.vti'))
    if not files: return None
    # 全局网格 1024x256 (Nj x Ni)
    grid=np.full((1024, 256), np.nan)
    for fp in files:
        # 从文件名提取块ID: rt2d_T{step}_B{bid}.vti
        m=_re.search(r'_B(\d+)\.vti$', fp)
        if not m: continue
        bid=int(m.group(1))
        d=parse_vti(fp)
        if d is None: continue
        arr=d.get(field_name)
        if arr is None: continue
        nx,ny=d['_nx'],d['_ny']  # 34, 66 (含ghost)
        arr2d=arr.reshape(ny,nx)
        # 块在全局网格中的位置
        bx, by = bid % 8, bid // 8
        gx_start, gy_start = bx * 32, by * 64
        # 提取内部数据 (跳过ghost cell: 行1..64, 列1..32)
        interior = arr2d[1:65, 1:33]  # 64x32
        grid[gy_start:gy_start+64, gx_start:gx_start+32] = interior
    return grid

def find_interface(phi_grid):
    """找phi=0.5的等高线位置, 返回每列的interface y坐标"""
    ny,nx=phi_grid.shape
    interface_y=np.full(nx, np.nan)
    for i in range(nx):
        col=phi_grid[:,i]
        for j in range(ny-1):
            if (col[j]-0.5)*(col[j+1]-0.5) <= 0 and col[j]!=col[j+1]:
                # 线性插值
                interface_y[i]=j + (0.5-col[j])/(col[j+1]-col[j])
                break
    return interface_y

def analyze_step(step):
    """分析指定时间步的对比"""
    step_int = int(step)
    print(f'\n{"="*60}')
    print(f'  T={step_int} (t*={step_int/36203:.3f})')
    print(f'{"="*60}')
    
    # Bom=41.18
    phi_mag=collect_field(step, '/home/dlf/myCode/FerroRTinstability/examples/rt2d', 'PHI')
    # Bom=0
    phi_nomag=collect_field(step, '/home/dlf/myCode/FerroRTinstability/examples/rt2d_bom0', 'PHI')
    
    if phi_mag is None or phi_nomag is None:
        print('  [数据不完整, 跳过]')
        return
    
    # 找界面位置
    iy_mag=find_interface(phi_mag)
    iy_nomag=find_interface(phi_nomag)
    
    # 计算spike(重相向下穿透)和bubble(轻相向上运动)位置
    # 界面初始在 y=512 (2L), spike向下 = min(interface_y), bubble向上 = max(interface_y)
    valid_mag=iy_mag[~np.isnan(iy_mag)]
    valid_nomag=iy_nomag[~np.isnan(iy_nomag)]
    
    if len(valid_mag)==0 or len(valid_nomag)==0:
        print('  [界面未找到, 跳过]')
        return
    
    spike_mag=valid_mag.min()
    bubble_mag=valid_mag.max()
    spike_nomag=valid_nomag.min()
    bubble_nomag=valid_nomag.max()
    
    print(f'  Bom=41.18: spike_y={spike_mag:.1f}  bubble_y={bubble_mag:.1f}  amp={bubble_mag-spike_mag:.1f}')
    print(f'  Bom=0:     spike_y={spike_nomag:.1f}  bubble_y={bubble_nomag:.1f}  amp={bubble_nomag-spike_nomag:.1f}')
    print(f'  差异: spike Δ={spike_mag-spike_nomag:+.1f}  bubble Δ={bubble_mag-bubble_nomag:+.1f}')
    
    # 磁场效果: 如果Bom=41.18的spike更低/bubble更高, 说明磁场促进了RT不稳定性
    if spike_mag < spike_nomag - 5:
        print(f'  >>> 磁场使spike穿透更深 (促进了RT不稳定性) <<<')
    elif spike_mag > spike_nomag + 5:
        print(f'  >>> 磁场使spike穿透较浅 (抑制了RT不稳定性) <<<')
    else:
        print(f'  >>> 磁场效果不明显 <<<')

# 分析所有可用的时间步
if __name__ == '__main__':
    steps=sys.argv[1:] if len(sys.argv)>1 else ['4000','8000','16000','32000','64000','80000']
    for s in steps:
        analyze_step(s)
