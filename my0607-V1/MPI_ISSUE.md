# MPI并行SegFault问题说明（GDB验证·最终版）

## 现象

`bubble2d.cpp` 在 256×512 网格配置下，多进程运行时发生段错误：

```
===================================================================================
=   BAD TERMINATION OF ONE OF YOUR APPLICATION PROCESSES
=   PID 3328935 RUNNING AT DENG
=   EXIT CODE: 9
=   CLEANING UP REMAINING PROCESSES
===================================================================================
YOUR APPLICATION TERMINATED WITH THE EXIT STRING: Segmentation fault (signal 11)
```

## 实测矩阵（exit 139 = SIGSEGV）

| Ni | Nj | Cell_Len | BlockCellLen | np | 结果 |
|----|----|----------|--------------|----|------|
| 256 | 512 | 1.0 | 32 | 1 | ✓ OK |
| 256 | 512 | 1.0 | 32 | 2 | ✗ exit=139 |
| 256 | 512 | 1.0 | 32 | 3 | ✗ exit=139 |
| 256 | 512 | 1.0 | 32 | 4 | ✗ exit=139 |
| 256 | 512 | 1.0 | 32 | 5 | ✗ exit=139 |
| 256 | 512 | 1.0 | 32 | 7 | ✗ exit=139 |
| 256 | 512 | 1.0 | 32 | 8 | ✗ exit=139 |
| 256 | 512 | 1.0 | 32 | 9 | ✓ OK（divide bug碰巧规避） |
| 200 | 400 | 1.0 | 32 | 4 | ✗ exit=139（GDB确认同一位置） |

**关键结论：整除性（np能否整除block数）不是崩溃条件。** np=3,5,7这些不能被整除的进程数也崩溃。崩溃的共同条件是 np>=2 —— 只要存在多块划分，就会有块在域内部有ghost面。

---

## GDB 确证：崩溃定位

```
$ mpiexec -n 4 gdb -batch -ex run -ex "bt full" -ex quit ./bubble2d.exe

Program received signal SIGSEGV, Segmentation fault.

#0  ff::FF2D::apply(cell) at src/ff/ff2d.hh:14
    -> T phi_i = cell.getNeighbor(i).template get<GenericRho>();
    i = 4            // D2Q9方向4 = 速度(0,-1) = 南方向

#1  tmp::SelectTask::execute(flag=2, cell) at src/utils/tmp.h:215
#2  tmp::TaskSelector::Execute(flag=2, cell) at src/utils/tmp.h:225
#3  BlockLattice::ApplyCellDynamics(flagarr) at src/data_struct/block_lattice.hh:68
    id = 2           // 局部cell id（扩展block坐标位置(2,0)，即南边界ghost层）
    voxNum = 33540   // 扩展block总cell数 = 258×130 (256×128 base, overlap=1)
```

### 证据摘要

| 证据 | 含义 |
|------|------|
| `#0 ff2d.hh:14` | 崩溃在ferrofluid计算核，**不在MPI通信代码** |
| `id=2, flag=BulkFlag(2)` | 扩展block南边界ghost cell被错误标记为体cell |
| `i=4 (南方向)` | `getNeighbor(4)` = `id + DeltaIndex[4]` = `2 + (-Nx)` = 负值 |
| `voxNum=33540` | 循环遍历了扩展block的**全部cell**（含ghost层） |

---

## 根因链：不是MPI通信bug，是计算核越界访问

### 步骤1：循环范围过大 — `src/data_struct/block_lattice.hh:53-71`

```cpp
template <typename CELLDYNAMICS, typename ArrayType>
void BlockLattice::ApplyCellDynamics(const ArrayType& flagarr) {
  Cell cell(0, *this);
  const std::size_t voxNum = this->getVoxNum();   // 行55: getN() = 扩展block的全部cell数
  for (std::size_t id = 0; id < voxNum; ++id) {   // 行62: 遍历 0..voxNum-1 = ALL cells含ghost
    cell.setId(id);                                // 行67
    CELLDYNAMICS::Execute(flagarr[id], cell);      // 行68: 对ghost层cell也执行任务
  }
}
```

`getVoxNum()` 在 `_VOX_ENABLED` 未定义时返回 `getN()` = 扩展block总数（base + overlap×2）。
对 `np=4` 的 256×128 block：extended N = 258×130 = **33540**。

### 步骤2：内部block的ghost cell被标记为BulkFlag — `examples/bubble2d/bubble2d.cpp:178-181`

```cpp
BlockFieldManager<FLAG, T, LatSet::d> FlagFM(Geo, VoidFlag);  // 行178: 默认VoidFlag=1
FlagFM.forEach(domain,                                          // 行179: 全局域(0,0)-(256,512)
  [&](FLAG& field, std::size_t id) { field.SetField(id, BulkFlag); }); // 行180: 域内→BulkFlag=2
```

`domain` = 全局物理域。对内部rank的block（如np=4的rank 1），扩展block南边界ghost cell物理中心 y=127.5 → 域内 → `BulkFlag=2` → 任务执行。

### 步骤3：内核无边界检查访问邻居 — `src/ff/ff2d.hh:13-14`

```cpp
void FF2D::apply(CELL& cell) {
  for (unsigned int i = 1; i < LatSet::q; ++i) {              // 行13: 遍历8个方向
    T phi_i = cell.getNeighbor(i).template get<GenericRho>(); // 行14: 无边界检查！
```

`getNeighbor()` = `src/data_struct/cell.h:146-148`：`Cell(Id + DeltaIndex[i], Lat)`。
D2Q9 DeltaIndex[4]（南）= `-Nx` = `-258`：
```
getNeighbor(4) = 2 + (-258) = -256
```
`Array[-256]` → **SIGSEGV**

---

## 全部 `ApplyCellDynamics` / `ApplyInnerCellDynamics` 函数梳理

代码库中共有 **3 个类**定义了此类函数，分布在 **2 个文件**中。

### 一、BlockLattice 层（3个实现）— `src/data_struct/block_lattice.hh`

#### 1. `ApplyCellDynamics(const ArrayType& flagarr)` — 行53-71 **【崩溃的直接来源】**

```cpp
const std::size_t voxNum = this->getVoxNum();       // 行55: getN() = 扩展block全部cell
for (std::size_t id = 0; id < voxNum; ++id) {       // 行62: 遍历0..N-1 = ALL cells含ghost
  CELLDYNAMICS::Execute(flagarr[id], cell);          // 行68: ghost边界cell也执行任务
}
```
- **问题**：遍历全部扩展cell，含ghost层边界cell
- **bubble2d调用**：行399（PFFFTask → FF2D, Laplacian, ChemPot）、行416（NS task）、行425（PF collision）

#### 2. `ApplyCellDynamics()` — 行75-92

```cpp
const std::size_t voxNum = this->getVoxNum();       // 行77: 同1号
for (std::size_t id = 0; id < voxNum; ++id) {       // 行84
  #ifdef _VOX_ENABLED
  cell.setId(map[id]);                                // 行86: ⚠️ 缺少 CELLDYNAMICS::apply(cell) !
  #else
  cell.setId(id);
  CELLDYNAMICS::apply(cell);                         // 行89
  #endif
}
```
- **问题A**：遍历全部扩展cell（同1号）
- **问题B（独立bug）**：`_VOX_ENABLED` 分支中 `apply()` 调用被遗漏——只 `setId`，不执行任务
- bubble2d未直接调用，通过 `BlockLatticeManager::ApplyCellDynamics()` 间接调用

#### 3. `ApplyCellDynamics(const Genericvector<std::size_t>& Idx)` — 行96-105

- **无问题**：索引由调用者显式指定，不自动遍历ghost层

---

### 二、BlockLatticeManager 层（6个重载）— `src/data_struct/block_lattice.hh:725-804`

全部委托，间接继承bug：

| 重载 | 行号 | 委托到 | bubble2d调用 |
|------|------|--------|-------------|
| `(count, BFM)` | 725-738 | `BlockLattice::1号(flagarr)` | — |
| `(BFM)` | 742-752 | 同上 | — |
| `(count)` | 756-766 | `BlockLattice::2号()` | — |
| `()` | 770-777 | 同上 | — |
| `(count, blockids)` | 781-792 | `BlockLattice::3号(Idx)` | — |
| `(blockids)` | 796-804 | 同上 | — |

---

### 三、BlockLatManagerCoupling 层（3个实现）— `src/data_struct/block_lattice.h:433-516` **【遗漏】**

此类用于跨 lattice 耦合计算（如 PF→NS 力传递）。代码内联在头文件中（模板类，非头文件模板实现）。

#### 1. `ApplyCellDynamics(count, BFM)` — 行433-462

```cpp
const std::size_t voxNum = blocklat0.getVoxNum();   // 行445: 同1号模式
for (std::size_t id = 0; id < voxNum; ++id) {       // 行449
  CELLDYNAMICS::Execute(flagArray[id], cell0, cell1);// 行453/458
}
```
- **问题**：遍历全部扩展cell（同 BlockLattice::1号）
- **bubble2d调用**：行409（STCoupling 表面张力）、行412（GravCoupling 重力）
- **任务不访问邻居**（FFSurfaceTension2D/FFGravityForce2D 只读写当前cell），所以不崩溃，但处理ghost cell是无效计算

#### 2. `ApplyCellDynamics(BFM)` — 行464-490

- **问题**：同上（`voxNum` 遍历），bubble2d未使用

#### 3. `ApplyCellDynamics()` — 行492-516

```cpp
const std::size_t voxNum = blocklat0.getVoxNum();   // 行501
for (std::size_t id = 0; id < voxNum; ++id) {       // 行505
  #ifdef _VOX_ENABLED
  cell0.setId(map[id]);
  cell1.setId(map[id]);
  #else
  cell0.setId(id);
  cell1.setId(id);
  #endif
  DYNAMICS::apply(cell0, cell1);                    // 行513: 在#endif外面，正确！
}
```
- **问题**：`voxNum` 遍历（同1号），但 `apply()` 不在 `#ifdef` 内，不存在遗漏bug
- bubble2d未使用

---

### 四、ApplyInnerCellDynamics —— 正确的实现（对照组）

`ApplyInnerCellDynamics` 是设计上正确的版本，它**只用 `overlap` 边界限制循环范围**：

#### BlockLattice::ApplyInnerCellDynamics — `block_lattice.hh:120-184`

```cpp
// 行120-150: 带flag版本
for (int j = this->getOverlap(); j < this->getNy() - this->getOverlap(); ++j) {
  for (int i = this->getOverlap(); i < this->getNx() - this->getOverlap(); ++i) {
    CELLDYNAMICS::Execute(flagarr[id], cell);
  }
}
```

#### BlockLatManagerCoupling::ApplyInnerCellDynamics — `block_lattice.h:519-618`

同样使用 `[overlap, N-overlap)` 边界，正确处理ghost层。

**结论**：`ApplyCellDynamics` 应该是 `ApplyInnerCellDynamics` 的"循环方式"，但错误地使用 `getVoxNum()` 遍历了全部 cell。

#### freeSurface.h 中的 ApplyInnerCellDynamics 调用

`src/lbm/freeSurface.h` 有 20 个调用点（406-954行），全部使用 `ApplyInnerCellDynamics`，不受此 bug 影响。

---

### CUDA 内核

#### `CuDevApplyCellDynamicsKernel` — `src/data_struct/block_lattice.h:183-188`

```cpp
std::size_t idx = blockIdx.x * blockDim.x + threadIdx.x;  // grid-stride，无N边界检查
CELLDYNAMICS::Execute(flagarr->operator[](idx), cell);
```
- 新版kernel无边界检查（`idx < N`的保护被移除），同样存在越界风险

### 六、受影响范围总结

| 示例/基准 | 文件 | 多block? | 受影响? | 说明 |
|-----------|------|---------|---------|------|
| cavity2d | cavity2d.cpp | 否 | ✓ 安全 (单block) | ghost全在域外→VoidFlag→跳过 |
| cavity3d | cavity3d.cpp | 否 | ✓ 安全 (单block) | 同上 |
| simpledrop2d | simpledrop2d.cpp | 否 | ✓ 安全 (单block) | 同上 |
| simpledrop2dMRT | simpledrop2d.cpp | 否 | ✓ 安全 (单block) | 同上 |
| simpledrop2dV2 | simpledrop2d.cpp | 否 | ✓ 安全 (单block) | 同上 |
| simpledrop3d | simpledrop3d.cpp | 否 | ✓ 安全 (单block) | 同上 |
| simpledrop3dMRT | simpledrop3d.cpp | 否 | ✓ 安全 (单block) | 同上 |
| simpledrop3dV2 | simpledrop3d.cpp | 否 | ✓ 安全 (单block) | 同上 |
| **bubble2d** | bubble2d.cpp | 是(MPI) | ✗ **崩溃** | ApplyCellDynamics + FF2D::getNeighbor |
| **cazsblock2d** | cazsblock2d.cpp | 是 | ✗ **受影响** | ApplyCellDynamics 遍历ghost |
| **cavblock2dmpi** | cavblock2dmpi.cpp | 是(MPI) | ✗ **受影响** | 同上 |
| **cavblock3dmpi** | cavblock3dmpi.cpp | 是(MPI) | ✗ **受影响** | 同上 |
| cavity2d_cu | cavity2d.cu | 否 | ✓ 安全 (CUDA) | 单block |
| cavity3d_cu | cavity3d.cu | 否 | ✓ 安全 (CUDA) | 单block |

### 七、函数影响矩阵

| 函数 | 文件:行号 | voxNum bug | 额外bug | 崩溃? | bubble2d使用 |
|------|----------|-----------|---------|-------|-------------|
| `BL::ApplyCellDynamics(flagarr)` | block_lattice.hh:53-71 | ✗ | — | **YES** | 行399,416,425 |
| `BL::ApplyCellDynamics()` | block_lattice.hh:75-92 | ✗ | `_VOX_ENABLED`缺apply | **YES** | 间接 |
| `BL::ApplyCellDynamics(Idx)` | block_lattice.hh:96-105 | ✓ | — | NO | 间接 |
| `BLM::ApplyCellDynamics(count,BFM)` | block_lattice.hh:725-738 | 继承1号 | — | 继承 | — |
| `BLM::ApplyCellDynamics(BFM)` | block_lattice.hh:742-752 | 继承1号 | — | 继承 | — |
| `BLM::ApplyCellDynamics(count)` | block_lattice.hh:756-766 | 继承2号 | 继承 | 继承 | — |
| `BLM::ApplyCellDynamics()` | block_lattice.hh:770-777 | 继承2号 | 继承 | 继承 | — |
| **`Coupling::ApplyCellDynamics(count,BFM)`** | block_lattice.h:433-462 | ✗ | — | 不崩溃(任务安全) | 行409,412 |
| **`Coupling::ApplyCellDynamics(BFM)`** | block_lattice.h:464-490 | ✗ | — | 不崩溃(任务安全) | — |
| **`Coupling::ApplyCellDynamics()`** | block_lattice.h:492-516 | ✗ | — | 不崩溃(任务安全) | — |
| **`BL::ApplyInnerCellDynamics(flagarr)`** | block_lattice.hh:120-150 | ✓** | — | NO | — |
| **`BL::ApplyInnerCellDynamics()`** | block_lattice.hh:154-184 | ✓** | — | NO | — |
| `CuDevApplyCellDynamicsKernel` | block_lattice.h:183-188 | ✗ | 无idx<N检查 | 可能 | cavity2d_cu |

> `BL` = BlockLattice, `BLM` = BlockLatticeManager, `Coupling` = BlockLatManagerCoupling
> ✓** = `ApplyInnerCellDynamics` 不使用 getVoxNum()，而是用 `[overlap, N-overlap)` 边界循环——**正确的实现方式**

---

## 全部代码错误位置汇总

### 崩溃直接原因（框架层）

| 文件 | 行号 | 函数/代码 | 问题 |
|------|------|----------|------|
| `src/data_struct/block_lattice.hh` | 55, 62 | `BL::ApplyCellDynamics(flagarr)` | `getVoxNum()` = getN(); `for(0..voxNum)` 遍历含ghost的全部cell |
| `src/data_struct/block_lattice.hh` | 68 | 同上 | `CELLDYNAMICS::Execute(flagarr[id], cell)` ghohst边界cell也执行 |
| `src/data_struct/block_lattice.hh` | 77, 84 | `BL::ApplyCellDynamics()` | 同1号; `_VOX_ENABLED`分支缺`apply()`调用（行86） |
| `src/data_struct/block_lattice.h` | 445, 449 | `Coupling::ApplyCellDynamics(count,BFM)` | `getVoxNum()` 遍历ghost，前inline在.h中 |
| `src/data_struct/block_lattice.h` | 474, 478 | `Coupling::ApplyCellDynamics(BFM)` | 同上 |
| `src/data_struct/block_lattice.h` | 501, 505 | `Coupling::ApplyCellDynamics()` | 同上 |
| `src/ff/ff2d.hh` | 14 | `FF2D::apply` | `cell.getNeighbor(i)` 无边界检查 |
| `src/ff/ff2d.hh` | 47 | `FFLaplacian2D::apply` | 同上 |
| `src/data_struct/cell.h` | 146-148 | `Cell::getNeighbor(i)` | `Id + DeltaIndex[i]`，无范围clamp |
| `src/data_struct/block_lattice.h` | 183-188 | `CuDevApplyCellDynamicsKernel` | 新版无 `idx < N` 检查

### 标记错误（bubble2d特定）

| 文件 | 行号 | 问题 |
|------|------|------|
| `examples/bubble2d/bubble2d.cpp` | 179-180 | `FlagFM.forEach(domain, ...)` 使用全局域 → 内部block ghost cell获BulkFlag |

### 额外发现的 `AABB::divide` bug（框架层）

| 文件 | 行号 | 问题 |
|------|------|------|
| `src/geometry/basic_geometry.hh` | 147 | `int rY = Ny - remainderY;` 余数分配逻辑**反向** |
| `src/geometry/basic_geometry.hh` | 127-169 | 完整函数体 |

```
// 当前(错): rY = Ny-remainder  →  只有1个block被扩大，其余不变
// 正确:    rY = remainder      →  remainder个block应获得partY+1
// 例: 512÷9: part=56, remainder=8
//   当前: 1个57, 8个56 → 505行 (≠512)
//   正确: 8个57, 1个56 → 512行 (=512)
```

---

## 解决方案

### 方案A（框架级修复，推荐）：限制循环范围为非ghost区域

修改 `src/data_struct/block_lattice.hh:53-71`（1号）和 `src/data_struct/block_lattice.hh:75-92`（2号）：

```cpp
// ---- 当前代码 (1号, 行53-71) ----
template <typename CELLDYNAMICS, typename ArrayType>
void BlockLattice::ApplyCellDynamics(const ArrayType& flagarr) {
  Cell cell(0, *this);
  const std::size_t voxNum = this->getVoxNum();                          // 行55
  #ifdef _VOX_ENABLED
  const VoxelMap& map = this->getBlock().getVoxelMap();
  #endif
  for (std::size_t id = 0; id < voxNum; ++id) {                          // 行62
    #ifdef _VOX_ENABLED
    cell.setId(map[id]);
    CELLDYNAMICS::Execute(flagarr[map[id]], cell);
    #else
    cell.setId(id);
    CELLDYNAMICS::Execute(flagarr[id], cell);                             // 行68
    #endif
  }
}

// ---- 修复后 (行53-71) ----
template <typename CELLDYNAMICS, typename ArrayType>
void BlockLattice::ApplyCellDynamics(const ArrayType& flagarr) {
  Cell cell(0, *this);
  #ifdef _VOX_ENABLED
  const VoxelMap& map = this->getBlock().getVoxelMap();
  const std::size_t voxNum = map.getVoxNum();
  for (std::size_t vid = 0; vid < voxNum; ++vid) {
    cell.setId(map[vid]);
    CELLDYNAMICS::Execute(flagarr[map[vid]], cell);
  }
  #else
  const auto& baseblock = this->getBlock().getBaseBlock();
  const int overlap = this->getBlock().getOverlap();
  for (int j = overlap; j < baseblock.getNy() - overlap; ++j) {
    for (int i = overlap; i < baseblock.getNx() - overlap; ++i) {
      std::size_t id = this->getIndex(Vector<int, LatSet::d>{i, j});
      cell.setId(id);
      CELLDYNAMICS::Execute(flagarr[id], cell);
    }
  }
  #endif
}
```

**同时修复2号函数（行75-92）的第二个bug**：`_VOX_ENABLED` 分支缺少 `apply()` 调用：

```cpp
// ---- 修复后 (行75-92) ----
template <typename CELLDYNAMICS>
void BlockLattice::ApplyCellDynamics() {
  Cell cell(0, *this);
  #ifdef _VOX_ENABLED
  const VoxelMap& map = this->getBlock().getVoxelMap();
  const std::size_t voxNum = map.getVoxNum();
  for (std::size_t vid = 0; vid < voxNum; ++vid) {
    cell.setId(map[vid]);
    CELLDYNAMICS::apply(cell);    // ← 当前缺失此行！
  }
  #else
  const auto& baseblock = this->getBlock().getBaseBlock();
  const int overlap = this->getBlock().getOverlap();
  for (int j = overlap; j < baseblock.getNy() - overlap; ++j) {
    for (int i = overlap; i < baseblock.getNx() - overlap; ++i) {
      std::size_t id = this->getIndex(Vector<int, LatSet::d>{i, j});
      cell.setId(id);
      CELLDYNAMICS::apply(cell);
    }
  }
  #endif
}
```

### 方案B：修复 `AABB::divide` 余数分配

修改 `src/geometry/basic_geometry.hh:147`：
```cpp
// 当前: int rY = Ny - remainderY;  // 行147
// 修复: int rY = remainderY;       // remainderY个block得 partY+1
```
**注意**：单独修此bug会消除np=9的"碰巧规避"，导致np=9也崩溃。必须先实施方案A。

### 方案C：修复CUDA内核边界检查

恢复 `src/data_struct/block_lattice.h:183-188` 中的 `idx < N` 边界检查。

### 临时绕开

- `mpiexec -n 1` — 单进程
- `mpiexec -n 9` — **不可靠，依赖divide bug**，域覆盖不完整（505/512行）

---

## 验证命令

```bash
# GDB 崩溃定位
mpiexec -n 4 gdb -batch -ex run -ex "bt full" -ex quit ./bubble2d.exe

# 编译时添加地址消毒器（需修改Makefile FLAGS）
# FLAGS += -fsanitize=address -fno-omit-frame-pointer
```
