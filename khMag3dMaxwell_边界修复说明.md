# khMag3dMaxwell 边界伪影修复说明

本备注记录 `examples/khMag3dMaxwell` 案例中“左右边界块状异常 / 周期缝伪影”的完整排查结论与解决方案。

---

## 1. 问题现象

- 云上正式网格（x∈[0,2L]、y∈[0,L]、z∈[0,4L]）在 `x=0` 与 `x=Ni-1` 附近出现明显的左右边界条带/块状异常。
- 异常在 `Bo_m=0`（无磁场）时同样存在，说明不是 Maxwell 磁力项本身造成的，而是 PF/NS 的周期边界数据错误。
- 之前仅修正周期 AABB 角柱、flag 设置顺序、角柱对角同步等，均不能消除该现象；云上重跑结果看起来与最初完全相同。

## 2. 根本原因

框架的 `FixedPeriodicBoundaryManager::Setup/SetupMPI` 在给周期 ghost 寻找 source 块时，使用：

```cpp
sblock.getSelfBlock().isInside(sourcePos)
```

并按块序号取第一个命中的块。

由于相邻块在 z 方向（块间 x/y 方向同理）共享 ghost 层，同一个 `sourcePos` 会同时落在两个 selfblock 内：
- 正确块的物理层；
- 相邻块的重叠 ghost 层。

按块序遍历时，管理器会先命中相邻块的 ghost 层，于是周期 ghost 的 PF/NS POP 被从**过期 ghost**复制，而不是从对侧物理列复制。随后 `Stream()` 把这些错误 POP 注入物理边界列，并在 `x=0 / x=Ni-1` 附近持续放大成左右边界伪影。

小网格实测证据（`PF_Per.Apply()` 后立即检查）：
- x 周期 ghost 与对侧物理列的 POP 间隙最大约 `7e-4`，且 `Apply()` 没有更新该 ghost；
- 该 ghost 的 `sourcePos` 同时命中 z 相邻块的 ghost 层和正确物理块，管理器选中了错误块。

## 3. 修复方案

只修改了 `examples/khMag3dMaxwell/khMag3dMaxwell.cpp`，未改动 `src/`。

在 PF/NS 每次碰撞后新增 `SyncPeriodicPops`，对周期 ghost 分布函数做一次显式、正确的全量同步：

1. 按 `GlobalBlockTable` 枚举所有 x/y 周期 ghost cell。
2. 将 ghost 坐标在 x/y 方向显式取模 wrap 到对侧**物理 cell**（z 不 wrap，z ghost 仍交给后面的壁面 ghost fix）。
3. 用 block 的 base 物理区间反查 source 块，绕开 `selfblock` 重叠导致的选块歧义。
4. 均匀块用 `(块长, min)` 编码做 O(1) 反查，避免每步 O(Nblock) 扫描。
5. MPI 通信按全局块条目对 `(ghostEntry, sourceEntry)` 分组：
   - source 在本地：打包物理 POP 发送；
   - ghost 在本地：接收并写回 ghost；
   - 非阻塞 iSend/iRecv + Waitall，两端使用一致的 pair tag，避免死锁。
6. 初始 setup 和每个时间步均调用，分别处理 PF 与 NS 两个 D3Q19 格子（fidx=0/100）。

要点：该同步放在 `PF_Per.Apply()` / `NS_Per.Apply()` 之后、`NormalFullCommunicate()` 与 `Stream()` 之前，因此保证进入 streaming 的所有周期 ghost POP 都严格等于对侧物理列。

## 4. 小网格验证

测试网格 `32×8×64`，`np=4`，`t=10`：

| 版本 | x=0 处 y 向伪差异 | y 周期缝最大差 |
|---|---|---|
| 修复前（仅 AABB + flag 顺序修正） | ~1.45e-1 | ~7e-2 |
| 本次全量周期 POP 同步 | 1.03e-5 | 2.2e-16 |

其他验证：
- `np=1/2/4` 均无死锁；
- `Bo_m=0` 与 `Bo_m=2547` 均通过；
- `64×16×128`、`np=4` 计时：修复前约 7.37s，修复后约 7.91s，额外开销约 7%。

说明：修复后 `x=0` 与 `x=Ni-1` 物理列的 PHI 仍可能有约 0.028 的差，这是初始界面
`h(x)=z0-A·cos(2πx/λ)` 在周期缝两侧本来就是两个相邻格子相位的差异，不是周期缝不连续；y 周期缝已严格到机器精度。

## 5. 使用方法

云上/集群上必须重新编译，不要使用旧的 `khMag3dMaxwell.exe`：

```bash
cd examples/khMag3dMaxwell
make -j
```

然后按原 ini 运行即可。`src/` 无需任何修改。

## 6. 相关结论

- 该问题属于框架周期边界管理器在“多块重叠 ghost 层处选 source 块歧义”的通用缺陷。
- 本修复在案例侧完成，不侵入框架；若后续其他多块周期案例出现同类边界条纹，可优先检查周期 ghost 是否被同步成对侧物理列，而不是相邻块的 ghost 层。

## 7. np=8/16 的 MPI Waitall 报错修复（最终结论）

在 `np=8` 和 `np=16` 下曾出现：

```text
Fatal error in PMPI_Waitall: Other MPI error
MPID_nem_lmt_shm_start_send ... unable to remove shared memory - unlink No such file or directory
```

完整定位（在本地复制 src 并对框架所有 `MPI_Waitall` 计数后确认）：

1. 报错里的 `count=176/286` 是**框架** `BlockFieldManager::AllMPINormalCommunicate`
   和 `BlockLatticeManager::MPINormalFullCommunicate` 的 `MPI_Waitall`
   request 数，不是案例代码的 Waitall。
2. 真正原因：MPICH CH3/nemesis 的 shared-memory LMT 传输在**双向大消息同时发生**时
   存在竞态——两端各自分配 shm region，又同时删除对方 lexicographically 较小的
   region，导致 `unlink` 时文件已被对端删掉（`MPIU_SHMW_Seg_detach` 报 ENOENT）。
   这和 request 数量本身无关；np=8 的块分布/消息大小最容易触发。
3. 案例代码本身还会叠加若干大消息，但不是报错的 Waitall。

最终修复（仍只改 `examples/khMag3dMaxwell/khMag3dMaxwell.cpp`）：

- 案例内所有自定义周期同步：
  - 发往同一 rank 的 pair 缓冲区合并成一个 MPI 消息，request 数降至最多 `nranks-1`；
  - 自己的 Waitall 使用显式 `MPI_Status` 数组；
  - 不再调用框架的 `NS_Per.Apply()` / `PF_Per.Apply()`（小网格 PHI 逐位一致，max diff = 0）。
- 在 `mpi().init()` 之前设置：
  ```cpp
  setenv("MPIR_CVAR_NOLOCAL", "1", 0);
  setenv("MPICH_NOLOCAL", "1", 0);
  setenv("MPIR_CVAR_CH3_NOLOCAL", "1", 0);   // MPICH <= 3.2
  ```
  这会让 MPICH 把全部进程按跨节点处理，**完全关闭 CH3/nemesis 的共享内存
  LMT 路径**，从而从根源上避开 `MPID_nem_lmt_shm_start_send` 竞态。
  该设置不改变计算结果，只改变 MPI 传输层选择。

验证：

- `32×8×64`、`np=4/8`：边界修复指标不变，PHI 与保留框架 Apply 的版本逐位一致；
- `512×8×512`、`np=8/16`：运行正常，无共享内存 LMT 错误；
- 本地开启 NOLOCAL 后，MPICH ch4 下大网格单步约慢 30%——这是避免该 MPI 内部
  竞态的代价；若集群 MPICH 升级到 4.x（ch4:ofi），可去掉这几行 setenv 再评估。

## 8. 重要提醒

- 集群上请务必重新 `make -j`，确认运行的是新生成的 `khMag3dMaxwell.exe`。
- 本修复不修改 `src/`，也不依赖特定 MPI 实现；聚合消息方案对 OpenMPI/MPICH 均适用。
