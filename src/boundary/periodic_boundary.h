#pragma once

#include "boundary/basic_boundary.h"
#include "lbm/moment.ur.h"

#ifdef MPI_ENABLED
#include <map>
#include <tuple>
#endif


template <typename T, typename LatSet>
class FixedPeriodicBoundary final : public AbstractBoundary {
 private:
  std::vector<std::size_t> VCells;
  std::vector<std::size_t> RCells;
  std::vector<std::size_t> BdCells;
  PopLattice<T, LatSet> &Lat;
  Geometry<T, LatSet::d> &Geo;
  GenericArray<std::uint8_t> &Field;
  std::uint8_t BdCellFlag;
  std::uint8_t voidFlag;
  AABB<T, LatSet::d> &Box;
  AABB<T, LatSet::d> &FromBox;
  std::string _name = "FixedPeriodic";

 public:
  FixedPeriodicBoundary(PopLattice<T, LatSet> &lat, AABB<T, LatSet::d> &box,
                        AABB<T, LatSet::d> &frombox, std::uint8_t cellflag,
                        std::uint8_t voidflag = std::uint8_t(1))
      : Lat(lat), Geo(lat.getGeo()), Field(lat.getGeo().getGeoFlagField().getField()),
        Box(box), FromBox(frombox), BdCellFlag(cellflag), voidFlag(voidflag) {
    Vector<T, LatSet::d> boxext = Box.getExtension();
    Vector<int, 3> SIZE(1);
    if constexpr (LatSet::d == 2) {
      SIZE[0] = int(boxext[0] / Geo.getVoxelSize()) + 2;
      SIZE[1] = int(boxext[1] / Geo.getVoxelSize()) + 2;
    } else if constexpr (LatSet::d == 3) {
      SIZE[0] = int(boxext[0] / Geo.getVoxelSize()) + 2;
      SIZE[1] = int(boxext[1] / Geo.getVoxelSize()) + 2;
      SIZE[2] = int(boxext[2] / Geo.getVoxelSize()) + 2;
    }
    int size = SIZE[0] * SIZE[1] * SIZE[2];
    VCells.reserve(size);
    RCells.reserve(size);
    BdCells.reserve(size);
    Geo.forEachVoxel(Box, BdCellFlag, [this](int id) { BdCells.push_back(id); });
    Vector<T, LatSet::d> dist = FromBox.getCenter() - Box.getCenter();
    for (int i = 0; i < LatSet::d; ++i) {
      if (dist[i] < T(-(1e-6))) {
        dist[i] -= Geo.getVoxelSize();
      } else if (dist[i] > T(1e-6)) {
        dist[i] += Geo.getVoxelSize();
      }
    }
    AABB<T, LatSet::d> extBox = Box.getExtended(Geo.getVoxelSize());
    Geo.forEachVoxel(extBox, voidFlag, [this, &dist](int id) {
      if (Geo.template hasNeighborFlag<LatSet>(id, BdCellFlag)) {
        VCells.push_back(id);
        int idr = Geo.findCellId(Geo.getVoxel(id) + dist);
        RCells.push_back(idr);
      }
    });
  }

  void Apply() override {
    int size = VCells.size();
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
    for (int i = 0; i < size; ++i) {
      BasicPopCell<T, LatSet> vcell(VCells[i], Lat);
      BasicPopCell<T, LatSet> rcell(RCells[i], Lat);
      for (int k = 1; k < LatSet::q; ++k) vcell[k] = rcell[k];
    }
  }
  void getinfo() override {
    std::cout << std::setw(18) << std::left << _name << std::setw(10) << std::left << BdCells.size()
              << std::endl;
  }
  void UpdateRho() override {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
    for (std::size_t id : BdCells) {
      BasicPopCell<T, LatSet> cell(id, Lat);
      moment::Rho<T, LatSet>::apply(cell, Lat.getRho(id));
    }
  }
  void UpdateU() override {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
    for (std::size_t id : BdCells) {
      BasicPopCell<T, LatSet> cell(id, Lat);
      moment::Velocity<T, LatSet>::apply(cell, Lat.getVelocity(id));
    }
  }
};


template <typename BLOCKLATTICEMANAGER, typename BLOCKFIELDMANAGER>
class FixedPeriodicBoundaryManager : public AbstractBlockBoundary {
 public:
  using BLOCKLATTICE = typename BLOCKLATTICEMANAGER::BLOCKLATTICE;
  using CELL = typename BLOCKLATTICE::CellType;
  using T = typename BLOCKLATTICE::FloatType;
  using LatSet = typename BLOCKLATTICE::LatticeSet;
  static constexpr unsigned int D = BLOCKLATTICE::LatticeSet::d;
  using ArrayType = typename BLOCKFIELDMANAGER::array_type;
  using GenericRho = typename CELL::GenericRho;

 private:
  struct CrossBlockPair {
    std::size_t ghostBlockIdx;
    std::size_t ghostId;
    std::size_t sourceBlockIdx;
    std::size_t sourceId;
  };

  struct BoundaryBoxPair {
    AABB<T, D> box0;
    AABB<T, D> box1;
    Vector<T, D> dist01;
    Vector<T, D> dist10;
  };

#ifdef MPI_ENABLED
  struct PeriodicMPISend {
    int targetRank;
    int sourceBlockId;
    std::size_t sourceBlockIdx;
    std::vector<std::size_t> sourceIds;
  };

  struct PeriodicMPIRecv {
    int sourceRank;
    int sourceBlockId;
    std::vector<std::size_t> ghostBlockIdxs;
    std::vector<std::size_t> ghostIds;
  };
#endif

  std::string _name;
  std::vector<CrossBlockPair> CrossBlockPairs;
  std::vector<BoundaryBoxPair> BoundaryBoxPairs;
  std::uint8_t BdCellFlag;
  std::uint8_t voidFlag;
  BLOCKLATTICEMANAGER &LatMan;
  BLOCKFIELDMANAGER &BlockFManager;

#ifdef MPI_ENABLED
  std::vector<PeriodicMPISend> MPISends;
  std::vector<PeriodicMPIRecv> MPIRecvs;
#endif

 public:
  FixedPeriodicBoundaryManager(std::string name, BLOCKLATTICEMANAGER &lat,
                               BLOCKFIELDMANAGER &BlockFM, std::uint8_t cellflag,
                               std::uint8_t voidflag = std::uint8_t(1))
      : _name(name), BdCellFlag(cellflag), voidFlag(voidflag), LatMan(lat),
        BlockFManager(BlockFM) {}

  void Setup(const AABB<T, D>& box0, NbrDirection dir0,
             const AABB<T, D>& box1, NbrDirection dir1) {
    Vector<T, D> dist01 = box1.getCenter() - box0.getCenter();
    for (unsigned int d = 0; d < D; ++d) {
      if (dist01[d] < T(-(1e-6))) dist01[d] += LatMan.getGeo().getVoxelSize();
      else if (dist01[d] > T(1e-6)) dist01[d] -= LatMan.getGeo().getVoxelSize();
    }
    Vector<T, D> dist10 = box0.getCenter() - box1.getCenter();
    for (unsigned int d = 0; d < D; ++d) {
      if (dist10[d] < T(-(1e-6))) dist10[d] += LatMan.getGeo().getVoxelSize();
      else if (dist10[d] > T(1e-6)) dist10[d] -= LatMan.getGeo().getVoxelSize();
    }

#ifndef MPI_ENABLED
    for (std::size_t i = 0; i < LatMan.getGeo().getBlockNum(); ++i) {
      auto& blockLat = LatMan.getBlockLat(i);
      auto& block = blockLat.getBlock();
      const auto& flagarr = BlockFManager.getBlockField(i).getField(0);

      auto& baseblock = block.getBaseBlock();
      auto& selfblock = block.getSelfBlock();

      auto processBox = [&](const AABB<T, D>& box, const Vector<T, D>& dist) {
        if (!isOverlapped(selfblock, box)) return;
        block.forEach(box, flagarr, BdCellFlag, [&](std::size_t id) {
          if (!baseblock.isInside(block.getLoc_t(id))) {
            Vector<T, D> ghostPos = block.getLoc_t(id);
            Vector<T, D> sourcePos = ghostPos + dist;

            std::size_t sourceBlockIdx = i;
            std::size_t sourceId;

            if (selfblock.isInside(sourcePos)) {
              sourceId = block.getIndex_t(sourcePos);
            } else {
              bool found = false;
              for (std::size_t bi = 0; bi < LatMan.getGeo().getBlockNum(); ++bi) {
                auto& otherBlock = LatMan.getGeo().getBlock(bi);
                if (otherBlock.getSelfBlock().isInside(sourcePos)) {
                  sourceBlockIdx = bi;
                  sourceId = otherBlock.getIndex_t(sourcePos);
                  found = true;
                  break;
                }
              }
              if (!found) return;
            }

            CrossBlockPairs.push_back({i, id, sourceBlockIdx, sourceId});
          }
        });
      };

      processBox(box0, dist01);
      processBox(box1, dist10);
    }

    MPI_RANK(0) {
      std::size_t sameBlockCount = 0;
      std::size_t crossBlockCount = 0;
      for (const auto& pair : CrossBlockPairs) {
        if (pair.ghostBlockIdx == pair.sourceBlockIdx) sameBlockCount++;
        else crossBlockCount++;
      }
      std::cout << "[FixedPeriodicBoundaryManager] Total pairs: " << CrossBlockPairs.size()
                << " (same-block: " << sameBlockCount << ", cross-block: " << crossBlockCount << ")" << std::endl;
    }
#endif

    BoundaryBoxPairs.push_back({box0, box1, dist01, dist10});
  }

#ifdef MPI_ENABLED
  template <typename GeoHelper>
  void SetupMPI(GeoHelper& geoHelper) {
    const auto& globalGeo = geoHelper.getBlockGeometry();
    int myRank = mpi().getRank();

    std::map<int, std::vector<std::tuple<std::size_t, std::size_t, Vector<T, D>>>> sendMap;
    std::map<int, std::vector<std::tuple<std::size_t, std::size_t, int, Vector<T, D>>>> recvMap;

    for (const auto& bp : BoundaryBoxPairs) {
      auto processDirection = [&](const AABB<T, D>& ghostBox, const Vector<T, D>& dist) {
        for (std::size_t bi = 0; bi < globalGeo.getBlockNum(); ++bi) {
          const auto& hblock = globalGeo.getBlock(bi);
          const auto& hbaseblock = hblock.getBaseBlock();
          const auto& hselfblock = hblock.getSelfBlock();
          int hblockRank = geoHelper.whichRank(hblock.getBlockId());

          if (!isOverlapped(hselfblock, ghostBox)) continue;

          Vector<int, D> idx_min, idx_max;
          hblock.getLocIdxRange(ghostBox, idx_min, idx_max);

          auto pairGhostCell = [&](const Vector<T, D>& pos,
                                   const Vector<T, D>& sourcePos) {
            std::size_t foundSbi = globalGeo.getBlockNum();
            int foundRank = -1;
            for (std::size_t sbi = 0; sbi < globalGeo.getBlockNum(); ++sbi) {
              const auto& sblock = globalGeo.getBlock(sbi);
              if (sblock.getSelfBlock().isInside(sourcePos)) {
                foundSbi = sbi;
                foundRank = geoHelper.whichRank(sblock.getBlockId());
                break;
              }
            }
            if (foundSbi >= globalGeo.getBlockNum()) return;
            const auto& sblock = globalGeo.getBlock(foundSbi);

            if (hblockRank == myRank && foundRank == myRank) {
              std::size_t ghostBlockIdx = LatMan.getGeo().findBlockIndex(hblock.getBlockId());
              std::size_t ghostId = LatMan.getBlockLat(ghostBlockIdx).getBlock().getIndex_t(pos);
              std::size_t sourceBlockIdx = LatMan.getGeo().findBlockIndex(sblock.getBlockId());
              std::size_t sourceId = LatMan.getBlockLat(sourceBlockIdx).getBlock().getIndex_t(sourcePos);
              CrossBlockPairs.push_back({ghostBlockIdx, ghostId, sourceBlockIdx, sourceId});
            } else if (hblockRank == myRank && foundRank != myRank) {
              std::size_t ghostBlockIdx = LatMan.getGeo().findBlockIndex(hblock.getBlockId());
              std::size_t ghostId = LatMan.getBlockLat(ghostBlockIdx).getBlock().getIndex_t(pos);
              int sourceBlockGlobalId = sblock.getBlockId();
              recvMap[foundRank].push_back({ghostBlockIdx, ghostId, sourceBlockGlobalId, pos});
            } else if (hblockRank != myRank && foundRank == myRank) {
              std::size_t sourceBlockIdx = LatMan.getGeo().findBlockIndex(sblock.getBlockId());
              std::size_t sourceId = LatMan.getBlockLat(sourceBlockIdx).getBlock().getIndex_t(sourcePos);
              sendMap[hblockRank].push_back({sourceBlockIdx, sourceId, pos});
            }
          };

          if constexpr (D == 2) {
            for (int j = idx_min[1]; j <= idx_max[1]; ++j) {
              for (int i = idx_min[0]; i <= idx_max[0]; ++i) {
                const Vector<T, 2> pt = hblock.getMinCenter() + (hblock.getVoxelSize() * Vector<T, 2>{T(i), T(j)});
                if (!ghostBox.isInside(pt)) continue;
                if (hbaseblock.isInside(pt)) continue;
                pairGhostCell(Vector<T, D>(pt), Vector<T, D>(pt) + dist);
              }
            }
          } else if constexpr (D == 3) {
            for (int k = idx_min[2]; k <= idx_max[2]; ++k) {
              for (int j = idx_min[1]; j <= idx_max[1]; ++j) {
                for (int i = idx_min[0]; i <= idx_max[0]; ++i) {
                  const Vector<T, 3> pt = hblock.getMinCenter() + (hblock.getVoxelSize() * Vector<T, 3>{T(i), T(j), T(k)});
                  if (!ghostBox.isInside(pt)) continue;
                  if (hbaseblock.isInside(pt)) continue;
                  pairGhostCell(Vector<T, D>(pt), Vector<T, D>(pt) + dist);
                }
              }
            }
          }
        }
      };

      processDirection(bp.box0, bp.dist01);
      processDirection(bp.box1, bp.dist10);
    }

    // IMPORTANT: Do NOT deduplicate send/recv entries by ghost position.
    // Different blocks can have ghost cells at the same physical position
    // (because SelfBlocks overlap at block boundaries).  Deduplication by
    // `pos` would drop one block's ghost cell, leaving it stuck (unpaired).
    // Deduplication by (blockIdx, cellId) would break MPI symmetry (send
    // count ≠ recv count).  The only correct option is no deduplication at
    // all — the double Setup() call produces identical duplicates on both
    // sides, so MPI send/recv sizes always match.

    for (auto& [rank, cells] : sendMap) {
      std::map<std::size_t, std::vector<std::size_t>> blockGroups;
      for (auto& [blockIdx, cellId, pos] : cells) {
        blockGroups[blockIdx].push_back(cellId);
      }
      for (auto& [blockIdx, cellIds] : blockGroups) {
        int sourceBlockGlobalId = static_cast<int>(LatMan.getGeo().getBlock(blockIdx).getBlockId());
        MPISends.push_back({rank, sourceBlockGlobalId, blockIdx, std::move(cellIds)});
      }
    }

    for (auto& [rank, cells] : recvMap) {
      std::map<int, std::vector<std::tuple<std::size_t, std::size_t>>> blockGroups;
      for (auto& [blockIdx, cellId, srcBlockId, pos] : cells) {
        blockGroups[srcBlockId].push_back({blockIdx, cellId});
      }
      for (auto& [srcBlockId, cellPairs] : blockGroups) {
        std::vector<std::size_t> blockIdxs;
        std::vector<std::size_t> ghostCellIds;
        for (auto& [blockIdx, cellId] : cellPairs) {
          blockIdxs.push_back(blockIdx);
          ghostCellIds.push_back(cellId);
        }
        MPIRecvs.push_back({rank, srcBlockId, std::move(blockIdxs), std::move(ghostCellIds)});
      }
    }
  }
#endif

  void Apply(std::int64_t count) override { Apply(); }

  void Apply() override {
    for (std::size_t j = 0; j < CrossBlockPairs.size(); ++j) {
      auto& pair = CrossBlockPairs[j];
      auto& ghostBlockLat = LatMan.getBlockLat(pair.ghostBlockIdx);
      auto& sourceBlockLat = LatMan.getBlockLat(pair.sourceBlockIdx);
      CELL ghostCell(pair.ghostId, ghostBlockLat);
      CELL sourceCell(pair.sourceId, sourceBlockLat);
      for (unsigned int k = 0; k < LatSet::q; ++k) {
        ghostCell[k] = sourceCell[k];
      }
      ghostCell.template get<GenericRho>() = sourceCell.template get<GenericRho>();
    }

#ifdef MPI_ENABLED
    if (MPISends.empty() && MPIRecvs.empty()) return;

    constexpr int PERIODIC_TAG_BASE = 9000;

    std::vector<std::vector<T>> SendBuffers(MPISends.size());
    std::vector<MPI_Request> SendRequests;

    for (std::size_t i = 0; i < MPISends.size(); ++i) {
      auto& send = MPISends[i];
      auto& buffer = SendBuffers[i];
      buffer.resize(send.sourceIds.size() * (LatSet::q + 1));
      auto& sourceBlockLat = LatMan.getBlockLat(send.sourceBlockIdx);
      std::size_t bufidx = 0;
      for (std::size_t id : send.sourceIds) {
        CELL sourceCell(id, sourceBlockLat);
        for (unsigned int k = 0; k < LatSet::q; ++k) {
          buffer[bufidx++] = sourceCell[k];
        }
        buffer[bufidx++] = sourceCell.template get<GenericRho>();
      }
      MPI_Request request;
      mpi().iSend(buffer.data(), static_cast<int>(buffer.size()), send.targetRank, &request,
                  PERIODIC_TAG_BASE + send.sourceBlockId);
      SendRequests.push_back(request);
    }

    std::vector<std::vector<T>> RecvBuffers(MPIRecvs.size());
    std::vector<MPI_Request> RecvRequests;

    for (std::size_t i = 0; i < MPIRecvs.size(); ++i) {
      auto& recv = MPIRecvs[i];
      auto& buffer = RecvBuffers[i];
      buffer.resize(recv.ghostIds.size() * (LatSet::q + 1));
      MPI_Request request;
      mpi().iRecv(buffer.data(), static_cast<int>(buffer.size()), recv.sourceRank, &request,
                  PERIODIC_TAG_BASE + recv.sourceBlockId);
      RecvRequests.push_back(request);
    }

    MPI_Waitall(static_cast<int>(SendRequests.size()), SendRequests.data(), MPI_STATUSES_IGNORE);

    for (std::size_t i = 0; i < MPIRecvs.size(); ++i) {
      MPI_Wait(&RecvRequests[i], MPI_STATUS_IGNORE);
      auto& recv = MPIRecvs[i];
      const auto& buffer = RecvBuffers[i];
      std::size_t bufidx = 0;
      for (std::size_t j = 0; j < recv.ghostIds.size(); ++j) {
        auto& ghostBlockLat = LatMan.getBlockLat(recv.ghostBlockIdxs[j]);
        CELL ghostCell(recv.ghostIds[j], ghostBlockLat);
        for (unsigned int k = 0; k < LatSet::q; ++k) {
          ghostCell[k] = buffer[bufidx++];
        }
        ghostCell.template get<GenericRho>() = buffer[bufidx++];
      }
    }
#endif
  }

  template <typename OtherBFM>
  void ApplyField(OtherBFM& fieldFM) {
    for (std::size_t j = 0; j < CrossBlockPairs.size(); ++j) {
      auto& pair = CrossBlockPairs[j];
      auto& ghostBF = fieldFM.getBlockField(static_cast<int>(pair.ghostBlockIdx));
      auto& sourceBF = fieldFM.getBlockField(static_cast<int>(pair.sourceBlockIdx));
      ghostBF.get(pair.ghostId) = sourceBF.get(pair.sourceId);
    }
#ifdef MPI_ENABLED
    if (MPISends.empty() && MPIRecvs.empty()) return;

    constexpr int PERIODIC_FIELD_TAG_BASE = 9100;

    std::vector<std::vector<T>> SendBuffers(MPISends.size());
    std::vector<MPI_Request> SendRequests;

    for (std::size_t i = 0; i < MPISends.size(); ++i) {
      auto& send = MPISends[i];
      auto& buffer = SendBuffers[i];
      buffer.resize(send.sourceIds.size());
      auto& sourceBF = fieldFM.getBlockField(static_cast<int>(send.sourceBlockIdx));
      std::size_t bufidx = 0;
      for (std::size_t id : send.sourceIds) {
        buffer[bufidx++] = sourceBF.get(id);
      }
      MPI_Request request;
      mpi().iSend(buffer.data(), static_cast<int>(buffer.size()), send.targetRank, &request,
                  PERIODIC_FIELD_TAG_BASE + send.sourceBlockId);
      SendRequests.push_back(request);
    }

    std::vector<std::vector<T>> RecvBuffers(MPIRecvs.size());
    std::vector<MPI_Request> RecvRequests;

    for (std::size_t i = 0; i < MPIRecvs.size(); ++i) {
      auto& recv = MPIRecvs[i];
      auto& buffer = RecvBuffers[i];
      buffer.resize(recv.ghostIds.size());
      MPI_Request request;
      mpi().iRecv(buffer.data(), static_cast<int>(buffer.size()), recv.sourceRank, &request,
                  PERIODIC_FIELD_TAG_BASE + recv.sourceBlockId);
      RecvRequests.push_back(request);
    }

    MPI_Waitall(static_cast<int>(SendRequests.size()), SendRequests.data(), MPI_STATUSES_IGNORE);

    for (std::size_t i = 0; i < MPIRecvs.size(); ++i) {
      MPI_Wait(&RecvRequests[i], MPI_STATUS_IGNORE);
      auto& recv = MPIRecvs[i];
      const auto& buffer = RecvBuffers[i];
      std::size_t bufidx = 0;
      for (std::size_t j = 0; j < recv.ghostIds.size(); ++j) {
        auto& ghostBF = fieldFM.getBlockField(static_cast<int>(recv.ghostBlockIdxs[j]));
        ghostBF.get(recv.ghostIds[j]) = buffer[bufidx++];
      }
    }
#endif
  }
};
