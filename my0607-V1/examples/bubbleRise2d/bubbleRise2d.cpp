// bubbleRise2d.cpp
// Bubble rising under buoyancy — phase field + Navier-Stokes (D2Q9)
// Two lattices sharing the same geometry:
//   PfLat: Allen-Cahn phase field (phi, grad, normal, chem potential)
//   NsLat: Navier-Stokes flow (rho, u, force)

#include "freelb.h"
#include "freelb.hh"
#include "pf/pf.h"

using T = FLOAT;
using LatSet = D2Q9<T>;

int Ni, Nj;
T Cell_Len, Interface_Width, Mobility;
int BlockCellLen, Thread_Num;
T Droplet_Radius, Droplet_CenterX, Droplet_CenterY;
T Kappa, Beta, Omega_phi;
int MaxStep, OutputStep;

void readParam() {
  iniReader param_reader("bubbleRise2d.ini");
  Thread_Num = param_reader.getValue<int>("parallel", "thread_num");
  Ni = param_reader.getValue<int>("Mesh", "Ni");
  Nj = param_reader.getValue<int>("Mesh", "Nj");
  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");
  BlockCellLen = param_reader.getValue<int>("Mesh", "BlockCellLen");
  Droplet_Radius = param_reader.getValue<T>("Droplet", "Radius");
  Droplet_CenterX = param_reader.getValue<T>("Droplet", "CenterX");
  Droplet_CenterY = param_reader.getValue<T>("Droplet", "CenterY");
  Interface_Width = param_reader.getValue<T>("Phase_Field", "Interface_Width");
  Mobility = param_reader.getValue<T>("Phase_Field", "Mobility");
  MaxStep = param_reader.getValue<int>("Simulation_Settings", "TotalStep");
  OutputStep = param_reader.getValue<int>("Simulation_Settings", "OutputStep");

  Kappa = T(3) * T(0.072) * Interface_Width * T(0.5);
  Beta  = T(12) * T(0.072) / Interface_Width;
  Omega_phi = T(1) / (T(0.5) + Mobility / LatSet::cs2);
  pf::Kappa<T> = Kappa;
  pf::Beta<T> = Beta;
  pf::Omega_phi<T> = Omega_phi;

  MPI_RANK(0) {
    std::cout << "Bubble Rise 2D | " << Ni << "x" << Nj
              << " | K=" << Kappa << " B=" << Beta << " w=" << Omega_phi << std::endl;
  }
}

int main(int argc, char* argv[]) {
  constexpr std::uint8_t VoidFlag = 1;
  constexpr std::uint8_t BulkFlag = 2;
  constexpr std::uint8_t BouncebackFlag = 8;

  mpi().init(&argc, &argv);
  readParam();

  // Converters
  BaseConverter<T> baseConv(LatSet::cs2);
  baseConv.Converter(T(1), T(1), T(1), T(Ni), T(1), T(0.1));
  PhaseFieldConverter<T> phiConv(baseConv);
  phiConv.Converter(Interface_Width, Mobility, T(0.072));

  // Geometry
  AABB<T, 2> domain(T(0), T(0), T(Ni)*Cell_Len, T(Nj)*Cell_Len);
  BlockGeometryHelper2D<T> gh(Ni, Nj, domain, Cell_Len, BlockCellLen);
  gh.CreateBlocks(1, mpi().getSize());
  gh.AdaptiveOptimization(mpi().getSize());
  gh.LoadBalancing(mpi().getSize());
  BlockGeometry2D<T> Geo(gh);

  // Flags
  BlockFieldManager<FLAG, T, LatSet::d> FlagFM(Geo, VoidFlag);
  FlagFM.forEach(domain, [&](FLAG& field, std::size_t id) { field.SetField(id, BulkFlag); });
  FlagFM.template SetupBoundary<LatSet>(domain, BulkFlag);

  // PF lattice
  using PfFields = TypePack<PHI<T>, POP<T,9>, VELOCITY<T,2>,
                            GRAD<T,2>, NORMAL<T,2>, INTERFACEWIDTH<T>>;
  using PfCell = Cell<T, LatSet, PfFields>;
  ValuePack PfInit(phiConv.getLatRhoInit(), T(0), Vector<T,2>{},
                   Vector<T,2>{}, Vector<T,2>{}, Interface_Width);
  BlockLatticeManager<T, LatSet, PfFields> PfLat(Geo, PfInit, phiConv);

  // NS lattice
  using NsFields = TypePack<RHO<T>, POP<T,9>, VELOCITY<T,2>, FORCE<T,2>>;
  using NsCell = Cell<T, LatSet, NsFields>;
  ValuePack NsInit(baseConv.getLatRhoInit(), T(0), Vector<T,2>{}, Vector<T,2>{});
  BlockLatticeManager<T, LatSet, NsFields> NsLat(Geo, NsInit, baseConv);

  // Init phase field (tanh bubble) - following simpledrop2d pattern
  {
    auto& pfPhi = PfLat.getField<PHI<T>>();
    Vector<T,2> c(Droplet_CenterX, Droplet_CenterY);
    for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
      const auto& block = Geo.getBlock(blockid);
      auto& blockPhi = pfPhi.getBlockField(blockid);
      int overlap = 0;
      const auto& proj = block.getProjection();
      for (int j = overlap; j < block.getNy() - overlap; ++j) {
        for (int i = overlap; i < block.getNx() - overlap; ++i) {
          std::size_t id = j * proj[1] + i;
          T dist = ((block.getMin() + Vector<T,2>(T(i), T(j)) * block.getVoxelSize()) - c).getnorm();
          blockPhi.get(id) = T(0.5)*(T(1) - std::tanh(T(2)*(dist - Droplet_Radius)/Interface_Width));
        }
      }
    }
  }
  // Init PF populations: f_k = w_k * phi (rest equilibrium) - following simpledrop2d pattern
  {
    auto& phiField = PfLat.getField<PHI<T>>();
    for (int block_idx = 0; block_idx < Geo.getBlockNum(); ++block_idx) {
      auto& blockLat = PfLat.getBlockLat(block_idx);
      auto& blockPhiField = phiField.getBlockField(block_idx);
      const auto& block = Geo.getBlock(block_idx);
      int overlap = 0;
      const auto& proj = block.getProjection();
      for (int j = overlap; j < block.getNy() - overlap; ++j) {
        for (int i = overlap; i < block.getNx() - overlap; ++i) {
          std::size_t id = j * proj[1] + i;
          PfCell cell(id, blockLat);
          T phi = blockPhiField.get(id);
          for (unsigned int k = 0; k < LatSet::q; ++k)
            cell[k] = latset::w<LatSet>(k) * phi;
        }
      }
    }
  }
  // Init NS velocity to zero
  NsLat.getField<VELOCITY<T,2>>().forEach(
    domain, FlagFM, BulkFlag,
    [&](auto& field, std::size_t id) { field.SetField(id, Vector<T,2>{}); });

  // Periodic top/bottom
  AABB<T,2> bb(Vector<T,2>(T(0), T(-1)*Cell_Len), Vector<T,2>(T(Ni)*Cell_Len, T(0)));
  AABB<T,2> tb(Vector<T,2>(T(0), T(Nj)*Cell_Len), Vector<T,2>(T(Ni)*Cell_Len, T(Nj+1)*Cell_Len));
  FixedPeriodicBoundaryManager<BlockLatticeManager<T,LatSet,PfFields>,
                               BlockFieldManager<FLAG,T,2>>
    PfPer("PfPer", PfLat, FlagFM, std::uint8_t(4), VoidFlag);
  PfPer.Setup(bb, NbrDirection::YN, tb, NbrDirection::YP);
  FixedPeriodicBoundaryManager<BlockLatticeManager<T,LatSet,NsFields>,
                               BlockFieldManager<FLAG,T,2>>
    NsPer("NsPer", NsLat, FlagFM, std::uint8_t(4), VoidFlag);
  NsPer.Setup(bb, NbrDirection::YN, tb, NbrDirection::YP);

  // Bounce-back left/right
  AABB<T,2> lw(Vector<T,2>(T(-1)*Cell_Len, T(0)), Vector<T,2>(T(0), T(Nj)*Cell_Len));
  AABB<T,2> rw(Vector<T,2>(T(Ni)*Cell_Len, T(0)), Vector<T,2>(T(Ni+1)*Cell_Len, T(Nj)*Cell_Len));
  FlagFM.forEach(lw, [&](FLAG& field, std::size_t id) { field.SetField(id, BouncebackFlag); });
  FlagFM.forEach(rw, [&](FLAG& field, std::size_t id) { field.SetField(id, BouncebackFlag); });

  using PfBB = BBLikeFixedBlockBdManager<bounceback::normal<PfCell>,
      BlockLatticeManager<T,LatSet,PfFields>, BlockFieldManager<FLAG,T,2>>;
  using NsBB = BBLikeFixedBlockBdManager<bounceback::anti_pressure<NsCell>,
      BlockLatticeManager<T,LatSet,NsFields>, BlockFieldManager<FLAG,T,2>>;
  BlockBoundaryManager PfBd(new PfBB("PfBB", PfLat, FlagFM, BouncebackFlag));
  BlockBoundaryManager NsBd(new NsBB("NsBB", NsLat, FlagFM, BouncebackFlag));

  // Tasks
  using PfGrad = tmp::TaskSelector<tmp::TupleWrapper<
      tmp::Key_TypePair<BulkFlag, pf::PFGradient2D<PfCell>>>, std::uint8_t, PfCell>;
  using PfNorm = tmp::TaskSelector<tmp::TupleWrapper<
      tmp::Key_TypePair<BulkFlag, pf::PFNormal2D<PfCell>>>, std::uint8_t, PfCell>;
  using PfChem = tmp::TaskSelector<tmp::TupleWrapper<
      tmp::Key_TypePair<BulkFlag, pf::PFChemPotential2D<PfCell>>>, std::uint8_t, PfCell>;
  using PfColl = tmp::TaskSelector<tmp::TupleWrapper<
      tmp::Key_TypePair<BulkFlag, collision::BGKSource<
          equilibrium::FirstOrder<PfCell>, NORMAL<T,2>, true>>>, std::uint8_t, PfCell>;
  using NsForce = tmp::TaskSelector<tmp::TupleWrapper<
      tmp::Key_TypePair<BulkFlag, force::Guo<FORCE<T,2>>>>, std::uint8_t, NsCell>;
  using NsColl = tmp::TaskSelector<tmp::TupleWrapper<
      tmp::Key_TypePair<BulkFlag, collision::BGKSource<
          equilibrium::SecondOrder<NsCell>, void, false>>>, std::uint8_t, NsCell>;

  // VTK
  vtmwriter::ScalarWriter PhiW("PHI", PfLat.getField<PHI<T>>());
  vtmwriter::VectorWriter VelW("VELOCITY", NsLat.getField<VELOCITY<T,2>>());
  vtmwriter::vtmWriter<T,2> MW("bubbleRise2d", Geo);
  MW.addWriterSet(PhiW, VelW);

  // Init comm
  PfPer.Apply(); NsPer.Apply();
  PfLat.NormalCommunicate(); NsLat.NormalCommunicate();

  Timer MainLoopTimer;
  T grav_acc = T(-9810) / baseConv.Conv_Acc;
  MW.WriteBinary(MainLoopTimer());

  // Main loop
  while (MainLoopTimer() < MaxStep) {
    ++MainLoopTimer;
    // PF: grad, normal, chem pot
    PfLat.template ApplyCellDynamics<PfGrad>(FlagFM);
    PfLat.getField<GRAD<T,2>>().Communicate();
    PfLat.template ApplyCellDynamics<PfNorm>(FlagFM);
    PfLat.getField<NORMAL<T,2>>().Communicate();
    PfLat.template ApplyCellDynamics<PfChem>(FlagFM);

    // Cross-coupling: density + forces
    {
      auto& pfPhi = PfLat.getField<PHI<T>>();
      auto& pfMu = PfLat.getField<INTERFACEWIDTH<T>>();
      auto& pfGr = PfLat.getField<GRAD<T,2>>();
      auto& nsRho = NsLat.getField<RHO<T>>();
      auto& nsFr = NsLat.getField<FORCE<T,2>>();
      T rl=T(1000), rg=T(1), dr=rl-rg;
      for (int ib = 0; ib < Geo.getBlockNum(); ++ib) {
        auto& bf_phi=pfPhi.getBlockField(ib);
        auto& bf_mu=pfMu.getBlockField(ib);
        auto& bf_gr=pfGr.getBlockField(ib);
        auto& bf_rho=nsRho.getBlockField(ib);
        auto&bf_fr=nsFr.getBlockField(ib);
        for (std::size_t id = Geo.getBlock(ib).getVoxNum(); id-- > 0;) {
          T phi = bf_phi.get(id);
          T rho = rg + dr * phi;
          bf_rho.get(id) = rho;
          T mu = bf_mu.get(id);
          auto gr = bf_gr.get(id);
          bf_fr.get(id) = Vector<T,2>(mu*gr[0], mu*gr[1] + (rho-rl)*grav_acc);
        }
      }
    }
    NsLat.getField<FORCE<T,2>>().Communicate();
    NsLat.getField<RHO<T>>().Communicate();

    // NS + PF collision
    NsLat.template ApplyCellDynamics<NsForce>(FlagFM);
    NsLat.template ApplyCellDynamics<NsColl>(FlagFM);
    PfLat.template ApplyCellDynamics<PfColl>(FlagFM);

    // BC + periodic
    PfBd.Apply(MainLoopTimer()); NsBd.Apply(MainLoopTimer());
    NsPer.Apply(); PfPer.Apply();

    // Stream + comm
    PfLat.Stream(); NsLat.Stream();
    PfLat.NormalCommunicate(); NsLat.NormalCommunicate();

    // Update phi from pop
    { auto& pfPhi = PfLat.getField<PHI<T>>();
      for (int ib = 0; ib < Geo.getBlockNum(); ++ib) {
        auto& blk = Geo.getBlock(ib);
        auto& blat = PfLat.getBlockLat(ib);
        auto& bf = pfPhi.getBlockField(ib);
        for (std::size_t id = blk.getVoxNum(); id-- > 0;) {
          PfCell cell(id, blat);
          T p = 0; for (unsigned int k = 0; k < LatSet::q; ++k) p += cell[k];
          bf.get(id) = p;
        }
      }
    }
    PfLat.getField<PHI<T>>().Communicate();

    // Copy NS velocity -> PF
    { auto& nv = NsLat.getField<VELOCITY<T,2>>();
      auto& pv = PfLat.getField<VELOCITY<T,2>>();
      for (int ib = 0; ib < Geo.getBlockNum(); ++ib) {
        auto& bN = nv.getBlockField(ib);
        auto& bP = pv.getBlockField(ib);
        for (std::size_t id = Geo.getBlock(ib).getVoxNum(); id-- > 0;)
          bP.get(id) = bN.get(id);
      }
    }

    if (MainLoopTimer() % OutputStep == 0) {
      if (mpi().getRank() == 0)
        std::cout << "Step " << MainLoopTimer() << "/" << MaxStep << std::endl;
      MW.WriteBinary(MainLoopTimer());
    }
  }
  MW.WriteBinary(MainLoopTimer());
  if (mpi().getRank() == 0) std::cout << "Done." << std::endl;
  return 0;
}
