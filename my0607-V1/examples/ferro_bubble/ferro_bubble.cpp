// ferro_bubble.cpp
// Ferrofluid bubble rise simulation — three-solver FDLBM
// Paper: Guo et al., Phys. Fluids 37, 022148 (2025), Section 3.3

#include "freelb.h"
#include "freelb.hh"

using T = FLOAT;
using LatSet9 = D2Q9<T>;
using LatSet5 = D2Q5<T>;

// ================================================================
// Global parameters (read from INI)
// ================================================================
int Ni, Nj, BlockCellLen, Thread_Num;
T Cell_Len;
T Drop_Radius, Drop_CX, Drop_CY;
T Rho_L, Rho_H, Eta_L, Eta_H, Mu_L, Mu_H, Mu0;
T Eo, Re, Gravity;
T Interface_Width, Mobility, Sigma;
T Beta, Kappa;
int MaxStep, OutputStep;

void readParam() {
  iniReader reader("ferro_bubble.ini");
  Thread_Num = reader.getValue<int>("parallel", "thread_num");
  Ni       = reader.getValue<int>("Mesh", "Ni");
  Nj       = reader.getValue<int>("Mesh", "Nj");
  Cell_Len = reader.getValue<T>("Mesh", "Cell_Len");
  BlockCellLen = reader.getValue<int>("Mesh", "BlockCellLen");
  Drop_Radius = reader.getValue<T>("Bubble", "Radius");
  Drop_CX = reader.getValue<T>("Bubble", "CenterX");
  Drop_CY = reader.getValue<T>("Bubble", "CenterY");
  Rho_L = reader.getValue<T>("Two_Phase", "rho_l");
  Rho_H = reader.getValue<T>("Two_Phase", "rho_h");
  Eta_L = reader.getValue<T>("Two_Phase", "eta_l");
  Eta_H = reader.getValue<T>("Two_Phase", "eta_h");
  Mu_L  = reader.getValue<T>("Two_Phase", "mu_l");
  Mu_H  = reader.getValue<T>("Two_Phase", "mu_h");
  Mu0   = reader.getValue<T>("Two_Phase", "mu0");
  Eo     = reader.getValue<T>("Two_Phase", "Eo");
  Re     = reader.getValue<T>("Two_Phase", "Re");
  Gravity = reader.getValue<T>("Two_Phase", "g");
  Interface_Width = reader.getValue<T>("Phase_Field", "Interface_Width");
  Mobility = reader.getValue<T>("Phase_Field", "Mobility");
  Sigma    = reader.getValue<T>("Phase_Field", "sigma");
  MaxStep    = reader.getValue<int>("Simulation_Settings", "TotalStep");
  OutputStep = reader.getValue<int>("Simulation_Settings", "OutputStep");

  Beta  = T{12} * Sigma / Interface_Width;
  Kappa = T{3} * Interface_Width * Sigma / T{2};

  MPI_RANK(0) {
    std::cout << "------ Ferrofluid Bubble Rise FDLBM ------\n";
    std::cout << "Grid: " << Ni << " x " << Nj << "\n";
    std::cout << "rho: " << Rho_L << "/" << Rho_H
              << "  eta: " << Eta_L << "/" << Eta_H
              << "  mu: " << Mu_L << "/" << Mu_H << "\n";
    std::cout << "Eo=" << Eo << " Re=" << Re
              << " W=" << Interface_Width << "\n";
    std::cout << "beta=" << Beta << " kappa=" << Kappa
              << " M=" << Mobility << "\n";
    std::cout << "------------------------------------------\n";
  }
}

// ================================================================
// Flags
// ================================================================
constexpr std::uint8_t VoidFlag = 1;
constexpr std::uint8_t BulkFlag = 2;
constexpr std::uint8_t WallFlag = 4;

int main(int argc, char* argv[]) {
  mpi().init(&argc, &argv);
  readParam();

  Printer::Print_BigBanner("Initializing Ferrofluid Bubble Rise FDLBM...");

  // ==============================================================
  // Geometry (shared across three solvers)
  // ==============================================================
  AABB<T, 2> domain(Vector<T, 2>(T{0}, T{0}),
                     Vector<T, 2>(T(Ni * Cell_Len), T(Nj * Cell_Len)));
  BlockGeometryHelper2D<T> GeoHelper(Ni, Nj, domain, Cell_Len, BlockCellLen);
  GeoHelper.CreateBlocks(1, mpi().getSize());
  GeoHelper.AdaptiveOptimization(mpi().getSize());
  GeoHelper.LoadBalancing(mpi().getSize());
  BlockGeometry2D<T> Geo(GeoHelper);

  // ==============================================================
  // Three TypePacks and Cell types
  // ==============================================================
  using PF_FIELDS = TypePack<
      FLAG, PHI<T>, POP<T, 9>, VELOCITY<T, 2>, GRAD<T, 2>,
      PHASENORMAL<T, 2>, PHASELAPLACIAN<T>, PHASECHEMPOTENTIAL<T>,
      INTERFACEWIDTH<T>, GBETA<T>>;
  using NS_FIELDS = TypePack<
      FLAG, POP<T, 9>, VELOCITY<T, 2>, RHO<T>, PRESSURE<T>, FORCE<T, 2>,
      VISCOSITY<T>, RHOL<T>, RHOH<T>, ETAL<T>, ETAH<T>>;
  using Mag_FIELDS = TypePack<
      FLAG, POP<T, 5>, PSI<T>, MAGH<T, 2>, MAGHSQ<T>,
      MAGPERMEABILITY<T>, FORCE<T, 2>, MUL<T>, MUH<T>>;

  using PF_CELL  = Cell<T, LatSet9, PF_FIELDS>;
  using NS_CELL  = Cell<T, LatSet9, NS_FIELDS>;
  using Mag_CELL = Cell<T, LatSet5, Mag_FIELDS>;

  // ==============================================================
  // Flag Managers
  // ==============================================================
  BlockFieldManager<FLAG, T, 2> FlagFM_PF(Geo, VoidFlag);
  BlockFieldManager<FLAG, T, 2> FlagFM_NS(Geo, VoidFlag);
  BlockFieldManager<FLAG, T, 2> FlagFM_Mag(Geo, VoidFlag);

  auto setBulk = [](FLAG& f, std::size_t id) { f.SetField(id, BulkFlag); };
  FlagFM_PF.forEach(domain, setBulk);
  FlagFM_NS.forEach(domain, setBulk);
  FlagFM_Mag.forEach(domain, setBulk);

  // Top/Bottom walls (bounce-back for NS)
  AABB<T, 2> topWall(Vector<T, 2>(T{0}, T(Nj * Cell_Len)),
                      Vector<T, 2>(T(Ni * Cell_Len), T((Nj + 1) * Cell_Len)));
  AABB<T, 2> botWall(Vector<T, 2>(T{0}, T(-Cell_Len)),
                      Vector<T, 2>(T(Ni * Cell_Len), T{0}));
  auto setWall = [](FLAG& f, std::size_t id) { f.SetField(id, WallFlag); };
  FlagFM_NS.forEach(topWall, setWall);
  FlagFM_NS.forEach(botWall, setWall);

  FlagFM_PF .template SetupBoundary<LatSet9>(domain, BulkFlag);
  FlagFM_NS .template SetupBoundary<LatSet9>(domain, BulkFlag);
  FlagFM_Mag.template SetupBoundary<LatSet5>(domain, BulkFlag);

  // ==============================================================
  // Lattice Managers
  // ==============================================================
  BaseConverter<T> Conv(LatSet9::cs2);
  Conv.Converter(T{1}, T{1}, T{1}, T(Ni), T{1}, T{1});

  // PF: FLAG,PHI,POP<9>,VEL,GRAD,PHASENORMAL,PHASELAPLACIAN,PHASECHEMPOT,W,GBETA
  ValuePack pf_init(BulkFlag, T{}, T{}, Vector<T,2>{}, Vector<T,2>{},
                    Vector<T,2>{}, T{}, T{}, Interface_Width, Beta);
  // NS: FLAG,POP<9>,VEL,RHO,PRESSURE,FORCE,VISCOSITY,RHOL,RHOH,ETAL,ETAH
  ValuePack ns_init(BulkFlag, T{}, Vector<T,2>{}, Rho_L, T{},
                    Vector<T,2>{}, Eta_L, Rho_L, Rho_H, Eta_L, Eta_H);
  // Mag: FLAG,POP<5>,PSI,MAGH,MAGHSQ,MAGPERMEABILITY,FORCE,MUL,MUH
  ValuePack mag_init(BulkFlag, T{}, T{}, Vector<T,2>{}, T{}, T{},
                     Vector<T,2>{}, Mu_L, Mu_H);

  BlockLatticeManager<T, LatSet9, PF_FIELDS>  PF_LatMan (Geo, pf_init,  Conv);
  BlockLatticeManager<T, LatSet9, NS_FIELDS>  NS_LatMan (Geo, ns_init,  Conv);
  BlockLatticeManager<T, LatSet5, Mag_FIELDS> Mag_LatMan(Geo, mag_init, Conv);

  // ==============================================================
  // Initialize PHI (tanh bubble profile)
  // ==============================================================
  Vector<T, 2> Center(Drop_CX, Drop_CY);
  auto& phiField = PF_LatMan.getField<PHI<T>>();
  for (int b = 0; b < Geo.getBlockNum(); ++b) {
    auto& blk  = Geo.getBlock(b);
    auto& bphi = phiField.getBlockField(b);
    T dx = blk.getVoxelSize();
    for (int j = 0; j < blk.getNy(); ++j) {
      for (int i = 0; i < blk.getNx(); ++i) {
        std::size_t id = j * blk.getProjection()[1] + i;
        T x = blk.getMin()[0] + T(i) * dx;
        T y = blk.getMin()[1] + T(j) * dx;
        T dist = std::sqrt((x - Center[0]) * (x - Center[0])
                         + (y - Center[1]) * (y - Center[1]));
        bphi.get(id) = T{0.5} * (T{1}
          - std::tanh(T{2} * (dist - Drop_Radius) / Interface_Width));
      }
    }
  }

  // ==============================================================
  // Initialize populations (rest state: f_i = w_i)
  // ==============================================================
  auto init_pops = [](auto& LatMan) {
    using LM = std::decay_t<decltype(LatMan)>;
    for (int b = 0; b < LatMan.getBlockLats().size(); ++b) {
      auto& lat = LatMan.getBlockLat(b);
      for (std::size_t id = 0; id < lat.getVoxNum(); ++id) {
        typename LM::CellType cell(id, lat);
        for (int k = 0; k < LM::LatticeSet::q; ++k)
          cell[k] = latset::w<typename LM::LatticeSet>(k);
      }
    }
  };
  init_pops(PF_LatMan);
  init_pops(NS_LatMan);
  init_pops(Mag_LatMan);

  // ==============================================================
  // Boundary Managers (periodic in x-direction)
  // ==============================================================
  using PF_BC  = FixedPeriodicBoundaryManager<
      BlockLatticeManager<T, LatSet9, PF_FIELDS>, BlockFieldManager<FLAG, T, 2>>;
  using NS_BC  = FixedPeriodicBoundaryManager<
      BlockLatticeManager<T, LatSet9, NS_FIELDS>, BlockFieldManager<FLAG, T, 2>>;
  using Mag_BC = FixedPeriodicBoundaryManager<
      BlockLatticeManager<T, LatSet5, Mag_FIELDS>, BlockFieldManager<FLAG, T, 2>>;

  PF_BC  pf_bc ("PF_BC",  PF_LatMan,  FlagFM_PF,  BulkFlag, VoidFlag);
  NS_BC  ns_bc ("NS_BC",  NS_LatMan,  FlagFM_NS,  BulkFlag, VoidFlag);
  Mag_BC mag_bc("Mag_BC", Mag_LatMan, FlagFM_Mag, BulkFlag, VoidFlag);

  pf_bc.Apply();  ns_bc.Apply();  mag_bc.Apply();
  PF_LatMan .NormalCommunicate();
  NS_LatMan .NormalCommunicate();
  Mag_LatMan.NormalCommunicate();


  Timer MainLoopTimer;
  Timer OutputTimer;

  // ==============================================================
  // VTK Writers
  // ==============================================================
  vtmo::ScalarWriter PhiW("phi",  PF_LatMan.getField<PHI<T>>());
  vtmo::VectorWriter VelW("u",    NS_LatMan.getField<VELOCITY<T, 2>>());
  vtmo::ScalarWriter MagW("psi",  Mag_LatMan.getField<PSI<T>>());
  vtmo::VectorWriter HW  ("H",    Mag_LatMan.getField<MAGH<T, 2>>());
  vtmo::ScalarWriter RhoW("rho",  NS_LatMan.getField<RHO<T>>());
  vtmo::vtmWriter<T, 2> MainWriter("ferro_bubble", Geo);
  MainWriter.addWriterSet(PhiW, VelW, MagW, HW, RhoW);
  MainWriter.WriteBinary(MainLoopTimer());

  // ==============================================================
  // Main loop — 10-step coupling schedule (Research.md Section 3.3)
  // ==============================================================

  Printer::Print_BigBanner("Starting Ferrofluid Bubble Rise Simulation...");

  while (MainLoopTimer() < MaxStep) {

    // ---- Step 1: PF pre-coupling ----
    PF_LatMan.template ApplyInnerCellDynamics<IsotropicGradient<PF_CELL>>();
    PF_LatMan.template ApplyInnerCellDynamics<ComputeNormal<PF_CELL>>();
    PF_LatMan.template ApplyInnerCellDynamics<IsotropicLaplacian<PF_CELL>>();
    PF_LatMan.template ApplyInnerCellDynamics<ChemicalPotential<PF_CELL>>();

    // ---- Step 2-3: PF → NS coupling ----
    for (int b = 0; b < Geo.getBlockNum(); ++b) {
      auto& pf = PF_LatMan.getBlockLat(b);
      auto& ns = NS_LatMan.getBlockLat(b);
      for (std::size_t id = 0; id < pf.getVoxNum(); ++id) {
        PF_CELL pfc(id, pf);  NS_CELL nsc(id, ns);
        PFtoNS_Properties<PF_CELL, NS_CELL>::apply(pfc, nsc);
        PFtoNS_Forces    <PF_CELL, NS_CELL>::apply(pfc, nsc);
      }
    }

    // ---- Step 4: PF → Mag coupling ----
    for (int b = 0; b < Geo.getBlockNum(); ++b) {
      auto& pf  = PF_LatMan.getBlockLat(b);
      auto& mag = Mag_LatMan.getBlockLat(b);
      for (std::size_t id = 0; id < pf.getVoxNum(); ++id) {
        PF_CELL pfc(id, pf);  Mag_CELL magc(id, mag);
        PFtoMag_Permeability<PF_CELL, Mag_CELL>::apply(pfc, magc);
      }
    }

    // ---- Step 5: Mag solve ----
    Mag_LatMan.NormalCommunicate();
    Mag_LatMan.template ApplyInnerCellDynamics<MagNeumannBC<Mag_CELL>>();
    Mag_LatMan.template ApplyInnerCellDynamics<MagMRT<Mag_CELL>>();
    Mag_LatMan.template ApplyInnerCellDynamics<moment::MagPotentialMoment<Mag_CELL>>();
    Mag_LatMan.Stream();
    Mag_LatMan.NormalCommunicate();
    Mag_LatMan.template ApplyInnerCellDynamics<MagGradient<Mag_CELL>>();
    Mag_LatMan.template ApplyInnerCellDynamics<MagHSq<Mag_CELL>>();
    Mag_LatMan.template ApplyInnerCellDynamics<MagForce<Mag_CELL>>();

    // ---- Step 6: Mag → NS coupling ----
    for (int b = 0; b < Geo.getBlockNum(); ++b) {
      auto& mag = Mag_LatMan.getBlockLat(b);
      auto& ns  = NS_LatMan.getBlockLat(b);
      for (std::size_t id = 0; id < mag.getVoxNum(); ++id) {
        Mag_CELL magc(id, mag);  NS_CELL nsc(id, ns);
        MagtoNS_Force<Mag_CELL, NS_CELL>::apply(magc, nsc);
      }
    }

    // ---- Step 7: NS solve ----
    NS_LatMan.NormalCommunicate();
    NS_LatMan.template ApplyInnerCellDynamics<NSMRT<NS_CELL>>();
    NS_LatMan.template ApplyInnerCellDynamics<moment::VelocityMoment<NS_CELL>>();
    NS_LatMan.template ApplyInnerCellDynamics<moment::PressureMoment<NS_CELL>>();
    NS_LatMan.Stream();

    // ---- Step 8: NS → PF coupling ----
    for (int b = 0; b < Geo.getBlockNum(); ++b) {
      auto& ns = NS_LatMan.getBlockLat(b);
      auto& pf = PF_LatMan.getBlockLat(b);
      for (std::size_t id = 0; id < ns.getVoxNum(); ++id) {
        NS_CELL nsc(id, ns);  PF_CELL pfc(id, pf);
        NStoPF_Velocity<NS_CELL, PF_CELL>::apply(nsc, pfc);
      }
    }

    // ---- Step 9: PF solve ----
    PF_LatMan.NormalCommunicate();
    PF_LatMan.template ApplyInnerCellDynamics<PhaseFieldMRT<PF_CELL>>();
    PF_LatMan.template ApplyInnerCellDynamics<moment::PhaseFieldMoment<PF_CELL>>();
    PF_LatMan.Stream();

    // ---- Step 10: Boundaries + Output ----
    pf_bc.Apply();  ns_bc.Apply();  mag_bc.Apply();
    PF_LatMan .NormalCommunicate();
    NS_LatMan .NormalCommunicate();
    Mag_LatMan.NormalCommunicate();

    ++MainLoopTimer;  ++OutputTimer;

    if (MainLoopTimer() % OutputStep == 0) {
      if (mpi().getRank() == 0) {
        std::cout << "Step " << MainLoopTimer() << " / " << MaxStep << "\n";
        OutputTimer.Print_InnerLoopPerformance(Geo.getTotalCellNum(), OutputStep);
        Printer::Endl();
      }
      MainWriter.WriteBinary(MainLoopTimer());
    }
  }

  MainWriter.WriteBinary(MainLoopTimer());
  MainLoopTimer.Print_MainLoopPerformance(Geo.getTotalCellNum());
  Printer::Print_BigBanner("Ferrofluid Bubble Rise Simulation Complete!");
  return 0;
}
