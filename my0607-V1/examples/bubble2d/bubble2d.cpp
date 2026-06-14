// bubble2d.cpp
// 2D bubble rising simulation - Phase field + NS flow coupling
// No magnetic field, no AMR, MPI parallel

#include "freelb.h"
#include "freelb.hh"
#include "ff/ff2d.h"
#include "lbm/force.h"

using T = FLOAT;
using LatSet = D2Q9<T>;

/*----------------------------------------------
            Simulation Parameters
-----------------------------------------------*/
int Ni;
int Nj;
T Cell_Len;
int BlockCellLen;
int Thread_Num;

// bubble parameters
T Bubble_Radius;
Vector<T, LatSet::d> Bubble_Center;

// phase field parameters
T Interface_Width;
T Mobility;
T Tau_phi;
T Omega_phi;
T Kappa;
T Beta;

// two-phase parameters
T rho_l;
T rho_h;
T eta_l;
T eta_h;
T sigma;
T gravity;
T Tau_ns;

// simulation
int MaxStep;
int OutputStep;

std::string work_dir;

void readParam() {
  iniReader param_reader("bubble2d.ini");
  work_dir = param_reader.getValue<std::string>("workdir", "workdir_");
  Thread_Num = param_reader.getValue<int>("parallel", "thread_num");

  Ni = param_reader.getValue<int>("Mesh", "Ni");
  Nj = param_reader.getValue<int>("Mesh", "Nj");
  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");
  BlockCellLen = param_reader.getValue<int>("Mesh", "BlockCellLen");

  Bubble_Radius = param_reader.getValue<T>("Bubble", "Radius");
  Bubble_Center[0] = param_reader.getValue<T>("Bubble", "CenterX");
  Bubble_Center[1] = param_reader.getValue<T>("Bubble", "CenterY");

  Interface_Width = param_reader.getValue<T>("Phase_Field", "Interface_Width");
  Mobility = param_reader.getValue<T>("Phase_Field", "Mobility");

  rho_l = param_reader.getValue<T>("Two_Phase", "rho_l");
  rho_h = param_reader.getValue<T>("Two_Phase", "rho_h");
  eta_l = param_reader.getValue<T>("Two_Phase", "eta_l");
  eta_h = param_reader.getValue<T>("Two_Phase", "eta_h");
  T Re = param_reader.getValue<T>("Two_Phase", "Re");
  T Eo = param_reader.getValue<T>("Two_Phase", "Eo");

  MaxStep = param_reader.getValue<int>("Simulation_Settings", "TotalStep");
  OutputStep = param_reader.getValue<int>("Simulation_Settings", "OutputStep");

  T D = T(2.0) * Bubble_Radius;

  // nu = eta / rho (using heavier fluid as reference)
  T nu = eta_h / rho_h;

  // U = Re * nu / D → g = U² / D
  T U_g = Re * nu / D;
  T g_abs = U_g * U_g / D;
  gravity = -g_abs;

  // sigma = Δρ * g * D² / Eo
  T DeltaRho = rho_h - rho_l;
  sigma = DeltaRho * g_abs * D * D / Eo;

  // derived from sigma: Eq.(14) β=12σ/W, κ=3Wσ/2
  Kappa = T(3.0) * Interface_Width * sigma * T(0.5);
  Beta = T(12.0) * sigma / Interface_Width;
  Tau_phi = T(0.5) + Mobility / LatSet::cs2;
  Omega_phi = T(1.0) / Tau_phi;

  // NS relaxation time  τ = 0.5 + ν/cs²
  Tau_ns = T(0.5) + nu / LatSet::cs2;

  MPI_RANK(0) {
    std::cout << "----------Bubble Rising Simulation----------\n";
    std::cout << "[Mesh]: " << Ni << "x" << Nj << "  BlockCellLen=" << BlockCellLen << "\n";
    std::cout << "[Bubble]: R=" << Bubble_Radius << " D=" << D
              << "  Center=(" << Bubble_Center[0] << "," << Bubble_Center[1] << ")\n";
    std::cout << "[TwoPhase]: rho_l=" << rho_l << " rho_h=" << rho_h
              << " eta_l=" << eta_l << " eta_h=" << eta_h
              << " Re=" << Re << " Eo=" << Eo << "\n";
    std::cout << "[Derived]: nu=" << nu << " U_g=" << U_g
              << " g=" << gravity << " sigma=" << sigma << "\n";
    std::cout << "[Phase]: W=" << Interface_Width << " M=" << Mobility
              << " beta=" << Beta << " kappa=" << Kappa
              << " tau_phi=" << Tau_phi << "\n";
    std::cout << "[Flow]: tau_ns=" << Tau_ns << "  omega=" << (T(1)/Tau_ns) << "\n";
    std::cout << "[Simulation]: MaxStep=" << MaxStep << "  OutputStep=" << OutputStep << "\n";
#ifdef _OPENMP
    std::cout << "[Parallel]: " << Thread_Num << " threads\n";
#endif
#ifdef MPI_ENABLED
    std::cout << "[Parallel]: " << mpi().getSize() << " MPI processes\n";
#endif
    std::cout << "--------------------------------------------\n";
  }

}

int main(int argc, char* argv[]) {
  constexpr std::uint8_t VoidFlag = std::uint8_t(1);
  constexpr std::uint8_t BulkFlag = std::uint8_t(2);
  constexpr std::uint8_t BouncebackFlag = std::uint8_t(4);

  mpi().init(&argc, &argv);
  MPI_DEBUG_WAIT

  Printer::Print_BigBanner(std::string("Initializing Bubble Rising..."));

  readParam();

  // ------------------ define converters ------------------
  BaseConverter<T> BaseConv(LatSet::cs2);
  BaseConv.SimplifiedConverterFromRT(Ni, T(0.01), Tau_ns);

  BaseConverter<T> PFBaseConv(LatSet::cs2);
  PFBaseConv.SimplifiedConverterFromRT(Ni, T(0.01), Tau_phi);

  UnitConvManager<T> ConvManager(&BaseConv);
  ConvManager.Check_and_Print();

  // ------------------ define geometry ------------------
  AABB<T, 2> domain(Vector<T, 2>(T(0), T(0)),
                    Vector<T, 2>(T(Ni * Cell_Len), T(Nj * Cell_Len)));

  BlockGeometryHelper2D<T> GeoHelper(Ni, Nj, domain, Cell_Len, BlockCellLen);
  GeoHelper.CreateBlocks(1,mpi().getSize());
  GeoHelper.AdaptiveOptimization(mpi().getSize());
  GeoHelper.LoadBalancing(mpi().getSize());

  BlockGeometry2D<T> Geo(GeoHelper);

  // ------------------ define flag field ------------------
  BlockFieldManager<FLAG, T, LatSet::d> FlagFM(Geo, VoidFlag);
  FlagFM.forEach(domain,
                 [&](FLAG& field, std::size_t id) { field.SetField(id, BulkFlag); });
  FlagFM.template SetupBoundary<LatSet>(domain, BouncebackFlag);

  // Write flag geometry for verification
  vtmo::ScalarWriter FlagWriter("flag", FlagFM);
  vtmo::vtmWriter<T, 2> GeoWriter("GeoFlag_Bubble", Geo, 1);
  GeoWriter.addWriterSet(FlagWriter);
  // GeoWriter.WriteBinary();

  // ------------------ define NS lattice ------------------
  using NSFIELDS = TypePack<RHO<T>, VELOCITY<T, 2>, POP<T, LatSet::q>,
                            FORCE<T, LatSet::d>, PHYSICAL_ETA<T>, OMEGA<T>,
                            PRESSURE<T>, CONSTFORCE<T, LatSet::d>>;
  ValuePack NSInitValues(BaseConv.getLatRhoInit(), Vector<T, 2>{T{0}, T{0}},
                         T{}, Vector<T, 2>{T{0}, T{0}}, T{}, BaseConv.getOMEGA(),
                         T{}, Vector<T, 2>{T{0}, gravity});
  using NSCELL = Cell<T, LatSet, NSFIELDS>;
  BlockLatticeManager<T, LatSet, NSFIELDS> NSLattice(Geo, NSInitValues, BaseConv);

  // ------------------ define PF lattice ------------------
  using PFFIELDS = TypePack<PHI<T>, POP<T, LatSet::q>, GRAD<T, LatSet::d>,
                            NORMAL<T, LatSet::d>, INTERFACEWIDTH<T>,
                            RHO_L<T>, RHO_H<T>, ETA_L<T>, ETA_H<T>,
                            BETA<T>, KAPPA<T>>;
  using PFFIELDREFS = TypePack<VELOCITY<T, LatSet::d>>;
  using PFFIELDPACK = TypePack<PFFIELDS, PFFIELDREFS>;
  ValuePack PFInitValues(T{}, T{}, Vector<T, 2>{T{0}, T{0}},
                         Vector<T, 2>{T{0}, T{0}}, Interface_Width,
                         rho_l, rho_h, eta_l, eta_h,
                         Beta, Kappa);
  using PFCELL = Cell<T, LatSet, ExtractFieldPack<PFFIELDPACK>::mergedpack>;
  BlockLatticeManager<T, LatSet, PFFIELDPACK> PFLattice(
    Geo, PFInitValues, PFBaseConv, &NSLattice.getField<VELOCITY<T, LatSet::d>>());


  // ---- Initialize PHI and POP fields ----
  // φ = 0 for light fluid (bubble interior), φ = 1 for heavy fluid (outside)
  // φ = 0.5 + 0.5*tanh(2*(r-R)/W) per Eq.(11)
  T R_phys = Bubble_Radius * Cell_Len;
  T xc_phys = Bubble_Center[0] * Cell_Len;
  T yc_phys = Bubble_Center[1] * Cell_Len;
  T W_phys = Interface_Width * Cell_Len;

  auto& phiField = PFLattice.getField<PHI<T>>();
  for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
    const auto& block = Geo.getBlock(blockid);
    const auto& proj = block.getProjection();
    auto& blockPhi = phiField.getBlockField(blockid);
    auto& blockLat = PFLattice.getBlockLat(blockid);
    T voxelSize = block.getVoxelSize();
    int overlap = 0;
    for (int j = overlap; j < block.getNy() - overlap; ++j) {
      for (int i = overlap; i < block.getNx() - overlap; ++i) {
        std::size_t id = j * proj[1] + i;
        T x = block.getMin()[0] + static_cast<T>(i) * voxelSize;
        T y = block.getMin()[1] + static_cast<T>(j) * voxelSize;
        T dx = x - xc_phys;
        T dy = y - yc_phys;
        T dist = std::sqrt(dx * dx + dy * dy);
        T phi = T(0.5) + T(0.5) * std::tanh(T(2.0) * (dist - R_phys) / W_phys);
        blockPhi.get(id) = phi;
      }
    }
  }

  for (int block_idx = 0; block_idx < Geo.getBlockNum(); ++block_idx) {
    auto& blockLat = PFLattice.getBlockLat(block_idx);
    auto& blockPhiField = phiField.getBlockField(block_idx);
    const auto& block = Geo.getBlock(block_idx);
    int overlap = 0;
    const auto& proj = block.getProjection();
    
    for (int j = overlap; j < block.getNy() - overlap; ++j){ 
      for(int i = overlap; i < block.getNx() - overlap; ++i){ 
        std::size_t id = j * proj[1] + i;
        PFCELL cell(id, blockLat);
        T phi = blockPhiField.get(id);
        
        for (unsigned int i = 0; i < LatSet::q; ++i) {
           T cu = 0;
           cell[i] = latset::w<LatSet>(i) * phi * (T(1) + LatSet::InvCs2 * cu);
           }
      }
    }
  }

  // Initialize PF interface width
  ff::BroadcastAllParams<T>(PFLattice,Interface_Width,rho_l, rho_h, eta_l, eta_h,gravity, Beta, Kappa);
  // Initialize NS populations to velocity-based equilibrium
  // geq[k] = w_k * (p/(rho*cs²) + e_k·u/cs² + ...)  with p=rho*cs², u=0  →  geq[k] = w_k
  Vector<T, 2> u_zero{T{0}, T{0}};
  T ns_rho_init = BaseConv.getLatRhoInit();
  T ns_p_init = ns_rho_init * LatSet::cs2;
  for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
    const auto& block = Geo.getBlock(blockid);
    const auto& proj = block.getProjection();
    auto& blockLat = NSLattice.getBlockLat(blockid);
    int overlap = 0;
    for (int j = overlap; j < block.getNy() - overlap; ++j) {
      for (int i = overlap; i < block.getNx() - overlap; ++i) {
        std::size_t id = j * proj[1] + i;
        T u2 = T{0};
        NSCELL cell(id, blockLat);
        for (unsigned int k = 0; k < LatSet::q; ++k) {
          T uc = u_zero * latset::c<LatSet>(k);
          cell[k] = latset::w<LatSet>(k) *
                    (ns_p_init * LatSet::InvCs2 / ns_rho_init
                     + LatSet::InvCs2 * uc
                     + uc * uc * T{0.5} * LatSet::InvCs4
                     - LatSet::InvCs2 * u2 * T{0.5});
        }
        cell.template get<PRESSURE<T>>() = ns_p_init;
      }
    }
  }



  // ------------------ define BCs ------------------
  // NS: bounceback on all walls
  BBLikeFixedBlockBdManager<bounceback::normal<NSCELL>,
                            BlockLatticeManager<T, LatSet, NSFIELDS>,
                            BlockFieldManager<FLAG, T, LatSet::d>>
    NS_BB("NS_BB", NSLattice, FlagFM, BouncebackFlag, VoidFlag);

  // PF: bounceback on all walls
  using PFBLKLAT = BlockLatticeManager<T, LatSet, PFFIELDPACK>;
  BBLikeFixedBlockBdManager<bounceback::normal<PFCELL>, PFBLKLAT,
                            BlockFieldManager<FLAG, T, LatSet::d>>
    PF_BB("PF_BB", PFLattice, FlagFM, BouncebackFlag, VoidFlag);

  // ------------------ define tasks ------------------
  using UpdatePhiTask = tmp::Key_TypePair<BulkFlag, moment::PhaseFieldPhi<PFCELL>>;
  using UpdatePhiNSTaskSelector = TaskSelector<std::uint8_t, PFCELL, UpdatePhiTask>;
  using UpdatePropTask = tmp::Key_TypePair<BulkFlag, moment::PropertyInterpolation<PFCELL, NSCELL>>;
  using InterpolationTaskSelector = CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, UpdatePropTask>;
  BlockLatManagerCoupling InterpolationCoupling(PFLattice, NSLattice);

  using UpdateGradTask = tmp::Key_TypePair<BulkFlag, ff::PhaseFieldGrad<PFCELL>>;
  using UpdateGradTaskSelector = TaskSelector<std::uint8_t, PFCELL, UpdateGradTask>;
  using UpdateNormalTask = tmp::Key_TypePair<BulkFlag, ff::PhaseFieldNormal<PFCELL>>;
  using UpdateNormalTaskSelector = TaskSelector<std::uint8_t, PFCELL, UpdateNormalTask>;

  using PhaseCollisionTask = tmp::Key_TypePair<BulkFlag, PhaseFieldBGK<PFCELL, equilibrium::SecondOrder<PFCELL>>>;
  using PhaseCollisionTaskSelector = TaskSelector<std::uint8_t, PFCELL, PhaseCollisionTask>;
  
  using NSMomentTask = tmp::Key_TypePair<BulkFlag, moment::VelocityMomenta<NSCELL>>;
  using NSMomentTaskSelector = TaskSelector<std::uint8_t, NSCELL, NSMomentTask>;

  using NSCollisionTask = tmp::Key_TypePair<BulkFlag, VelocityBGK<NSCELL,equilibrium::NSFieldEquilibrium<NSCELL>>>;
  using NSCollisionTaskSelector = TaskSelector<std::uint8_t, NSCELL, NSCollisionTask>;
  
  using NSGravityTask = tmp::Key_TypePair<BulkFlag, force::Gravity<PFCELL, NSCELL>>;
  using NSSurfaceTensionTask = tmp::Key_TypePair<BulkFlag, force::SurfaceTension<PFCELL, NSCELL>>;
  using NSForceTaskSelector = CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, NSGravityTask, NSSurfaceTensionTask>;
  BlockLatManagerCoupling NSForceCoupling(PFLattice, NSLattice);
  
  using NSVelPressTask = tmp::Key_TypePair<BulkFlag, moment::VelocityMomenta<NSCELL>>;
  using NSVelPressTaskSelector = TaskSelector<std::uint8_t, NSCELL, NSVelPressTask>;

  // ------------------ writers ------------------

  // ------------------ writers ------------------
  vtmo::ScalarWriter PHIWriter("PHI", PFLattice.getField<PHI<T>>());
  vtmo::VectorWriter GRADWriter("GRAD", PFLattice.getField<GRAD<T, 2>>());
  vtmo::VectorWriter NormalWriter("NORMAL", PFLattice.getField<NORMAL<T, 2>>());
  vtmo::VectorWriter VecWriter("Velocity", NSLattice.getField<VELOCITY<T, 2>>());
  vtmo::ScalarWriter RhoWriter("Rho", NSLattice.getField<RHO<T>>());
  vtmo::VectorWriter ForceWriter("Force", NSLattice.getField<FORCE<T, 2>>());

  vtmo::vtmWriter<T, 2> MainWriter("bubble2d", Geo);
  MainWriter.addWriterSet(PHIWriter, GRADWriter, NormalWriter,
                          VecWriter, RhoWriter, ForceWriter);

  // ------------------ timer ------------------
  Timer MainLoopTimer;
  Timer OutputTimer;

  PFLattice.NormalCommunicate();
  NSLattice.NormalCommunicate();
  // MainWriter.WriteBinary(MainLoopTimer());

  Printer::Print_BigBanner(std::string("Start Calculation..."));

  // TODO: Main simulation loop will be re-implemented after LBM/FF refactoring

  while (MainLoopTimer() < MaxStep) {
    
    PFLattice.template ApplyCellDynamics<UpdatePhiNSTaskSelector>(FlagFM);
    PFLattice.template getField<PHI<T>>().Communicate();
    
    InterpolationCoupling.template ApplyCellDynamics<InterpolationTaskSelector>(FlagFM);
    NSLattice.template getField<RHO<T>>().Communicate();
    NSLattice.template getField<OMEGA<T>>().Communicate();
    
    PFLattice.template ApplyCellDynamics<UpdateGradTaskSelector>(FlagFM);
    PFLattice.template getField<GRAD<T, 2>>().Communicate();
    PFLattice.template ApplyCellDynamics<UpdateNormalTaskSelector>(FlagFM);
    PFLattice.template getField<NORMAL<T, 2>>().Communicate();
    
    PFLattice.template ApplyCellDynamics<PhaseCollisionTaskSelector>(FlagFM);
    PF_BB.Apply();
    PFLattice.Stream();
    PFLattice.NormalCommunicate();

    // Clear total force and accumulate from current phi/velocity fields
    NSLattice.template getField<FORCE<T, 2>>().InitValue(Vector<T, 2>{T{0}, T{0}});
    NSForceCoupling.template ApplyCellDynamics<NSForceTaskSelector>(FlagFM);
    NSLattice.template getField<FORCE<T, 2>>().Communicate();

    // NS collision uses FORCE from current time step
    NSLattice.template ApplyCellDynamics<NSCollisionTaskSelector>(FlagFM);
    NS_BB.Apply();
    NSLattice.Stream();
    NSLattice.NormalCommunicate();

    NSLattice.template ApplyCellDynamics<NSVelPressTaskSelector>(FlagFM);
    NSLattice.template getField<VELOCITY<T, 2>>().Communicate();
    NSLattice.template getField<PRESSURE<T>>().Communicate();

    ++MainLoopTimer;
    ++OutputTimer;
    if (MainLoopTimer() % OutputStep == 0) {
      OutputTimer.Print_InnerLoopPerformance(Geo.getTotalCellNum(), OutputStep);
      Printer::Endl();
      MainWriter.WriteBinary(MainLoopTimer());
    }
  }
  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  return 0;
}
