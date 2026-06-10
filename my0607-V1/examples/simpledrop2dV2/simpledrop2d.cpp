// simple_droplet2d.cpp
// Simple 2D droplet moving with constant velocity - no fluid dynamics, just advection

#include "freelb.h"
#include "freelb.hh"
#include "ff/ff2d.h"

using T = FLOAT;
using LatSet = D2Q9<T>;  

// 全局变量
Vector<T, LatSet::d> Droplet_Velocity_Global;
T Beta_Global;
T Kappa_Global;
T Mobility_Global;
T Omega_phi_Global;
std::string TestCase_Global;
T U0_Zalesak_Global;
T SlotWidth_Global;
T U0_Vortex_Global;
int N_Vortex_Global;
T T_Vortex_Global;

/*----------------------------------------------
                Simulation Parameters
-----------------------------------------------*/
// geometry
int Ni, Nj;
T Cell_Len;
int BlockCellLen;
int Thread_Num;
// droplets parameters
T Droplet_Radius;  // lattice units
Vector<T, LatSet::d> Droplet_Center;  // lattice units
Vector<T, LatSet::d> Droplet_Velocity; // lattice units/time step
T Droplet_Velocity_Mag;

// phase field parameters
T Mobility;          // mobility M
T Interface_Width;   // interface width ξ
// numerical parameters
T epsilon = T(1e-6); // prevent division by zero
// Simulation settings
int MaxStep;
int OutputStep;

void readParam() {
  iniReader param_reader("simpledrop2d.ini");
  // parallel
  Thread_Num = param_reader.getValue<int>("parallel", "thread_num");

  Ni = param_reader.getValue<int>("Mesh", "Ni");
  Nj = param_reader.getValue<int>("Mesh", "Nj");
  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");
  BlockCellLen = param_reader.getValue<int>("Mesh", "BlockCellLen");
  
  // droplet parameters (in lattice units)
  Droplet_Radius = param_reader.getValue<T>("Droplet", "Radius");
  Droplet_Center[0] = param_reader.getValue<T>("Droplet", "CenterX");
  Droplet_Center[1] = param_reader.getValue<T>("Droplet", "CenterY");
  Droplet_Velocity[0] = param_reader.getValue<T>("Droplet", "VelocityX");
  Droplet_Velocity[1] = param_reader.getValue<T>("Droplet", "VelocityY");
  Droplet_Velocity_Mag = std::sqrt(Droplet_Velocity[0]*Droplet_Velocity[0] +
                                  Droplet_Velocity[1]*Droplet_Velocity[1]);
  // phase field parameters
  Interface_Width = param_reader.getValue<T>("Phase_Field", "Interface_Width");
  Mobility = param_reader.getValue<T>("Phase_Field", "Mobility");
  // Beta, Kappa, Tau_phi, Omega_phi will be computed by Converter in main()
  // Simulation settings
  MaxStep = param_reader.getValue<int>("Simulation_Settings", "TotalStep");
  OutputStep = param_reader.getValue<int>("Simulation_Settings", "OutputStep");
  TestCase_Global = param_reader.getValue<std::string>("Simulation_Settings", "TestCase");
  U0_Zalesak_Global = param_reader.getValue<T>("Zalesak", "U0");
  SlotWidth_Global = param_reader.getValue<T>("Zalesak", "SlotWidth");
  U0_Vortex_Global = param_reader.getValue<T>("Vortex", "U0");
  N_Vortex_Global = param_reader.getValue<int>("Vortex", "n");
  T_Vortex_Global = param_reader.getValue<T>("Vortex", "T");
  MPI_RANK(0) {
    std::cout << "------------Phase Field 2D Droplet Simulation:-------------\n" << std::endl;
    std::cout << "[Phase Field Parameters]:" << std::endl;
    std::cout << "  Interface Width:  " << Interface_Width << " lu" << std::endl;
    std::cout << "  Mobility:         " << Mobility << std::endl;
    std::cout << "  (tau/omega/beta/kappa from Converter in main())" << std::endl;
    std::cout << "[Simulation_Settings]:\n"
              << "  TotalStep:        " << MaxStep << "\n"
              << "  OutputStep:       " << OutputStep << "\n"
              << "  Grid:             " << Ni << " × " << Nj << "\n"
#ifdef _OPENMP
              << "  Running on " << Thread_Num << " threads\n"
#endif
#ifdef MPI_ENABLED
               << "  Running on " << mpi().getSize() << " processes\n"
#endif
              << "----------------------------------------------" << std::endl;
  }
}

int main(int argc, char* argv[]) {
  constexpr std::uint8_t VoidFlag = std::uint8_t(1);
  constexpr std::uint8_t BulkFlag = std::uint8_t(2);
  constexpr std::uint8_t PeriodicFlag = std::uint8_t(4);
  constexpr std::uint8_t BouncebackFlag = std::uint8_t(8);
  bool mpiEnable = false;
  int expected_blocks = 1;

  mpi().init(&argc, &argv);
  if (mpi().getRank() == 0) {
      std::cout << "Total processes: " << mpi().getSize() << std::endl;
      std::cout << "Expected blocks per process: " << expected_blocks << std::endl;
  }

  MPI_DEBUG_WAIT

  Printer::Print_BigBanner(std::string("Initializing Simple 2D Droplet Simulation..."));

  readParam();
  
  Droplet_Velocity_Global = Droplet_Velocity;
  Mobility_Global = Mobility;

  T sigma_lat = T(0.072);
  BaseConverter<T> BaseConv(LatSet::cs2);
  BaseConv.ConvertFromTimeStep(T(1), T(1), T(1), T(Ni), T(1), T(0.1));
  PhaseFieldConverter<T> PhaseConv(BaseConv, LatSet::cs2);
  PhaseConv.fromLattice(Interface_Width, Mobility, sigma_lat);

  Beta_Global = PhaseConv.getLatticeBeta();
  Kappa_Global = PhaseConv.getLatticeKappa();
  Omega_phi_Global = PhaseConv.getOMEGA();

  UnitConvManager<T> ConvManager(&BaseConv, &PhaseConv);
  ConvManager.Check_and_Print();

  MPI_RANK(0) {
    std::cout << "  [Converter]: tau_phi=" << PhaseConv.getLattice_RT()
              << "  omega=" << PhaseConv.getOMEGA()
              << "  beta=" << PhaseConv.getLatticeBeta()
              << "  kappa=" << PhaseConv.getLatticeKappa() << std::endl;
  }

  bool isTaylorGreen = false;//(TestCase_Global == "taylor_green");
  bool isZalesak = false;//(TestCase_Global == "zalesak");
  bool isVortex = true;
  
  // ------------------ define geometry -----------------
  AABB<T, 2> domain(Vector<T, 2>(T(0), T(0)), 
                    Vector<T, 2>(T(Ni * Cell_Len), T(Nj * Cell_Len)));

  AABB<T, 2> left_box(Vector<T, 2>(T(-Cell_Len), T(0)), 
                      Vector<T, 2>(T(0), T(Nj * Cell_Len)));
  AABB<T, 2> right_box(Vector<T, 2>(T((Ni) * Cell_Len), T(0)), 
                      Vector<T, 2>(T((Ni+1) * Cell_Len), T(Nj * Cell_Len)));
  AABB<T, 2> front_box(Vector<T, 2>(T(0), T(-Cell_Len)), 
                      Vector<T, 2>(T(Ni * Cell_Len), T(0)));
  AABB<T, 2> back_box(Vector<T, 2>(T(0), T((Nj) * Cell_Len)), 
                      Vector<T, 2>(T(Ni * Cell_Len), T((Nj+1) * Cell_Len)));
  AABB<T, 2> bottom_left_box(Vector<T, 2>(T(-Cell_Len), T(-Cell_Len)), 
                             Vector<T, 2>(T(0), T(0)));
  AABB<T, 2> bottom_right_box(Vector<T, 2>(T(Ni * Cell_Len), T(-Cell_Len)), 
                              Vector<T, 2>(T((Ni+1) * Cell_Len), T(0)));
  AABB<T, 2> top_left_box(Vector<T, 2>(T(-Cell_Len), T(Nj * Cell_Len)), 
                          Vector<T, 2>(T(0), T((Nj+1) * Cell_Len)));
  AABB<T, 2> top_right_box(Vector<T, 2>(T(Ni * Cell_Len), T(Nj * Cell_Len)), 
                           Vector<T, 2>(T((Ni+1) * Cell_Len), T((Nj+1) * Cell_Len)));

  BlockGeometryHelper2D<T> GeoHelper(Ni, Nj, domain, Cell_Len, BlockCellLen);
  GeoHelper.CreateBlocks(1, mpi().getSize());
  GeoHelper.AdaptiveOptimization(mpi().getSize());
  GeoHelper.LoadBalancing(mpi().getSize());

  BlockGeometry2D<T> Geo(GeoHelper);

  // ------------------ define flag field ------------------
  BlockFieldManager<FLAG, T, LatSet::d> FlagFM(Geo, VoidFlag);

  FlagFM.forEach(
    domain, [&](FLAG& field, std::size_t id) { field.SetField(id, BulkFlag); });

  FlagFM.forEach(
    left_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(
    right_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(
    front_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(
    back_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(
    bottom_left_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(
    bottom_right_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(
    top_left_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(
    top_right_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });

  FlagFM.template SetupBoundary<LatSet>(domain, BulkFlag);
  
  vtmwriter::ScalarWriter FlagWriter("flag", FlagFM);
  vtmwriter::vtmWriter<T, 2> GeoWriter("GeoFlag", Geo);
  GeoWriter.addWriterSet(FlagWriter);
  GeoWriter.WriteBinary();

  // ------------------ define lattice ------------------
  using FIELDS = TypePack<PHI<T>, POP<T, LatSet::q>, VELOCITY<T, LatSet::d>, GRAD<T, LatSet::d>, NORMAL<T, LatSet::d>, INTERFACEWIDTH<T>>;
  using CELL = Cell<T, LatSet, FIELDS>;
  
  ValuePack InitValues(T{}, T{}, T{}, Vector<T, LatSet::d>{T(0), T(0)}, Vector<T, LatSet::d>{T(0), T(0)}, Interface_Width);
  
  BlockLatticeManager<T, LatSet, FIELDS> DropletLattice(Geo, InitValues, PhaseConv);

  // VELOCITY 和分布函数将在初始化循环中按具体测试用例同步设置
  auto& phiField = DropletLattice.getField<PHI<T>>();
  for(int blockid = 0; blockid < Geo.getBlockNum(); ++blockid){
    const auto& block = Geo.getBlock(blockid);
    const auto& proj = block.getProjection();
    auto& blockPhi = phiField.getBlockField(blockid);
    T voxelSize = block.getVoxelSize();
    int overlap = 0;
    for (int j = overlap; j < block.getNy() - overlap; ++j) {
      for (int i = overlap; i < block.getNx() - overlap; ++i) {
        std::size_t id = j * proj[1] + i;
        T x = block.getMin()[0] + static_cast<T>(i) * voxelSize;
        T y = block.getMin()[1] + static_cast<T>(j) * voxelSize;

        T phi;
        if (isVortex) {
          T dist = std::sqrt((x - Droplet_Center[0])*(x - Droplet_Center[0]) + 
                            (y - Droplet_Center[1])*(y - Droplet_Center[1]));
          phi = T(0.5) * (T(1) - std::tanh(T(2) * (dist - Droplet_Radius) / Interface_Width));
        } else if (isZalesak) {
          T d = static_cast<T>(Ni) * Cell_Len;
          T cx = Droplet_Center[0];
          T cy = Droplet_Center[1];
          T R = Droplet_Radius;
          T half_slot = SlotWidth_Global * T(0.5);
          T dist = std::sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy));
          T phi_disk = T(0.5) * (T(1) - std::tanh(T(2) * (dist - R) / Interface_Width));
          T notch_x = T(0.5) * (T(1) - std::tanh(T(2) * (std::abs(x - cx) - half_slot) / Interface_Width));
          T notch_y = T(0.5) * (T(1) + std::tanh(T(2) * (y - (cy + R / T(3))) / Interface_Width))
                     * T(0.5) * (T(1) + std::tanh(T(2) * (cy + R - y) / Interface_Width));
          phi = phi_disk * (T(1) - notch_x * notch_y);
        } else {
          T dist = std::sqrt((x - Droplet_Center[0])*(x - Droplet_Center[0]) + 
                            (y - Droplet_Center[1])*(y - Droplet_Center[1]));
          phi = T(0.5) * (T(1) - std::tanh(2 * (dist - Droplet_Radius) / Interface_Width));
        }
        blockPhi.get(id) = phi;
      }
    }
  }
  
  // 同时初始化速度场和分布函数，确保 f_i 与局部速度 u 一致
  auto& velocityField = DropletLattice.getField<VELOCITY<T, LatSet::d>>();
  for (int block_idx = 0; block_idx < Geo.getBlockNum(); ++block_idx) {
    auto& blockLat = DropletLattice.getBlockLat(block_idx);
    auto& blockPhiField = phiField.getBlockField(block_idx);
    auto& blockVelocity = velocityField.getBlockField(block_idx);
    const auto& block = Geo.getBlock(block_idx);
    int overlap = 0;
    const auto& proj = block.getProjection();
    T voxelSize = block.getVoxelSize();
    T d = static_cast<T>(Ni) * Cell_Len;
    T n_pi = static_cast<T>(N_Vortex_Global) * M_PI;

    for (int j = overlap; j < block.getNy() - overlap; ++j){
      for(int i = overlap; i < block.getNx() - overlap; ++i){
        std::size_t id = j * proj[1] + i;
        CELL cell(id, blockLat);
        T phi = blockPhiField.get(id);

        Vector<T, LatSet::d> vel;
        if (isVortex) {
          T x = block.getMin()[0] + static_cast<T>(i) * voxelSize;
          T y = block.getMin()[1] + static_cast<T>(j) * voxelSize;
          T xi = x / d + T(0.5);
          T eta = y / d + T(0.5);
          T ux = -U0_Vortex_Global * std::sin(n_pi * xi) * std::sin(n_pi * eta);
          T uy = -U0_Vortex_Global * std::cos(n_pi * xi) * std::cos(n_pi * eta);
          vel = Vector<T, 2>(ux, uy);
        } else {
          vel = Droplet_Velocity_Global;
        }
        blockVelocity.get(id) = vel;

        T usqr = vel.getnorm2();
        for (unsigned int k = 0; k < LatSet::q; ++k) {
           T cu = latset::c<LatSet>(k) * vel;
           cell[k] = latset::w<LatSet>(k) * phi *
             (T(1) + LatSet::InvCs2 * cu + cu * cu * T(0.5) * LatSet::InvCs4 -
              LatSet::InvCs2 * usqr * T(0.5));
           }
      }
    }
  }

  FixedPeriodicBoundaryManager<BlockLatticeManager<T, LatSet, FIELDS>, BlockFieldManager<FLAG, T, 2>> 
    Drop_Periodic("Drop_Periodic", DropletLattice, FlagFM, PeriodicFlag, VoidFlag);
  Drop_Periodic.Setup(left_box, NbrDirection::XN, right_box, NbrDirection::XP);
  Drop_Periodic.Setup(front_box, NbrDirection::YN, back_box, NbrDirection::YP);
  Drop_Periodic.Setup(bottom_left_box, NbrDirection::XN, top_right_box, NbrDirection::XP);
  Drop_Periodic.Setup(bottom_right_box, NbrDirection::XP, top_left_box, NbrDirection::XN);
#ifdef MPI_ENABLED
  Drop_Periodic.SetupMPI(GeoHelper);
#endif
  
  using NormalTask = tmp::Key_TypePair<BulkFlag, ff::FF2D<CELL>>;
  using NormalTaskCollection = tmp::TupleWrapper<NormalTask>;
  using DropletNormalTask = tmp::TaskSelector<NormalTaskCollection, std::uint8_t, CELL>;

  using BulkTask = tmp::Key_TypePair<BulkFlag, 
                                     collision::BGKSource<
                                                 equilibrium::SecondOrder<CELL>,
                                                 NORMAL<T, LatSet::d>,
                                                  true
                                                        >
                                    >;
  using CollisionTaskCollection = tmp::TupleWrapper<BulkTask>;
  using DropletCollisionTask = tmp::TaskSelector<CollisionTaskCollection, std::uint8_t, CELL>;

  Timer MainLoopTimer;
  Timer OutputTimer;

  vtmo::ScalarWriter PHIWriter("PHI", DropletLattice.getField<PHI<T>>());
  vtmo::VectorWriter GRADWriter("GRAD", DropletLattice.getField<GRAD<T, 2>>());
  vtmo::VectorWriter VELOCITYWriter("VELOCITY", DropletLattice.getField<VELOCITY<T, LatSet::d>>());
  vtmo::vtmWriter<T, 2> MainWriter("simpledrop2d", Geo);
  MainWriter.addWriterSet(PHIWriter, GRADWriter, VELOCITYWriter);
  MainWriter.WriteBinary(MainLoopTimer());

  Drop_Periodic.Apply();
  DropletLattice.NormalCommunicate();

  Printer::Print_BigBanner(std::string("Start Simple 2D Droplet Simulation..."));

  while (MainLoopTimer() < MaxStep) {
    if(isTaylorGreen){
      // 计算当前物理时间，使用Period参数控制变化周期
      T Period = 16000;
      T U0_Global = 0.03;
      T t = static_cast<T>(MainLoopTimer());
      T cos_pi_t = std::cos(M_PI * t / Period);
      auto& velocityField = DropletLattice.getField<VELOCITY<T, LatSet::d>>();

      // 直接在主循环中实现2D泰勒-格林涡速度场计算
      // u_x = 2*U0*sin²(πx)*sin(2πy)*cos(πt)
      // u_y = -U0*sin²(πy)*sin(2πx)*cos(πt)
      for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
        const auto& block = Geo.getBlock(blockid);
        auto& blockVelocity = velocityField.getBlockField(blockid);
        T voxelSize = block.getVoxelSize();
        const auto& proj = block.getProjection();

        // 遍历所有单元格（包括幽灵层）
        for (int j = 0; j < block.getNy(); ++j) {
          for (int i = 0; i < block.getNx(); ++i) {
            std::size_t id = j * proj[1] + i;

            // 计算物理坐标（归一化到 [0, 1] 范围）
            T x = (block.getMin()[0] + static_cast<T>(i) * voxelSize) / (Ni * Cell_Len);
            T y = (block.getMin()[1] + static_cast<T>(j) * voxelSize) / (Nj * Cell_Len);
            
            // 2D泰勒-格林涡速度场公式
            T sin_pi_x = std::sin(M_PI * x);
            T sin_pi_y = std::sin(M_PI * y);
            T sin_2pi_x = std::sin(2 * M_PI * x);
            T sin_2pi_y = std::sin(2 * M_PI * y);

            T ux = 1 * U0_Global * sin_pi_x * sin_pi_x * sin_2pi_y * cos_pi_t;
            T uy = -U0_Global * sin_pi_y * sin_pi_y * sin_2pi_x * cos_pi_t;

            // 速度限制：确保速度不超过最大允许值（LBM稳定性要求）
            // 声速 cs = 1/sqrt(3) ≈ 0.577，通常速度应小于 0.1-0.2
            const T max_velocity = T(0.2);
            T velocity_mag = std::sqrt(ux * ux + uy * uy);
            if (velocity_mag > max_velocity) {
              T scale = max_velocity / velocity_mag;
              ux *= scale;
              uy *= scale;
            }

            // blockVelocity.get(id) = Vector<T, 2>(ux, uy);
            blockVelocity.get(id) = Vector<T, 2>(ux, uy);
          }
        }
      }
    } else if(isVortex) {
      T d = static_cast<T>(Ni) * Cell_Len;
      T t = static_cast<T>(MainLoopTimer());
      T cos_pi_t_T = std::cos(M_PI * t / T_Vortex_Global);
      T n_pi = static_cast<T>(N_Vortex_Global) * M_PI;
      auto& velocityField = DropletLattice.getField<VELOCITY<T, LatSet::d>>();

      for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
        const auto& block = Geo.getBlock(blockid);
        auto& blockVelocity = velocityField.getBlockField(blockid);
        T voxelSize = block.getVoxelSize();
        const auto& proj = block.getProjection();

        for (int j = 0; j < block.getNy(); ++j) {
          for (int i = 0; i < block.getNx(); ++i) {
            std::size_t id = j * proj[1] + i;
            T x = block.getMin()[0] + static_cast<T>(i) * voxelSize;
            T y = block.getMin()[1] + static_cast<T>(j) * voxelSize;

            T xi = x / d + 0.5;
            T eta = y / d + 0.5;

            T ux = -U0_Vortex_Global * std::sin(n_pi * xi) * std::sin(n_pi * eta) * cos_pi_t_T;
            T uy = -U0_Vortex_Global * std::cos(n_pi * xi) * std::cos(n_pi * eta) * cos_pi_t_T;

            blockVelocity.get(id) = Vector<T, 2>(ux, uy);
          }
        }
      }
    } else if(isZalesak) {
      T d = static_cast<T>(Ni) * Cell_Len;
      T omega = U0_Zalesak_Global * M_PI / d;
      auto& velocityField = DropletLattice.getField<VELOCITY<T, LatSet::d>>();

      for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
        const auto& block = Geo.getBlock(blockid);
        auto& blockVelocity = velocityField.getBlockField(blockid);
        T voxelSize = block.getVoxelSize();
        const auto& proj = block.getProjection();

        for (int j = 0; j < block.getNy(); ++j) {
          for (int i = 0; i < block.getNx(); ++i) {
            std::size_t id = j * proj[1] + i;
            T x = block.getMin()[0] + static_cast<T>(i) * voxelSize;
            T y = block.getMin()[1] + static_cast<T>(j) * voxelSize;

            T ux = -omega * (y - T(0.5) * d);
            T uy =  omega * (x - T(0.5) * d);

            blockVelocity.get(id) = Vector<T, 2>(ux, uy);
          }
        }
      }
    }
    DropletLattice.getField<VELOCITY<T, 2>>().Communicate();

    for (int block_idx = 0; block_idx < Geo.getBlockNum(); ++block_idx) {
        auto& blockLat = DropletLattice.getBlockLat(block_idx);
        auto& blockPhiField = phiField.getBlockField(block_idx);
        const auto& block = Geo.getBlock(block_idx);
        int overlap = block.getOverlap();
        const auto& proj = block.getProjection();
        
        for (int j = overlap; j < block.getNy() - overlap; ++j) {
            for (int i = overlap; i < block.getNx() - overlap; ++i) {
                std::size_t id = j * proj[1] + i;
                CELL cell(id, blockLat);
                T phi_new = T(0);
                for (unsigned int k = 0; k < LatSet::q; ++k) {
                    phi_new += cell[k];
                }
                blockPhiField.get(id) = phi_new;
            }
        }
    }

    DropletLattice.getField<PHI<T>>().Communicate();

    DropletLattice.template ApplyCellDynamics<DropletNormalTask>(FlagFM);
    DropletLattice.getField<NORMAL<T, LatSet::d>>().Communicate();
    DropletLattice.template ApplyCellDynamics<DropletCollisionTask>(FlagFM);

    Drop_Periodic.Apply();

    DropletLattice.Stream();

    DropletLattice.NormalCommunicate();
    
    ++MainLoopTimer;
    ++OutputTimer;

    if (MainLoopTimer() % OutputStep == 0) {
      DropletLattice.getField<GRAD<T, LatSet::d>>().Communicate();
      DropletLattice.getField<VELOCITY<T, 2>>().Communicate();
      DropletLattice.getField<PHI<T>>().Communicate();
      OutputTimer.Print_InnerLoopPerformance(Geo.getTotalCellNum(), OutputStep);
      Printer::Endl();
      MainWriter.WriteBinary(MainLoopTimer());
    }
  }

  Printer::Print_BigBanner(std::string("Simple 2D Droplet Simulation Complete!"));
  MainWriter.WriteBinary(MainLoopTimer());
  MainLoopTimer.Print_MainLoopPerformance(Geo.getTotalCellNum());
  Printer::Print("Total PhysTime", BaseConv.getPhysTime(MainLoopTimer()));
  Printer::Endl();

  return 0;
}
