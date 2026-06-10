// simple_droplet3d.cpp
// Simple 3D droplet moving with constant velocity - no fluid dynamics, just advection

#include "freelb.h"
#include "freelb.hh"
#include "ff/ff2d.h"

using T = FLOAT;
using LatSet = D3Q19<T>;  

// 全局变量
Vector<T, LatSet::d> Droplet_Velocity_Global;
T Beta_Global;
T Kappa_Global;
T Mobility_Global;
T Omega_phi_Global;

/*----------------------------------------------
                Simulation Parameters
-----------------------------------------------*/
// geometry
int Nx, Ny, Nz;
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
T Interface_Width;   // interface width W
// numerical parameters
T epsilon = T(1e-6); // prevent division by zero
// Simulation settings
int MaxStep;
int OutputStep;

void readParam() {
  iniReader param_reader("simpledrop3d.ini");
  // parallel
  Thread_Num = param_reader.getValue<int>("parallel", "thread_num");

  Nx = param_reader.getValue<int>("Mesh", "Nx");
  Ny = param_reader.getValue<int>("Mesh", "Ny");
  Nz = param_reader.getValue<int>("Mesh", "Nz");
  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");
  BlockCellLen = param_reader.getValue<int>("Mesh", "BlockCellLen");
  
  // droplet parameters (in lattice units)
  Droplet_Radius = param_reader.getValue<T>("Droplet", "Radius");
  Droplet_Center[0] = param_reader.getValue<T>("Droplet", "CenterX");
  Droplet_Center[1] = param_reader.getValue<T>("Droplet", "CenterY");
  Droplet_Center[2] = param_reader.getValue<T>("Droplet", "CenterZ");
  Droplet_Velocity[0] = param_reader.getValue<T>("Droplet", "VelocityX");
  Droplet_Velocity[1] = param_reader.getValue<T>("Droplet", "VelocityY");
  Droplet_Velocity[2] = param_reader.getValue<T>("Droplet", "VelocityZ");
  Droplet_Velocity_Mag = std::sqrt(Droplet_Velocity[0]*Droplet_Velocity[0] +
                                  Droplet_Velocity[1]*Droplet_Velocity[1] +
                                  Droplet_Velocity[2]*Droplet_Velocity[2]);
  // phase field parameters
  Interface_Width = param_reader.getValue<T>("Phase_Field", "Interface_Width");
  Mobility = param_reader.getValue<T>("Phase_Field", "Mobility");
  // Beta, Kappa, Tau_phi, Omega_phi will be computed by Converter in main()
  // Simulation settings
  MaxStep = param_reader.getValue<int>("Simulation_Settings", "TotalStep");
  OutputStep = param_reader.getValue<int>("Simulation_Settings", "OutputStep");
  MPI_RANK(0) {
    std::cout << "------------Phase Field 3D Droplet Simulation:-------------\n" << std::endl;
    std::cout << "[Phase Field Parameters]:" << std::endl;
    std::cout << "  Interface Width:  " << Interface_Width << " lu" << std::endl;
    std::cout << "  Mobility:         " << Mobility << std::endl;
    std::cout << "  (tau/omega/beta/kappa from Converter in main())" << std::endl;
    std::cout << "[Simulation_Settings]:\n"
              << "  TotalStep:        " << MaxStep << "\n"
              << "  OutputStep:       " << OutputStep << "\n"
              << "  Grid:             " << Nx << " × " << Ny << " × " << Nz << "\n"
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
  int expected_blocks = 1; // 期望的每个进程的块数

  mpi().init(&argc, &argv);
  // 检查进程数和几何划分是否匹配
  if (mpi().getRank() == 0) {
      std::cout << "Total processes: " << mpi().getSize() << std::endl;
      std::cout << "Expected blocks per process: " << expected_blocks << std::endl;
  }

  MPI_DEBUG_WAIT

  Printer::Print_BigBanner(std::string("Initializing Simple 3D Droplet Simulation..."));

  readParam();
  
  // 设置全局变量
  Droplet_Velocity_Global = Droplet_Velocity;
  Mobility_Global = Mobility;

  T sigma_lat = T(0.072);
  BaseConverter<T> BaseConv(LatSet::cs2);
  BaseConv.ConvertFromTimeStep(T(1), T(1), T(1), T(Nx), T(1), T(0.1));
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

  // ------------------ define geometry -----------------
  // 使用周期性边界条件
  AABB<T, 3> domain(Vector<T, 3>(T(0), T(0), T(0)), 
                    Vector<T, 3>(T(Nx * Cell_Len), T(Ny * Cell_Len), T(Nz * Cell_Len)));

  // 在创建周期性边界管理器后，需要正确设置边界盒子
  AABB<T, 3> left_box(Vector<T, 3>(T(-Cell_Len), T(0), T(0)), 
                      Vector<T, 3>(T(0), T(Ny * Cell_Len), T(Nz * Cell_Len)));
  AABB<T, 3> right_box(Vector<T, 3>(T((Nx) * Cell_Len), T(0), T(0)), 
                      Vector<T, 3>(T((Nx+1) * Cell_Len), T(Ny * Cell_Len), T(Nz * Cell_Len)));
  AABB<T, 3> front_box(Vector<T, 3>(T(0), T(-Cell_Len), T(0)), 
                      Vector<T, 3>(T(Nx * Cell_Len), T(0), T(Nz * Cell_Len)));
  AABB<T, 3> back_box(Vector<T, 3>(T(0), T((Ny) * Cell_Len), T(0)), 
                      Vector<T, 3>(T(Nx * Cell_Len), T((Ny+1) * Cell_Len), T(Nz * Cell_Len)));
  AABB<T, 3> bottom_box(Vector<T, 3>(T(0), T(0), T(-Cell_Len)), 
                        Vector<T, 3>(T(Nx * Cell_Len), T(Ny * Cell_Len), T(0)));
  AABB<T, 3> top_box(Vector<T, 3>(T(0), T(0), T((Nz) * Cell_Len)), 
                    Vector<T, 3>(T(Nx * Cell_Len), T(Ny * Cell_Len), T((Nz+1) * Cell_Len)));

  BlockGeometryHelper3D<T> GeoHelper(Nx, Ny, Nz, domain, Cell_Len, BlockCellLen);
  GeoHelper.CreateBlocks(1, mpi().getSize(), 1);
  GeoHelper.AdaptiveOptimization(mpi().getSize());
  GeoHelper.LoadBalancing(mpi().getSize());

  BlockGeometry3D<T> Geo(GeoHelper);

  // ------------------ define flag field ------------------
  BlockFieldManager<FLAG, T, LatSet::d> FlagFM(Geo, VoidFlag);

  FlagFM.forEach(
    domain, [&](FLAG& field, std::size_t id) { field.SetField(id, BulkFlag); });

  // 设置周期性边界（所有边界都是周期性的）
  FlagFM.template SetupBoundary<LatSet>(domain, PeriodicFlag);

  FlagFM.forEach(
    left_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(
    right_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(
    front_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(
    back_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(
    bottom_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(
    top_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  
  vtmwriter::ScalarWriter FlagWriter("flag", FlagFM);
  vtmwriter::vtmWriter<T, 3> GeoWriter("GeoFlag", Geo);
  GeoWriter.addWriterSet(FlagWriter);
  GeoWriter.WriteBinary();

  // ------------------ define lattice ------------------
  using FIELDS = TypePack<PHI<T>, POP<T, LatSet::q>, VELOCITY<T, LatSet::d>, GRAD<T, LatSet::d>, NORMAL<T, LatSet::d>, INTERFACEWIDTH<T>>;
  using CELL = Cell<T, LatSet, FIELDS>;
  
  ValuePack InitValues(T{}, T{}, T{}, Vector<T, LatSet::d>{T(0), T(0), T(0)}, Vector<T, LatSet::d>{T(0), T(0), T(0)}, T{});
  
  // 创建格子
  BlockLatticeManager<T, LatSet, FIELDS> DropletLattice(Geo, InitValues, PhaseConv);

  DropletLattice.getField<VELOCITY<T, LatSet::d>>().forEach(
    domain, FlagFM, BulkFlag,
    [&](auto& field, std::size_t id) { field.SetField(id, Droplet_Velocity_Global); });

  // 初始化相场（球形液滴）- 在每个进程中分别处理各自负责的块
  auto& phiField = DropletLattice.getField<PHI<T>>();
  for(int blockid = 0; blockid < Geo.getBlockNum(); ++blockid){
    const auto& block = Geo.getBlock(blockid);
    const auto& proj = block.getProjection();
    auto& blockPhi = phiField.getBlockField(blockid);
    T voxelSize = block.getVoxelSize();
    int overlap = block.getOverlap();
    for (int k = overlap; k < block.getNz() - overlap; ++k) {
      for (int j = overlap; j < block.getNy() - overlap; ++j) {
        for (int i = overlap; i < block.getNx() - overlap; ++i) {
          std::size_t id = k * proj[2] + j * proj[1] + i;
          T x = block.getMin()[0] + static_cast<T>(i) * voxelSize;
          T y = block.getMin()[1] + static_cast<T>(j) * voxelSize;
          T z = block.getMin()[2] + static_cast<T>(k) * voxelSize;

          T dist = std::sqrt((x - Droplet_Center[0])*(x - Droplet_Center[0]) + 
                            (y - Droplet_Center[1])*(y - Droplet_Center[1]) +
                            (z - Droplet_Center[2])*(z - Droplet_Center[2]));
          
          // tanh界面
          T phi = T(0.5) * (T(1) - std::tanh(2 * (dist - Droplet_Radius) / Interface_Width));
          blockPhi.get(id) = phi;
        }
      }
    }
  }
  

  // 初始化分布函数为平衡态
  for (int block_idx = 0; block_idx < Geo.getBlockNum(); ++block_idx) {
    auto& blockLat = DropletLattice.getBlockLat(block_idx);
    auto& blockPhiField = phiField.getBlockField(block_idx);
    const auto& block = Geo.getBlock(block_idx);
    int overlap = 0; //block.getOverlap();
    const auto& proj = block.getProjection();
    
    for (int k = overlap; k < block.getNz() - overlap; ++k) { 
      for(int j = overlap; j < block.getNy() - overlap; ++j){
        for(int i = overlap; i < block.getNx() - overlap; ++i){ 
          std::size_t id = k * proj[2] + j * proj[1] + i;
          CELL cell(id, blockLat);
          T phi = blockPhiField.get(id);
          
          // 计算平衡分布函数
          // std::array<T, LatSet::q> geq;
          for (unsigned int i = 0; i < LatSet::q; ++i) {
             T cu = latset::c<LatSet>(i) * Droplet_Velocity_Global;
             cell[i] = latset::w<LatSet>(i) * phi * (T(1) + LatSet::InvCs2 * cu);
             }
        }
      }
    }
  }


  FixedPeriodicBoundaryManager<BlockLatticeManager<T, LatSet, FIELDS>, BlockFieldManager<FLAG, T, 3>> 
    Drop_Periodic("Drop_Periodic", DropletLattice, FlagFM, PeriodicFlag, VoidFlag);
  Drop_Periodic.Setup(left_box, NbrDirection::XN, right_box, NbrDirection::XP);
  Drop_Periodic.Setup(front_box, NbrDirection::YN, back_box, NbrDirection::YP);
  Drop_Periodic.Setup(bottom_box, NbrDirection::ZN, top_box, NbrDirection::ZP);
  
  // ------------------ define task/ dynamics ------------------

  using NormalTask = tmp::Key_TypePair<BulkFlag, ff::FF3D<CELL>>;
  using NormalTaskCollection = tmp::TupleWrapper<NormalTask>;
  using DropletNormalTask = tmp::TaskSelector<NormalTaskCollection, std::uint8_t, CELL>;

  using BulkTask = tmp::Key_TypePair<BulkFlag, 
                                     collision::BGKSource<
                                                 equilibrium::FirstOrder<CELL>,
                                                 NORMAL<T, LatSet::d>,
                                                  true
                                                        >
                                    >;
  using CollisionTaskCollection = tmp::TupleWrapper<BulkTask>;
  using DropletCollisionTask = tmp::TaskSelector<CollisionTaskCollection, std::uint8_t, CELL>;

  
  // count and timer
  Timer MainLoopTimer;
  Timer OutputTimer;

  // writer
  vtmo::ScalarWriter PHIWriter("PHI", DropletLattice.getField<PHI<T>>());
  vtmo::VectorWriter GRADWriter("GRAD", DropletLattice.getField<GRAD<T, 3>>());
  vtmo::vtmWriter<T, 3> MainWriter("simpledrop3d", Geo);
  MainWriter.addWriterSet(PHIWriter, GRADWriter);
  MainWriter.WriteBinary();

  Printer::Print_BigBanner(std::string("Start Simple 3D Droplet Simulation..."));

  while (MainLoopTimer() < MaxStep) {
    for (int block_idx = 0; block_idx < Geo.getBlockNum(); ++block_idx) {
      auto& blockLat = DropletLattice.getBlockLat(block_idx);
      auto& blockPhiField = phiField.getBlockField(block_idx);
      for (std::size_t id = 0; id < blockLat.getVoxNum(); ++id) {
        CELL cell(id, blockLat);
        T phi_new = T(0);
        for (unsigned int i = 0; i < LatSet::q; ++i) {
          phi_new += cell[i];
        }
        blockPhiField.get(id) = phi_new;
      }
    }

    DropletLattice.getField<PHI<T>>().Communicate();
    DropletLattice.template ApplyCellDynamics<DropletNormalTask>(FlagFM);
    DropletLattice.getField<NORMAL<T, LatSet::d>>().Communicate();
    DropletLattice.template ApplyCellDynamics<DropletCollisionTask>(FlagFM);
    
    DropletLattice.Stream();
    
    // 应用周期性边界条件
    Drop_Periodic.Apply();
    
    // 通信更新重叠区域
    DropletLattice.NormalCommunicate();
  
    ++MainLoopTimer;
    ++OutputTimer;

    if (MainLoopTimer() % OutputStep == 0) {
      DropletLattice.getField<GRAD<T, LatSet::d>>().Communicate();
      DropletLattice.getField<VELOCITY<T, 3>>().Communicate();
      DropletLattice.getField<PHI<T>>().Communicate();
      OutputTimer.Print_InnerLoopPerformance(Geo.getTotalCellNum(), OutputStep);
      Printer::Endl();
      MainWriter.WriteBinary(MainLoopTimer());

    }
  }

  Printer::Print_BigBanner(std::string("Simple 3D Droplet Simulation Complete!"));
  MainWriter.WriteBinary(MainLoopTimer());
  MainLoopTimer.Print_MainLoopPerformance(Geo.getTotalCellNum());
  Printer::Print("Total PhysTime", BaseConv.getPhysTime(MainLoopTimer()));
  Printer::Endl();



  return 0;
}