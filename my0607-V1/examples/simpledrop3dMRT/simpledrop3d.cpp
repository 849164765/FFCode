#include "freelb.h"
#include "freelb.hh"
#include "ff/ff2d.h"
#include "lbm/collisionMRT.h"

using T = FLOAT;
using LatSet = D3Q19<T>;

Vector<T, LatSet::d> Droplet_Velocity_Global;
T Beta_Global;
T Kappa_Global;
T Mobility_Global;
T Omega_phi_Global;
std::string TestCase_Global;
T U0_Global;
T SlotWidth_Global;
T SlotHeight_Global;
T Tf_Global;
T T_Global;

bool isSlottedSphere;
bool isVortex;
bool isDeformation;
bool isShear;
bool isPeriodic;

int Nx, Ny, Nz;
T Cell_Len;
int BlockCellLen;
int Thread_Num;
T Droplet_Radius;
Vector<T, LatSet::d> Droplet_Center;
Vector<T, LatSet::d> Droplet_Velocity;
T Droplet_Velocity_Mag;

T Mobility;
T Sigma;
T Interface_Width;
T epsilon = T(1e-6);
int MaxStep;
int OutputStep;

void readParam() {
  iniReader param_reader("simpledrop3d.ini");
  Thread_Num = param_reader.getValue<int>("parallel", "thread_num");

  Nx = param_reader.getValue<int>("Mesh", "Nx");
  Ny = param_reader.getValue<int>("Mesh", "Ny");
  Nz = param_reader.getValue<int>("Mesh", "Nz");
  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");
  BlockCellLen = param_reader.getValue<int>("Mesh", "BlockCellLen");

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

  Interface_Width = param_reader.getValue<T>("Phase_Field", "Interface_Width");
  Mobility = param_reader.getValue<T>("Phase_Field", "Mobility");
  Sigma = param_reader.getValue<T>("Phase_Field", "Sigma");
  // Beta, Kappa, Tau_phi, Omega_phi will be computed by Converter in main()

  OutputStep = param_reader.getValue<int>("Simulation_Settings", "OutputStep");
  TestCase_Global = param_reader.getValue<std::string>("Simulation_Settings", "TestCase");

  isSlottedSphere = false;
  isVortex = true;
  isDeformation = false;
  isShear = false;
  isPeriodic = true;

  if (isSlottedSphere) {
    U0_Global = param_reader.getValue<T>("SlottedSphere", "U0");
    SlotWidth_Global = param_reader.getValue<T>("SlottedSphere", "SlotWidth");
    SlotHeight_Global = param_reader.getValue<T>("SlottedSphere", "SlotHeight");
    Tf_Global = T(1);
    T_Global = static_cast<T>(Nx) / U0_Global;
    MaxStep = static_cast<int>(Tf_Global * T_Global);
  } else if (isVortex) {
    U0_Global = param_reader.getValue<T>("Vortex", "U0");
    T_Global = param_reader.getValue<T>("Vortex", "T");
    Tf_Global = param_reader.getValue<T>("Vortex", "Tf");
    MaxStep = static_cast<int>(Tf_Global * T_Global);
  } else if (isDeformation) {
    U0_Global = param_reader.getValue<T>("Deformation", "U0");
    T_Global = param_reader.getValue<T>("Deformation", "T");
    Tf_Global = param_reader.getValue<T>("Deformation", "Tf");
    MaxStep = static_cast<int>(Tf_Global * T_Global);
  } else if (isShear) {
    U0_Global = param_reader.getValue<T>("Shear", "U0");
    T_Global = param_reader.getValue<T>("Shear", "T");
    Tf_Global = param_reader.getValue<T>("Shear", "Tf");
    MaxStep = static_cast<int>(Tf_Global * T_Global);
  } else {
    U0_Global = T(0.01);
    Tf_Global = T(1);
    T_Global = static_cast<T>(Nx) / U0_Global;
    MaxStep = static_cast<int>(Tf_Global * T_Global);
  }

  MPI_RANK(0) {
    std::cout << "------------Phase Field 3D Droplet Simulation (MRT):-------------\n" << std::endl;
    std::cout << "[Phase Field Parameters]:" << std::endl;
    std::cout << "  Interface Width:  " << Interface_Width << " lu" << std::endl;
    std::cout << "  Mobility:         " << Mobility << std::endl;
    std::cout << "  Sigma:            " << Sigma << std::endl;
    std::cout << "  (tau/omega/beta/kappa from Converter in main())" << std::endl;
    std::cout << "[Simulation_Settings]:\n"
              << "  TotalStep:        " << MaxStep << "\n"
              << "  OutputStep:       " << OutputStep << "\n"
              << "  Grid:             " << Nx << " x " << Ny << " x " << Nz << "\n"
              << "  U0:               " << U0_Global << "\n"
              << "  TestCase:         " << TestCase_Global << "\n"
              << "  T (half period):  " << T_Global << "\n"
              << "  Tf (half periods):" << Tf_Global << "\n"
              << "  Full period:      " << T(2) * T_Global << "\n"
              << "  MaxStep (auto):   " << MaxStep << "\n"
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

  mpi().init(&argc, &argv);

  MPI_DEBUG_WAIT

  Printer::Print_BigBanner(std::string("Initializing Simple 3D Droplet Simulation (MRT)..."));

  readParam();

  Droplet_Velocity_Global = Droplet_Velocity;
  Mobility_Global = Mobility;

  BaseConverter<T> BaseConv(LatSet::cs2);
  BaseConv.ConvertFromTimeStep(T(1), T(1), T(1), T(Nx), T(1), T(0.1));
  PhaseFieldConverter<T> PhaseConv(BaseConv, LatSet::cs2);
  PhaseConv.fromLattice(Interface_Width, Mobility, Sigma);

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

  AABB<T, 3> domain(Vector<T, 3>(T(0), T(0), T(0)),
                    Vector<T, 3>(T(Nx * Cell_Len), T(Ny * Cell_Len), T(Nz * Cell_Len)));

  AABB<T, 3> left_box(Vector<T, 3>(T(-Cell_Len), T(0), T(0)),
                      Vector<T, 3>(T(0), T(Ny * Cell_Len), T(Nz * Cell_Len)));
  AABB<T, 3> right_box(Vector<T, 3>(T(Nx * Cell_Len), T(0), T(0)),
                       Vector<T, 3>(T((Nx+1) * Cell_Len), T(Ny * Cell_Len), T(Nz * Cell_Len)));
  AABB<T, 3> front_box(Vector<T, 3>(T(0), T(-Cell_Len), T(0)),
                       Vector<T, 3>(T(Nx * Cell_Len), T(0), T(Nz * Cell_Len)));
  AABB<T, 3> back_box(Vector<T, 3>(T(0), T(Ny * Cell_Len), T(0)),
                      Vector<T, 3>(T(Nx * Cell_Len), T((Ny+1) * Cell_Len), T(Nz * Cell_Len)));
  AABB<T, 3> bottom_box(Vector<T, 3>(T(0), T(0), T(-Cell_Len)),
                        Vector<T, 3>(T(Nx * Cell_Len), T(Ny * Cell_Len), T(0)));
  AABB<T, 3> top_box(Vector<T, 3>(T(0), T(0), T(Nz * Cell_Len)),
                     Vector<T, 3>(T(Nx * Cell_Len), T(Ny * Cell_Len), T((Nz+1) * Cell_Len)));

  AABB<T, 3> left_front_box(Vector<T, 3>(T(-Cell_Len), T(-Cell_Len), T(0)),
                            Vector<T, 3>(T(0), T(0), T(Nz * Cell_Len)));
  AABB<T, 3> right_front_box(Vector<T, 3>(T(Nx * Cell_Len), T(-Cell_Len), T(0)),
                             Vector<T, 3>(T((Nx+1) * Cell_Len), T(0), T(Nz * Cell_Len)));
  AABB<T, 3> left_back_box(Vector<T, 3>(T(-Cell_Len), T(Ny * Cell_Len), T(0)),
                           Vector<T, 3>(T(0), T((Ny+1) * Cell_Len), T(Nz * Cell_Len)));
  AABB<T, 3> right_back_box(Vector<T, 3>(T(Nx * Cell_Len), T(Ny * Cell_Len), T(0)),
                            Vector<T, 3>(T((Nx+1) * Cell_Len), T((Ny+1) * Cell_Len), T(Nz * Cell_Len)));
  AABB<T, 3> left_bottom_box(Vector<T, 3>(T(-Cell_Len), T(0), T(-Cell_Len)),
                             Vector<T, 3>(T(0), T(Ny * Cell_Len), T(0)));
  AABB<T, 3> right_bottom_box(Vector<T, 3>(T(Nx * Cell_Len), T(0), T(-Cell_Len)),
                              Vector<T, 3>(T((Nx+1) * Cell_Len), T(Ny * Cell_Len), T(0)));
  AABB<T, 3> left_top_box(Vector<T, 3>(T(-Cell_Len), T(0), T(Nz * Cell_Len)),
                          Vector<T, 3>(T(0), T(Ny * Cell_Len), T((Nz+1) * Cell_Len)));
  AABB<T, 3> right_top_box(Vector<T, 3>(T(Nx * Cell_Len), T(0), T(Nz * Cell_Len)),
                           Vector<T, 3>(T((Nx+1) * Cell_Len), T(Ny * Cell_Len), T((Nz+1) * Cell_Len)));
  AABB<T, 3> front_bottom_box(Vector<T, 3>(T(0), T(-Cell_Len), T(-Cell_Len)),
                              Vector<T, 3>(T(Nx * Cell_Len), T(0), T(0)));
  AABB<T, 3> back_bottom_box(Vector<T, 3>(T(0), T(Ny * Cell_Len), T(-Cell_Len)),
                             Vector<T, 3>(T(Nx * Cell_Len), T((Ny+1) * Cell_Len), T(0)));
  AABB<T, 3> front_top_box(Vector<T, 3>(T(0), T(-Cell_Len), T(Nz * Cell_Len)),
                           Vector<T, 3>(T(Nx * Cell_Len), T(0), T((Nz+1) * Cell_Len)));
  AABB<T, 3> back_top_box(Vector<T, 3>(T(0), T(Ny * Cell_Len), T(Nz * Cell_Len)),
                          Vector<T, 3>(T(Nx * Cell_Len), T((Ny+1) * Cell_Len), T((Nz+1) * Cell_Len)));

  AABB<T, 3> left_front_bottom_box(Vector<T, 3>(T(-Cell_Len), T(-Cell_Len), T(-Cell_Len)),
                                   Vector<T, 3>(T(0), T(0), T(0)));
  AABB<T, 3> right_front_bottom_box(Vector<T, 3>(T(Nx * Cell_Len), T(-Cell_Len), T(-Cell_Len)),
                                    Vector<T, 3>(T((Nx+1) * Cell_Len), T(0), T(0)));
  AABB<T, 3> left_back_bottom_box(Vector<T, 3>(T(-Cell_Len), T(Ny * Cell_Len), T(-Cell_Len)),
                                  Vector<T, 3>(T(0), T((Ny+1) * Cell_Len), T(0)));
  AABB<T, 3> right_back_bottom_box(Vector<T, 3>(T(Nx * Cell_Len), T(Ny * Cell_Len), T(-Cell_Len)),
                                   Vector<T, 3>(T((Nx+1) * Cell_Len), T((Ny+1) * Cell_Len), T(0)));
  AABB<T, 3> left_front_top_box(Vector<T, 3>(T(-Cell_Len), T(-Cell_Len), T(Nz * Cell_Len)),
                                Vector<T, 3>(T(0), T(0), T((Nz+1) * Cell_Len)));
  AABB<T, 3> right_front_top_box(Vector<T, 3>(T(Nx * Cell_Len), T(-Cell_Len), T(Nz * Cell_Len)),
                                 Vector<T, 3>(T((Nx+1) * Cell_Len), T(0), T((Nz+1) * Cell_Len)));
  AABB<T, 3> left_back_top_box(Vector<T, 3>(T(-Cell_Len), T(Ny * Cell_Len), T(Nz * Cell_Len)),
                               Vector<T, 3>(T(0), T((Ny+1) * Cell_Len), T((Nz+1) * Cell_Len)));
  AABB<T, 3> right_back_top_box(Vector<T, 3>(T(Nx * Cell_Len), T(Ny * Cell_Len), T(Nz * Cell_Len)),
                                Vector<T, 3>(T((Nx+1) * Cell_Len), T((Ny+1) * Cell_Len), T((Nz+1) * Cell_Len)));

  BlockGeometryHelper3D<T> GeoHelper(Nx, Ny, Nz, domain, Cell_Len, BlockCellLen);
  GeoHelper.CreateBlocks(1, mpi().getSize(), 1);
  GeoHelper.AdaptiveOptimization(mpi().getSize());
  GeoHelper.LoadBalancing(mpi().getSize());

  BlockGeometry3D<T> Geo(GeoHelper);

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
    bottom_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(
    top_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });

  FlagFM.forEach(
    left_front_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(
    right_front_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(
    left_back_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(
    right_back_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(
    left_bottom_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(
    right_bottom_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(
    left_top_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(
    right_top_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(
    front_bottom_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(
    back_bottom_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(
    front_top_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(
    back_top_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });

  FlagFM.forEach(
    left_front_bottom_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(
    right_front_bottom_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(
    left_back_bottom_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(
    right_back_bottom_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(
    left_front_top_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(
    right_front_top_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(
    left_back_top_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(
    right_back_top_box, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });

  FlagFM.template SetupBoundary<LatSet>(domain, PeriodicFlag);

  using FIELDS = TypePack<PHI<T>, POP<T, LatSet::q>, VELOCITY<T, LatSet::d>, GRAD<T, LatSet::d>, NORMAL<T, LatSet::d>, INTERFACEWIDTH<T>>;
  using CELL = Cell<T, LatSet, FIELDS>;

  ValuePack InitValues(T{}, T{}, T{}, Vector<T, LatSet::d>{T(0), T(0), T(0)}, Vector<T, LatSet::d>{T(0), T(0), T(0)}, Interface_Width);

  BlockLatticeManager<T, LatSet, FIELDS> DropletLattice(Geo, InitValues, PhaseConv);

  auto& phiField = DropletLattice.getField<PHI<T>>();
  T d = static_cast<T>(Nx) * Cell_Len;

  for(int blockid = 0; blockid < Geo.getBlockNum(); ++blockid){
    const auto& block = Geo.getBlock(blockid);
    const auto& proj = block.getProjection();
    auto& blockPhi = phiField.getBlockField(blockid);
    T voxelSize = block.getVoxelSize();
    int overlap = 0;
    for (int k = overlap; k < block.getNz() - overlap; ++k) {
      for (int j = overlap; j < block.getNy() - overlap; ++j) {
        for (int i = overlap; i < block.getNx() - overlap; ++i) {
          std::size_t id = k * proj[2] + j * proj[1] + i;
          T x = block.getMin()[0] + static_cast<T>(i) * voxelSize;
          T y = block.getMin()[1] + static_cast<T>(j) * voxelSize;
          T z = block.getMin()[2] + static_cast<T>(k) * voxelSize;

          T phi;
          if (isSlottedSphere) {
            T half_slot = SlotWidth_Global * T(0.5);
            T slot_bottom = Droplet_Center[1] + Droplet_Radius - SlotHeight_Global;
            T dist = std::sqrt((x - Droplet_Center[0]) * (x - Droplet_Center[0]) + (y - Droplet_Center[1]) * (y - Droplet_Center[1]) + (z - Droplet_Center[2]) * (z - Droplet_Center[2]));
            T phi_sphere = T(0.5) * (T(1) - std::tanh(T(2) * (dist - Droplet_Radius) / Interface_Width));
            T notch_x = T(0.5) * (T(1) - std::tanh(T(2) * (std::abs(x - Droplet_Center[0]) - half_slot) / Interface_Width));
            T notch_y = T(0.5) * (T(1) + std::tanh(T(2) * (y - slot_bottom) / Interface_Width));
            phi = phi_sphere * (T(1) - notch_x * notch_y);
          } else {
            T dist = std::sqrt((x - Droplet_Center[0])*(x - Droplet_Center[0]) +
                              (y - Droplet_Center[1])*(y - Droplet_Center[1]) +
                              (z - Droplet_Center[2])*(z - Droplet_Center[2]));
            phi = T(0.5) * (T(1) - std::tanh(T(2) * (dist - Droplet_Radius) / Interface_Width));
          }
          blockPhi.get(id) = phi;
        }
      }
    }
  }

  BlockFieldManager<PHI<T>, T, LatSet::d> phiInitialFM(phiField);

  auto& velocityField = DropletLattice.getField<VELOCITY<T, LatSet::d>>();
  for (int block_idx = 0; block_idx < Geo.getBlockNum(); ++block_idx) {
    auto& blockLat = DropletLattice.getBlockLat(block_idx);
    auto& blockPhiField = phiField.getBlockField(block_idx);
    auto& blockVelocity = velocityField.getBlockField(block_idx);
    const auto& block = Geo.getBlock(block_idx);
    int overlap = 0;
    const auto& proj = block.getProjection();
    T voxelSize = block.getVoxelSize();

    for (int k = overlap; k < block.getNz() - overlap; ++k) {
      for (int j = overlap; j < block.getNy() - overlap; ++j) {
        for (int i = overlap; i < block.getNx() - overlap; ++i) {
          std::size_t id = k * proj[2] + j * proj[1] + i;
          CELL cell(id, blockLat);
          T phi = blockPhiField.get(id);

          Vector<T, LatSet::d> vel;
          if (isSlottedSphere) {
            T x_norm = (block.getMin()[0] + static_cast<T>(i) * voxelSize) / d;
            T y_norm = (block.getMin()[1] + static_cast<T>(j) * voxelSize) / d;
            T ux = T(2) * M_PI * U0_Global * (T(0.5) - y_norm);
            T uy = T(2) * M_PI * U0_Global * (x_norm - T(0.5));
            T uz = T(0);
            vel = Vector<T, 3>(ux, uy, uz);
          } else if (isVortex) {
            T x_norm = (block.getMin()[0] + static_cast<T>(i) * voxelSize) / d;
            T y_norm = (block.getMin()[1] + static_cast<T>(j) * voxelSize) / d;
            T z_norm = (block.getMin()[2] + static_cast<T>(k) * voxelSize) / d;
            T sin2_pi_x = std::pow(std::sin(M_PI * x_norm), T(2));
            T sin2_pi_y = std::pow(std::sin(M_PI * y_norm), T(2));
            T sin2_pi_z = std::pow(std::sin(M_PI * z_norm), T(2));
            T sin_2pi_x = std::sin(T(2) * M_PI * x_norm);
            T sin_2pi_y = std::sin(T(2) * M_PI * y_norm);
            T sin_2pi_z = std::sin(T(2) * M_PI * z_norm);
            T ux = T(2) * U0_Global * sin2_pi_x * sin_2pi_y * sin_2pi_z;
            T uy = -U0_Global * sin2_pi_y * sin_2pi_z * sin_2pi_x;
            T uz = -U0_Global * sin2_pi_z * sin_2pi_x * sin_2pi_y;
            vel = Vector<T, 3>(ux, uy, uz);
          } else if (isDeformation) {
            T xh = (block.getMin()[0] + static_cast<T>(i) * voxelSize) / d - T(0.5);
            T yh = (block.getMin()[1] + static_cast<T>(j) * voxelSize) / d - T(0.5);
            T zh = (block.getMin()[2] + static_cast<T>(k) * voxelSize) / d - T(0.5);
            T ux = (U0_Global / T(2)) * (std::sin(T(4)*M_PI*xh)*std::sin(T(4)*M_PI*yh) + std::cos(T(4)*M_PI*zh)*std::cos(T(4)*M_PI*xh));
            T uy = (U0_Global / T(2)) * (std::sin(T(4)*M_PI*yh)*std::sin(T(4)*M_PI*zh) + std::cos(T(4)*M_PI*xh)*std::cos(T(4)*M_PI*yh));
            T uz = (U0_Global / T(2)) * (std::sin(T(4)*M_PI*zh)*std::sin(T(4)*M_PI*xh) + std::cos(T(4)*M_PI*yh)*std::cos(T(4)*M_PI*zh));
            vel = Vector<T, 3>(ux, uy, uz);
          } else if (isShear) {
            T xh = (block.getMin()[0] + static_cast<T>(i) * voxelSize) / d - T(0.5);
            T yh = (block.getMin()[1] + static_cast<T>(j) * voxelSize) / d - T(0.5);
            T zh = (block.getMin()[2] + static_cast<T>(k) * voxelSize) / d - T(0.5);
            T ux = M_PI * U0_Global * std::cos(M_PI * xh) * (std::sin(M_PI * zh) - std::sin(M_PI * yh));
            T uy = M_PI * U0_Global * std::cos(M_PI * yh) * (std::sin(M_PI * xh) - std::sin(M_PI * zh));
            T uz = M_PI * U0_Global * std::cos(M_PI * zh) * (std::sin(M_PI * yh) - std::sin(M_PI * xh));
            vel = Vector<T, 3>(ux, uy, uz);
          } else {
            vel = Droplet_Velocity_Global;
          }
          blockVelocity.get(id) = vel;

          T usqr = vel.getnorm2();
          for (unsigned int q = 0; q < LatSet::q; ++q) {
             T cu = latset::c<LatSet>(q) * vel;
             cell[q] = latset::w<LatSet>(q) * phi *
               (T(1) + LatSet::InvCs2 * cu + cu * cu * T(0.5) * LatSet::InvCs4 -
                LatSet::InvCs2 * usqr * T(0.5));
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

  Drop_Periodic.Setup(left_front_box, static_cast<NbrDirection>(NbrDirection::XN | NbrDirection::YN), right_back_box, static_cast<NbrDirection>(NbrDirection::XP | NbrDirection::YP));
  Drop_Periodic.Setup(right_front_box, static_cast<NbrDirection>(NbrDirection::XP | NbrDirection::YN), left_back_box, static_cast<NbrDirection>(NbrDirection::XN | NbrDirection::YP));
  Drop_Periodic.Setup(left_bottom_box, static_cast<NbrDirection>(NbrDirection::XN | NbrDirection::ZN), right_top_box, static_cast<NbrDirection>(NbrDirection::XP | NbrDirection::ZP));
  Drop_Periodic.Setup(right_bottom_box, static_cast<NbrDirection>(NbrDirection::XP | NbrDirection::ZN), left_top_box, static_cast<NbrDirection>(NbrDirection::XN | NbrDirection::ZP));
  Drop_Periodic.Setup(front_bottom_box, static_cast<NbrDirection>(NbrDirection::YN | NbrDirection::ZN), back_top_box, static_cast<NbrDirection>(NbrDirection::YP | NbrDirection::ZP));
  Drop_Periodic.Setup(back_bottom_box, static_cast<NbrDirection>(NbrDirection::YP | NbrDirection::ZN), front_top_box, static_cast<NbrDirection>(NbrDirection::YN | NbrDirection::ZP));

  Drop_Periodic.Setup(left_front_bottom_box, static_cast<NbrDirection>(NbrDirection::XN | NbrDirection::YN | NbrDirection::ZN), right_back_top_box, static_cast<NbrDirection>(NbrDirection::XP | NbrDirection::YP | NbrDirection::ZP));
  Drop_Periodic.Setup(right_front_bottom_box, static_cast<NbrDirection>(NbrDirection::XP | NbrDirection::YN | NbrDirection::ZN), left_back_top_box, static_cast<NbrDirection>(NbrDirection::XN | NbrDirection::YP | NbrDirection::ZP));
  Drop_Periodic.Setup(left_back_bottom_box, static_cast<NbrDirection>(NbrDirection::XN | NbrDirection::YP | NbrDirection::ZN), right_front_top_box, static_cast<NbrDirection>(NbrDirection::XP | NbrDirection::YN | NbrDirection::ZP));
  Drop_Periodic.Setup(right_back_bottom_box, static_cast<NbrDirection>(NbrDirection::XP | NbrDirection::YP | NbrDirection::ZN), left_front_top_box, static_cast<NbrDirection>(NbrDirection::XN | NbrDirection::YN | NbrDirection::ZP));

#ifdef MPI_ENABLED
  Drop_Periodic.SetupMPI(GeoHelper);
#endif

  using NormalTask = tmp::Key_TypePair<BulkFlag, ff::FF3D<CELL>>;
  using NormalTaskCollection = tmp::TupleWrapper<NormalTask>;
  using DropletNormalTask = tmp::TaskSelector<NormalTaskCollection, std::uint8_t, CELL>;

  using BulkTask = tmp::Key_TypePair<BulkFlag,
                                     collision::MRTSource<
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
  vtmo::VectorWriter GRADWriter("GRAD", DropletLattice.getField<GRAD<T, 3>>());
  vtmo::VectorWriter VELOCITYWriter("VELOCITY", DropletLattice.getField<VELOCITY<T, LatSet::d>>());
  vtmo::vtmWriter<T, 3> MainWriter("simpledrop3d", Geo);
  MainWriter.addWriterSet(PHIWriter, GRADWriter, VELOCITYWriter);
  MainWriter.WriteBinary(MainLoopTimer());

  Drop_Periodic.Apply();
  DropletLattice.NormalCommunicate();

  Printer::Print_BigBanner(std::string("Start Simple 3D Droplet Simulation (MRT)..."));

  while (MainLoopTimer() < MaxStep) {
    T t = static_cast<T>(MainLoopTimer());
    T cos_pi_t_over_T = std::cos(M_PI * t / T_Global);

    if (isSlottedSphere) {
      for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
        const auto& block = Geo.getBlock(blockid);
        auto& blockVelocity = velocityField.getBlockField(blockid);
        T voxelSize = block.getVoxelSize();
        const auto& proj = block.getProjection();

        for (int k = 0; k < block.getNz(); ++k) {
          for (int j = 0; j < block.getNy(); ++j) {
            for (int i = 0; i < block.getNx(); ++i) {
              std::size_t id = k * proj[2] + j * proj[1] + i;
              T x_norm = (block.getMin()[0] + static_cast<T>(i) * voxelSize) / d;
              T y_norm = (block.getMin()[1] + static_cast<T>(j) * voxelSize) / d;
              T ux = T(2) * M_PI * U0_Global * (T(0.5) - y_norm);
              T uy = T(2) * M_PI * U0_Global * (x_norm - T(0.5));
              T uz = T(0);
              blockVelocity.get(id) = Vector<T, 3>(ux, uy, uz);
            }
          }
        }
      }
    } else if (isVortex) {
      for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
        const auto& block = Geo.getBlock(blockid);
        auto& blockVelocity = velocityField.getBlockField(blockid);
        T voxelSize = block.getVoxelSize();
        const auto& proj = block.getProjection();

        for (int k = 0; k < block.getNz(); ++k) {
          for (int j = 0; j < block.getNy(); ++j) {
            for (int i = 0; i < block.getNx(); ++i) {
              std::size_t id = k * proj[2] + j * proj[1] + i;
              T x_norm = (block.getMin()[0] + static_cast<T>(i) * voxelSize) / d;
              T y_norm = (block.getMin()[1] + static_cast<T>(j) * voxelSize) / d;
              T z_norm = (block.getMin()[2] + static_cast<T>(k) * voxelSize) / d;
              T sin2_pi_x = std::pow(std::sin(M_PI * x_norm), T(2));
              T sin2_pi_y = std::pow(std::sin(M_PI * y_norm), T(2));
              T sin2_pi_z = std::pow(std::sin(M_PI * z_norm), T(2));
              T sin_2pi_x = std::sin(T(2) * M_PI * x_norm);
              T sin_2pi_y = std::sin(T(2) * M_PI * y_norm);
              T sin_2pi_z = std::sin(T(2) * M_PI * z_norm);
              T ux = T(2) * U0_Global * sin2_pi_x * sin_2pi_y * sin_2pi_z * cos_pi_t_over_T;
              T uy = -U0_Global * sin2_pi_y * sin_2pi_z * sin_2pi_x * cos_pi_t_over_T;
              T uz = -U0_Global * sin2_pi_z * sin_2pi_x * sin_2pi_y * cos_pi_t_over_T;
              blockVelocity.get(id) = Vector<T, 3>(ux, uy, uz);
            }
          }
        }
      }
    } else if (isDeformation) {
      for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
        const auto& block = Geo.getBlock(blockid);
        auto& blockVelocity = velocityField.getBlockField(blockid);
        T voxelSize = block.getVoxelSize();
        const auto& proj = block.getProjection();

        for (int k = 0; k < block.getNz(); ++k) {
          for (int j = 0; j < block.getNy(); ++j) {
            for (int i = 0; i < block.getNx(); ++i) {
              std::size_t id = k * proj[2] + j * proj[1] + i;
              T xh = (block.getMin()[0] + static_cast<T>(i) * voxelSize) / d - T(0.5);
              T yh = (block.getMin()[1] + static_cast<T>(j) * voxelSize) / d - T(0.5);
              T zh = (block.getMin()[2] + static_cast<T>(k) * voxelSize) / d - T(0.5);
              T ux = (U0_Global / T(2)) * (std::sin(T(4)*M_PI*xh)*std::sin(T(4)*M_PI*yh) + std::cos(T(4)*M_PI*zh)*std::cos(T(4)*M_PI*xh)) * cos_pi_t_over_T;
              T uy = (U0_Global / T(2)) * (std::sin(T(4)*M_PI*yh)*std::sin(T(4)*M_PI*zh) + std::cos(T(4)*M_PI*xh)*std::cos(T(4)*M_PI*yh)) * cos_pi_t_over_T;
              T uz = (U0_Global / T(2)) * (std::sin(T(4)*M_PI*zh)*std::sin(T(4)*M_PI*xh) + std::cos(T(4)*M_PI*yh)*std::cos(T(4)*M_PI*zh)) * cos_pi_t_over_T;
              blockVelocity.get(id) = Vector<T, 3>(ux, uy, uz);
            }
          }
        }
      }
    } else if (isShear) {
      for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
        const auto& block = Geo.getBlock(blockid);
        auto& blockVelocity = velocityField.getBlockField(blockid);
        T voxelSize = block.getVoxelSize();
        const auto& proj = block.getProjection();

        for (int k = 0; k < block.getNz(); ++k) {
          for (int j = 0; j < block.getNy(); ++j) {
            for (int i = 0; i < block.getNx(); ++i) {
              std::size_t id = k * proj[2] + j * proj[1] + i;
              T xh = (block.getMin()[0] + static_cast<T>(i) * voxelSize) / d - T(0.5);
              T yh = (block.getMin()[1] + static_cast<T>(j) * voxelSize) / d - T(0.5);
              T zh = (block.getMin()[2] + static_cast<T>(k) * voxelSize) / d - T(0.5);
              T ux = M_PI * U0_Global * std::cos(M_PI * xh) * (std::sin(M_PI * zh) - std::sin(M_PI * yh)) * cos_pi_t_over_T;
              T uy = M_PI * U0_Global * std::cos(M_PI * yh) * (std::sin(M_PI * xh) - std::sin(M_PI * zh)) * cos_pi_t_over_T;
              T uz = M_PI * U0_Global * std::cos(M_PI * zh) * (std::sin(M_PI * yh) - std::sin(M_PI * xh)) * cos_pi_t_over_T;
              blockVelocity.get(id) = Vector<T, 3>(ux, uy, uz);
            }
          }
        }
      }
    }
    DropletLattice.getField<VELOCITY<T, 3>>().Communicate();

    for (int block_idx = 0; block_idx < Geo.getBlockNum(); ++block_idx) {
        auto& blockLat = DropletLattice.getBlockLat(block_idx);
        auto& blockPhiField = phiField.getBlockField(block_idx);
        for (std::size_t id = 0; id < blockLat.getVoxNum(); ++id) {
          CELL cell(id, blockLat);
          T phi_new = T(0);
          for (unsigned int q = 0; q < LatSet::q; ++q) {
              phi_new += cell[q];
          }
          blockPhiField.get(id) = phi_new;
        }
    }

    DropletLattice.getField<PHI<T>>().Communicate();

    DropletLattice.template ApplyCellDynamics<DropletNormalTask>(FlagFM);
    DropletLattice.getField<NORMAL<T, LatSet::d>>().Communicate();
    DropletLattice.template ApplyCellDynamics<DropletCollisionTask>(FlagFM);

    DropletLattice.Stream();

    Drop_Periodic.Apply();

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

  Printer::Print_BigBanner(std::string("Simple 3D Droplet Simulation (MRT) Complete!"));
  MainWriter.WriteBinary(MainLoopTimer());
  MainLoopTimer.Print_MainLoopPerformance(Geo.getTotalCellNum());
  Printer::Print("Total PhysTime", BaseConv.getPhysTime(MainLoopTimer()));
  Printer::Endl();

  T l2_sum = T(0);
  T L0 = static_cast<T>(Nx);
  for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
    auto& blockPhi = phiField.getBlockField(blockid);
    auto& blockPhiInit = phiInitialFM.getBlockField(blockid);
    const auto& block = Geo.getBlock(blockid);
    int overlap = block.getOverlap();
    const auto& proj = block.getProjection();
    for (int k = overlap; k < block.getNz() - overlap; ++k) {
      for (int j = overlap; j < block.getNy() - overlap; ++j) {
        for (int i = overlap; i < block.getNx() - overlap; ++i) {
          std::size_t id = k * proj[2] + j * proj[1] + i;
          T diff = blockPhi.get(id) - blockPhiInit.get(id);
          l2_sum += diff * diff;
        }
      }
    }
  }
  T l2_norm = std::sqrt(l2_sum);
  T l2_error = l2_norm / (L0 * L0);
  IF_MPI_RANK(0) {
    std::cout << "============ L2 Error ============" << std::endl;
    std::cout << "  L2 norm:           " << l2_norm << std::endl;
    std::cout << "  L0^2:              " << L0 * L0 << std::endl;
    std::cout << "  L2 error (paper):  " << l2_error << std::endl;
    std::cout << "==================================" << std::endl;
  }

  return 0;
}
