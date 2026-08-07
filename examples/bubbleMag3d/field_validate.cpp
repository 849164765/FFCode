// field_validate.cpp — 磁场解析解验证: 静止球体在匀强磁场中的磁标势/磁场
//
// 与 bubbleMag3d 相同的网格/参数/ψ求解器, 但主循环只做 MF 求解 (球体冻结):
//   0a 系数更新 -> 0b 壁面钉扎 -> 子迭代(碰撞+迁移+宏量) -> H 计算 -> 同步
// 输出 PHI/PSI/HX/HY/HZ/HMAG, 供 validate_field.py 与解析解对比:
//
// 解析解 (μ_in=1 球体在 μ_out=9 铁磁流体, 均匀场 H0 沿 +z):
//   球内 (r<R):  H = 3μ_out/(μ_in+2μ_out)·H0 ẑ = 1.42105·H0 ẑ (均匀)
//   球外 (r>R):  ψ = -H0·z + C·H0·R³·z/r³,  C = (μ_in-μ_out)/(μ_in+2μ_out) = -8/19
//                H_z(r轴) = H0·(1 - 16/19·(R/r)³)   (极轴上方场减弱)
//                H_r(赤道) = H0·(1 + 8/19·(R/r)³)    (赤道场增强)
//   界面边界条件: H 切向连续, μH 法向连续 (解析解自动满足)
#include "freelb.h"
#include "freelb.hh"
#include "ff/ff2d.h"

using T = FLOAT;
using LatSet = D3Q19<T>;
using MFLatSet = D3Q7<T>;
using namespace mfield;

// ---- Simulation Parameters ----
int Ni, Nj, Nk;
T Cell_Len;
int BlockCellLen, Thread_Num;

// bubble
T Bubble_Radius;
Vector<T, 3> Bubble_Center;

// phase field
T Interface_Width, Mobility, Tau_phi, Omega_phi, Kappa, Beta;

// two-phase
T rho_l, rho_h, eta_l, eta_h, sigma, gravity, Eo, Re, U_g, Tau_ns;

// magnetic field
T chi_l, chi_h, mu_l, mu_h, H0, Bom, DeltaRho;

int PsiSolver_Iter; T PsiSolver_K;

int MaxStep, OutputStep;
std::string work_dir;

void readParam() {
  iniReader r("bubbleMag3d.ini");
  work_dir = r.getValue<std::string>("workdir","workdir_");
  Thread_Num = r.getValue<int>("parallel","thread_num");
  Ni = r.getValue<int>("Mesh","Ni"); Nj = r.getValue<int>("Mesh","Nj");
  Nk = r.getValue<int>("Mesh","Nk");
  Cell_Len = r.getValue<T>("Mesh","Cell_Len");
  BlockCellLen = r.getValue<int>("Mesh","BlockCellLen");
  Bubble_Radius = r.getValue<T>("Bubble","Radius");
  Bubble_Center[0]=r.getValue<T>("Bubble","CenterX");
  Bubble_Center[1]=r.getValue<T>("Bubble","CenterY");
  Bubble_Center[2]=r.getValue<T>("Bubble","CenterZ");
  Interface_Width=r.getValue<T>("Phase_Field","Interface_Width");
  Mobility=r.getValue<T>("Phase_Field","Mobility");
  rho_l=r.getValue<T>("Two_Phase","rho_l");
  rho_h=r.getValue<T>("Two_Phase","rho_h");
  Eo=r.getValue<T>("Two_Phase","Eo");
  Re=r.getValue<T>("Two_Phase","Re");
  U_g=r.getValue<T>("Two_Phase","U_g");
  MaxStep=r.getValue<int>("Simulation_Settings","TotalStep");
  OutputStep=r.getValue<int>("Simulation_Settings","OutputStep");
  // Magnetic
  chi_l=r.getValue<T>("Magnetic_Field","chi_l");
  chi_h=r.getValue<T>("Magnetic_Field","chi_h");
  mu_l=r.getValue<T>("Magnetic_Field","mu_l");
  mu_h=r.getValue<T>("Magnetic_Field","mu_h");
  Bom=r.getValue<T>("Magnetic_Field","Bom");
  PsiSolver_Iter=r.getValue<int>("Magnetic_Field","PsiSolver_Iter");
  PsiSolver_K=r.getValue<T>("Magnetic_Field","PsiSolver_K");

  eta_l=T(0.0568)/T(100.0); eta_h=T(0.0568);
  T D_bubble = T{2} * Bubble_Radius;
  T gy_sqrt = Re * eta_h / (rho_h * std::pow(D_bubble, T{1.5}));
  gravity = gy_sqrt * gy_sqrt;
  sigma = gravity * rho_h * D_bubble * D_bubble / Eo;
  H0 = std::sqrt(T(2.0) * Bom * sigma / D_bubble);
  Beta=T(12.0)*sigma/Interface_Width;
  Kappa=T(3.0)*Interface_Width*sigma*T{0.5};
  Tau_phi=T(3.0)*Mobility+T(0.5); Omega_phi=T(1.0)/Tau_phi;
  Tau_ns=T(0.5)+eta_h/rho_h/LatSet::cs2;
  DeltaRho=rho_h-rho_l;

  MPI_RANK(0){
    printf("---- Magnetic Field Validation (static sphere in uniform H) ----\n");
    printf("Mesh: %dx%dx%d  R=%.1f center=(%.0f,%.0f,%.0f)  W=%.1f\n",
           Ni,Nj,Nk,Bubble_Radius,Bubble_Center[0],Bubble_Center[1],Bubble_Center[2],Interface_Width);
    printf("mu: in=%.1f out=%.1f  H0=%.5f  PsiSolver: K=%.3f iter=%d\n",
           mu_l,mu_h,H0,PsiSolver_K,PsiSolver_Iter);
    printf("解析解: H_in=%.5f (%.5f·H0), C=%.5f\n",
           T{3}*mu_h/(mu_l+T{2}*mu_h)*H0, T{3}*mu_h/(mu_l+T{2}*mu_h),
           (mu_l-mu_h)/(mu_l+T{2}*mu_h));
    printf("------------------------------------------\n");
  }
}

int main(int argc, char* argv[]) {
  constexpr std::uint8_t VoidFlag=1,BulkFlag=2,BouncebackFlag=4,PeriodicFlag=8;
  mpi().init(&argc,&argv); MPI_DEBUG_WAIT
  Printer::Print_BigBanner(std::string("Magnetic Field Validation..."));
  readParam();

  // -- converters --
  BaseConverter<T> BaseConv(LatSet::cs2);
  BaseConv.SimplifiedConverterFromRT(Ni,T(0.01),Tau_ns);
  BaseConverter<T> PFBaseConv(LatSet::cs2);
  PFBaseConv.SimplifiedConverterFromRT(Ni,T(0.01),Tau_phi);
  BaseConverter<T> MFBaseConv(MFLatSet::cs2);
  MFBaseConv.SimplifiedConverterFromRT(Ni,T(0.01),T(1.0));
  UnitConvManager<T> ConvManager(&BaseConv); ConvManager.Check_and_Print();

  // -- geometry --
  AABB<T,3> domain({0,0,0},{T(Ni*Cell_Len),T(Nj*Cell_Len),T(Nk*Cell_Len)});
  AABB<T,3> left  ({T(-Cell_Len),0,0},{0,T(Nj*Cell_Len),T(Nk*Cell_Len)});
  AABB<T,3> right ({T(Ni*Cell_Len),0,0},{T((Ni+1)*Cell_Len),T(Nj*Cell_Len),T(Nk*Cell_Len)});
  AABB<T,3> front ({0,T(-Cell_Len),0},{T(Ni*Cell_Len),0,T(Nk*Cell_Len)});
  AABB<T,3> back  ({0,T(Nj*Cell_Len),0},{T(Ni*Cell_Len),T((Nj+1)*Cell_Len),T(Nk*Cell_Len)});
  BlockGeometryHelper3D<T> GeoHelper(Ni,Nj,Nk,domain,Cell_Len,BlockCellLen);
  GeoHelper.CreateBlocks(2,2,4);
  GeoHelper.AdaptiveOptimization(mpi().getSize());
  GeoHelper.LoadBalancing(mpi().getSize());
  BlockGeometry3D<T> Geo(GeoHelper);

  // ---- 全局 block 表 (跨 rank 周期 ghost 同步, 同 bubbleMag3d) ----
  std::vector<double> GlobalBlockTable;
  {
    const int nranks = mpi().getSize();
    const T xLeft = T{0}, xRight = T(Ni) * Cell_Len;
    const T yFront = T{0}, yBack = T(Nj) * Cell_Len;
    std::vector<double> local;
    for (int b = 0; b < Geo.getBlockNum(); ++b) {
      const auto& bk = Geo.getBlock(b);
      local.push_back(double(mpi().getRank()));
      local.push_back(bk.getMin()[0]); local.push_back(bk.getMin()[1]); local.push_back(bk.getMin()[2]);
      local.push_back(bk.getNx()); local.push_back(bk.getNy()); local.push_back(bk.getNz());
      local.push_back(bk.getMin()[0] < xLeft  + Cell_Len * T{0.5} ? 1.0 : 0.0);
      local.push_back(bk.getMax()[0] > xRight - Cell_Len * T{0.5} ? 1.0 : 0.0);
      local.push_back(bk.getMin()[1] < yFront + Cell_Len * T{0.5} ? 1.0 : 0.0);
      local.push_back(bk.getMax()[1] > yBack  - Cell_Len * T{0.5} ? 1.0 : 0.0);
    }
    std::vector<int> counts(nranks, 0);
    {
      double tmp[1], dummy[1];
      if (mpi().getRank() == 0) {
        for (int r = 0; r < nranks; ++r) {
          tmp[0] = double(local.size());
          mpi().sendRecv(tmp, dummy, 1, r, r, 1000);
          counts[r] = int(dummy[0]);
        }
      } else {
        tmp[0] = double(local.size());
        mpi().sendRecv(tmp, dummy, 1, 0, 0, 1000);
      }
    }
    mpi().bCast(counts.data(), nranks, 0);
    int total = 0;
    for (int c : counts) total += c;
    std::vector<double> all(total, 0.0);
    {
      std::vector<double> dummy;
      if (mpi().getRank() == 0) {
        std::copy(local.begin(), local.end(), all.begin());
        int off = counts[0];
        for (int r = 1; r < nranks; ++r) {
          dummy.assign(counts[r], 0.0);
          mpi().sendRecv(dummy.data(), all.data() + off, counts[r], r, r, 1001);
          off += counts[r];
        }
      } else {
        dummy.assign(counts[mpi().getRank()], 0.0);
        mpi().sendRecv(local.data(), dummy.data(), int(local.size()), 0, 0, 1001);
      }
    }
    mpi().bCast(all.data(), total, 0);
    GlobalBlockTable = all;
  }

  // -- flag --
  BlockFieldManager<FLAG,T,3> FlagFM(Geo,VoidFlag);
  FlagFM.forEach(domain,[&](FLAG&f,std::size_t id){f.SetField(id,BulkFlag);});
  FlagFM.forEach(left,[&](FLAG&f,std::size_t id){f.SetField(id,PeriodicFlag);});
  FlagFM.forEach(right,[&](FLAG&f,std::size_t id){f.SetField(id,PeriodicFlag);});
  FlagFM.forEach(front,[&](FLAG&f,std::size_t id){f.SetField(id,PeriodicFlag);});
  FlagFM.forEach(back,[&](FLAG&f,std::size_t id){f.SetField(id,PeriodicFlag);});
  FlagFM.template SetupBoundary<LatSet>(domain,BouncebackFlag);

  // -- NS lattice (仅作为 PF 的速度参考场, 不参与演化) --
  using NSFIELDS=TypePack<DENSITY<T>,VELOCITY<T,3>,POP<T,LatSet::q>,FORCE<T,3>,OMEGA<T>,PRESSURE<T>>;
  T omega_ns=T{1}/Tau_ns;
  ValuePack NSI(T{1},Vector<T,3>{0,0,0},T{},Vector<T,3>{0,0,0},omega_ns,T{});
  using NSCELL=Cell<T,LatSet,NSFIELDS>;
  BlockLatticeManager<T,LatSet,NSFIELDS> NSLattice(Geo,NSI,BaseConv);

  // -- PF lattice (仅存储冻结的 φ 球体) --
  using PFFIELDS=TypePack<PHI<T>,POP<T,LatSet::q>,GRAD<T,3>,NORMAL<T,3>,INTERFACEWIDTH<T>,
    ff::LAPLACIAN<T>,ff::CHEMICALPOTENTIAL<T>,
    ff::GRAVITY<T>,ff::BETA<T>,ff::KAPPA<T>,
    ff::RHO_L<T>,ff::RHO_H<T>,ff::ETA_L<T>,ff::ETA_H<T>,ff::DELTARHO<T>>;
  using PFREF=TypePack<VELOCITY<T,3>>;
  using PFPACK=TypePack<PFFIELDS,PFREF>;
  ValuePack PFI(T{},T{},Vector<T,3>{0,0,0},Vector<T,3>{0,0,0},Interface_Width,
    T{},T{},gravity,Beta,Kappa,rho_l,rho_h,eta_l,eta_h,DeltaRho);
  using PFCELL=Cell<T,LatSet,ExtractFieldPack<PFPACK>::mergedpack>;
  BlockLatticeManager<T,LatSet,PFPACK> PFLattice(Geo,PFI,PFBaseConv,
    &NSLattice.getField<VELOCITY<T,3>>());

  ff::BroadcastAllParams<T>(PFLattice,rho_l,rho_h,eta_l,eta_h,gravity,Beta,Kappa);
  PFLattice.template getField<ff::DELTARHO<T>>().InitValue(DeltaRho);

  // -- MF lattice (D3Q7) --
  using MFFIELDS=TypePack<PSI<T>,OMEGA_PSI<T>,MU_PERCELL<T>,CHI_PERCELL<T>,
    HX<T>,HY<T>,HZ<T>,HMAG<T>,POP<T,MFLatSet::q>,
    MU_L<T>,MU_H<T>,CHI_L<T>,CHI_H<T>,H_0<T>,PSI_K<T>>;
  using MFREF=TypePack<PHI<T>>;
  using MFPACK=TypePack<MFFIELDS,MFREF>;
  ValuePack MFI(T{},T{1.0},T{mu_l},T{chi_l},T{},T{},T{},T{},T{},T{},
    mu_l,mu_h,chi_l,chi_h,H0,PsiSolver_K);
  using MFCELL=Cell<T,MFLatSet,ExtractFieldPack<MFPACK>::mergedpack>;
  BlockLatticeManager<T,MFLatSet,MFPACK> MFLattice(Geo,MFI,MFBaseConv,
    &PFLattice.getField<PHI<T>>());
  BroadcastAllMFParams3D<T>(MFLattice,mu_l,mu_h,chi_l,chi_h,H0,PsiSolver_K);
  MFLattice.getField<OMEGA_PSI<T>>().InitValue(T{1.0});

  // -- init phi (冻结的 tanh 球体) --
  T R_phys=Bubble_Radius*Cell_Len, xc=Bubble_Center[0]*Cell_Len;
  T yc=Bubble_Center[1]*Cell_Len, zc=Bubble_Center[2]*Cell_Len;
  T W_phys=Interface_Width*Cell_Len;
  auto& phiField=PFLattice.getField<PHI<T>>();
  for(int b=0;b<Geo.getBlockNum();++b){
    const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
    auto& bPhi=phiField.getBlockField(b);
    T vs=bk.getVoxelSize(),mx=bk.getMin()[0],my=bk.getMin()[1],mz=bk.getMin()[2];
    int ov=0;
    for(int k=ov;k<bk.getNz()-ov;++k){
      T z=mz+T(k)*vs, dz=z-zc;
      for(int j=ov;j<bk.getNy()-ov;++j){
        T y=my+T(j)*vs, dy=y-yc;
        for(int i=ov;i<bk.getNx()-ov;++i){
          T x=mx+T(i)*vs, dx=x-xc;
          T dist=std::sqrt(dx*dx+dy*dy+dz*dz);
          T phi=T{0.5}+T{0.5}*std::tanh(T{2.0}*(dist-R_phys)/W_phys);
          bPhi.get(k*pr[2]+j*pr[1]+i)=phi;
        }
      }
    }
  }

  // init MF: psi = -H0*z (uniform field along +z), g_k = w_k*psi
  {
    auto& psiF=MFLattice.getField<PSI<T>>();
    for(int b=0;b<Geo.getBlockNum();++b){
      auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
      auto& bPsi=psiF.getBlockField(b);
      T vs=bk.getVoxelSize(),mz=bk.getMin()[2];
      for(int k=0;k<bk.getNz();++k){
        T z=mz+T(k)*vs, psi=-H0*z;
        for(int j=0;j<bk.getNy();++j){
          for(int i=0;i<bk.getNx();++i){
            std::size_t id=k*pr[2]+j*pr[1]+i; MFCELL c(id,bl);
            for(unsigned kk=0;kk<MFLatSet::q;++kk) c[kk]=latset::w<MFLatSet>(kk)*psi;
            bPsi.get(id)=psi;
          }
        }
      }
    }
  }

  // -- BCs: 仅 MF 周期边界 (x/y) --
  using LM_MF=BlockLatticeManager<T,MFLatSet,MFPACK>;
  using FM=BlockFieldManager<FLAG,T,3>;

  FixedPeriodicBoundaryManager<LM_MF,FM> MF_Per("MF_Per",MFLattice,FlagFM,PeriodicFlag,VoidFlag);
  MF_Per.Setup(left,NbrDirection::XN,right,NbrDirection::XP);
  MF_Per.Setup(right,NbrDirection::XP,left,NbrDirection::XN);
  MF_Per.Setup(front,NbrDirection::YN,back,NbrDirection::YP);
  MF_Per.Setup(back,NbrDirection::YP,front,NbrDirection::YN);
#ifdef MPI_ENABLED
  MF_Per.SetupMPI(GeoHelper);
#endif

  // MF coupling: PF→MF (coeff update)
  using MCT=tmp::Key_TypePair<BulkFlag,MFUpdateCoeffs3D<PFCELL,MFCELL>>;
  using MCSel=CoupledTaskSelector<std::uint8_t,PFCELL,MFCELL,MCT>;
  BlockLatManagerCoupling MCC(PFLattice,MFLattice);

  // Writers
  vtmo::ScalarWriter PW("PHI",PFLattice.getField<PHI<T>>());
  vtmo::ScalarWriter PS("PSI",MFLattice.getField<PSI<T>>());
  vtmo::ScalarWriter Hxw("HX",MFLattice.getField<HX<T>>());
  vtmo::ScalarWriter Hyw("HY",MFLattice.getField<HY<T>>());
  vtmo::ScalarWriter Hzw("HZ",MFLattice.getField<HZ<T>>());
  vtmo::ScalarWriter Hmw("HMAG",MFLattice.getField<HMAG<T>>());
  vtmo::vtmWriter<T,3> MW("fieldcheck",Geo);
  MW.addWriterSet(PW,PS,Hxw,Hyw,Hzw,Hmw);

  // ===== initial setup =====
  PFLattice.NormalFullCommunicate(); NSLattice.NormalFullCommunicate(); MFLattice.NormalFullCommunicate();

  Printer::Print_BigBanner(std::string("Start MF Solve (sphere frozen)..."));
  Timer t; Timer ot;
  T H_global=T(Nk)*Cell_Len;

  while(t()<MaxStep){
    // 0a: Update per-cell mu, chi, omega_psi from phi (冻结球体, 每步刷新)
    MCC.ApplyInnerCellDynamics<MCSel>(t(),FlagFM);
    CommunicateOMEGAPSI<T>(MFLattice);

    // 0b: Set wall psi = -H0*z  and wall pops = feq(psi_bc)
    {
      auto& psiF=MFLattice.getField<PSI<T>>();
      const T nwall=Cell_Len*T{3.0};
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bPsi=psiF.getBlockField(b);
        int nx=bk.getNx(),ny=bk.getNy(),nz=bk.getNz(),ov=bk.getOverlap();
        T minZ=bk.getMin()[2],vs=bk.getVoxelSize();
        for(int kk=0;kk<nz;++kk){
          T z=minZ+T(kk-ov)*vs;
          if((z<=nwall&&z>=-nwall)||(z<=H_global+nwall&&z>=H_global-nwall)){
            T psi_w=-H0*z;
            for(int jj=0;jj<ny;++jj){
              for(int ii=0;ii<nx;++ii){
                std::size_t id=kk*pr[2]+jj*pr[1]+ii; MFCELL c(id,bl);
                for(unsigned k2=0;k2<MFLatSet::q;++k2) c[k2]=latset::w<MFLatSet>(k2)*psi_w;
                bPsi.get(id)=psi_w;
              }
            }
          }
        }
      }
    }

    // SyncMFPeriodicGhosts: 跨 rank 周期 ghost 同步 (同 bubbleMag3d)
    auto SyncMFPeriodicGhosts = [&](auto& field, int fidx) {
      const T xLeft = T{0};
      const T xRight = T(Ni) * Cell_Len;
      const T yFront = T{0};
      const T yBack = T(Nj) * Cell_Len;
      const int nEntry = int(GlobalBlockTable.size() / 11);
      const int myRank = mpi().getRank();
      for (int b = 0; b < Geo.getBlockNum(); ++b) {
        const auto& bk = Geo.getBlock(b);
        const bool atL = bk.getMin()[0] < xLeft + Cell_Len * T{0.5};
        const bool atR = bk.getMax()[0] > xRight - Cell_Len * T{0.5};
        const bool atF = bk.getMin()[1] < yFront + Cell_Len * T{0.5};
        const bool atB = bk.getMax()[1] > yBack - Cell_Len * T{0.5};
        if (!atL && !atR && !atF && !atB) continue;
        const auto& pr = bk.getProjection();
        const int nx = bk.getNx(), ny = bk.getNy(), nz = bk.getNz();
        auto& f1 = field.getBlockField(b);
        if (atL || atR) {
          int eL = -1, eP = -1;
          for (int e = 0; e < nEntry; ++e) {
            const double* en = &GlobalBlockTable[e * 11];
            if (en[2] != bk.getMin()[1] || en[3] != bk.getMin()[2]) continue;
            if (atL && en[8] > T{0.5}) eP = e;
            if (atR && en[7] > T{0.5}) eP = e;
            if (int(en[0]) == myRank && en[1] == bk.getMin()[0] &&
                en[2] == bk.getMin()[1] && en[3] == bk.getMin()[2])
              eL = e;
          }
          if (eP >= 0 && eL >= 0) {
            const double* enP = &GlobalBlockTable[eP * 11];
            const int pRank = int(enP[0]);
            const int plane = ny * nz;
            if (pRank == myRank) {
              int pb = -1;
              for (int bb = 0; bb < Geo.getBlockNum(); ++bb) {
                const auto& bkk = Geo.getBlock(bb);
                if (bkk.getMin()[0] == enP[1] && bkk.getMin()[1] == enP[2] &&
                    bkk.getMin()[2] == enP[3]) { pb = bb; break; }
              }
              if (pb >= 0) {
                auto& f2 = field.getBlockField(pb);
                const int pnx2 = Geo.getBlock(pb).getNx();
                if (atL) for (int kk = 0; kk < nz; ++kk)
                  for (int jj = 0; jj < ny; ++jj)
                    f1.get(kk * pr[2] + jj * pr[1] + 0) = f2.get(kk * pr[2] + jj * pr[1] + pnx2 - 2);
                if (atR) for (int kk = 0; kk < nz; ++kk)
                  for (int jj = 0; jj < ny; ++jj)
                    f1.get(kk * pr[2] + jj * pr[1] + nx - 1) = f2.get(kk * pr[2] + jj * pr[1] + 1);
              }
            } else {
              const int tag = fidx * 1000 + std::min(eL, eP) * 2 + 0;
              const int srcI = atL ? 1 : (nx - 2);
              const int dstI = atL ? 0 : (nx - 1);
              std::vector<T> sbuf(plane), rbuf(plane);
              int kk = 0;
              for (int k = 0; k < nz; ++k)
                for (int j = 0; j < ny; ++j)
                  sbuf[kk++] = f1.get(k * pr[2] + j * pr[1] + srcI);
              mpi().sendRecv(sbuf.data(), rbuf.data(), plane, pRank, pRank, tag);
              kk = 0;
              for (int k = 0; k < nz; ++k)
                for (int j = 0; j < ny; ++j)
                  f1.get(k * pr[2] + j * pr[1] + dstI) = rbuf[kk++];
            }
          }
        }
        if (atF || atB) {
          int eL = -1, eP = -1;
          for (int e = 0; e < nEntry; ++e) {
            const double* en = &GlobalBlockTable[e * 11];
            if (en[1] != bk.getMin()[0] || en[3] != bk.getMin()[2]) continue;
            if (atF && en[10] > T{0.5}) eP = e;
            if (atB && en[9]  > T{0.5}) eP = e;
            if (int(en[0]) == myRank && en[1] == bk.getMin()[0] &&
                en[2] == bk.getMin()[1] && en[3] == bk.getMin()[2])
              eL = e;
          }
          if (eP >= 0 && eL >= 0) {
            const double* enP = &GlobalBlockTable[eP * 11];
            const int pRank = int(enP[0]);
            const int plane = nx * nz;
            if (pRank == myRank) {
              int pb = -1;
              for (int bb = 0; bb < Geo.getBlockNum(); ++bb) {
                const auto& bkk = Geo.getBlock(bb);
                if (bkk.getMin()[0] == enP[1] && bkk.getMin()[1] == enP[2] &&
                    bkk.getMin()[2] == enP[3]) { pb = bb; break; }
              }
              if (pb >= 0) {
                auto& f2 = field.getBlockField(pb);
                const int pny2 = Geo.getBlock(pb).getNy();
                if (atF) for (int kk = 0; kk < nz; ++kk)
                  for (int ii = 0; ii < nx; ++ii)
                    f1.get(kk * pr[2] + 0 * pr[1] + ii) = f2.get(kk * pr[2] + (pny2 - 2) * pr[1] + ii);
                if (atB) for (int kk = 0; kk < nz; ++kk)
                  for (int ii = 0; ii < nx; ++ii)
                    f1.get(kk * pr[2] + (ny - 1) * pr[1] + ii) = f2.get(kk * pr[2] + 1 * pr[1] + ii);
              }
            } else {
              const int tag = fidx * 1000 + std::min(eL, eP) * 2 + 1;
              const int srcJ = atF ? 1 : (ny - 2);
              const int dstJ = atF ? 0 : (ny - 1);
              std::vector<T> sbuf(plane), rbuf(plane);
              int kk = 0;
              for (int k = 0; k < nz; ++k)
                for (int i = 0; i < nx; ++i)
                  sbuf[kk++] = f1.get(k * pr[2] + srcJ * pr[1] + i);
              mpi().sendRecv(sbuf.data(), rbuf.data(), plane, pRank, pRank, tag);
              kk = 0;
              for (int k = 0; k < nz; ++k)
                for (int i = 0; i < nx; ++i)
                  f1.get(k * pr[2] + dstJ * pr[1] + i) = rbuf[kk++];
            }
          }
        }
      }
    };

    // 0c-0e+: ψ 求解器子迭代 (同 bubbleMag3d)
    for(int sub=0;sub<PsiSolver_Iter;++sub){
      // 0c: MF collision
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        int ov=bk.getOverlap();
        for(int k=ov;k<bk.getNz()-ov;++k){
          for(int j=ov;j<bk.getNy()-ov;++j){
            for(int i=ov;i<bk.getNx()-ov;++i){
              MFCELL c(k*pr[2]+j*pr[1]+i,bl);
              collision::MRTDiffusion<MFCELL,OMEGA_PSI<T>>::apply(c);
            }
          }
        }
      }
      MF_Per.Apply();
      // 0d: Stream
      MFLattice.NormalFullCommunicate();
      MFLattice.Stream();
      MFLattice.NormalFullCommunicate();
      // 0e: PSI = Σg_i
      {
        auto& psiF=MFLattice.getField<PSI<T>>();
        for(int b=0;b<Geo.getBlockNum();++b){
          auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
          auto& bPsi=psiF.getBlockField(b); int ov=bk.getOverlap();
          for(int k=ov;k<bk.getNz()-ov;++k){
            for(int j=ov;j<bk.getNy()-ov;++j){
              for(int i=ov;i<bk.getNx()-ov;++i){
                std::size_t id=k*pr[2]+j*pr[1]+i; MFCELL c(id,bl);
                T psi=0; for(unsigned k2=0;k2<MFLatSet::q;++k2)psi+=c[k2];
                bPsi.get(id)=psi;
              }
            }
          }
        }
      }
      CommunicatePSI<T>(MFLattice);
      // 0e+: re-pin walls
      {
        auto& psiF=MFLattice.getField<PSI<T>>();
        const T nwall=Cell_Len*T{3.0};
        for(int b=0;b<Geo.getBlockNum();++b){
          auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
          auto& bPsi=psiF.getBlockField(b);
          int nx=bk.getNx(),ny=bk.getNy(),nz=bk.getNz(),ov=bk.getOverlap();
          T minZ=bk.getMin()[2],vs=bk.getVoxelSize();
          for(int kk=0;kk<nz;++kk){
            T z=minZ+T(kk-ov)*vs;
            if((z<=nwall&&z>=-nwall)||(z<=H_global+nwall&&z>=H_global-nwall)){
              T psi_w=-H0*z;
              for(int jj=0;jj<ny;++jj){
                for(int ii=0;ii<nx;++ii){
                  std::size_t id=kk*pr[2]+jj*pr[1]+ii; MFCELL c(id,bl);
                  for(unsigned k2=0;k2<MFLatSet::q;++k2) c[k2]=latset::w<MFLatSet>(k2)*psi_w;
                  bPsi.get(id)=psi_w;
                }
              }
            }
          }
        }
      }
    }
    SyncMFPeriodicGhosts(MFLattice.getField<PSI<T>>(), 0);

    // 0f: Compute H = -∇ψ
    for(int b=0;b<Geo.getBlockNum();++b){
      auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
      int ov=bk.getOverlap();
      for(int k=ov;k<bk.getNz()-ov;++k){
        for(int j=ov;j<bk.getNy()-ov;++j){
          for(int i=ov;i<bk.getNx()-ov;++i){
            MFCELL c(k*pr[2]+j*pr[1]+i,bl);
            MFComputeH3D<MFCELL>::apply(c);
          }
        }
      }
    }
    CommunicateAllMFFields3D<T>(MFLattice);
    SyncMFPeriodicGhosts(MFLattice.getField<HX<T>>(), 1);
    SyncMFPeriodicGhosts(MFLattice.getField<HY<T>>(), 2);
    SyncMFPeriodicGhosts(MFLattice.getField<HZ<T>>(), 3);
    SyncMFPeriodicGhosts(MFLattice.getField<HMAG<T>>(), 4);

    ++t; ++ot;
    if(t()%OutputStep==0){
      ot.Print_InnerLoopPerformance(Geo.getTotalCellNum(),OutputStep);
      Printer::Endl();
      MW.WriteBinary(t());
    }
  }

  Printer::Print_BigBanner(std::string("Field Validation Complete!"));
  t.Print_MainLoopPerformance(Geo.getTotalCellNum());
  return 0;
}
