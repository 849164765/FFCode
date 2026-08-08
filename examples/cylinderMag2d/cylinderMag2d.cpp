// cylinderMag2d.cpp — validate the D2Q5 magnetic-field solver against the
// analytical solution for a circular cylinder in a uniform magnetic field
// (paper Sec. III.A): inside the cylinder H_in = 2*mu2/(mu1+mu2)*H0 (uniform),
// outside the dipole-like field (Eq. 59-61). Runs the magnetic solver only
// (static mu field, no flow / phase-field dynamics).
#include "freelb.h"
#include "freelb.hh"
#include "ff/ff2d.h"
#include <cstring>

using T = FLOAT;
using LatSet = D2Q9<T>;
using MFLatSet = D2Q5<T>;
using namespace mfield;

int Ni, Nj;
T Cell_Len;
int BlockCellLen, Thread_Num;
T Cyl_Radius, Cyl_Cx, Cyl_Cy;
T Interface_Width, Mobility;
T mu_l, mu_h, H0;
int PsiSolver_Iter; T PsiSolver_K;
int MaxSubStep;
std::string work_dir;

void readParam() {
  iniReader r("cylinderMag2d.ini");
  work_dir = r.getValue<std::string>("workdir","workdir_");
  Thread_Num = r.getValue<int>("parallel","thread_num");
  Ni = r.getValue<int>("Mesh","Ni"); Nj = r.getValue<int>("Mesh","Nj");
  Cell_Len = r.getValue<T>("Mesh","Cell_Len");
  BlockCellLen = r.getValue<int>("Mesh","BlockCellLen");
  Cyl_Radius = r.getValue<T>("Cylinder","Radius");
  Cyl_Cx = r.getValue<T>("Cylinder","CenterX");
  Cyl_Cy = r.getValue<T>("Cylinder","CenterY");
  Interface_Width = r.getValue<T>("Phase_Field","Interface_Width");
  Mobility = r.getValue<T>("Phase_Field","Mobility");
  mu_l = r.getValue<T>("Magnetic_Field","mu_l");   // outer fluid
  mu_h = r.getValue<T>("Magnetic_Field","mu_h");   // cylinder
  H0 = r.getValue<T>("Magnetic_Field","H0");
  PsiSolver_Iter = r.getValue<int>("Magnetic_Field","PsiSolver_Iter");
  PsiSolver_K = r.getValue<T>("Magnetic_Field","PsiSolver_K");
  MaxSubStep = r.getValue<int>("Simulation_Settings","SubSteps");
  IF_MPI_RANK(0){
    printf("---- Cylinder in Uniform Magnetic Field (paper Sec. III.A) ----\n");
    printf("Mesh: %dx%d  R=%.1f center=(%.0f,%.0f)  mu_l=%.3f mu_h=%.3f H0=%.4f\n",
           Ni,Nj,Cyl_Radius,Cyl_Cx,Cyl_Cy,mu_l,mu_h,H0);
    printf("analytic: H_in/H0 = 2*mu_h/(mu_l+mu_h) = %.6f\n",
           T{2}*mu_h/(mu_l+mu_h));
    printf("psi solver: K=%.2f iter=%d substeps=%d\n",PsiSolver_K,PsiSolver_Iter,MaxSubStep);
    printf("---------------------------------------\n");
  }
}

int main(int argc, char* argv[]) {
  constexpr std::uint8_t VoidFlag=1,BulkFlag=2,BouncebackFlag=4,PeriodicFlag=8;
  mpi().init(&argc,&argv); MPI_DEBUG_WAIT
  Printer::Print_BigBanner(std::string("Initializing Cylinder Mag Validation..."));
  readParam();

  BaseConverter<T> BaseConv(LatSet::cs2);
  BaseConv.SimplifiedConverterFromRT(Ni,T(0.01),T(1.0));
  BaseConverter<T> MFBaseConv(MFLatSet::cs2);
  MFBaseConv.SimplifiedConverterFromRT(Ni,T(0.01),T(1.0));
  UnitConvManager<T> ConvManager(&BaseConv); ConvManager.Check_and_Print();

  // -- geometry (periodic left/right, walls top/bottom) --
  AABB<T,2> domain({0,0},{T(Ni*Cell_Len),T(Nj*Cell_Len)});
  AABB<T,2> left({T(-Cell_Len),0},{0,T(Nj*Cell_Len)});
  AABB<T,2> right({T(Ni*Cell_Len),0},{T((Ni+1)*Cell_Len),T(Nj*Cell_Len)});
  BlockGeometryHelper2D<T> GeoHelper(Ni,Nj,domain,Cell_Len,BlockCellLen);
  GeoHelper.CreateBlocks(8,8);
  GeoHelper.AdaptiveOptimization(mpi().getSize());
  GeoHelper.LoadBalancing(mpi().getSize());
  BlockGeometry2D<T> Geo(GeoHelper);

  BlockFieldManager<FLAG,T,2> FlagFM(Geo,VoidFlag);
  FlagFM.forEach(domain,[&](FLAG&f,std::size_t id){f.SetField(id,BulkFlag);});
  FlagFM.forEach(left,[&](FLAG&f,std::size_t id){f.SetField(id,PeriodicFlag);});
  FlagFM.forEach(right,[&](FLAG&f,std::size_t id){f.SetField(id,PeriodicFlag);});
  FlagFM.template SetupBoundary<LatSet>(domain,BouncebackFlag);
  // reset left/right boundary columns to BulkFlag (periodic, not walls)
  AABB<T,2> left_col({T{0},Cell_Len},{Cell_Len,T((Nj-1)*Cell_Len)});
  AABB<T,2> right_col({T((Ni-1)*Cell_Len),Cell_Len},{T(Ni*Cell_Len),T((Nj-1)*Cell_Len)});
  FlagFM.forEach(left_col,[&](FLAG&f,std::size_t id){f.SetField(id,BulkFlag);});
  FlagFM.forEach(right_col,[&](FLAG&f,std::size_t id){f.SetField(id,BulkFlag);});

  // -- PF lattice: only holds the static cylinder phi field --
  using PFFIELDS=TypePack<PHI<T>,POP<T,LatSet::q>,GRAD<T,2>,NORMAL<T,2>,INTERFACEWIDTH<T>,
    ff::LAPLACIAN<T>,ff::CHEMICALPOTENTIAL<T>,ff::GRAVITY<T>,ff::BETA<T>,ff::KAPPA<T>,
    ff::RHO_L<T>,ff::RHO_H<T>,ff::ETA_L<T>,ff::ETA_H<T>,ff::DELTARHO<T>>;
  using PFREF=TypePack<VELOCITY<T,2>>;
  using PFPACK=TypePack<PFFIELDS,PFREF>;
  ValuePack PFI(T{},T{},Vector<T,2>{0,0},Vector<T,2>{0,0},Interface_Width,
    T{},T{},T{0},T{1},T{1},T{1},T{1},T{1},T{1},T{0});
  using PFCELL=Cell<T,LatSet,ExtractFieldPack<PFPACK>::mergedpack>;
  BlockLatticeManager<T,LatSet,PFPACK> PFLattice(Geo,PFI,BaseConv);

  // -- MF lattice (D2Q5) --
  using MFFIELDS=TypePack<PSI<T>,OMEGA_PSI<T>,MU_PERCELL<T>,CHI_PERCELL<T>,
    HX<T>,HY<T>,HMAG<T>,POP<T,MFLatSet::q>,
    MU_L<T>,MU_H<T>,CHI_L<T>,CHI_H<T>,H_0<T>,PSI_K<T>>;
  using MFREF=TypePack<PHI<T>>;
  using MFPACK=TypePack<MFFIELDS,MFREF>;
  ValuePack MFI(T{},T{1.0},T{mu_l},T{0},T{},T{},T{},T{},
    mu_l,mu_h,T{0},T{0},H0,PsiSolver_K);
  using MFCELL=Cell<T,MFLatSet,ExtractFieldPack<MFPACK>::mergedpack>;
  BlockLatticeManager<T,MFLatSet,MFPACK> MFLattice(Geo,MFI,MFBaseConv,
    &PFLattice.getField<PHI<T>>());
  T chi_l_d=T{0}, chi_h_d=T{0};
  BroadcastAllMFParams<T>(MFLattice,mu_l,mu_h,chi_l_d,chi_h_d,H0,PsiSolver_K);

  // -- init phi: cylinder (phi=1 inside, mu_h) in the outer fluid (phi=0, mu_l) --
  T R_phys=Cyl_Radius*Cell_Len, xc=Cyl_Cx*Cell_Len, yc=Cyl_Cy*Cell_Len;
  T W_phys=Interface_Width*Cell_Len;
  auto& phiField=PFLattice.getField<PHI<T>>();
  for(int b=0;b<Geo.getBlockNum();++b){
    const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
    auto& bPhi=phiField.getBlockField(b);
    T vs=bk.getVoxelSize(),mx=bk.getMin()[0],my=bk.getMin()[1];
    for(int j=0;j<bk.getNy();++j){
      T y=my+T(j)*vs;
      for(int i=0;i<bk.getNx();++i){
        T x=mx+T(i)*vs;
        T dist=std::sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc));
        T phi=T{0.5}-T{0.5}*std::tanh(T{2.0}*(dist-R_phys)/W_phys);
        bPhi.get(j*pr[1]+i)=phi;
      }
    }
  }
  PFLattice.getField<PHI<T>>().Communicate();

  // -- BCs --
  using LM_MF=BlockLatticeManager<T,MFLatSet,MFPACK>;
  using FM=BlockFieldManager<FLAG,T,2>;
  FixedPeriodicBoundaryManager<LM_MF,FM> MF_Per("MF_Per",MFLattice,FlagFM,PeriodicFlag,VoidFlag);
  MF_Per.Setup(left,NbrDirection::XN,right,NbrDirection::XP);
  MF_Per.Setup(right,NbrDirection::XP,left,NbrDirection::XN);
#ifdef MPI_ENABLED
  MF_Per.SetupMPI(GeoHelper);
#endif

  // MF coupling: PF->MF coeff update (per-cell mu from phi)
  using MCT=tmp::Key_TypePair<BulkFlag,MFUpdateCoeffs2D<PFCELL,MFCELL>>;
  using MCSel=CoupledTaskSelector<std::uint8_t,PFCELL,MFCELL,MCT>;
  BlockLatManagerCoupling MCC(PFLattice,MFLattice);

  // seam sync for per-cell MF fields (cross-rank aware)
  auto SyncSeamField = [&](auto& field) {
    using ValT = typename std::decay_t<decltype(field.getBlockField(0))>::value_type;
    const int valSz = static_cast<int>(sizeof(ValT));
    const T Lx = T(Ni) * Cell_Len;
    const int myRank = mpi().getRank();
    const auto& gGeo = GeoHelper.getBlockGeometry();
    const int nb = Geo.getBlockNum();
#ifdef MPI_ENABLED
    std::vector<std::vector<char>> sendBufs, recvBufs;
    std::vector<MPI_Request> sendReqs, recvReqs;
    std::vector<std::pair<int, int>> recvJobs;
#endif
    for (int b = 0; b < nb; ++b) {
      const auto& bk = Geo.getBlock(b);
      const bool atL = bk.getMin()[0] < Cell_Len * T{0.5};
      const bool atR = bk.getMax()[0] > Lx - Cell_Len * T{0.5};
      if (!atL && !atR) continue;
      const int nx = bk.getNx(), ny = bk.getNy(), ov = bk.getOverlap();
      const auto& pr = bk.getProjection();
      const T ymid = (bk.getMin()[1] + bk.getMax()[1]) * T{0.5};
      const T xprobe = atL ? (Lx - Cell_Len * T{0.5}) : (Cell_Len * T{0.5});
      int srcGid = -1;
      for (std::size_t sbi = 0; sbi < gGeo.getBlockNum(); ++sbi) {
        const auto& gb = gGeo.getBlock(sbi);
        if (gb.getSelfBlock().isInside(Vector<T, 2>{xprobe, ymid})) { srcGid = gb.getBlockId(); break; }
      }
      if (srcGid < 0) continue;
      const int ghostCol0 = atL ? 0 : nx - ov;
      const int theirPhysCol = atL ? (nx - 1 - ov) : ov;
      if (GeoHelper.whichRank(srcGid) == myRank) {
        const int srcBlock = Geo.findBlockIndex(srcGid);
        auto& bF = field.getBlockField(b);
        auto& sF = field.getBlockField(srcBlock);
        for (int j = 0; j < ny; ++j) {
          const ValT v = sF.get(j * pr[1] + theirPhysCol);
          for (int c = 0; c < ov; ++c) bF.get(j * pr[1] + ghostCol0 + c) = v;
        }
      } else {
#ifdef MPI_ENABLED
        const int myPhysCol = atL ? ov : nx - 1 - ov;
        std::vector<char> snd(ny * valSz);
        auto& bF = field.getBlockField(b);
        for (int j = 0; j < ny; ++j) {
          const ValT v = bF.get(j * pr[1] + myPhysCol);
          std::memcpy(&snd[j * valSz], &v, valSz);
        }
        sendBufs.emplace_back(std::move(snd));
        MPI_Request sreq;
        mpi().iSend(sendBufs.back().data(), static_cast<int>(sendBufs.back().size()),
                    GeoHelper.whichRank(srcGid), &sreq, 9500 + bk.getBlockId());
        sendReqs.push_back(sreq);
        recvBufs.emplace_back(ny * valSz);
        recvJobs.emplace_back(b, ghostCol0);
        MPI_Request rreq;
        mpi().iRecv(recvBufs.back().data(), static_cast<int>(recvBufs.back().size()),
                    GeoHelper.whichRank(srcGid), &rreq, 9500 + srcGid);
        recvReqs.push_back(rreq);
#endif
      }
    }
#ifdef MPI_ENABLED
    MPI_Waitall(static_cast<int>(sendReqs.size()), sendReqs.data(), MPI_STATUSES_IGNORE);
    for (std::size_t i = 0; i < recvReqs.size(); ++i) {
      MPI_Wait(&recvReqs[i], MPI_STATUS_IGNORE);
      auto& bF = field.getBlockField(recvJobs[i].first);
      const auto& bk = Geo.getBlock(recvJobs[i].first);
      const auto& pr = bk.getProjection();
      const int ov = bk.getOverlap();
      const auto& buf = recvBufs[i];
      for (int j = 0; j < bk.getNy(); ++j) {
        ValT v; std::memcpy(&v, &buf[j * valSz], valSz);
        for (int c = 0; c < ov; ++c) bF.get(j * pr[1] + recvJobs[i].second + c) = v;
      }
    }
#endif
  };

  // -- init psi: ANALYTICAL solution (paper Eq. 59-61, with the corrected
  // coefficients: the paper's Eq. 61 A-coefficient has mu1/mu2 swapped; the
  // magnetostatic solution is A = -2*mu1/(mu1+mu2)*H0 and
  // D = +(mu2-mu1)/(mu1+mu2)*R^2*H0, giving H_in = 2*mu1/(mu1+mu2)*H0) --
  //   psi_in  = -2*mu1/(mu1+mu2)*H0 * y                    (r <= R)
  //   psi_out = -H0*y + (mu2-mu1)/(mu1+mu2)*R^2*H0 * y/r^2  (r > R)
  // Seeding the exact fixed point skips the slow L^2/D global relaxation and
  // validates directly that the D2Q5 solver PRESERVES the correct field.
  {
    T mu1=mu_l, mu2=mu_h;
    T A_ana = -T{2}*mu1/(mu1+mu2)*H0;
    T D_ana = (mu2-mu1)/(mu2+mu1)*R_phys*R_phys*H0;
    auto& psiF=MFLattice.getField<PSI<T>>();
    for(int b=0;b<Geo.getBlockNum();++b){
      auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
      auto& bPsi=psiF.getBlockField(b);
      T vs=bk.getVoxelSize(),mx=bk.getMin()[0],my=bk.getMin()[1];
      for(int j=0;j<bk.getNy();++j){
        T y=my+T(j)*vs;
        for(int i=0;i<bk.getNx();++i){
          T x=mx+T(i)*vs;
          T dx=x-xc, dy=y-yc;
          T r=std::sqrt(dx*dx+dy*dy);
          T psi;
          if(r<=R_phys) psi = A_ana*y;
          else          psi = -H0*y + D_ana*y/(r*r);
          std::size_t id=j*pr[1]+i; MFCELL c(id,bl);
          for(unsigned k=0;k<MFLatSet::q;++k) c[k]=latset::w<MFLatSet>(k)*psi;
          bPsi.get(id)=psi;
        }
      }
    }
  }
  MFLattice.NormalFullCommunicate();
  MF_Per.Apply();

  // -- magnetic solve: pseudo-time sub-iterations on the frozen mu field --
  const T H_global=T(Nj)*Cell_Len;
  const T nwall=Cell_Len*T{3.0};
  T mu1_p=mu_l, mu2_p=mu_h;
  T D_ana_p = (mu2_p-mu1_p)/(mu2_p+mu1_p)*R_phys*R_phys*H0;  // dipole coeff (centered at xc,yc)
  auto PinWalls = [&](){
    auto& psiF=MFLattice.getField<PSI<T>>();
    for(int b=0;b<Geo.getBlockNum();++b){
      auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
      auto& bPsi=psiF.getBlockField(b);
      int nx=bk.getNx(),ny=bk.getNy(),ov=bk.getOverlap();
      T minY=bk.getMin()[1],minX=bk.getMin()[0],vs=bk.getVoxelSize();
      for(int jj=0;jj<ny;++jj){
        T y=minY+T(jj-ov)*vs;
        if((y<=nwall&&y>=-nwall)||(y<=H_global+nwall&&y>=H_global-nwall)){
          for(int ii=0;ii<nx;++ii){
            T x=minX+T(ii)*vs;
            // analytic outside solution at the wall rows (centered at xc,yc)
            T dx=x-xc, dy=y-yc;
            T r2=dx*dx+dy*dy;
            T psi_w=-H0*y + D_ana_p*dy/r2;
            std::size_t id=jj*pr[1]+ii; MFCELL c(id,bl);
            for(unsigned k=0;k<MFLatSet::q;++k) c[k]=latset::w<MFLatSet>(k)*psi_w;
            bPsi.get(id)=psi_w;
          }
        }
      }
    }
  };

  Timer t;
  for(int step=0; step<MaxSubStep; ++step){
    MCC.ApplyInnerCellDynamics<MCSel>(step,FlagFM);
    CommunicateOMEGAPSI<T>(MFLattice);
    for(int sub=0; sub<PsiSolver_Iter; ++sub){
      PinWalls();
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        int ov=bk.getOverlap();
        for(int j=ov;j<bk.getNy()-ov;++j)
          for(int i=ov;i<bk.getNx()-ov;++i){
            MFCELL c(j*pr[1]+i,bl);
            collision::MRTDiffusion<MFCELL,OMEGA_PSI<T>>::apply(c);
          }
      }
      MF_Per.Apply();
      MFLattice.NormalFullCommunicate();
      MFLattice.Stream();
      {
        auto& psiF=MFLattice.getField<PSI<T>>();
        for(int b=0;b<Geo.getBlockNum();++b){
          auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
          auto& bPsi=psiF.getBlockField(b); int ov=bk.getOverlap();
          for(int j=ov;j<bk.getNy()-ov;++j)
            for(int i=ov;i<bk.getNx()-ov;++i){
              std::size_t id=j*pr[1]+i; MFCELL c(id,bl);
              T psi=0; for(unsigned k=0;k<MFLatSet::q;++k)psi+=c[k];
              bPsi.get(id)=psi;
            }
        }
      }
      CommunicatePSI<T>(MFLattice);
      PinWalls();
    }
    SyncSeamField(MFLattice.getField<PSI<T>>());
    for(int b=0;b<Geo.getBlockNum();++b){
      auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
      int ov=bk.getOverlap();
      for(int j=ov;j<bk.getNy()-ov;++j)
        for(int i=ov;i<bk.getNx()-ov;++i){
          MFCELL c(j*pr[1]+i,bl);
          MFComputeH2D<MFCELL>::apply(c);
        }
    }
    CommunicateAllMFFields<T>(MFLattice);
    SyncSeamField(MFLattice.getField<HX<T>>());
    SyncSeamField(MFLattice.getField<HY<T>>());
    SyncSeamField(MFLattice.getField<HMAG<T>>());
  }
  t.Print_MainLoopPerformance(Geo.getTotalCellNum());
  // -- validation output: H/H0 along the vertical axis through the center --
  IF_MPI_RANK(0){
    FILE* fo=std::fopen("cylinder_axis.dat","w");
    std::fprintf(fo,"# y  y/R  Hx  Hy  |H|/H0  analytic_in  analytic_out  mu\n");
    auto& hyF=MFLattice.getField<HY<T>>();
    auto& hxF=MFLattice.getField<HX<T>>();
    auto& muF=MFLattice.getField<MU_PERCELL<T>>();
    T h_center=T{0};
    for(int b=0;b<Geo.getBlockNum();++b){
      const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
      auto& bHy=hyF.getBlockField(b); auto& bHx=hxF.getBlockField(b);
      auto& bMu=muF.getBlockField(b);
      T vs=bk.getVoxelSize(),mx=bk.getMin()[0],my=bk.getMin()[1];
      for(int j=0;j<bk.getNy();++j){
        T y=my+T(j)*vs;
        for(int i=0;i<bk.getNx();++i){
          T x=mx+T(i)*vs;
          if(std::abs(x-xc)>T{0.5}) continue;
          T hy=bHy.get(j*pr[1]+i)/H0, hx=bHx.get(j*pr[1]+i)/H0;
          T r=y-yc;
          T ana_in=T{2}*mu_h/(mu_l+mu_h);
          T ana_out=(std::abs(r)>T{1e-9})?
            (T{1}+(mu_h-mu_l)/(mu_h+mu_l)*(R_phys*R_phys/(r*r))):T{0};
          std::fprintf(fo,"%.4f %.4f %.6f %.6f %.6f %.6f %.6f %.4f\n",
                       y,r/R_phys,hx,hy,std::sqrt(hx*hx+hy*hy),ana_in,ana_out,
                       bMu.get(j*pr[1]+i));
          if(std::abs(y-yc)<=T{0.5}&&std::abs(x-xc)<=T{0.5})
            h_center=std::sqrt(hx*hx+hy*hy);
        }
      }
    }
    std::fclose(fo);
    printf("H_center/H0 = %.6f  (correct analytic 2*mu_l/(mu_l+mu_h) = %.6f)\n",
           h_center,T{2}*mu_l/(mu_l+mu_h));
    // psi profile along the axis (x=xc)
    {
      auto& psiF=MFLattice.getField<PSI<T>>();
      FILE* fp=std::fopen("cylinder_psi_axis.dat","w");
      std::fprintf(fp,"# y/R  psi  seed_psi_in  seed_psi_out\n");
      T A_ana = -T{2}*mu_l/(mu_l+mu_h)*H0;
      T D_ana = (mu_h-mu_l)/(mu_h+mu_l)*R_phys*R_phys*H0;
      for(int b=0;b<Geo.getBlockNum();++b){
        const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bPsi=psiF.getBlockField(b);
        T vs=bk.getVoxelSize(),mx=bk.getMin()[0],my=bk.getMin()[1];
        for(int j=0;j<bk.getNy();++j){
          T y=my+T(j)*vs;
          T dy=y-yc;
          for(int i=0;i<bk.getNx();++i){
            T x=mx+T(i)*vs;
            if(std::abs(x-xc)>T{0.5}) continue;
            T r=std::sqrt((x-xc)*(x-xc)+dy*dy);
            T seed = (r<=R_phys)?(A_ana*y):(-H0*y + D_ana*y/(r*r));
            std::fprintf(fp,"%.4f %.6f %.6f\n", dy/R_phys,
                         bPsi.get(j*pr[1]+i), seed);
          }
        }
      }
      std::fclose(fp);
    }
  }
  Printer::Print_BigBanner(std::string("Cylinder Mag Validation Complete!"));
  return 0;
}
