// rosensweig2d.cpp — 2D Rosensweig instability (Phase field + NS + Magnetic)
#include "freelb.h"
#include "freelb.hh"
#include "ff/ff2d.h"

using T = FLOAT;
using LatSet = D2Q9<T>;
using MFLatSet = D2Q5<T>;
using namespace mfield;

// ---- Simulation Parameters ----
int Ni, Nj;
T Cell_Len;
int BlockCellLen, Thread_Num;

// phase field
T Interface_Width, Mobility, Tau_phi, Omega_phi, Kappa, Beta;

// two-phase
T rho_l, rho_h, eta_l, eta_h, sigma, gravity, Tau_ns, DeltaRho;

// interface geometry
T Interface_Y, Perturb_Amplitude, Perturb_Wavelength;

// magnetic field
T chi_l, chi_h, mu_l, mu_h, H0, H0_phys_kAm;

// psi-solver: boosted relaxation tau = 0.5 + PsiSolver_K*mu (K=0.5 -> smooth
// fixed point, ~100 sub-iterations; K=3.0 = 1/cs^2 reproduces original solve)
int PsiSolver_Iter; T PsiSolver_K;

int MaxStep, OutputStep;
std::string work_dir;

void readParam() {
  iniReader r("rosensweig2d.ini");
  work_dir = r.getValue<std::string>("workdir","workdir_");
  Thread_Num = r.getValue<int>("parallel","thread_num");
  Ni = r.getValue<int>("Mesh","Ni"); Nj = r.getValue<int>("Mesh","Nj");
  Cell_Len = r.getValue<T>("Mesh","Cell_Len");
  BlockCellLen = r.getValue<int>("Mesh","BlockCellLen");
  Interface_Width=r.getValue<T>("Phase_Field","Interface_Width");
  Mobility=r.getValue<T>("Phase_Field","Mobility");
  rho_l=r.getValue<T>("Two_Phase","rho_l");
  rho_h=r.getValue<T>("Two_Phase","rho_h");
  eta_l=r.getValue<T>("Two_Phase","eta_l");
  eta_h=r.getValue<T>("Two_Phase","eta_h");
  sigma=r.getValue<T>("Two_Phase","sigma");
  gravity=r.getValue<T>("Two_Phase","gravity");
  Interface_Y=r.getValue<T>("Interface","Interface_Y");
  Perturb_Amplitude=r.getValue<T>("Interface","Perturb_Amplitude");
  Perturb_Wavelength=r.getValue<T>("Interface","Perturb_Wavelength");
  MaxStep=r.getValue<int>("Simulation_Settings","TotalStep");
  OutputStep=r.getValue<int>("Simulation_Settings","OutputStep");
  // Magnetic
  chi_l=r.getValue<T>("Magnetic_Field","chi_l");
  chi_h=r.getValue<T>("Magnetic_Field","chi_h");
  mu_l=r.getValue<T>("Magnetic_Field","mu_l");
  mu_h=r.getValue<T>("Magnetic_Field","mu_h");
  H0_phys_kAm=r.getValue<T>("Magnetic_Field","H0_phys_kAm");
  PsiSolver_Iter=r.getValue<int>("Magnetic_Field","PsiSolver_Iter");
  PsiSolver_K=r.getValue<T>("Magnetic_Field","PsiSolver_K");

  // Derived parameters
  DeltaRho = rho_h - rho_l;
  Beta = T(12.0) * sigma / Interface_Width;
  Kappa = T(3.0) * Interface_Width * sigma * T(0.5);
  Tau_phi = T(3.0) * Mobility + T(0.5); Omega_phi = T(1.0) / Tau_phi;
  Tau_ns = T(0.5) + eta_h / rho_h / LatSet::cs2;

  // H0 conversion: physical (kA/m) -> lattice via Cowley-Rosensweig critical field
  // H_c_lattice = sqrt(2 * (1/mu_r+1)/((1/mu_r-1)^2) * sqrt(sigma*g*DeltaRho))
  // H_c_phys = 4.7 kA/m
  // H0_lattice = H0_phys * (H_c_lattice / 4.7)
  T mu_r = mu_h / mu_l;  // relative permeability
  T inv_mu_r = T(1.0) / mu_r;
  T cr_coeff = (inv_mu_r + T(1.0)) / ((inv_mu_r - T(1.0)) * (inv_mu_r - T(1.0)));
  T sg_drho = std::sqrt(sigma * gravity * DeltaRho);
  T Hc_lattice = std::sqrt(T(2.0) * cr_coeff * sg_drho);
  T Hc_phys = T(4.7);  // kA/m
  H0 = H0_phys_kAm * (Hc_lattice / Hc_phys);

  MPI_RANK(0){
    printf("---- Rosensweig Instability (Ferrofluid) ----\n");
    printf("Mesh: %dx%d  BlockCellLen=%d\n",Ni,Nj,BlockCellLen);
    printf("rho: l=%.4f h=%.4f  eta: l=%.4f h=%.4f  sigma=%.5e g=%.5e\n",rho_l,rho_h,eta_l,eta_h,sigma,gravity);
    printf("Interface: Y=%.1f W=%.1f M=%.3f tau_phi=%.3f\n",Interface_Y,Interface_Width,Mobility,Tau_phi);
    printf("Magnetic: chi=(%.1f,%.1f) mu=(%.1f,%.1f) mu_r=%.2f\n",chi_l,chi_h,mu_l,mu_h,mu_r);
    printf("H0: phys=%.2f kA/m -> lattice=%.5f (Hc_lattice=%.5f, Hc_phys=4.7)\n",H0_phys_kAm,H0,Hc_lattice);
    printf("PsiSolver: K=%.3f iter=%d\n",PsiSolver_K,PsiSolver_Iter);
    printf("Perturb: A=%.2f lambda=%.2f\n",Perturb_Amplitude,Perturb_Wavelength);
    printf("----------------------------------------------\n");
  }
}

int main(int argc, char* argv[]) {
  constexpr std::uint8_t VoidFlag=1,BulkFlag=2,BouncebackFlag=4,PeriodicFlag=8;
  mpi().init(&argc,&argv); MPI_DEBUG_WAIT
  Printer::Print_BigBanner(std::string("Initializing Rosensweig Instability..."));
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
  AABB<T,2> domain({0,0},{T(Ni*Cell_Len),T(Nj*Cell_Len)});
  AABB<T,2> left({T(-Cell_Len),0},{0,T(Nj*Cell_Len)});
  AABB<T,2> right({T(Ni*Cell_Len),0},{T((Ni+1)*Cell_Len),T(Nj*Cell_Len)});
  BlockGeometryHelper2D<T> GeoHelper(Ni,Nj,domain,Cell_Len,BlockCellLen);
  GeoHelper.CreateBlocks(3,3);
  GeoHelper.AdaptiveOptimization(mpi().getSize());
  GeoHelper.LoadBalancing(mpi().getSize());
  BlockGeometry2D<T> Geo(GeoHelper);

  // -- flag --
  BlockFieldManager<FLAG,T,2> FlagFM(Geo,VoidFlag);
  FlagFM.forEach(domain,[&](FLAG&f,std::size_t id){f.SetField(id,BulkFlag);});
  FlagFM.forEach(left,[&](FLAG&f,std::size_t id){f.SetField(id,PeriodicFlag);});
  FlagFM.forEach(right,[&](FLAG&f,std::size_t id){f.SetField(id,PeriodicFlag);});
  FlagFM.template SetupBoundary<LatSet>(domain,BouncebackFlag);
  // 修复: SetupBoundary 会把 x 向邻居为周期 ghost 的接缝列 (x=0.5 / x=251.5)
  // 误标为 BouncebackFlag —— 其 x 邻居虽位于 domain AABB 之外, 却是周期 ghost
  // 而非实体边界。将这些接缝列恢复为 BulkFlag, 使 RoC/碰撞/力/宏量照常作用于
  // 接缝, 消除左右边界的"壁面爬升"伪影。上下墙行 (j=ov / j=ny-1-ov) 保持
  // BouncebackFlag (真实 no-slip 墙)。
  {
    // 仅跳过真实墙行: j=ov (y=0.5, 仅含 y=0 的块) 和 j=ny-1-ov (y=83.5,
    // 仅含 y=H_global 的块); 其它块同列行是块间内部边界, 不是墙
    const T Hglobal = T(Nj) * Cell_Len;
    for (int b = 0; b < Geo.getBlockNum(); ++b) {
      const auto& bk = Geo.getBlock(b); const auto& pr = bk.getProjection();
      auto& bF = FlagFM.getBlockField(b).getField(0);
      const int nx = bk.getNx(), ny = bk.getNy(), ov = bk.getOverlap();
      T my = bk.getMin()[1], My = bk.getMax()[1];
      const bool hasBottomWall = my < Cell_Len * T{1.5};
      const bool hasTopWall = My > Hglobal - Cell_Len * T{1.5};
      for (int j = ov; j < ny - ov; ++j) {
        if (hasBottomWall && j == ov) continue;
        if (hasTopWall && j == ny - 1 - ov) continue;
        for (int i = ov; i < nx - ov; ++i) {
          std::size_t id = j * pr[1] + i;
          if (bF[id] == BouncebackFlag &&
              (bF[j * pr[1] + i - 1] == PeriodicFlag ||
               bF[j * pr[1] + i + 1] == PeriodicFlag)) {
            bF.set(id, BulkFlag);
          }
        }
      }
    }
  }

  // -- NS lattice --
  using NSFIELDS=TypePack<DENSITY<T>,VELOCITY<T,2>,POP<T,LatSet::q>,FORCE<T,2>,OMEGA<T>,PRESSURE<T>>;
  T omega_ns=T{1}/Tau_ns;
  ValuePack NSI(T{1},Vector<T,2>{0,0},T{},Vector<T,2>{0,0},omega_ns,T{});
  using NSCELL=Cell<T,LatSet,NSFIELDS>;
  BlockLatticeManager<T,LatSet,NSFIELDS> NSLattice(Geo,NSI,BaseConv);

  // -- PF lattice --
  using PFFIELDS=TypePack<PHI<T>,POP<T,LatSet::q>,GRAD<T,2>,NORMAL<T,2>,INTERFACEWIDTH<T>,
    ff::LAPLACIAN<T>,ff::CHEMICALPOTENTIAL<T>,
    ff::GRAVITY<T>,ff::BETA<T>,ff::KAPPA<T>,
    ff::RHO_L<T>,ff::RHO_H<T>,ff::ETA_L<T>,ff::ETA_H<T>,ff::DELTARHO<T>>;
  using PFREF=TypePack<VELOCITY<T,2>>;
  using PFPACK=TypePack<PFFIELDS,PFREF>;
  ValuePack PFI(T{},T{},Vector<T,2>{0,0},Vector<T,2>{0,0},Interface_Width,
    T{},T{},gravity,Beta,Kappa,rho_l,rho_h,eta_l,eta_h,DeltaRho);
  using PFCELL=Cell<T,LatSet,ExtractFieldPack<PFPACK>::mergedpack>;
  BlockLatticeManager<T,LatSet,PFPACK> PFLattice(Geo,PFI,PFBaseConv,
    &NSLattice.getField<VELOCITY<T,2>>());

  ff::BroadcastAllParams<T>(PFLattice,rho_l,rho_h,eta_l,eta_h,gravity,Beta,Kappa);
  PFLattice.template getField<ff::DELTARHO<T>>().InitValue(DeltaRho);

  // -- MF lattice (D2Q5) --
  using MFFIELDS=TypePack<PSI<T>,OMEGA_PSI<T>,MU_PERCELL<T>,CHI_PERCELL<T>,
    HX<T>,HY<T>,HMAG<T>,POP<T,MFLatSet::q>,
    MU_L<T>,MU_H<T>,CHI_L<T>,CHI_H<T>,H_0<T>,PSI_K<T>>;
  using MFREF=TypePack<PHI<T>>;
  using MFPACK=TypePack<MFFIELDS,MFREF>;
  ValuePack MFI(T{},T{1.0},T{mu_l},T{chi_l},T{},T{},T{},T{},
    mu_l,mu_h,chi_l,chi_h,H0,PsiSolver_K);
  using MFCELL=Cell<T,MFLatSet,ExtractFieldPack<MFPACK>::mergedpack>;
  BlockLatticeManager<T,MFLatSet,MFPACK> MFLattice(Geo,MFI,MFBaseConv,
    &PFLattice.getField<PHI<T>>());
  BroadcastAllMFParams<T>(MFLattice,mu_l,mu_h,chi_l,chi_h,H0,PsiSolver_K);
  MFLattice.getField<OMEGA_PSI<T>>().InitValue(T{1.0});

  // -- init phi (horizontal layer with cosine perturbation) --
  // ferrofluid (phi=1) below interface, organic solvent (phi=0) above
  // interface at y = Interface_Y + Perturb_Amplitude * cos(2*pi*x/lambda)
  T W_phys = Interface_Width * Cell_Len;
  auto& phiField = PFLattice.getField<PHI<T>>();
  for(int b=0;b<Geo.getBlockNum();++b){
    const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
    auto& bPhi=phiField.getBlockField(b); auto& bl=PFLattice.getBlockLat(b);
    T vs=bk.getVoxelSize(),mx=bk.getMin()[0],my=bk.getMin()[1];
    int ov=0;
    for(int j=ov;j<bk.getNy()-ov;++j){
      T y=my+T(j)*vs;
      for(int i=ov;i<bk.getNx()-ov;++i){
        T x=mx+T(i)*vs;
        // perturbed interface position
        T y_int = Interface_Y * Cell_Len
                + Perturb_Amplitude * Cell_Len
                * std::cos(T(2.0)*M_PI*x/(Perturb_Wavelength*Cell_Len));
        // phi=1 below interface (ferrofluid), phi=0 above (organic solvent)
        // phi = 0.5 - 0.5*tanh(2*(y-y_int)/W)  -> phi=1 when y<<y_int, phi=0 when y>>y_int
        T phi = T{0.5} - T{0.5}*std::tanh(T{2.0}*(y-y_int)/W_phys);
        bPhi.get(j*pr[1]+i)=phi;
      }
    }
  }
  // init PF pops
  for(int b=0;b<Geo.getBlockNum();++b){
    auto& bl=PFLattice.getBlockLat(b); auto& bPhi=phiField.getBlockField(b);
    const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
    int ov=0;
    for(int j=ov;j<bk.getNy()-ov;++j)
      for(int i=ov;i<bk.getNx()-ov;++i){
        std::size_t id=j*pr[1]+i; PFCELL c(id,bl); T phi=bPhi.get(id);
        for(unsigned k=0;k<LatSet::q;++k) c[k]=latset::w<LatSet>(k)*phi*(T{1}+LatSet::InvCs2*T{0});
      }
  }
  PFLattice.getField<INTERFACEWIDTH<T>>().InitValue(Interface_Width);

  // init NS pops (p=0, u=0)
  Vector<T,2> uz{0,0}; T pz=0;
  for(int b=0;b<Geo.getBlockNum();++b){
    auto& bl=NSLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
    int ov=0;
    for(int j=ov;j<bk.getNy()-ov;++j)
      for(int i=ov;i<bk.getNx()-ov;++i){
        std::size_t id=j*pr[1]+i; NSCELL c(id,bl);
        for(unsigned k=0;k<LatSet::q;++k){
          T uc=uz*latset::c<LatSet>(k);
          c[k]=latset::w<LatSet>(k)*(pz+LatSet::InvCs2*uc+uc*uc*T{0.5}*LatSet::InvCs4-LatSet::InvCs2*T{0});
        }
      }
  }

  // init MF: psi = -H0*y (uniform field), g_k = w_k*psi
  {
    auto& psiF=MFLattice.getField<PSI<T>>();
    for(int b=0;b<Geo.getBlockNum();++b){
      auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
      auto& bPsi=psiF.getBlockField(b);
      T vs=bk.getVoxelSize(),my=bk.getMin()[1];
      for(int j=0;j<bk.getNy();++j){
        T y=my+T(j)*vs, psi=-H0*y;
        for(int i=0;i<bk.getNx();++i){
          std::size_t id=j*pr[1]+i; MFCELL c(id,bl);
          for(unsigned k=0;k<MFLatSet::q;++k) c[k]=latset::w<MFLatSet>(k)*psi;
          bPsi.get(id)=psi;
        }
      }
    }
  }

  // -- BCs --
  using LM_NS=BlockLatticeManager<T,LatSet,NSFIELDS>;
  using LM_PF=BlockLatticeManager<T,LatSet,PFPACK>;
  using LM_MF=BlockLatticeManager<T,MFLatSet,MFPACK>;
  using FM=BlockFieldManager<FLAG,T,2>;

  BBLikeFixedBlockBdManager<bounceback::normal<NSCELL>,LM_NS,FM>
    NS_BB("NS_BB",NSLattice,FlagFM,BouncebackFlag,VoidFlag);
  BBLikeFixedBlockBdManager<bounceback::normal<PFCELL>,LM_PF,FM>
    PF_BB("PF_BB",PFLattice,FlagFM,BouncebackFlag,VoidFlag);

  FixedPeriodicBoundaryManager<LM_NS,FM> NS_Per("NS_Per",NSLattice,FlagFM,PeriodicFlag,VoidFlag);
  NS_Per.Setup(left,NbrDirection::XN,right,NbrDirection::XP);
  NS_Per.Setup(right,NbrDirection::XP,left,NbrDirection::XN);
  FixedPeriodicBoundaryManager<LM_PF,FM> PF_Per("PF_Per",PFLattice,FlagFM,PeriodicFlag,VoidFlag);
  PF_Per.Setup(left,NbrDirection::XN,right,NbrDirection::XP);
  PF_Per.Setup(right,NbrDirection::XP,left,NbrDirection::XN);
  FixedPeriodicBoundaryManager<LM_MF,FM> MF_Per("MF_Per",MFLattice,FlagFM,PeriodicFlag,VoidFlag);
  MF_Per.Setup(left,NbrDirection::XN,right,NbrDirection::XP);
  MF_Per.Setup(right,NbrDirection::XP,left,NbrDirection::XN);
#ifdef MPI_ENABLED
  NS_Per.SetupMPI(GeoHelper); PF_Per.SetupMPI(GeoHelper); MF_Per.SetupMPI(GeoHelper);
#endif

  // -- MF wrapped-ghost field sync plans (MPI-aware) --
  // The x-wrapped periodic ghosts (x=-0.5 / x=252.5) of the per-cell MF fields
  // (PSI/HX/HY/HMAG) must mirror the opposite edge's physical seam column.
  // MF_Per.Apply() copies only pops + GenericRho, and GenericRho for MFCELL
  // resolves to PHI (see FindGenericRhoType in src/utils/tmp.h), so the wrapped
  // ghost columns of these fields are never exchanged by the framework. The
  // former per-call sync scanned only local blocks and became a no-op on
  // multi-rank runs (edge blocks live on different ranks), leaving the ghosts
  // frozen at their T0 values while the physical seam column drifts
  // (demagnetization) — MFComputeH2D / MFMagneticForce2D then read stale
  // neighbors at the seam → spurious HX → the wall climbing seen on the
  // 128-rank cloud run. Build the partner plans once: same-rank partners are
  // copied directly, cross-rank partners are exchanged with non-blocking sends
  // keyed by the sender's global block id (same convention as
  // FixedPeriodicBoundaryManager).
  struct MFGhostSyncPlan {
    std::size_t ghostBlockIdx; // local edge block holding the wrapped ghost
    std::size_t srcBlockIdx;   // local partner (same rank)
    int ghostCol;              // wrapped-ghost column (0 or nx-1)
    int seamCol;               // MY physical seam column sent to the partner (1 or nx-2)
    int srcCol;                // partner's opposite-edge column copied into ghostCol (nx-2 or 1)
    int ny;
    bool crossRank;
    int partnerRank, sendTag, recvTag; // crossRank only
  };
  std::vector<MFGhostSyncPlan> MFGhostSyncPlans;
  {
    const T xLeft = T{0};
    const T xRight = T(Ni) * Cell_Len;
    auto edgeOf = [&](const auto& bk, bool& atL, bool& atR) {
      atL = bk.getMin()[0] < xLeft + Cell_Len * T{0.5};
      atR = bk.getMax()[0] > xRight - Cell_Len * T{0.5};
    };
    for (std::size_t b = 0; b < Geo.getBlockNum(); ++b) {
      const auto& bk = Geo.getBlock(b);
      bool atL, atR; edgeOf(bk, atL, atR);
      if (!atL && !atR) continue;
      MFGhostSyncPlan pl;
      pl.ghostBlockIdx = b;
      pl.ghostCol = atL ? 0 : bk.getNx() - 1;
      pl.seamCol = atL ? 1 : bk.getNx() - 2;   // my column that wraps to the partner
      pl.srcCol = atL ? bk.getNx() - 2 : 1;    // partner's column that wraps to me
      pl.ny = bk.getNy();
      pl.crossRank = false;
      // same-rank partner first (serial and same-rank edge pairs)
      std::size_t partner = Geo.getBlockNum();
      for (std::size_t bp = 0; bp < Geo.getBlockNum(); ++bp) {
        const auto& bk2 = Geo.getBlock(bp);
        if (bk2.getMin()[1] != bk.getMin()[1]) continue;
        bool pL, pR; edgeOf(bk2, pL, pR);
        if (!((atL && pR) || (atR && pL) || (atL && atR && bp == b))) continue;
        partner = bp; break;
      }
      if (partner < Geo.getBlockNum()) {
        pl.srcBlockIdx = partner;
        MFGhostSyncPlans.push_back(pl);
      } else {
#ifdef MPI_ENABLED
        // remote partner: locate it in the global geometry (same y-row,
        // opposite edge) and exchange by rank
        const auto& globalGeo = GeoHelper.getBlockGeometry();
        int partnerId = -1;
        for (std::size_t g = 0; g < globalGeo.getBlockNum(); ++g) {
          const auto& gb = globalGeo.getBlock(g);
          if (gb.getMin()[1] != bk.getMin()[1]) continue;
          bool gL, gR; edgeOf(gb, gL, gR);
          if (!((atL && gR) || (atR && gL))) continue;
          partnerId = gb.getBlockId(); break;
        }
        if (partnerId < 0) continue; // geometry invariant broken — skip
        pl.crossRank = true;
        pl.partnerRank = GeoHelper.whichRank(partnerId);
        pl.sendTag = bk.getBlockId(); // my sends carry my block id
        pl.recvTag = partnerId;       // I expect the partner's block id
        MFGhostSyncPlans.push_back(pl);
#endif
      }
    }
  }

  // Execute the MFGhostSyncPlans for one per-cell field: direct column copies
  // for same-rank partners, non-blocking column exchange for cross-rank ones.
  auto SyncMFPeriodicGhosts = [&](auto& field) {
#ifdef MPI_ENABLED
    constexpr int MF_SYNC_TAG_BASE = 9500; // distinct from PERIODIC_TAG_BASE
    std::vector<std::vector<T>> sendBufs, recvBufs;
    std::vector<MPI_Request> sendReqs, recvReqs;
    std::vector<std::size_t> recvBlock;
    std::vector<int> recvCol, recvNy;
#endif
    for (const auto& pl : MFGhostSyncPlans) {
      auto& fg = field.getBlockField(pl.ghostBlockIdx);
      const auto& pr = Geo.getBlock(pl.ghostBlockIdx).getProjection();
      if (!pl.crossRank) {
        auto& fs = field.getBlockField(pl.srcBlockIdx);
        for (int j = 0; j < pl.ny; ++j)
          fg.get(j * pr[1] + pl.ghostCol) = fs.get(j * pr[1] + pl.srcCol);
      }
#ifdef MPI_ENABLED
      else {
        // pack my physical seam column and post the exchange
        sendBufs.emplace_back(pl.ny);
        auto& sb = sendBufs.back();
        for (int j = 0; j < pl.ny; ++j) sb[j] = fg.get(j * pr[1] + pl.seamCol);
        MPI_Request rq;
        mpi().iSend(sb.data(), pl.ny, pl.partnerRank, &rq, MF_SYNC_TAG_BASE + pl.sendTag);
        sendReqs.push_back(rq);
        recvBufs.emplace_back(pl.ny);
        recvBlock.push_back(pl.ghostBlockIdx);
        recvCol.push_back(pl.ghostCol);
        recvNy.push_back(pl.ny);
        mpi().iRecv(recvBufs.back().data(), pl.ny, pl.partnerRank, &rq, MF_SYNC_TAG_BASE + pl.recvTag);
        recvReqs.push_back(rq);
      }
#endif
    }
#ifdef MPI_ENABLED
    if (!sendReqs.empty()) {
      MPI_Waitall(static_cast<int>(sendReqs.size()), sendReqs.data(), MPI_STATUSES_IGNORE);
      for (std::size_t i = 0; i < recvReqs.size(); ++i) {
        MPI_Wait(&recvReqs[i], MPI_STATUS_IGNORE);
        auto& fg = field.getBlockField(recvBlock[i]);
        const auto& pr = Geo.getBlock(recvBlock[i]).getProjection();
        const auto& rb = recvBufs[i];
        for (int j = 0; j < recvNy[i]; ++j) fg.get(j * pr[1] + recvCol[i]) = rb[j];
      }
    }
#endif
  };

  // -- PF tasks --
  using PFNT=tmp::Key_TypePair<BulkFlag,ff::FF2D<PFCELL>>;
  using PFLT=tmp::Key_TypePair<BulkFlag,ff::FFLaplacian2D<PFCELL>>;
  using PFCT=tmp::Key_TypePair<BulkFlag,ff::FFChemPotential2D<PFCELL>>;
  using PFSelN=TaskSelector<std::uint8_t,PFCELL,PFNT>;
  using PFSelL=TaskSelector<std::uint8_t,PFCELL,PFLT>;
  using PFSelC=TaskSelector<std::uint8_t,PFCELL,PFCT>;
  using PFColT=tmp::Key_TypePair<BulkFlag,
    collision::MRTSource<equilibrium::FirstOrder<PFCELL>,NORMAL<T,2>,true,true>>;
  using PFPerT=tmp::Key_TypePair<PeriodicFlag,collision::PeriodicBoundary<PFCELL>>;
  using PFAll=tmp::TupleWrapper<PFColT,PFPerT>;
  using PFSel=tmp::TaskSelector<PFAll,std::uint8_t,PFCELL>;

  // -- NS tasks --
  using NSMT=tmp::Key_TypePair<BulkFlag,collision::MRTForce<NSCELL,FORCE<T,2>>>;
  using NSPT=tmp::Key_TypePair<PeriodicFlag,collision::PeriodicBoundary<NSCELL>>;
  using NSAll=tmp::TupleWrapper<NSMT,NSPT>;
  using NSSel=tmp::TaskSelector<NSAll,std::uint8_t,NSCELL>;

  // -- Coupling tasks --
  using STT=tmp::Key_TypePair<BulkFlag,ff::FFSurfaceTension2D<PFCELL,NSCELL>>;
  BlockLatManagerCoupling STC(PFLattice,NSLattice);
  using GrT=tmp::Key_TypePair<BulkFlag,ff::FFGravityForce2D<PFCELL,NSCELL>>;
  BlockLatManagerCoupling GrC(PFLattice,NSLattice);
  using PrT=tmp::Key_TypePair<BulkFlag,ff::FFPreForce2D<PFCELL,NSCELL>>;
  BlockLatManagerCoupling PrC(PFLattice,NSLattice);
  using ViT=tmp::Key_TypePair<BulkFlag,ff::FFViscoForce2D<PFCELL,NSCELL>>;
  BlockLatManagerCoupling ViC(PFLattice,NSLattice);
  using RoT=tmp::Key_TypePair<BulkFlag,ff::FFRhoOmegaUpdate2D<PFCELL,NSCELL>>;
  BlockLatManagerCoupling RoC(PFLattice,NSLattice);

  // MF coupling: PF→MF (coeff update)
  using MCT=tmp::Key_TypePair<BulkFlag,MFUpdateCoeffs2D<PFCELL,MFCELL>>;
  using MCSel=CoupledTaskSelector<std::uint8_t,PFCELL,MFCELL,MCT>;
  BlockLatManagerCoupling MCC(PFLattice,MFLattice);

  // Writers
  vtmo::ScalarWriter PW("PHI",PFLattice.getField<PHI<T>>());
  vtmo::ScalarWriter PS("PSI",MFLattice.getField<PSI<T>>());
  vtmo::VectorWriter VW("Velocity",NSLattice.getField<VELOCITY<T,2>>());
  vtmo::ScalarWriter Dw("Density",NSLattice.getField<DENSITY<T>>());
  vtmo::VectorWriter Fw("Force",NSLattice.getField<FORCE<T,2>>());
  vtmo::ScalarWriter Hxw("HX",MFLattice.getField<HX<T>>());
  vtmo::ScalarWriter Hyw("HY",MFLattice.getField<HY<T>>());
  vtmo::vtmWriter<T,2> MW("rosensweig2d",Geo);
  MW.addWriterSet(PW,PS,VW,Dw,Fw,Hxw,Hyw);

  // ===== initial setup =====
  PFLattice.NormalFullCommunicate(); NSLattice.NormalFullCommunicate(); MFLattice.NormalFullCommunicate();
  NS_Per.Apply(); PF_Per.Apply();

  PFLattice.template ApplyInnerCellDynamics<PFSelN>(FlagFM);
  PFLattice.template ApplyInnerCellDynamics<PFSelL>(FlagFM);
  PFLattice.template ApplyInnerCellDynamics<PFSelC>(FlagFM);
  PFLattice.getField<NORMAL<T,2>>().Communicate();
  PFLattice.getField<GRAD<T,2>>().Communicate();
  ff::CommunicateAllSelfFields<T>(PFLattice);
  RoC.ApplyInnerCellDynamics<CoupledTaskSelector<std::uint8_t,PFCELL,NSCELL,RoT>>(0,FlagFM);
  MW.WriteBinary(0);

  Printer::Print_BigBanner(std::string("Start Calculation..."));
  Timer t; Timer ot;
  T H_global=T(Nj)*Cell_Len;

  while(t()<MaxStep){
    // ===== Phase 0: Magnetic field solve =====
    // 0a: Update per-cell mu, chi, omega_psi from phi
    MCC.ApplyInnerCellDynamics<MCSel>(t(),FlagFM);
    CommunicateOMEGAPSI<T>(MFLattice);

    // 0b: Set wall psi = -H0*y  and wall pops = feq(psi_bc)
    // NOTE: lattice row jj holds physical y = minY + (jj-ov)*vs (ov halo
    // rows below the block). Pinning by absolute y (halo + wall rows) gives
    // the exact linear profile, so H == H0 everywhere in the far field.
    {
      auto& psiF=MFLattice.getField<PSI<T>>();
      const T nwall=Cell_Len*T{3.0};
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bPsi=psiF.getBlockField(b);
        int nx=bk.getNx(),ny=bk.getNy(),ov=bk.getOverlap();
        T minY=bk.getMin()[1],vs=bk.getVoxelSize();
        for(int jj=0;jj<ny;++jj){
          T y=minY+T(jj-ov)*vs;
          if((y<=nwall&&y>=-nwall)||(y<=H_global+nwall&&y>=H_global-nwall)){
            T psi_w=-H0*y;
            for(int ii=0;ii<nx;++ii){
              std::size_t id=jj*pr[1]+ii; MFCELL c(id,bl);
              for(unsigned k=0;k<MFLatSet::q;++k) c[k]=latset::w<MFLatSet>(k)*psi_w;
              bPsi.get(id)=psi_w;
            }
          }
        }
      }
    }

    // 0c-0e+: psi-solver sub-iterations (PsiSolver_Iter repeats of
    // collision+stream+macro). Each sub-iteration is one D2Q5 BGK diffusion
    // step with the boosted relaxation omega_psi = 1/(0.5+K*mu); ~100
    // iterations converge the solve to its (smooth) fixed point, removing the
    // grid-scale |H| oscillation that caused the faceted bubble outline.
    for(int sub=0;sub<PsiSolver_Iter;++sub){
      // 0c: MF collision (direct iteration)
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        int ov=bk.getOverlap();
        for(int j=ov;j<bk.getNy()-ov;++j){
          for(int i=ov;i<bk.getNx()-ov;++i){
            MFCELL c(j*pr[1]+i,bl);
            collision::MRTDiffusion<MFCELL,OMEGA_PSI<T>>::apply(c);
          }
        }
      }
      MF_Per.Apply();

      // 0d: Stream
      MFLattice.NormalFullCommunicate();
      MFLattice.Stream();
      MFLattice.NormalFullCommunicate();

      // 0e: PSI = Σg_i (interior only, walls keep psi_bc)
      {
        auto& psiF=MFLattice.getField<PSI<T>>();
        for(int b=0;b<Geo.getBlockNum();++b){
          auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
          auto& bPsi=psiF.getBlockField(b); int ov=bk.getOverlap();
          for(int j=ov;j<bk.getNy()-ov;++j){
            for(int i=ov;i<bk.getNx()-ov;++i){
              std::size_t id=j*pr[1]+i; MFCELL c(id,bl);
              T psi=0; for(unsigned k=0;k<MFLatSet::q;++k)psi+=c[k];
              bPsi.get(id)=psi;
            }
          }
        }
      }
      CommunicatePSI<T>(MFLattice);

      // 0e+: re-pin walls (0e overwrote them with the drifted Σg) so the
      // stored field is exactly -H0*y at the walls/halo before H is computed.
      {
        auto& psiF=MFLattice.getField<PSI<T>>();
        const T nwall=Cell_Len*T{3.0};
        for(int b=0;b<Geo.getBlockNum();++b){
          auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
          auto& bPsi=psiF.getBlockField(b);
          int nx=bk.getNx(),ny=bk.getNy(),ov=bk.getOverlap();
          T minY=bk.getMin()[1],vs=bk.getVoxelSize();
          for(int jj=0;jj<ny;++jj){
            T y=minY+T(jj-ov)*vs;
            if((y<=nwall&&y>=-nwall)||(y<=H_global+nwall&&y>=H_global-nwall)){
              T psi_w=-H0*y;
              for(int ii=0;ii<nx;++ii){
                std::size_t id=jj*pr[1]+ii; MFCELL c(id,bl);
                for(unsigned k=0;k<MFLatSet::q;++k) c[k]=latset::w<MFLatSet>(k)*psi_w;
                bPsi.get(id)=psi_w;
              }
            }
          }
        }
      }
    }
    // One sync per step is enough: the sub-iterations never read neighbor PSI
    // (collision uses own pops, 0e uses own cells), and the H ghosts are only
    // read by A4 below, after the sync that follows 0f.
    SyncMFPeriodicGhosts(MFLattice.getField<PSI<T>>());

    // 0f: Compute H = -∇ψ
    for(int b=0;b<Geo.getBlockNum();++b){
      auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
      int ov=bk.getOverlap();
      for(int j=ov;j<bk.getNy()-ov;++j){
        for(int i=ov;i<bk.getNx()-ov;++i){
          MFCELL c(j*pr[1]+i,bl);
          MFComputeH2D<MFCELL>::apply(c);
        }
      }
    }
    CommunicateAllMFFields<T>(MFLattice);
    SyncMFPeriodicGhosts(MFLattice.getField<HX<T>>());
    SyncMFPeriodicGhosts(MFLattice.getField<HY<T>>());
    SyncMFPeriodicGhosts(MFLattice.getField<HMAG<T>>());

    // ===== Phase A: Force setup =====
    RoC.ApplyInnerCellDynamics<CoupledTaskSelector<std::uint8_t,PFCELL,NSCELL,RoT>>(t(),FlagFM);
    NSLattice.getField<FORCE<T,2>>().InitValue(Vector<T,2>{0,0});
    STC.ApplyInnerCellDynamics<CoupledTaskSelector<std::uint8_t,PFCELL,NSCELL,STT>>(t(),FlagFM);

    // A4: Magnetic force (if H0>0)
    // Paper Eq. (8): F_m = (μ₀χ/2)∇(|H|²) = μ₀χ|H|∇|H| (μ₀≡1, Kelvin form).
    // The full interfacial stress (χ_h/2)∫φ·d|H|² is delivered ONLY when the
    // μ-jump is co-located with the φ-midpoint — i.e. the interface is sharp
    // (Interface_Width ≈ 1-2 cells, as in the paper's AMR finest level).
    if(H0>T{0}){
      // Force-free wall bands: the wall ψ-pins and the D2Q5 ψ-solver's ~1.4%
      // far-field slope mismatch leave a spurious |H| spike (|H| ~ 4-6·H0,
      // worst at the x-corners) within ~6 rows of each wall. The Kelvin force
      // there (χ|H|∇|H|) is a numerical artifact, so it is zeroed in the
      // bands |y|<=12 and |y-H_global|<=12. (Physically the far-field Kelvin
      // force is ~0 anyway: |H| is uniform, ∇|H|≈0.)
      const T band=T{12.0};
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& pf_bl=PFLattice.getBlockLat(b); auto& mf_bl=MFLattice.getBlockLat(b);
        auto& ns_bl=NSLattice.getBlockLat(b);
        const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        int ov=bk.getOverlap();
        T minY=bk.getMin()[1],vs=bk.getVoxelSize();
        for(int j=ov;j<bk.getNy()-ov;++j){
          T y=minY+T(j-ov)*vs;
          if(y<=band||y>=H_global-band) continue;
          for(int i=ov;i<bk.getNx()-ov;++i){
            std::size_t id=j*pr[1]+i;
            PFCELL pf(id,pf_bl); MFCELL mf(id,mf_bl); NSCELL ns(id,ns_bl);
            MFMagneticForce2D<PFCELL,MFCELL,NSCELL>::apply(pf,mf,ns);
          }
        }
      }
    }

    GrC.ApplyInnerCellDynamics<CoupledTaskSelector<std::uint8_t,PFCELL,NSCELL,GrT>>(t(),FlagFM);
    PrC.ApplyInnerCellDynamics<CoupledTaskSelector<std::uint8_t,PFCELL,NSCELL,PrT>>(t(),FlagFM);
    NSLattice.getField<FORCE<T,2>>().Communicate();

    // ===== Phase B-C: PF + NS collision =====
    PFLattice.template ApplyInnerCellDynamics<PFSel>(FlagFM);
    PF_Per.Apply(); PFLattice.NormalFullCommunicate();

    ViC.ApplyInnerCellDynamics<CoupledTaskSelector<std::uint8_t,PFCELL,NSCELL,ViT>>(t(),FlagFM);
    NSLattice.template ApplyInnerCellDynamics<NSSel>(FlagFM);
    NS_Per.Apply(); NSLattice.NormalFullCommunicate();

    // ===== Phase D: Streaming =====
    PF_BB.Apply(t());
    { // PF Y ghost pop fix
      T Hg=T(Nj)*Cell_Len;
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& bl=PFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        int nx=bk.getNx(),ny=bk.getNy(),ov=bk.getOverlap();
        T my=bk.getMin()[1],My=bk.getMax()[1];
        if(my<Cell_Len*T{1.5})
          for(int i=0;i<nx;++i){PFCELL g(0*pr[1]+i,bl),w(ov*pr[1]+i,bl); for(unsigned k=0;k<LatSet::q;++k)g[k]=w[k];}
        if(My>Hg-Cell_Len*T{1.5})
          for(int i=0;i<nx;++i){PFCELL g((ny-1)*pr[1]+i,bl),w((ny-1-ov)*pr[1]+i,bl); for(unsigned k=0;k<LatSet::q;++k)g[k]=w[k];}
      }
    }
    PFLattice.Stream(); PFLattice.NormalFullCommunicate();

    NS_BB.Apply(t());
    { // NS Y ghost pop fix
      T Hg=T(Nj)*Cell_Len;
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& bl=NSLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        int nx=bk.getNx(),ny=bk.getNy(),ov=bk.getOverlap();
        T my=bk.getMin()[1],My=bk.getMax()[1];
        if(my<Cell_Len*T{1.5})
          for(int i=0;i<nx;++i){NSCELL g(0*pr[1]+i,bl),w(ov*pr[1]+i,bl); for(unsigned k=0;k<LatSet::q;++k)g[k]=w[k];}
        if(My>Hg-Cell_Len*T{1.5})
          for(int i=0;i<nx;++i){NSCELL g((ny-1)*pr[1]+i,bl),w((ny-1-ov)*pr[1]+i,bl); for(unsigned k=0;k<LatSet::q;++k)g[k]=w[k];}
      }
    }
    NSLattice.Stream(); NSLattice.NormalFullCommunicate();

    // ===== Phase E: Macro update =====
    // phi = sum(g_i)
    {
      auto& pF=PFLattice.getField<PHI<T>>();
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& bl=PFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bP=pF.getBlockField(b); int ov=bk.getOverlap();
        for(int j=ov;j<bk.getNy()-ov;++j)
          for(int i=ov;i<bk.getNx()-ov;++i){
            std::size_t id=j*pr[1]+i; PFCELL c(id,bl);
            T pn=0; for(unsigned k=0;k<LatSet::q;++k)pn+=c[k];
            if(pn<T{0})pn=T{0}; if(pn>T{1})pn=T{1}; bP.get(id)=pn;
          }
      }
    }
    PFLattice.getField<PHI<T>>().Communicate();

    // wall phi BC: bottom=1 (ferrofluid), top=0 (organic solvent)
    {
      auto& pF=PFLattice.getField<PHI<T>>();
      for(int b=0;b<Geo.getBlockNum();++b){
        const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bP=pF.getBlockField(b); int nx=bk.getNx(),ny=bk.getNy(),ov=bk.getOverlap();
        T my=bk.getMin()[1],My=bk.getMax()[1];
        if(my<Cell_Len*T{1.5}) for(int i=0;i<nx;++i)bP.get(ov*pr[1]+i)=T{1};
        if(My>H_global-Cell_Len*T{1.5}){int jj=ny-1-ov;for(int i=0;i<nx;++i)bP.get(jj*pr[1]+i)=T{0};}
      }
    }
    PFLattice.getField<PHI<T>>().Communicate();

    // wall phi ghost fix: bottom=1, top=0
    {
      auto& pF=PFLattice.getField<PHI<T>>();
      for(int b=0;b<Geo.getBlockNum();++b){
        const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bP=pF.getBlockField(b); int nx=bk.getNx(),ny=bk.getNy(),ov=bk.getOverlap();
        T my=bk.getMin()[1],My=bk.getMax()[1];
        if(my<Cell_Len*T{1.5}) for(int j=0;j<ov;++j) for(int i=0;i<nx;++i) bP.get(j*pr[1]+i)=T{1};
        if(My>H_global-Cell_Len*T{1.5}) for(int j=ny-ov;j<ny;++j) for(int i=0;i<nx;++i) bP.get(j*pr[1]+i)=T{0};
      }
    }

    // gradients
    PFLattice.template ApplyInnerCellDynamics<PFSelN>(FlagFM);
    PFLattice.template ApplyInnerCellDynamics<PFSelL>(FlagFM);
    PFLattice.template ApplyInnerCellDynamics<PFSelC>(FlagFM);
    PFLattice.getField<NORMAL<T,2>>().Communicate();
    PFLattice.getField<GRAD<T,2>>().Communicate();
    ff::CommunicateAllSelfFields<T>(PFLattice);

    // wall grad/chempot fix (PF)
    {
      auto& gF=PFLattice.getField<GRAD<T,2>>();
      for(int b=0;b<Geo.getBlockNum();++b){
        const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bG=gF.getBlockField(b); int nx=bk.getNx(),ov=bk.getOverlap();
        T my=bk.getMin()[1],My=bk.getMax()[1];
        if(my<Cell_Len*T{1.5}) for(int i=ov;i<nx-ov;++i) bG.get(ov*pr[1]+i)[1]=bG.get((ov+1)*pr[1]+i)[1];
        if(My>H_global-Cell_Len*T{1.5}){int jj=bk.getNy()-1-ov;for(int i=ov;i<nx-ov;++i)bG.get(jj*pr[1]+i)[1]=bG.get((jj-1)*pr[1]+i)[1];}
      }
    }
    {
      auto& cF=PFLattice.getField<ff::CHEMICALPOTENTIAL<T>>();
      for(int b=0;b<Geo.getBlockNum();++b){
        const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bC=cF.getBlockField(b); int nx=bk.getNx(),ov=bk.getOverlap();
        T my=bk.getMin()[1],My=bk.getMax()[1];
        if(my<Cell_Len*T{1.5}) for(int i=ov;i<nx-ov;++i) bC.get(ov*pr[1]+i)=(T{4}*bC.get((ov+1)*pr[1]+i)-bC.get((ov+2)*pr[1]+i))/T{3};
        if(My>H_global-Cell_Len*T{1.5}){int jj=bk.getNy()-1-ov;for(int i=ov;i<nx-ov;++i)bC.get(jj*pr[1]+i)=(T{4}*bC.get((jj-1)*pr[1]+i)-bC.get((jj-2)*pr[1]+i))/T{3};}
      }
    }
    ff::CommunicateAllSelfFields<T>(PFLattice);

    // NS macro
    {
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& bl=NSLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& rF=NSLattice.getField<DENSITY<T>>(); auto& pF=NSLattice.getField<PRESSURE<T>>();
        auto& vF=NSLattice.getField<VELOCITY<T,2>>(); auto& fF=NSLattice.getField<FORCE<T,2>>();
        auto& bR=rF.getBlockField(b); auto& bP=pF.getBlockField(b);
        auto& bV=vF.getBlockField(b); auto& bF=fF.getBlockField(b);
        int ov=bk.getOverlap();
        for(int j=ov;j<bk.getNy()-ov;++j)
          for(int i=ov;i<bk.getNx()-ov;++i){
            std::size_t id=j*pr[1]+i; NSCELL c(id,bl);
            T p=0,ux=0,uy=0;
            for(unsigned k=0;k<LatSet::q;++k){p+=c[k];ux+=latset::c<LatSet>(k)[0]*c[k];uy+=latset::c<LatSet>(k)[1]*c[k];}
            T rho=bR.get(id); auto F=bF.get(id);
            bP.get(id)=p; bV.get(id)=Vector<T,2>{ux+T{0.5}*F[0]/rho,uy+T{0.5}*F[1]/rho};
          }
      }
    }
    NSLattice.getField<VELOCITY<T,2>>().Communicate();
    NSLattice.getField<PRESSURE<T>>().Communicate();
    NSLattice.getField<DENSITY<T>>().Communicate();
    NSLattice.getField<OMEGA<T>>().Communicate();

    ++t; ++ot;
    if(t()%OutputStep==0){
      ot.Print_InnerLoopPerformance(Geo.getTotalCellNum(),OutputStep);
      Printer::Endl();
      MW.WriteBinary(t());
    }
  }

  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  t.Print_MainLoopPerformance(Geo.getTotalCellNum());
  return 0;
}
