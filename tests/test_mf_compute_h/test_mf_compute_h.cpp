// test_mf_compute_h.cpp — 1000-step test for MFComputeH2D (H = -∇ψ)
// Sets ψ = -H₀·y, verifies H_x=0, H_y=H₀, |H|=H₀
#include "freelb.h"
#include "freelb.hh"

using T = double;
using MFLatSet = D2Q5<T>;
using namespace mfield;

constexpr int N = 8;
constexpr T Cell_Len = 1.0;
constexpr int BlockCellLen = 8;
constexpr T H0 = T{0.1};
constexpr T tol = 1e-12;

int main(int argc, char** argv) {
  mpi().init(&argc, &argv);
  constexpr std::uint8_t VoidFlag = 1, BulkFlag = 2;

  AABB<T,2> domain({0,0},{T(N*Cell_Len),T(N*Cell_Len)});
  BlockGeometryHelper2D<T> Gh(N,N,domain,Cell_Len,BlockCellLen);
  Gh.CreateBlocks(1,1); Gh.AdaptiveOptimization(1); Gh.LoadBalancing(1);
  BlockGeometry2D<T> Geo(Gh);

  BlockFieldManager<FLAG,T,2> FlagFM(Geo,VoidFlag);
  FlagFM.forEach(domain,[&](FLAG&f,std::size_t i){f.SetField(i,BulkFlag);});

  BaseConverter<T> Conv(MFLatSet::cs2);
  Conv.SimplifiedConverterFromRT(N,T(0.01),T{1.0});

  using MFF = TypePack<PSI<T>,POP<T,MFLatSet::q>,HX<T>,HY<T>,HMAG<T>>;
  ValuePack Init(T{},T{},T{},T{},T{});
  using MFCELL = Cell<T,MFLatSet,MFF>;
  BlockLatticeManager<T,MFLatSet,MFF> MFLattice(Geo,Init,Conv);

  // Set ψ = -H₀·y everywhere
  T H_global = T(N)*Cell_Len;
  {
    auto& psiF = MFLattice.getField<PSI<T>>();
    for(int b=0;b<Geo.getBlockNum();++b){
      auto& bl = MFLattice.getBlockLat(b);
      const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
      auto& bPsi=psiF.getBlockField(b);
      T vs=bk.getVoxelSize(),minY=bk.getMin()[1];
      for(int j=0;j<bk.getNy();++j){
        T y=minY+T(j)*vs, psi=-H0*y;
        for(int i=0;i<bk.getNx();++i)
          bPsi.get(j*pr[1]+i)=psi;
      }
    }
  }

  MFLattice.getField<PSI<T>>().Communicate();

  // Run 1000 calls
  for(int step=0;step<1000;++step){
    for(int b=0;b<Geo.getBlockNum();++b){
      auto& bl=MFLattice.getBlockLat(b);
      const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
      int ov=bk.getOverlap();
      for(int j=ov;j<bk.getNy()-ov;++j){
        for(int i=ov;i<bk.getNx()-ov;++i){
          MFCELL cell(j*pr[1]+i,bl);
          MFComputeH2D<MFCELL>::apply(cell);
        }
      }
    }
  }

  // Verify at center cell
  T max_err=0;
  {
    const auto& bk=Geo.getBlock(0); const auto& pr=bk.getProjection();
    auto& hxF=MFLattice.getField<HX<T>>();
    auto& hyF=MFLattice.getField<HY<T>>();
    auto& hmF=MFLattice.getField<HMAG<T>>();
    auto& bHx=hxF.getBlockField(0); auto& bHy=hyF.getBlockField(0); auto& bHm=hmF.getBlockField(0);
    int ov=bk.getOverlap();
    for(int j=ov;j<bk.getNy()-ov;++j){
      for(int i=ov;i<bk.getNx()-ov;++i){
        std::size_t id=j*pr[1]+i;
        T e_hx=std::abs(bHx.get(id)-T{0});
        T e_hy=std::abs(bHy.get(id)-H0);
        T e_hm=std::abs(bHm.get(id)-H0);
        if(e_hx>max_err)max_err=e_hx;
        if(e_hy>max_err)max_err=e_hy;
        if(e_hm>max_err)max_err=e_hm;
      }
    }
    // Print center cell
    int mid_j=(bk.getNy()-ov+ov)/2, mid_i=(bk.getNx()-ov+ov)/2;
    std::size_t id=mid_j*pr[1]+mid_i;
    printf("Center: Hx=%.12f Hy=%.12f |H|=%.12f (expected Hx=0, Hy=%.1f, |H|=%.1f)\n",
           bHx.get(id),bHy.get(id),bHm.get(id),H0,H0);
  }

  bool pass=max_err<tol;
  printf("Step 1000: max_err=%.2e %s\n",max_err,pass?"PASS":"FAIL");
  return pass?0:1;
}
