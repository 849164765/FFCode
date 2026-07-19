// test_mf_force.cpp — 1000-step test for MFMagneticForce2D
// Uniform H field → zero spurious force validation
#include "freelb.h"
#include "freelb.hh"

using T = double;
using MFLatSet = D2Q5<T>;
using NSLatSet = D2Q9<T>;
using namespace mfield;

constexpr int N = 8, BlockCellLen = 8;
constexpr T Cell_Len = 1.0;
constexpr T H0 = T{0.1};
constexpr T tol = 1e-12;

int main(int argc, char** argv) {
  mpi().init(&argc, &argv);
  constexpr std::uint8_t VoidFlag=1,BulkFlag=2;

  AABB<T,2> domain({0,0},{T(N*Cell_Len),T(N*Cell_Len)});
  BlockGeometryHelper2D<T> Gh(N,N,domain,Cell_Len,BlockCellLen);
  Gh.CreateBlocks(1,1); Gh.AdaptiveOptimization(1); Gh.LoadBalancing(1);
  BlockGeometry2D<T> Geo(Gh);

  BlockFieldManager<FLAG,T,2> FlagFM(Geo,VoidFlag);
  FlagFM.forEach(domain,[&](FLAG&f,std::size_t i){f.SetField(i,BulkFlag);});

  // PF lattice: PHI only (phi=0.5)
  BaseConverter<T> PFConv(NSLatSet::cs2); PFConv.SimplifiedConverterFromRT(N,T(0.01),T{1.0});
  using PFF = TypePack<PHI<T>>; ValuePack PFI(T{0.5});
  using PFCELL = Cell<T,NSLatSet,PFF>;
  BlockLatticeManager<T,NSLatSet,PFF> PFLattice(Geo,PFI,PFConv);

  // MF lattice: HMAG, HX, HY, CHI_L, CHI_H (uniform |H|=H0)
  BaseConverter<T> MFConv(MFLatSet::cs2); MFConv.SimplifiedConverterFromRT(N,T(0.01),T{1.0});
  using MFF = TypePack<HMAG<T>,HX<T>,HY<T>,CHI_L<T>,CHI_H<T>>;
  ValuePack MFI(H0,T{0},H0,T{0},T{1}); // HMAG=0.1, HX=0, HY=0.1, chi_l=0, chi_h=1
  using MFCELL = Cell<T,MFLatSet,MFF>;
  BlockLatticeManager<T,MFLatSet,MFF> MFLattice(Geo,MFI,MFConv);

  // NS lattice: FORCE<2> only
  BaseConverter<T> NSConv(NSLatSet::cs2); NSConv.SimplifiedConverterFromRT(N,T(0.01),T{1.0});
  using NSFF = TypePack<FORCE<T,2>>; ValuePack NSI(Vector<T,2>{0,0});
  using NSCELL = Cell<T,NSLatSet,NSFF>;
  BlockLatticeManager<T,NSLatSet,NSFF> NSLattice(Geo,NSI,NSConv);

  // init MF block-level chi params (already set via Init)
  MFLattice.getField<CHI_L<T>>().InitValue(T{0});
  MFLattice.getField<CHI_H<T>>().InitValue(T{1});

  // run 1000 steps: compute F_mag for each cell
  for(int step=0;step<1000;++step){
    NSLattice.getField<FORCE<T,2>>().InitValue(Vector<T,2>{0,0});
    for(int b=0;b<Geo.getBlockNum();++b){
      auto& pf_bl = PFLattice.getBlockLat(b);
      auto& mf_bl = MFLattice.getBlockLat(b);
      auto& ns_bl = NSLattice.getBlockLat(b);
      const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
      int ov=bk.getOverlap();
      for(int j=ov;j<bk.getNy()-ov;++j){
        for(int i=ov;i<bk.getNx()-ov;++i){
          std::size_t id=j*pr[1]+i;
          PFCELL pf_cell(id,pf_bl);
          MFCELL mf_cell(id,mf_bl);
          NSCELL ns_cell(id,ns_bl);
          MFMagneticForce2D<PFCELL,MFCELL,NSCELL>::apply(pf_cell,mf_cell,ns_cell);
        }
      }
    }
  }

  // verify max|F_mag| < tol
  T max_err=0;
  {
    auto& fF=NSLattice.getField<FORCE<T,2>>();
    for(int b=0;b<Geo.getBlockNum();++b){
      auto& bF=fF.getBlockField(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
      int ov=bk.getOverlap();
      for(int j=ov;j<bk.getNy()-ov;++j){
        for(int i=ov;i<bk.getNx()-ov;++i){
          auto F=bF.get(j*pr[1]+i);
          T e=std::abs(F[0])+std::abs(F[1]);
          if(e>max_err)max_err=e;
        }
      }
    }
  }
  bool pass=max_err<tol;
  printf("Step 1000: max|F_mag|=%.2e (expected 0) %s\n",max_err,pass?"PASS":"FAIL");
  return pass?0:1;
}
