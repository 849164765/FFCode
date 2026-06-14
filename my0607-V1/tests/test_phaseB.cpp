// test_phaseB.cpp — Phase B: h-LBE collision & normal tests
// Validates H_Collision (Eq.6) and UpdatePhaseField (Eq.12)

#include "freelb.h"
#include "freelb.hh"
#include "lbm/phase_gradient.ur.h"
#include "lbm/h_collision.ur.h"

using T = double;
using LatSet = D2Q9<T>;
constexpr std::uint8_t Void = 1, Bulk = 2;


// ===========================================================================
// Test 1: ComputeNormalFromGradient
// ===========================================================================
static int testNormal() {
  constexpr int Mn=8;
  AABB<T,2> dm(Vector<T,2>(T{0},T{0}), Vector<T,2>(T(Mn),T(Mn)));
  BlockGeometryHelper2D<T> gh(Mn,Mn,dm,T{1},Mn);
  gh.CreateBlocks(1,mpi().getSize()); gh.AdaptiveOptimization(mpi().getSize());
  gh.LoadBalancing(mpi().getSize());
  BlockGeometry2D<T> geo(gh);
  BlockFieldManager<FLAG,T,2> fm(geo,Void);
  fm.forEach(dm,[&](FLAG& f,std::size_t id){f.SetField(id,Bulk);});
  fm.template SetupBoundary<LatSet>(dm,Bulk);
  BaseConverter<T> cv(LatSet::cs2); cv.SimplifiedConverterFromRT(Mn,T{1},T{1.0});
  using PF=TypePack<PHI<T>,GRAD<T,2>,NORMAL<T,2>>;
  ValuePack ivs(T{},Vector<T,2>{T{0},T{0}},Vector<T,2>{T{0},T{0}});
  BlockLatticeManager<T,LatSet,PF> Lat(geo,ivs,cv);

  // φ(x,y)=x² via Cell
  auto& bl=Lat.getBlockLat(0); const auto& b=bl.getBlock(); int ol=b.getOverlap();
  const auto& pr=b.getProjection(); T vs=b.getVoxelSize(), bx=b.getMin()[0];
  using MC=Cell<T,LatSet,PF>;
  for(int j=0;j<b.getNy();++j) for(int i=0;i<b.getNx();++i)
    { MC c(j*pr[1]+i,bl); T x=bx+T(i)*vs; c.template get<PHI<T>>()=x*x; }
  Lat.template ApplyInnerCellDynamics<phase_gradient::ComputeGradientPhi<MC>>();
  Lat.template ApplyInnerCellDynamics<phase_gradient::ComputeNormalFromGradient<MC>>();

  int pass=0,fail=0; T eps=T{1e-10};
  for(int j=ol;j<b.getNy()-ol;++j) for(int i=ol;i<b.getNx()-ol;++i) {
    MC c(j*pr[1]+i,bl); const auto& g=c.template get<GRAD<T,2>>();
    const auto& n=c.template get<NORMAL<T,2>>(); T gm=g.getnorm();
    Vector<T,2> en{T{0},T{0}};
    if(gm>=T{0.005}) for(int d=0;d<2;++d) en[d]=g[d]/gm;
    bool ok=true; for(int d=0;d<2;++d) if(std::abs(n[d]-en[d])>eps) ok=false;
    ok?pass++:fail++;
  }
  if(mpi().isMainProcessor())
    std::cout<<"[normal] "<<pass<<" passed, "<<fail<<" failed"<<std::endl;
  return fail;
}


// ===========================================================================
// Test 2: Equilibrium invariance — h at equilibrium stays there after collision
// ===========================================================================
static int testCollisionEq() {
  constexpr int Mn=8;
  AABB<T,2> dm(Vector<T,2>(T{0},T{0}), Vector<T,2>(T(Mn),T(Mn)));
  BlockGeometryHelper2D<T> gh(Mn,Mn,dm,T{1},Mn);
  gh.CreateBlocks(1,mpi().getSize()); gh.AdaptiveOptimization(mpi().getSize());
  gh.LoadBalancing(mpi().getSize());
  BlockGeometry2D<T> geo(gh);
  BlockFieldManager<FLAG,T,2> fm(geo,Void);
  fm.forEach(dm,[&](FLAG& f,std::size_t id){f.SetField(id,Bulk);});
  fm.template SetupBoundary<LatSet>(dm,Bulk);
  BaseConverter<T> cv(LatSet::cs2); cv.SimplifiedConverterFromRT(Mn,T{1},T{1.0});
  using PF=TypePack<PHI<T>,GRAD<T,2>,NORMAL<T,2>,VELOCITY<T,2>,INTERFACEWIDTH<T>,POP<T,LatSet::q>>;
  ValuePack ivs(T{},Vector<T,2>{T{0},T{0}},Vector<T,2>{T{0},T{0}},Vector<T,2>{T{0},T{0}},T{4},T{});
  BlockLatticeManager<T,LatSet,PF> Lat(geo,ivs,cv);

  auto& bl=Lat.getBlockLat(0); const auto& b=bl.getBlock();
  int ol=b.getOverlap(); const auto& p=b.getProjection();

  // Init all fields via Cell
  using MC=Cell<T,LatSet,PF>;
  for(int j=0;j<b.getNy();++j) for(int i=0;i<b.getNx();++i) {
    MC c(j*p[1]+i,bl);
    c.template get<PHI<T>>()=T{0.5};
    c.template get<VELOCITY<T,2>>()=Vector<T,2>{T{0},T{0}};
    c.template get<NORMAL<T,2>>()=Vector<T,2>{T{0},T{0}};
    c.template get<INTERFACEWIDTH<T>>()=T{4};
  }

  // Init h to equilibrium (u=0, phi=0.5 → h_eq_α = w_α * 0.5)
  MC cell(ol*p[1]+ol,bl);
  for(unsigned int k=0;k<LatSet::q;++k) cell[k]=latset::w<LatSet>(k)*T{0.5};
  T old_h[LatSet::q];
  for(unsigned int k=0;k<LatSet::q;++k) old_h[k]=cell[k];

  // Collision: n=0 → zero forcing, equilibrium unchanged
  phase_field::H_Collision<MC>::apply(cell);

  int pass=0,fail=0; T eps=T{1e-10};
  for(unsigned int k=0;k<LatSet::q;++k) {
    T d=std::abs(cell[k]-old_h[k]);
    d<=eps?pass++:fail++;
    if(d>eps) std::cerr<<"  dir["<<k<<"] before="<<old_h[k]<<" after="<<cell[k]<<std::endl;
  }
  if(mpi().isMainProcessor())
    std::cout<<"[h-collision-eq] "<<pass<<" passed, "<<fail<<" failed"<<std::endl;
  return fail;
}


// ===========================================================================
// Test 3: φ = Σ h_α (Eq.12)
// ===========================================================================
static int testPhiUpdate() {
  constexpr int Mn=8;
  AABB<T,2> dm(Vector<T,2>(T{0},T{0}), Vector<T,2>(T(Mn),T(Mn)));
  BlockGeometryHelper2D<T> gh(Mn,Mn,dm,T{1},Mn);
  gh.CreateBlocks(1,mpi().getSize()); gh.AdaptiveOptimization(mpi().getSize());
  gh.LoadBalancing(mpi().getSize());
  BlockGeometry2D<T> geo(gh);
  BlockFieldManager<FLAG,T,2> fm(geo,Void);
  fm.forEach(dm,[&](FLAG& f,std::size_t id){f.SetField(id,Bulk);});
  fm.template SetupBoundary<LatSet>(dm,Bulk);
  BaseConverter<T> cv(LatSet::cs2); cv.SimplifiedConverterFromRT(Mn,T{1},T{1.0});
  using PF=TypePack<PHI<T>,POP<T,LatSet::q>>;
  ValuePack ivs(T{},T{});
  BlockLatticeManager<T,LatSet,PF> Lat(geo,ivs,cv);
  auto& bl=Lat.getBlockLat(0); const auto& b=bl.getBlock();
  int ol=b.getOverlap(); const auto& p=b.getProjection();

  using MC=Cell<T,LatSet,PF>;
  MC cell(ol*p[1]+ol,bl);
  for(unsigned int k=0;k<LatSet::q;++k) cell[k]=latset::w<LatSet>(k)*T{0.5};
  phase_field::UpdatePhaseField<MC>::apply(cell);
  T phi=cell.template get<PHI<T>>();
  bool ok=std::abs(phi-T{0.5})<T{1e-14};
  if(mpi().isMainProcessor()) std::cout<<"[phi-update] phi="<<phi<<" "<<(ok?"PASS":"FAIL")<<std::endl;
  return ok?0:1;
}


int main(int argc,char** argv){
  mpi().init(&argc,&argv);
  int total=0;
  total+=testNormal();
  total+=testCollisionEq();
  total+=testPhiUpdate();
  if(mpi().isMainProcessor()){
    if(total==0) std::cout<<"\n=== Phase B: ALL TESTS PASSED ==="<<std::endl;
    else std::cout<<"\n=== Phase B: "<<total<<" FAILURES ==="<<std::endl;
  }
  return total;
}
