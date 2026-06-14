// test_phaseC.cpp — Phase C: g-LBE hydrodynamics
// Validates C1 (density), C5 (chemical potential), C13 (pressure/velocity),
// G_Collision equilibrium invariance

#include "freelb.h"
#include "freelb.hh"
#include "lbm/phase_c.ur.h"

using T = double;
using LatSet = D2Q9<T>;
constexpr uint8_t Void=1, Bulk=2;

static auto makeGeo(int Mn=8) {
  T cl=T{1};
  AABB<T,2> dm(Vector<T,2>(T{0},T{0}), Vector<T,2>(T(Mn*cl),T(Mn*cl)));
  BlockGeometryHelper2D<T> gh(Mn,Mn,dm,cl,Mn);
  gh.CreateBlocks(1,mpi().getSize()); gh.AdaptiveOptimization(mpi().getSize());
  gh.LoadBalancing(mpi().getSize());
  BlockGeometry2D<T> g(gh);
  BlockFieldManager<FLAG,T,2> fm(g,Void);
  fm.forEach(dm,[&](FLAG& f,std::size_t id){f.SetField(id,Bulk);});
  fm.template SetupBoundary<LatSet>(dm,Bulk);
  return std::make_pair(std::move(g), std::move(fm));
}


// ===========================================================================
// Test 1: ComputeDensityFromPhase (C1)
// φ=0.3, ρ_L=1, ρ_H=10 → ρ = 1 + 0.3*9 = 3.7
// ===========================================================================
static int testDensity() {
  auto [geo,fm]=makeGeo();
  BaseConverter<T> cv(LatSet::cs2); cv.SimplifiedConverterFromRT(8,T{1},T{1.0});
  using FI=TypePack<PHI<T>, RHO<T>, ff::RHO_L<T>, ff::RHO_H<T>>;
  ValuePack ivs(T{},T{},T{1},T{10});
  BlockLatticeManager<T,LatSet,FI> Lat(geo,ivs,cv);
  using MC=Cell<T,LatSet,FI>;
  auto& bl=Lat.getBlockLat(0); int ol=bl.getBlock().getOverlap();
  const auto& pr=bl.getBlock().getProjection();
  MC cell(ol*pr[1]+ol,bl);
  Lat.template getField<ff::RHO_L<T>>().InitValue(T{1});
  Lat.template getField<ff::RHO_H<T>>().InitValue(T{10});
  cell.template get<PHI<T>>()=T{0.3};
  phase_field::ComputeDensityFromPhase<MC>::apply(cell);
  T rho=cell.template get<RHO<T>>(), exp=T{3.7};
  bool ok=std::abs(rho-exp)<T{1e-14};
  if(mpi().isMainProcessor())
    std::cout<<"[density] phi=0.3 rho="<<rho<<" exp=3.7 "<<(ok?"PASS":"FAIL")<<std::endl;
  return ok?0:1;
}


// ===========================================================================
// Test 2: ComputeChemicalPotential (C5)
// φ=0.5, β=0.25, κ=1, ∇²φ=2 → μ_φ = 4*0.25*0.5*(-0.5)*0 - 2 = -2
// ===========================================================================
static int testChemPot() {
  auto [geo,fm]=makeGeo();
  BaseConverter<T> cv(LatSet::cs2); cv.SimplifiedConverterFromRT(8,T{1},T{1.0});
  using FI=TypePack<PHI<T>, ff::LAPLACIAN<T>, ff::CHEMICALPOTENTIAL<T>, ff::BETA<T>, ff::KAPPA<T>>;
  ValuePack ivs(T{},T{},T{},T{0.25},T{1});
  BlockLatticeManager<T,LatSet,FI> Lat(geo,ivs,cv);
  using MC=Cell<T,LatSet,FI>;
  auto& bl=Lat.getBlockLat(0); int ol=bl.getBlock().getOverlap();
  const auto& pr=bl.getBlock().getProjection();
  MC cell(ol*pr[1]+ol,bl);
  Lat.template getField<ff::BETA<T>>().InitValue(T{0.25});
  Lat.template getField<ff::KAPPA<T>>().InitValue(T{1});
  cell.template get<PHI<T>>()=T{0.5};
  cell.template get<ff::LAPLACIAN<T>>()=T{2};
  phase_field::ComputeChemicalPotential<MC>::apply(cell);
  T mu=cell.template get<ff::CHEMICALPOTENTIAL<T> >(), exp=T{-2};
  bool ok=std::abs(mu-exp)<T{1e-14};
  if(mpi().isMainProcessor())
    std::cout<<"[chemPot] phi=0.5 mu="<<mu<<" exp=-2 "<<(ok?"PASS":"FAIL")<<std::endl;
  return ok?0:1;
}


// ===========================================================================
// Test 3: C13 UpdatePressureVelocity — p*=Σg, u=Σg*e+F/(2ρ)
// g_α = w_α*1.0, F=0, ρ=1 → p*=1, u=(0,0)
// ===========================================================================
static int testMacroUpdate() {
  auto [geo,fm]=makeGeo();
  BaseConverter<T> cv(LatSet::cs2); cv.SimplifiedConverterFromRT(8,T{1},T{1.0});
  using FI=TypePack<RHO<T>, VELOCITY<T,2>, FORCE<T,2>, POP<T,LatSet::q>>;
  ValuePack ivs(T{},Vector<T,2>{T{0},T{0}},Vector<T,2>{T{0},T{0}},T{});
  BlockLatticeManager<T,LatSet,FI> Lat(geo,ivs,cv);
  using MC=Cell<T,LatSet,FI>;
  auto& bl=Lat.getBlockLat(0); int ol=bl.getBlock().getOverlap();
  const auto& pr=bl.getBlock().getProjection();
  MC cell(ol*pr[1]+ol,bl);
  for(unsigned int k=0;k<LatSet::q;++k) cell[k]=latset::w<LatSet>(k);
  cell.template get<FORCE<T,2>>()=Vector<T,2>{T{0},T{0}};
  cell.template get<RHO<T>>()=T{1};
  phase_field::UpdatePressureVelocity<MC>::apply(cell);
  T ps=cell.template get<RHO<T>>();
  Vector<T,2> u=cell.template get<VELOCITY<T,2>>();
  bool ok=std::abs(ps-T{1})<T{1e-14};
  for(int d=0;d<2;++d) if(std::abs(u[d])>T{1e-14}) ok=false;
  if(mpi().isMainProcessor())
    std::cout<<"[macro] p*="<<ps<<" u=("<<u[0]<<","<<u[1]<<") "<<(ok?"PASS":"FAIL")<<std::endl;
  return ok?0:1;
}


// ===========================================================================
// Test 4: G_Collision equilibrium invariance
// g at equilibrium, zero force → collision preserves it
// ===========================================================================
static int testGCollisionEq() {
  auto [geo,fm]=makeGeo();
  BaseConverter<T> cv(LatSet::cs2); cv.SimplifiedConverterFromRT(8,T{1},T{1.0});
  using FI=TypePack<RHO<T>, VELOCITY<T,2>, FORCE<T,2>, OMEGA<T>, POP<T,LatSet::q>>;
  ValuePack ivs(T{},Vector<T,2>{T{0},T{0}},Vector<T,2>{T{0},T{0}},T{0.5},T{});
  BlockLatticeManager<T,LatSet,FI> Lat(geo,ivs,cv);
  using MC=Cell<T,LatSet,FI>;
  auto& bl=Lat.getBlockLat(0); int ol=bl.getBlock().getOverlap();
  const auto& pr=bl.getBlock().getProjection();
  MC cell(ol*pr[1]+ol,bl);
  T ps=T{0.5}; Vector<T,2> uu{T{0},T{0}}; Vector<T,2> ff{T{0},T{0}};
  cell.template get<RHO<T>>()=ps;
  cell.template get<VELOCITY<T,2>>()=uu;
  cell.template get<FORCE<T,2>>()=ff;
  cell.template get<OMEGA<T>>()=T{0.5};
  // Init g_eq = p*·w_α + (Γ_α - w_α)   Eq.17
  T u2=uu[0]*uu[0]+uu[1]*uu[1];
  for(unsigned int k=0;k<LatSet::q;++k){
    T cu=latset::c<LatSet>(k)[0]*uu[0]+latset::c<LatSet>(k)[1]*uu[1];
    T Gamma=latset::w<LatSet>(k)*(T{1}+cu*LatSet::InvCs2+cu*cu*LatSet::InvCs4*T{0.5}-u2*LatSet::InvCs2*T{0.5});
    cell[k]=ps*latset::w<LatSet>(k)+(Gamma-latset::w<LatSet>(k));
  }
  T old[LatSet::q]; for(unsigned int k=0;k<LatSet::q;++k) old[k]=cell[k];
  phase_field::G_Collision<MC>::apply(cell);
  int pass=0,fail=0; T eps=T{1e-12};
  for(unsigned int k=0;k<LatSet::q;++k){
    if(std::abs(cell[k]-old[k])<=eps) pass++; else fail++;
  }
  if(mpi().isMainProcessor())
    std::cout<<"[g-collision-eq] "<<pass<<" passed, "<<fail<<" failed"<<std::endl;
  return fail;
}

int main(int argc,char** argv){
  mpi().init(&argc,&argv);
  int total=0;
  total+=testDensity();
  total+=testChemPot();
  total+=testMacroUpdate();
  total+=testGCollisionEq();
  if(mpi().isMainProcessor()){
    if(total==0) std::cout<<"\n=== Phase C: ALL TESTS PASSED ==="<<std::endl;
    else std::cout<<"\n=== Phase C: "<<total<<" FAILURES ==="<<std::endl;
  }
  return total;
}
