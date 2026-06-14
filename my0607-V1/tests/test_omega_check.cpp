#include "freelb.h"
#include "freelb.hh"
#include "lbm/h_collision.ur.h"
using T=double; using LatSet=D2Q9<T>;
int main(int argc,char** argv){
  mpi().init(&argc,&argv);
  constexpr int Mn=8;
  AABB<T,2> dm(Vector<T,2>(T{0},T{0}),Vector<T,2>(T(Mn),T(Mn)));
  BlockGeometryHelper2D<T> gh(Mn,Mn,dm,T{1},Mn);
  gh.CreateBlocks(1,mpi().getSize());gh.AdaptiveOptimization(mpi().getSize());gh.LoadBalancing(mpi().getSize());
  BlockGeometry2D<T> geo(gh);
  BlockFieldManager<FLAG,T,2> fm(geo,1);fm.forEach(dm,[&](FLAG& f,std::size_t id){f.SetField(id,2);});fm.template SetupBoundary<LatSet>(dm,2);
  BaseConverter<T> cv(LatSet::cs2);cv.Converter(T{1},T{1},T{1},T(Mn),T{1},T{0.1});
  using PF=TypePack<PHI<T>,GRAD<T,2>,NORMAL<T,2>,VELOCITY<T,2>,INTERFACEWIDTH<T>,POP<T,LatSet::q>>;
  ValuePack ivs(T{},Vector<T,2>{T{0},T{0}},Vector<T,2>{T{0},T{0}},Vector<T,2>{T{0},T{0}},T{4},T{});
  BlockLatticeManager<T,LatSet,PF> Lat(geo,ivs,cv);
  auto& bl=Lat.getBlockLat(0);const auto& b=bl.getBlock();int ol=b.getOverlap();
  using MC=Cell<T,LatSet,PF>;
  // Setup cell
  MC cell(ol*b.getProjection()[1]+ol,bl);
  cell.template get<PHI<T>>()=T{0.5};
  cell.template get<VELOCITY<T,2>>()=Vector<T,2>{T{0},T{0}};
  cell.template get<NORMAL<T,2>>()=Vector<T,2>{T{0},T{0}};
  cell.template get<INTERFACEWIDTH<T>>()=T{4};
  // Check omega
  std::cout<<"Omega="<<cell.getOmega()<<" _Omega="<<cell.get_Omega()<<" fOmega="<<cell.getfOmega()<<std::endl;
  std::cout<<"Sum="<<cell.getOmega()+cell.get_Omega()<<std::endl;
  // Init h
  for(unsigned int k=0;k<LatSet::q;++k) cell[k]=latset::w<LatSet>(k)*T{0.5};
  T h_before[LatSet::q];for(unsigned int k=0;k<LatSet::q;++k)h_before[k]=cell[k];
  // H_Collision
  phase_field::H_Collision<MC>::apply(cell);
  // Print differences
  for(unsigned int k=0;k<LatSet::q;++k){
    T h_eq=cell.template get<PHI<T>>()*latset::w<LatSet>(k);
    T expected=cell.getOmega()*h_eq + cell.get_Omega()*h_before[k];
    std::cout<<"  k="<<k<<" w="<<latset::w<LatSet>(k)
             <<" before="<<h_before[k]<<" after="<<cell[k]
             <<" h_eq="<<h_eq<<" exp="<<expected<<" diff="<<(cell[k]-expected)<<std::endl;
  }
  return 0;
}
