#include "freelb.h"
#include "freelb.hh"
using T = double;
using LatSet = D2Q9<T>;
int main(int argc, char** argv) {
  mpi().init(&argc, &argv);
  constexpr int Mn=8; T cl=T{1};
  AABB<T,2> dm(Vector<T,2>(T{0},T{0}), Vector<T,2>(T(Mn*cl),T(Mn*cl)));
  BlockGeometryHelper2D<T> gh(Mn,Mn,dm,cl,Mn);
  gh.CreateBlocks(1,mpi().getSize()); gh.AdaptiveOptimization(mpi().getSize()); gh.LoadBalancing(mpi().getSize());
  BlockGeometry2D<T> geo(gh);
  BlockFieldManager<FLAG,T,2> fm(geo,1);
  fm.forEach(dm,[&](FLAG& f,std::size_t id){f.SetField(id,2);});
  fm.template SetupBoundary<LatSet>(dm,2);
  BaseConverter<T> cv(LatSet::cs2); cv.Converter(T{1},T{1},T{1},T(Mn),T{1},T{0.1});
  using PF=TypePack<PHI<T>,GRAD<T,2>,NORMAL<T,2>,VELOCITY<T,2>,INTERFACEWIDTH<T>,POP<T,LatSet::q>>;
  ValuePack ivs(T{},Vector<T,2>{T{0},T{0}},Vector<T,2>{T{0},T{0}},Vector<T,2>{T{0},T{0}},T{4},T{});
  BlockLatticeManager<T,LatSet,PF> Lat(geo,ivs,cv);
  auto& bl=Lat.getBlockLat(0); const auto& b=bl.getBlock(); int ol=b.getOverlap(); const auto& p=b.getProjection();
  std::cout<<"Lattice ready. Nx="<<b.getNx()<<std::endl;
  using MC=Cell<T,LatSet,PF>;
  // Use Cell-based access
  for(int j=0;j<b.getNy();++j) for(int i=0;i<b.getNx();++i) {
    MC cell(j*p[1]+i,bl);
    cell.template get<PHI<T>>()=T{0.5};
    cell.template get<VELOCITY<T,2>>()=Vector<T,2>{T{0},T{0}};
    cell.template get<NORMAL<T,2>>()=Vector<T,2>{T{0},T{0}};
    cell.template get<INTERFACEWIDTH<T>>()=T{4};
  }
  std::cout<<"Init via Cell OK"<<std::endl;
  // Init h to equilibrium
  MC cell(ol*p[1]+ol,bl);
  for(unsigned int k=0;k<LatSet::q;++k) cell[k]=latset::w<LatSet>(k)*T{0.5};
  std::cout<<"POP init OK"<<std::endl;
  // Collision test
  phase_field::H_Collision<MC>::apply(cell);
  std::cout<<"Collision OK"<<std::endl;
  return 0;
}
