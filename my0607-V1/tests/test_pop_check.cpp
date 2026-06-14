#include "freelb.h"
#include "freelb.hh"
using T = double;
using LatSet = D2Q9<T>;
int main(int argc, char** argv) {
  mpi().init(&argc, &argv);
  constexpr int Mn=8; T cl=T{1};
  AABB<T,2> dm(Vector<T,2>(T{0},T{0}), Vector<T,2>(T(Mn*cl),T(Mn*cl)));
  BlockGeometryHelper2D<T> gh(Mn,Mn,dm,cl,Mn);
  gh.CreateBlocks(1,mpi().getSize());
  gh.AdaptiveOptimization(mpi().getSize());
  gh.LoadBalancing(mpi().getSize());
  BlockGeometry2D<T> geo(gh);
  BlockFieldManager<FLAG,T,2> fm(geo,1);
  fm.forEach(dm,[&](FLAG& f,std::size_t id){f.SetField(id,2);});
  fm.template SetupBoundary<LatSet>(dm,2);
  BaseConverter<T> cv(LatSet::cs2);
  cv.Converter(T{1},T{1},T{1},T(Mn),T{1},T{0.1});
  std::cout << "=== PHI only ===" << std::endl;
  { using PF=TypePack<PHI<T>>; ValuePack ivs(T{}); BlockLatticeManager<T,LatSet,PF> Lat(geo,ivs,cv); std::cout << "OK" << std::endl; }
  std::cout << "=== PHI+GRAD+NORMAL ===" << std::endl;
  { using PF=TypePack<PHI<T>,GRAD<T,2>,NORMAL<T,2>>; ValuePack ivs(T{},Vector<T,2>{T{0},T{0}},Vector<T,2>{T{0},T{0}}); BlockLatticeManager<T,LatSet,PF> Lat(geo,ivs,cv); std::cout << "OK" << std::endl; }
  std::cout << "=== PHI+POP ===" << std::endl;
  { using PF=TypePack<PHI<T>,POP<T,LatSet::q>>; ValuePack ivs(T{},T{}); BlockLatticeManager<T,LatSet,PF> Lat(geo,ivs,cv); std::cout << "OK" << std::endl; }
  std::cout << "=== ALL ===" << std::endl;
  { using PF=TypePack<PHI<T>,GRAD<T,2>,VELOCITY<T,2>,POP<T,LatSet::q>>; ValuePack ivs(T{},Vector<T,2>{T{0},T{0}},Vector<T,2>{T{0},T{0}},T{}); BlockLatticeManager<T,LatSet,PF> Lat(geo,ivs,cv); std::cout << "OK" << std::endl; }
  std::cout << "DONE" << std::endl;
  return 0;
}
