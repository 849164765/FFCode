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
  std::cout << "=== PHI+GRAD+NORMAL+VEL+IWIDTH+POP ===" << std::endl;
  { using PF = TypePack<PHI<T>,GRAD<T,2>,NORMAL<T,2>,VELOCITY<T,2>,INTERFACEWIDTH<T>,POP<T,LatSet::q>>;
    ValuePack ivs(T{},Vector<T,2>{T{0},T{0}},Vector<T,2>{T{0},T{0}},Vector<T,2>{T{0},T{0}},T{4},T{});
    BlockLatticeManager<T,LatSet,PF> Lat(geo,ivs,cv);
    auto& bl = Lat.getBlockLat(0);
    std::cout << "  getBlockLat OK" << std::endl;
    const auto& b = bl.getBlock();
    std::cout << "  getBlock OK Nx=" << b.getNx() << " overlap=" << b.getOverlap() << std::endl;
    int ol=b.getOverlap();
    const auto& p=b.getProjection();
    for (int j=0;j<b.getNy();++j)
      for (int i=0;i<b.getNx();++i) {
        std::size_t id=j*p[1]+i;
        Lat.getField<PHI<T>>().getBlockField(0).get(id)=T{0.5};
        Lat.getField<VELOCITY<T,2>>().getBlockField(0).get(0,id)=Vector<T,2>{T{0},T{0}};
        Lat.getField<NORMAL<T,2>>().getBlockField(0).get(0,id)=Vector<T,2>{T{0},T{0}};
        Lat.getField<INTERFACEWIDTH<T>>().getBlockField(0).get(id)=T{4};
      }
    std::cout << "  Init fields OK" << std::endl;
    using MC = Cell<T,LatSet,PF>;
    int ci=ol, cj=ol;
    std::size_t cid=cj*p[1]+ci;
    MC cell(cid,bl);
    std::cout << "  Cell created OK" << std::endl;
    for (unsigned int k=0;k<LatSet::q;++k) cell[k]=latset::w<LatSet>(k)*T{0.5};
    std::cout << "  POP init OK" << std::endl;
  }
  std::cout << "DONE" << std::endl;
  return 0;
}
