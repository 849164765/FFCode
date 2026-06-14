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

  auto doTest = [&](const char* name){
    std::cout << "=== " << name << " ===" << std::endl;
  };
  
  doTest("1: PHI+VEL+IWIDTH+POP");
  { using PF = TypePack<PHI<T>,VELOCITY<T,2>,INTERFACEWIDTH<T>,POP<T,LatSet::q>>;
    ValuePack ivs(T{},Vector<T,2>{T{0},T{0}},T{4},T{});
    BlockLatticeManager<T,LatSet,PF> Lat(geo,ivs,cv);
    auto& bl=Lat.getBlockLat(0); const auto& b=bl.getBlock(); const auto& p=b.getProjection();
    std::cout << "  fill fields..." << std::endl;
    for(int j=0;j<b.getNy();++j) for(int i=0;i<b.getNx();++i) {
      std::size_t id=j*p[1]+i;
      auto& bf=Lat.getField<INTERFACEWIDTH<T>>().getBlockField(0); bf.get(id)=T{4};
    }
    std::cout << "  OK" << std::endl;
  }
  
  doTest("2: PHI+GRAD+NORMAL+POP");
  { using PF = TypePack<PHI<T>,GRAD<T,2>,NORMAL<T,2>,POP<T,LatSet::q>>;
    ValuePack ivs(T{},Vector<T,2>{T{0},T{0}},Vector<T,2>{T{0},T{0}},T{});
    BlockLatticeManager<T,LatSet,PF> Lat(geo,ivs,cv);
    auto& bl=Lat.getBlockLat(0); const auto& b=bl.getBlock(); const auto& p=b.getProjection();
    std::cout << "  fill fields..." << std::endl;
    for(int j=0;j<b.getNy();++j) for(int i=0;i<b.getNx();++i) {
      std::size_t id=j*p[1]+i;
      Lat.getField<NORMAL<T,2>>().getBlockField(0).get(0,id)=Vector<T,2>{T{0},T{0}};
    }
    std::cout << "  OK" << std::endl;
  }
  
  std::cout << "DONE" << std::endl;
  return 0;
}
