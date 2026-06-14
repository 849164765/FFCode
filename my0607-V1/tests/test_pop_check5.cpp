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
  { std::cout<<"=== NORMAL only ==="<<std::endl;
    using PF=TypePack<NORMAL<T,2>>;
    ValuePack ivs(Vector<T,2>{T{0},T{0}});
    BlockLatticeManager<T,LatSet,PF> Lat(geo,ivs,cv);
    auto& bl=Lat.getBlockLat(0); const auto& b=bl.getBlock(); const auto& p=b.getProjection();
    std::cout<<"  Nx="<<b.getNx()<<" Ny="<<b.getNy()<<std::endl;
    for(int j=0;j<b.getNy();++j) for(int i=0;i<b.getNx();++i)
      Lat.getField<NORMAL<T,2>>().getBlockField(0).get(0,j*p[1]+i)=Vector<T,2>{T{0},T{0}};
    std::cout<<"  NORMAL write OK"<<std::endl;
  }
  { std::cout<<"=== POP only ==="<<std::endl;
    using PF=TypePack<POP<T,LatSet::q>>;
    ValuePack ivs(T{});
    BlockLatticeManager<T,LatSet,PF> Lat(geo,ivs,cv);
    auto& bl=Lat.getBlockLat(0); const auto& b=bl.getBlock();
    std::cout<<"  Nx="<<b.getNx()<<" Ny="<<b.getNy()<<std::endl;
    using MC=Cell<T,LatSet,PF>;
    int ol=b.getOverlap(); const auto& p=b.getProjection();
    MC cell(ol*p[1]+ol,bl);
    std::cout<<"  cell created, writing pops..."<<std::endl;
    for(unsigned int k=0;k<LatSet::q;++k) cell[k]=T{0.5};
    std::cout<<"  POP write OK"<<std::endl;
  }
  std::cout<<"DONE"<<std::endl;
  return 0;
}
