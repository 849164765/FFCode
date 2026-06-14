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
  { std::cout<<"=== NORMAL test ==="<<std::endl;
    using PF=TypePack<NORMAL<T,2>>;
    ValuePack ivs(Vector<T,2>{T{0},T{0}});
    BlockLatticeManager<T,LatSet,PF> Lat(geo,ivs,cv);
    auto& bl=Lat.getBlockLat(0); const auto& b=bl.getBlock();
    std::cout<<"  Nx="<<b.getNx()<<" Ny="<<b.getNy()<<std::endl;
    std::cout<<"  getField..."<<std::endl;
    auto& nf = Lat.template getField<NORMAL<T,2>>();
    std::cout<<"  getBlockField(0)..."<<std::endl;
    auto& bf = nf.getBlockField(0);
    std::cout<<"  getFieldPtr..."<<std::endl;
    auto* ptr = bf.getField(0).getdataPtr();
    std::cout<<"  ptr="<<(void*)ptr<<" size="<<bf.getField(0).getSize()<<std::endl;
    std::cout<<"  writing..."<<std::endl;
    const auto& p=b.getProjection();
    for(int j=0;j<b.getNy();++j) for(int i=0;i<b.getNx();++i) {
      std::size_t id=j*p[1]+i;
      bf.get(0,id)=Vector<T,2>{T{0},T{0}}; // or: bf.get(0,id)[0]=0; bf.get(0,id)[1]=0;
    }
    std::cout<<"  OK"<<std::endl;
  }
  std::cout<<"DONE"<<std::endl;
  return 0;
}
