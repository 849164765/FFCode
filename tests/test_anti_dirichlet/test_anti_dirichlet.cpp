// test_anti_dirichlet.cpp — 1000-step Dirichlet BC test (D2Q5)
// Uses "full-way" equilibrium BC at walls: set wall populations to feq at psi_bc
// Verifies 1D diffusion ψ(y)=y/H
#include "freelb.h"
#include "freelb.hh"

using T = double;
using MFLatSet = D2Q5<T>;
using namespace mfield;

constexpr int N = 16;
constexpr T Cell_Len = 1.0;
constexpr int BlockCellLen = 16;
constexpr T omega_test = T{1.0};
constexpr int TotalSteps = 1000;
constexpr T tol = 0.06;  // full-way BC has O(Δx) error ≈ 1/N

int main(int argc, char** argv) {
  mpi().init(&argc, &argv);

  constexpr std::uint8_t VoidFlag = std::uint8_t(1);
  constexpr std::uint8_t BulkFlag = std::uint8_t(2);

  AABB<T,2> domain(Vector<T,2>{0,0}, Vector<T,2>{T(N*Cell_Len),T(N*Cell_Len)});
  BlockGeometryHelper2D<T> GeoHelper(N,N,domain,Cell_Len,BlockCellLen);
  GeoHelper.CreateBlocks(1,1);
  GeoHelper.AdaptiveOptimization(1);
  GeoHelper.LoadBalancing(1);
  BlockGeometry2D<T> Geo(GeoHelper);

  BlockFieldManager<FLAG,T,2> FlagFM(Geo,VoidFlag);
  FlagFM.forEach(domain,[&](FLAG&f,std::size_t id){f.SetField(id,BulkFlag);});

  BaseConverter<T> MFBaseConv(MFLatSet::cs2);
  MFBaseConv.SimplifiedConverterFromRT(N,T(0.01),T{1.0}/omega_test);

  using MFF = TypePack<PSI<T>,OMEGA_PSI<T>,POP<T,MFLatSet::q>>;
  ValuePack MFInit(T{0.5},omega_test,T{});
  using MFCELL = Cell<T,MFLatSet,MFF>;
  BlockLatticeManager<T,MFLatSet,MFF> MFLattice(Geo,MFInit,MFBaseConv);

  // init: psi=0.5, g_k=w_k*psi everywhere
  MFLattice.getField<OMEGA_PSI<T>>().InitValue(omega_test);
  {
    auto& psiF = MFLattice.getField<PSI<T>>();
    for(int b=0;b<Geo.getBlockNum();++b){
      auto& bl = MFLattice.getBlockLat(b);
      const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
      auto& bPsi=psiF.getBlockField(b);
      for(int j=0;j<bk.getNy();++j){
        for(int i=0;i<bk.getNx();++i){
          std::size_t id=j*pr[1]+i; MFCELL cell(id,bl);
          for(unsigned k=0;k<MFLatSet::q;++k) cell[k]=latset::w<MFLatSet>(k)*T{0.5};
          bPsi.get(id)=T{0.5};
        }
      }
    }
  }

  MFLattice.NormalFullCommunicate();

  Timer timer;
  T H_global = T(N)*Cell_Len;
  for(int step=0;step<TotalSteps;++step){
    // a) Full-way wall BC: set wall cell populations = feq(psi_bc)
    {
      auto& psiF = MFLattice.getField<PSI<T>>();
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& bl = MFLattice.getBlockLat(b);
        const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bPsi=psiF.getBlockField(b);
        int nx=bk.getNx(),ny=bk.getNy(),ov=bk.getOverlap();
        T minY=bk.getMin()[1],maxY=bk.getMax()[1],vs=bk.getVoxelSize();

        auto apply_wall_row = [&](int jj, T psi_bc) {
          for(int ii=0;ii<nx;++ii){
            std::size_t id=jj*pr[1]+ii;
            MFCELL cell(id,bl);
            for(unsigned k=0;k<MFLatSet::q;++k)
              cell[k]=latset::w<MFLatSet>(k)*psi_bc;
            bPsi.get(id)=psi_bc;
          }
        };

        if(minY<Cell_Len*T{1.5})
          for(int jj=0;jj<=ov;++jj) apply_wall_row(jj, T{0.0});
        if(maxY>H_global-Cell_Len*T{1.5})
          for(int jj=ny-ov-1;jj<ny;++jj) apply_wall_row(jj, T{1.0});

        // left wall: psi = y/H
        for(int jj=0;jj<ny;++jj){
          T y=minY+T(jj)*vs, ev=T{1.0}*y/H_global;
          for(int ii=0;ii<=ov;++ii){
            std::size_t id=jj*pr[1]+ii; MFCELL cell(id,bl);
            for(unsigned k=0;k<MFLatSet::q;++k) cell[k]=latset::w<MFLatSet>(k)*ev;
            bPsi.get(id)=ev;
          }
        }
        // right wall
        for(int jj=0;jj<ny;++jj){
          T y=minY+T(jj)*vs, ev=T{1.0}*y/H_global;
          for(int ii=nx-ov-1;ii<nx;++ii){
            std::size_t id=jj*pr[1]+ii; MFCELL cell(id,bl);
            for(unsigned k=0;k<MFLatSet::q;++k) cell[k]=latset::w<MFLatSet>(k)*ev;
            bPsi.get(id)=ev;
          }
        }
      }
    }

    // b) Collision on interior cells
    for(int b=0;b<Geo.getBlockNum();++b){
      auto& bl=MFLattice.getBlockLat(b);
      const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
      int ov=bk.getOverlap();
      for(int j=ov;j<bk.getNy()-ov;++j){
        for(int i=ov;i<bk.getNx()-ov;++i){
          MFCELL cell(j*pr[1]+i,bl);
          collision::MRTDiffusion<MFCELL,OMEGA_PSI<T>>::apply(cell);
        }
      }
    }

    // c) Stream
    MFLattice.NormalFullCommunicate();
    MFLattice.Stream();
    MFLattice.NormalFullCommunicate();

    // d) psi = sum(g_i) for interior
    {
      auto& psiF=MFLattice.getField<PSI<T>>();
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& bl=MFLattice.getBlockLat(b);
        const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bPsi=psiF.getBlockField(b);
        int ov=bk.getOverlap();
        for(int j=ov;j<bk.getNy()-ov;++j){
          for(int i=ov;i<bk.getNx()-ov;++i){
            std::size_t id=j*pr[1]+i; MFCELL cell(id,bl);
            T psi=T{0};
            for(unsigned k=0;k<MFLatSet::q;++k)psi+=cell[k];
            bPsi.get(id)=psi;
          }
        }
      }
    }
    ++timer;
  }

  // verify
  T max_err=0; int count=0;
  {
    auto& psiF=MFLattice.getField<PSI<T>>();
    const auto& bk=Geo.getBlock(0); const auto& pr=bk.getProjection();
    auto& bPsi=psiF.getBlockField(0);
    int ov=bk.getOverlap();
    T vs=bk.getVoxelSize(),minY=bk.getMin()[1];
    printf("psi profile (center column):\n");
    for(int j=ov;j<bk.getNy()-ov;++j){
      T y=minY+T(j)*vs,ev=T{1.0}*y/H_global;
      T val=bPsi.get(j*pr[1]+ov+(bk.getNx())/2-ov);
      T err=std::abs(val-ev);
      printf("  y=%.1f  psi=%.6f  expected=%.6f  err=%.6f\n",y,val,ev,err);
      if(err>max_err)max_err=err; ++count;
    }
  }

  bool pass=max_err<tol;
  printf("Step %d: max|psi-y/H|=%.6f %s\n",TotalSteps,max_err,pass?"PASS":"FAIL");
  return pass?0:1;
}
