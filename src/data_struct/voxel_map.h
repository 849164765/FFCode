// voxel_map.h
// map between valid voxel and cell of AABB
// this may be useful for complex geometry

#pragma once

#include <cstddef>
#include "utils/util.h"

#ifdef _VOX_ENABLED

class VoxelMap {
 private:
  // voxel number
  std::size_t VoxNum;

  // voxel to AABB cell: size = VoxNum
  std::size_t* Vox_to_AABB;

 public:
  VoxelMap() : VoxNum(0), Vox_to_AABB(nullptr) {}
  VoxelMap(std::size_t voxnum) : 
    VoxNum(voxnum), Vox_to_AABB(new std::size_t[voxnum]{}) {}

  ~VoxelMap(){ if (Vox_to_AABB) delete[] Vox_to_AABB; }

  void Init(std::size_t voxnum) {
    VoxNum = voxnum;
    Vox_to_AABB = new std::size_t[VoxNum]{};
  }

  // get idx without boundary check
  std::size_t operator[](std::size_t VoxId) const { return Vox_to_AABB[VoxId]; }

  std::size_t* getVox_to_AABB() { return Vox_to_AABB; }

  std::size_t getVoxNum() const { return VoxNum; }

};

#endif