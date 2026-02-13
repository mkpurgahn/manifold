// Copyright 2021 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// Properties computation for the C11 Manifold port.

#include "manifold_impl.h"

double manifold_impl_get_volume(const ManifoldImpl *impl) {
  double vol = 0.0;
  size_t numTri = manifold_impl_num_tri(impl);
  for (size_t tri = 0; tri < numTri; tri++) {
    ManifoldVec3 a = impl->vertPos.data[impl->halfedge.data[3 * tri].startVert];
    ManifoldVec3 b = impl->vertPos.data[impl->halfedge.data[3 * tri + 1].startVert];
    ManifoldVec3 c = impl->vertPos.data[impl->halfedge.data[3 * tri + 2].startVert];
    vol += vec3_dot(a, vec3_cross(b, c));
  }
  return vol / 6.0;
}

double manifold_impl_get_surface_area(const ManifoldImpl *impl) {
  double area = 0.0;
  size_t numTri = manifold_impl_num_tri(impl);
  for (size_t tri = 0; tri < numTri; tri++) {
    ManifoldVec3 a = impl->vertPos.data[impl->halfedge.data[3 * tri].startVert];
    ManifoldVec3 b = impl->vertPos.data[impl->halfedge.data[3 * tri + 1].startVert];
    ManifoldVec3 c = impl->vertPos.data[impl->halfedge.data[3 * tri + 2].startVert];
    ManifoldVec3 cross = vec3_cross(vec3_sub(b, a), vec3_sub(c, a));
    area += vec3_length(cross);
  }
  return area / 2.0;
}
