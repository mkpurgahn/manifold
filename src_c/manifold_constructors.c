// Copyright 2021 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// Primitive constructors for the C11 Manifold port.

#include "manifold_impl.h"

void manifold_impl_shape(ManifoldImpl *impl, int shape, ManifoldMat3x4 m) {
  manifold_impl_init(impl);

  ManifoldVecVec3 vertPos = {0};
  ManifoldVecIVec3 triVerts = {0};

  switch (shape) {
    case 0: { // Tetrahedron
      vec_vec3_push(&vertPos, manifold_vec3(-1, -1, 1));
      vec_vec3_push(&vertPos, manifold_vec3(-1, 1, -1));
      vec_vec3_push(&vertPos, manifold_vec3(1, -1, -1));
      vec_vec3_push(&vertPos, manifold_vec3(1, 1, 1));
      vec_ivec3_push(&triVerts, manifold_ivec3(2, 0, 1));
      vec_ivec3_push(&triVerts, manifold_ivec3(0, 3, 1));
      vec_ivec3_push(&triVerts, manifold_ivec3(2, 3, 0));
      vec_ivec3_push(&triVerts, manifold_ivec3(3, 2, 1));
      break;
    }
    case 1: { // Cube
      vec_vec3_push(&vertPos, manifold_vec3(0, 0, 0));
      vec_vec3_push(&vertPos, manifold_vec3(0, 0, 1));
      vec_vec3_push(&vertPos, manifold_vec3(0, 1, 0));
      vec_vec3_push(&vertPos, manifold_vec3(0, 1, 1));
      vec_vec3_push(&vertPos, manifold_vec3(1, 0, 0));
      vec_vec3_push(&vertPos, manifold_vec3(1, 0, 1));
      vec_vec3_push(&vertPos, manifold_vec3(1, 1, 0));
      vec_vec3_push(&vertPos, manifold_vec3(1, 1, 1));
      vec_ivec3_push(&triVerts, manifold_ivec3(1, 0, 4));
      vec_ivec3_push(&triVerts, manifold_ivec3(2, 4, 0));
      vec_ivec3_push(&triVerts, manifold_ivec3(1, 3, 0));
      vec_ivec3_push(&triVerts, manifold_ivec3(3, 1, 5));
      vec_ivec3_push(&triVerts, manifold_ivec3(3, 2, 0));
      vec_ivec3_push(&triVerts, manifold_ivec3(3, 7, 2));
      vec_ivec3_push(&triVerts, manifold_ivec3(5, 4, 6));
      vec_ivec3_push(&triVerts, manifold_ivec3(5, 1, 4));
      vec_ivec3_push(&triVerts, manifold_ivec3(6, 4, 2));
      vec_ivec3_push(&triVerts, manifold_ivec3(7, 6, 2));
      vec_ivec3_push(&triVerts, manifold_ivec3(7, 3, 5));
      vec_ivec3_push(&triVerts, manifold_ivec3(7, 5, 6));
      break;
    }
    case 2: { // Octahedron
      vec_vec3_push(&vertPos, manifold_vec3(1, 0, 0));
      vec_vec3_push(&vertPos, manifold_vec3(-1, 0, 0));
      vec_vec3_push(&vertPos, manifold_vec3(0, 1, 0));
      vec_vec3_push(&vertPos, manifold_vec3(0, -1, 0));
      vec_vec3_push(&vertPos, manifold_vec3(0, 0, 1));
      vec_vec3_push(&vertPos, manifold_vec3(0, 0, -1));
      vec_ivec3_push(&triVerts, manifold_ivec3(0, 2, 4));
      vec_ivec3_push(&triVerts, manifold_ivec3(1, 5, 3));
      vec_ivec3_push(&triVerts, manifold_ivec3(2, 1, 4));
      vec_ivec3_push(&triVerts, manifold_ivec3(3, 5, 0));
      vec_ivec3_push(&triVerts, manifold_ivec3(1, 3, 4));
      vec_ivec3_push(&triVerts, manifold_ivec3(0, 5, 2));
      vec_ivec3_push(&triVerts, manifold_ivec3(3, 0, 4));
      vec_ivec3_push(&triVerts, manifold_ivec3(2, 5, 1));
      break;
    }
  }

  // Apply transform to vertices
  impl->vertPos = vec_vec3_create_n(vertPos.len);
  for (size_t i = 0; i < vertPos.len; i++) {
    impl->vertPos.data[i] = mat3x4_transform_point(m, vertPos.data[i]);
  }

  ManifoldVecIVec3 emptyTriVert = {0};
  manifold_impl_create_halfedges(impl, &triVerts, &emptyTriVert);
  manifold_impl_initialize_original(impl);
  manifold_impl_calculate_bbox(impl);
  manifold_impl_set_epsilon(impl, -1.0, false);
  manifold_impl_sort_geometry(impl);
  manifold_impl_set_normals_and_coplanar(impl);

  vec_vec3_free(&vertPos);
  vec_ivec3_free(&triVerts);
  vec_ivec3_free(&emptyTriVert);
}

void manifold_impl_tetrahedron(ManifoldImpl *impl) {
  manifold_impl_shape(impl, 0, mat3x4_identity());
}

void manifold_impl_cube(ManifoldImpl *impl, ManifoldMat3x4 m) {
  manifold_impl_shape(impl, 1, m);
}

void manifold_impl_octahedron(ManifoldImpl *impl, ManifoldMat3x4 m) {
  manifold_impl_shape(impl, 2, m);
}
