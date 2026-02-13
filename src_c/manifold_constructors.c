// Copyright 2021 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// Primitive constructors for the C11 Manifold port.

#include "manifold_impl.h"
#include "manifold_polygon.h"
#include <math.h>

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

static inline double cosd(double deg) { return cos(deg * MANIFOLD_PI / 180.0); }
static inline double sind(double deg) { return sin(deg * MANIFOLD_PI / 180.0); }

void manifold_impl_extrude(ManifoldImpl *impl,
                           const ManifoldVec2 *polyVerts,
                           const int *polySizes, int nPolys,
                           double height, int nDivisions,
                           double twistDegrees, ManifoldVec2 scaleTop) {
  manifold_impl_init(impl);
  if (nPolys <= 0 || height <= 0.0) return;

  if (scaleTop.x < 0.0) scaleTop.x = 0.0;
  if (scaleTop.y < 0.0) scaleTop.y = 0.0;
  nDivisions++;

  bool isCone = (scaleTop.x == 0.0 && scaleTop.y == 0.0);

  // Count total cross-section vertices
  int nCrossSection = 0;
  for (int i = 0; i < nPolys; i++) nCrossSection += polySizes[i];

  ManifoldVecVec3 vertPos = {0};
  ManifoldVecIVec3 triVerts = {0};

  // Add bottom layer vertices
  int offset = 0;
  for (int p = 0; p < nPolys; p++) {
    for (int v = 0; v < polySizes[p]; v++) {
      vec_vec3_push(&vertPos, manifold_vec3(polyVerts[offset + v].x,
                                             polyVerts[offset + v].y, 0.0));
    }
    offset += polySizes[p];
  }

  // Add intermediate and top layers
  for (int i = 1; i <= nDivisions; i++) {
    double alpha = (double)i / (double)nDivisions;
    double phi = alpha * twistDegrees;
    double sx = 1.0 + alpha * (scaleTop.x - 1.0);
    double sy = 1.0 + alpha * (scaleTop.y - 1.0);
    double cosP = cosd(phi), sinP = sind(phi);

    offset = 0;
    int j = 0;
    for (int p = 0; p < nPolys; p++) {
      for (int v = 0; v < polySizes[p]; v++) {
        size_t thisOfs = (size_t)(offset + nCrossSection * i);
        size_t thisVert = (size_t)v + thisOfs;
        size_t lastVert = (v == 0 ? (size_t)polySizes[p] : (size_t)v) - 1 + thisOfs;

        if (i == nDivisions && isCone) {
          vec_ivec3_push(&triVerts, manifold_ivec3(
              nCrossSection * i + j,
              (int)(lastVert - (size_t)nCrossSection),
              (int)(thisVert - (size_t)nCrossSection)));
        } else {
          double px = polyVerts[offset + v].x;
          double py = polyVerts[offset + v].y;
          double tx = sx * px * cosP - sy * py * sinP;
          double ty = sx * px * sinP + sy * py * cosP;
          vec_vec3_push(&vertPos, manifold_vec3(tx, ty, height * alpha));
          vec_ivec3_push(&triVerts, manifold_ivec3(
              (int)thisVert, (int)lastVert, (int)(thisVert - (size_t)nCrossSection)));
          vec_ivec3_push(&triVerts, manifold_ivec3(
              (int)lastVert, (int)(lastVert - (size_t)nCrossSection),
              (int)(thisVert - (size_t)nCrossSection)));
        }
      }
      j++;
      offset += polySizes[p];
    }
  }

  // Add cone apex vertices if needed
  if (isCone) {
    for (int j = 0; j < nPolys; j++) {
      vec_vec3_push(&vertPos, manifold_vec3(0.0, 0.0, height));
    }
  }

  // Triangulate bottom face (using first polygon only for simplicity)
  // Use ear-clipping for each polygon
  offset = 0;
  for (int p = 0; p < nPolys; p++) {
    ManifoldVecIVec3 bottom = manifold_triangulate_polygon(
        polyVerts + offset, NULL, (size_t)polySizes[p]);
    for (size_t t = 0; t < bottom.len; t++) {
      // Bottom: reverse winding
      vec_ivec3_push(&triVerts, manifold_ivec3(
          bottom.data[t].x + offset,
          bottom.data[t].z + offset,
          bottom.data[t].y + offset));
      // Top: normal winding
      if (!isCone) {
        vec_ivec3_push(&triVerts, manifold_ivec3(
            bottom.data[t].x + offset + nCrossSection * nDivisions,
            bottom.data[t].y + offset + nCrossSection * nDivisions,
            bottom.data[t].z + offset + nCrossSection * nDivisions));
      }
    }
    vec_ivec3_free(&bottom);
    offset += polySizes[p];
  }

  impl->vertPos = vec_vec3_create_n(vertPos.len);
  for (size_t i = 0; i < vertPos.len; i++) {
    impl->vertPos.data[i] = vertPos.data[i];
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

// Forward declaration
int manifold_get_circular_segments(double radius);

void manifold_impl_revolve(ManifoldImpl *impl,
                           const ManifoldVec2 *polyVerts,
                           const int *polySizes, int nPolys,
                           int circularSegments, double revolveDegrees) {
  manifold_impl_init(impl);
  if (nPolys <= 0) return;

  double radius = 0;
  int offset = 0;
  for (int p = 0; p < nPolys; p++) {
    for (int v = 0; v < polySizes[p]; v++) {
      if (polyVerts[offset + v].x > radius)
        radius = polyVerts[offset + v].x;
    }
    offset += polySizes[p];
  }
  if (radius <= 0.0) return;

  if (revolveDegrees > 360.0) revolveDegrees = 360.0;
  bool isFullRevolution = (revolveDegrees == 360.0);

  int nDivisions = circularSegments > 2 ? circularSegments :
                   (int)(manifold_get_circular_segments(radius) * revolveDegrees / 360.0);
  if (nDivisions < 3) nDivisions = 3;

  int nSlices = isFullRevolution ? nDivisions : nDivisions + 1;
  double dPhi = revolveDegrees / nDivisions;

  ManifoldVecVec3 vertPos = {0};
  ManifoldVecIVec3 triVerts = {0};

  offset = 0;
  for (int p = 0; p < nPolys; p++) {
    int nPoly = polySizes[p];
    int nPosVerts = 0, nAxisVerts = 0;
    for (int v = 0; v < nPoly; v++) {
      if (polyVerts[offset + v].x > 0) nPosVerts++;
      else nAxisVerts++;
    }

    for (int polyVert = 0; polyVert < nPoly; polyVert++) {
      int startPosIndex = (int)vertPos.len;
      ManifoldVec2 curr = polyVerts[offset + polyVert];
      ManifoldVec2 prev = polyVerts[offset + (polyVert == 0 ? nPoly - 1 : polyVert - 1)];

      int prevStartPosIndex = startPosIndex +
          (polyVert == 0 ? nAxisVerts + (nSlices * nPosVerts) : 0) +
          (prev.x == 0.0 ? -1 : -nSlices);

      for (int slice = 0; slice < nSlices; slice++) {
        double phi = slice * dPhi;
        if (slice == 0 || curr.x > 0) {
          vec_vec3_push(&vertPos, manifold_vec3(
              curr.x * cosd(phi), curr.x * sind(phi), curr.y));
        }

        if (isFullRevolution || slice > 0) {
          int lastSlice = (slice == 0 ? nDivisions : slice) - 1;
          if (curr.x > 0.0) {
            int psi = (prev.x == 0.0) ? prevStartPosIndex : prevStartPosIndex + lastSlice;
            vec_ivec3_push(&triVerts, manifold_ivec3(
                startPosIndex + slice, startPosIndex + lastSlice, psi));
          }
          if (prev.x > 0.0) {
            int csi = (curr.x == 0.0) ? startPosIndex : startPosIndex + slice;
            vec_ivec3_push(&triVerts, manifold_ivec3(
                prevStartPosIndex + lastSlice, prevStartPosIndex + slice, csi));
          }
        }
      }
    }
    offset += nPoly;
  }

  impl->vertPos = vec_vec3_create_n(vertPos.len);
  for (size_t i = 0; i < vertPos.len; i++) {
    impl->vertPos.data[i] = vertPos.data[i];
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
