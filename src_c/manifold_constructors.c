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

  // Triangulate bottom and top faces with hole support
  // Use bridge-based hole elimination for proper hole handling
  ManifoldVecIVec3 bottom = manifold_triangulate_with_holes(
      polyVerts, polySizes, nPolys, 0);
  for (size_t t = 0; t < bottom.len; t++) {
    // Bottom: reverse winding
    vec_ivec3_push(&triVerts, manifold_ivec3(
        bottom.data[t].x,
        bottom.data[t].z,
        bottom.data[t].y));
    // Top: normal winding
    if (!isCone) {
      vec_ivec3_push(&triVerts, manifold_ivec3(
          bottom.data[t].x + nCrossSection * nDivisions,
          bottom.data[t].y + nCrossSection * nDivisions,
          bottom.data[t].z + nCrossSection * nDivisions));
    }
  }
  vec_ivec3_free(&bottom);

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

  // Clip polygons to x >= 0 (matching C++ behavior)
  // First pass: count total output vertices
  ManifoldVec2 *clippedVerts = NULL;
  int *clippedSizes = NULL;
  int nClippedPolys = 0;
  size_t totalClipped = 0;
  size_t clippedCap = 0;

  int polyOffset = 0;
  for (int p = 0; p < nPolys; p++) {
    int nPoly = polySizes[p];
    const ManifoldVec2 *poly = polyVerts + polyOffset;

    // Find first vertex with x >= 0
    int startIdx = -1;
    for (int i = 0; i < nPoly; i++) {
      if (poly[i].x >= 0) { startIdx = i; break; }
    }
    if (startIdx < 0) { polyOffset += nPoly; continue; }

    // Clip this polygon
    size_t clipStart = totalClipped;
    int i = startIdx;
    do {
      if (poly[i].x >= 0) {
        // Ensure capacity
        if (totalClipped >= clippedCap) {
          clippedCap = clippedCap == 0 ? 64 : clippedCap * 2;
          clippedVerts = (ManifoldVec2*)realloc(clippedVerts,
              clippedCap * sizeof(ManifoldVec2));
        }
        clippedVerts[totalClipped++] = poly[i];
      }
      int next = (i + 1) % nPoly;
      // If crossing the x=0 boundary, add intersection point
      if ((poly[next].x < 0) != (poly[i].x < 0)) {
        double y = poly[next].y - poly[next].x *
                   (poly[i].y - poly[next].y) /
                   (poly[i].x - poly[next].x);
        if (totalClipped >= clippedCap) {
          clippedCap = clippedCap == 0 ? 64 : clippedCap * 2;
          clippedVerts = (ManifoldVec2*)realloc(clippedVerts,
              clippedCap * sizeof(ManifoldVec2));
        }
        clippedVerts[totalClipped++] = (ManifoldVec2){0.0, y};
      }
      i = next;
    } while (i != startIdx);

    int clipLen = (int)(totalClipped - clipStart);
    if (clipLen >= 3) {
      nClippedPolys++;
      clippedSizes = (int*)realloc(clippedSizes, nClippedPolys * sizeof(int));
      clippedSizes[nClippedPolys - 1] = clipLen;
    } else {
      totalClipped = clipStart; // revert
    }
    polyOffset += nPoly;
  }

  if (nClippedPolys == 0) {
    free(clippedVerts);
    free(clippedSizes);
    return;
  }

  double radius = 0;
  for (size_t i = 0; i < totalClipped; i++) {
    if (clippedVerts[i].x > radius) radius = clippedVerts[i].x;
  }
  if (radius <= 0.0) {
    free(clippedVerts);
    free(clippedSizes);
    return;
  }

  if (revolveDegrees > 360.0) revolveDegrees = 360.0;
  bool isFullRevolution = (revolveDegrees == 360.0);

  int nDivisions = circularSegments > 2 ? circularSegments :
                   (int)(manifold_get_circular_segments(radius) * revolveDegrees / 360.0);
  if (nDivisions < 3) nDivisions = 3;

  int nSlices = isFullRevolution ? nDivisions : nDivisions + 1;
  double dPhi = revolveDegrees / nDivisions;

  ManifoldVecVec3 vertPos = {0};
  ManifoldVecIVec3 triVerts = {0};

  // For partial revolution: track start/end positions
  int *startPoses = NULL;
  int *endPoses = NULL;
  int poseCount = 0;
  if (!isFullRevolution) {
    int totalPolyVerts = 0;
    for (int p = 0; p < nClippedPolys; p++) totalPolyVerts += clippedSizes[p];
    startPoses = (int*)malloc(totalPolyVerts * sizeof(int));
    endPoses = (int*)malloc(totalPolyVerts * sizeof(int));
  }

  int clipOffset = 0;
  for (int p = 0; p < nClippedPolys; p++) {
    int nPoly = clippedSizes[p];
    const ManifoldVec2 *poly = clippedVerts + clipOffset;

    int nPosVerts = 0, nAxisVerts = 0;
    for (int v = 0; v < nPoly; v++) {
      if (poly[v].x > 0) nPosVerts++;
      else nAxisVerts++;
    }

    for (int polyVert = 0; polyVert < nPoly; polyVert++) {
      int startPosIndex = (int)vertPos.len;

      if (!isFullRevolution) startPoses[poseCount] = startPosIndex;

      ManifoldVec2 curr = poly[polyVert];
      ManifoldVec2 prev = poly[polyVert == 0 ? nPoly - 1 : polyVert - 1];

      int prevStartPosIndex =
          startPosIndex +
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
            int psi = (prev.x == 0.0) ? prevStartPosIndex
                                       : prevStartPosIndex + lastSlice;
            vec_ivec3_push(&triVerts, manifold_ivec3(
                startPosIndex + slice, startPosIndex + lastSlice, psi));
          }
          if (prev.x > 0.0) {
            int csi = (curr.x == 0.0) ? startPosIndex
                                       : startPosIndex + slice;
            vec_ivec3_push(&triVerts, manifold_ivec3(
                prevStartPosIndex + lastSlice, prevStartPosIndex + slice, csi));
          }
        }
      }

      if (!isFullRevolution) endPoses[poseCount] = (int)vertPos.len - 1;
      poseCount++;
    }
    clipOffset += nPoly;
  }

  // Add front and back cap triangles if not a full revolution
  if (!isFullRevolution) {
    // Triangulate the clipped cross-section
    ManifoldVecIVec3 capTris = manifold_triangulate_with_holes(
        clippedVerts, clippedSizes, nClippedPolys, 0);
    for (size_t t = 0; t < capTris.len; t++) {
      // Front cap
      vec_ivec3_push(&triVerts, manifold_ivec3(
          startPoses[capTris.data[t].x],
          startPoses[capTris.data[t].y],
          startPoses[capTris.data[t].z]));
      // Back cap (reversed winding)
      vec_ivec3_push(&triVerts, manifold_ivec3(
          endPoses[capTris.data[t].z],
          endPoses[capTris.data[t].y],
          endPoses[capTris.data[t].x]));
    }
    vec_ivec3_free(&capTris);
  }

  free(startPoses);
  free(endPoses);

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
  free(clippedVerts);
  free(clippedSizes);
}
