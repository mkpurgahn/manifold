// Copyright 2021 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// Polygon triangulation for the C11 Manifold port.
// Implements ear-clipping algorithm.

#ifndef MANIFOLD_POLYGON_H
#define MANIFOLD_POLYGON_H

#include "manifold_vec.h"
#include "manifold_vec_math.h"

// A single polygon contour: array of vec2 points
typedef struct {
  ManifoldVec2 *data;
  size_t len;
  size_t cap;
} ManifoldSimplePolygon;

// A polygon with holes: array of contours
typedef struct {
  ManifoldSimplePolygon *data;
  size_t len;
  size_t cap;
} ManifoldPolygons;

// Indexed polygon vertex
typedef struct {
  ManifoldVec2 pos;
  int idx;
} ManifoldPolyVert;

typedef struct {
  ManifoldPolyVert *data;
  size_t len;
  size_t cap;
} ManifoldSimplePolygonIdx;

// Helpers
static inline void simple_polygon_push(ManifoldSimplePolygon *p, ManifoldVec2 v) {
  if (p->len >= p->cap) {
    p->cap = p->cap == 0 ? 16 : p->cap * 2;
    p->data = (ManifoldVec2 *)realloc(p->data, p->cap * sizeof(ManifoldVec2));
  }
  p->data[p->len++] = v;
}

static inline void simple_polygon_free(ManifoldSimplePolygon *p) {
  free(p->data); p->data = NULL; p->len = 0; p->cap = 0;
}

static inline void simple_polygon_idx_push(ManifoldSimplePolygonIdx *p,
                                            ManifoldPolyVert v) {
  if (p->len >= p->cap) {
    p->cap = p->cap == 0 ? 16 : p->cap * 2;
    p->data = (ManifoldPolyVert *)realloc(p->data, p->cap * sizeof(ManifoldPolyVert));
  }
  p->data[p->len++] = v;
}

static inline void simple_polygon_idx_free(ManifoldSimplePolygonIdx *p) {
  free(p->data); p->data = NULL; p->len = 0; p->cap = 0;
}

// Check if point p is inside triangle (a,b,c) using cross products
static inline bool point_in_triangle_2d(ManifoldVec2 p, ManifoldVec2 a,
                                         ManifoldVec2 b, ManifoldVec2 c) {
  double d1 = vec2_cross(vec2_sub(b, a), vec2_sub(p, a));
  double d2 = vec2_cross(vec2_sub(c, b), vec2_sub(p, b));
  double d3 = vec2_cross(vec2_sub(a, c), vec2_sub(p, c));
  bool has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
  bool has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);
  return !(has_neg && has_pos);
}

// Ear-clipping triangulation of a simple polygon
// Returns triangles as ManifoldVecIVec3
static inline ManifoldVecIVec3 manifold_triangulate_polygon(
    const ManifoldVec2 *verts, const int *indices, size_t n) {
  ManifoldVecIVec3 tris = {0};
  if (n < 3) return tris;

  // Make a working copy of indices
  int *idx = (int *)malloc(n * sizeof(int));
  ManifoldVec2 *pts = (ManifoldVec2 *)malloc(n * sizeof(ManifoldVec2));
  for (size_t i = 0; i < n; i++) {
    idx[i] = indices ? indices[i] : (int)i;
    pts[i] = verts[i];
  }

  size_t remaining = n;
  while (remaining > 2) {
    bool found = false;
    for (size_t i = 0; i < remaining; i++) {
      size_t prev = (i + remaining - 1) % remaining;
      size_t next = (i + 1) % remaining;

      ManifoldVec2 a = pts[prev], b = pts[i], c = pts[next];
      double cross = vec2_cross(vec2_sub(b, a), vec2_sub(c, a));

      // Only consider CCW ears
      if (cross <= 0) continue;

      // Check no other vertex is inside this ear
      bool ear = true;
      for (size_t j = 0; j < remaining; j++) {
        if (j == prev || j == i || j == next) continue;
        if (point_in_triangle_2d(pts[j], a, b, c)) {
          ear = false;
          break;
        }
      }

      if (ear) {
        vec_ivec3_push(&tris, manifold_ivec3(idx[prev], idx[i], idx[next]));
        // Remove vertex i
        for (size_t k = i; k < remaining - 1; k++) {
          idx[k] = idx[k + 1];
          pts[k] = pts[k + 1];
        }
        remaining--;
        found = true;
        break;
      }
    }
    if (!found) {
      // Degenerate polygon - force triangulation
      if (remaining >= 3) {
        vec_ivec3_push(&tris, manifold_ivec3(idx[0], idx[1], idx[2]));
        for (size_t k = 1; k < remaining - 1; k++) {
          idx[k] = idx[k + 1];
          pts[k] = pts[k + 1];
        }
        remaining--;
      } else {
        break;
      }
    }
  }

  free(idx);
  free(pts);
  return tris;
}

// Simple fan triangulation for convex polygons
static inline ManifoldVecIVec3 manifold_fan_triangulate(int startIdx, int n) {
  ManifoldVecIVec3 tris = {0};
  for (int i = 1; i < n - 1; i++) {
    vec_ivec3_push(&tris, manifold_ivec3(startIdx, startIdx + i, startIdx + i + 1));
  }
  return tris;
}

// Triangulate a polygon with holes using bridge-based hole elimination.
// polyVerts: flat array of all polygon vertices (outer first, then holes)
// polySizes: array of vertex counts per polygon
// nPolys: number of polygons (1 outer + n-1 holes)
// baseIndex: starting vertex index offset for the output triangles
// Returns triangles as ManifoldVecIVec3
static inline ManifoldVecIVec3 manifold_triangulate_with_holes(
    const ManifoldVec2 *polyVerts, const int *polySizes, int nPolys,
    int baseIndex) {
  if (nPolys <= 1) {
    // No holes, just triangulate the single polygon
    return manifold_triangulate_polygon(polyVerts, NULL,
        nPolys == 1 ? (size_t)polySizes[0] : 0);
  }

  // Build combined polygon by bridging holes into the outer contour
  // Strategy: for each hole, find its rightmost vertex, then find the closest
  // visible vertex on the outer polygon, and insert a bridge (duplicate vertices)

  int outerN = polySizes[0];
  // Start with outer polygon vertices
  size_t totalVerts = 0;
  for (int i = 0; i < nPolys; i++) totalVerts += (size_t)polySizes[i];

  // Working arrays for the merged polygon
  ManifoldVec2 *merged = (ManifoldVec2 *)malloc(
      (totalVerts + 2 * (size_t)(nPolys - 1)) * sizeof(ManifoldVec2));
  int *mergedIdx = (int *)malloc(
      (totalVerts + 2 * (size_t)(nPolys - 1)) * sizeof(int));
  size_t mergedLen = (size_t)outerN;

  // Copy outer polygon
  for (int i = 0; i < outerN; i++) {
    merged[i] = polyVerts[i];
    mergedIdx[i] = baseIndex + i;
  }

  // Process each hole
  int holeOffset = outerN;
  for (int h = 1; h < nPolys; h++) {
    int holeN = polySizes[h];
    const ManifoldVec2 *hole = polyVerts + holeOffset;

    // Find the rightmost vertex in the hole
    int rightmost = 0;
    for (int i = 1; i < holeN; i++) {
      if (hole[i].x > hole[rightmost].x ||
          (hole[i].x == hole[rightmost].x && hole[i].y > hole[rightmost].y)) {
        rightmost = i;
      }
    }

    // Find the closest visible vertex in the merged polygon
    ManifoldVec2 holeVert = hole[rightmost];
    int bestOuter = 0;
    double bestDist = 1e308;
    for (size_t i = 0; i < mergedLen; i++) {
      double dx = merged[i].x - holeVert.x;
      double dy = merged[i].y - holeVert.y;
      double dist = dx * dx + dy * dy;
      if (dist < bestDist) {
        bestDist = dist;
        bestOuter = (int)i;
      }
    }

    // Insert bridge: outer[bestOuter] -> hole[rightmost] -> ... -> hole[rightmost] -> outer[bestOuter]
    size_t newLen = mergedLen + (size_t)holeN + 2;
    ManifoldVec2 *newMerged = (ManifoldVec2 *)malloc(newLen * sizeof(ManifoldVec2));
    int *newIdx = (int *)malloc(newLen * sizeof(int));
    size_t k = 0;

    // Copy outer up to and including bestOuter
    for (int i = 0; i <= bestOuter; i++) {
      newMerged[k] = merged[i];
      newIdx[k] = mergedIdx[i];
      k++;
    }

    // Insert hole starting at rightmost, going around
    for (int i = 0; i < holeN; i++) {
      int hi = (rightmost + i) % holeN;
      newMerged[k] = hole[hi];
      newIdx[k] = baseIndex + holeOffset + hi;
      k++;
    }

    // Close bridge back: duplicate hole[rightmost] and outer[bestOuter]
    newMerged[k] = hole[rightmost];
    newIdx[k] = baseIndex + holeOffset + rightmost;
    k++;
    newMerged[k] = merged[bestOuter];
    newIdx[k] = mergedIdx[bestOuter];
    k++;

    // Copy rest of outer
    for (size_t i = (size_t)(bestOuter + 1); i < mergedLen; i++) {
      newMerged[k] = merged[i];
      newIdx[k] = mergedIdx[i];
      k++;
    }

    free(merged);
    free(mergedIdx);
    merged = newMerged;
    mergedIdx = newIdx;
    mergedLen = k;

    holeOffset += holeN;
  }

  // Triangulate the merged polygon
  ManifoldVecIVec3 tris = manifold_triangulate_polygon(merged, mergedIdx, mergedLen);

  free(merged);
  free(mergedIdx);
  return tris;
}

#endif // MANIFOLD_POLYGON_H
