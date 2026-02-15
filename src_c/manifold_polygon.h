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

// Check if two 2D points are coincident
static inline bool vec2_coincident(ManifoldVec2 a, ManifoldVec2 b) {
  return a.x == b.x && a.y == b.y;
}

// Ear-clipping triangulation of a simple polygon (supports bridged holes)
// Returns triangles as ManifoldVecIVec3 with indices from `indices` array
// (or flat 0..n-1 if indices is NULL)
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
  int maxIter = (int)(n * n);  // safety limit
  while (remaining > 2 && maxIter-- > 0) {
    bool found = false;
    for (size_t i = 0; i < remaining; i++) {
      size_t prev = (i + remaining - 1) % remaining;
      size_t next = (i + 1) % remaining;

      ManifoldVec2 a = pts[prev], b = pts[i], c = pts[next];
      double cross = vec2_cross(vec2_sub(b, a), vec2_sub(c, a));

      // Only consider CCW ears (positive cross product)
      if (cross < -1e-10) continue;

      // Check no other vertex is strictly inside this ear
      bool ear = true;
      for (size_t j = 0; j < remaining; j++) {
        if (j == prev || j == i || j == next) continue;
        if (vec2_coincident(pts[j], a) || vec2_coincident(pts[j], b) ||
            vec2_coincident(pts[j], c))
          continue;
        // Use strict containment: boundary points don't block ears
        double e1 = vec2_cross(vec2_sub(b, a), vec2_sub(pts[j], a));
        double e2 = vec2_cross(vec2_sub(c, b), vec2_sub(pts[j], b));
        double e3 = vec2_cross(vec2_sub(a, c), vec2_sub(pts[j], c));
        bool strictly_inside = (e1 > 0 && e2 > 0 && e3 > 0) ||
                                (e1 < 0 && e2 < 0 && e3 < 0);
        if (strictly_inside) {
          ear = false;
          break;
        }
      }

      if (ear) {
        // Always generate the triangle (even if zero-area)
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
      // Force-clip: remove any vertex to make progress
      bool forced = false;
      for (size_t i = 0; i < remaining; i++) {
        size_t prev = (i + remaining - 1) % remaining;
        size_t next = (i + 1) % remaining;
        vec_ivec3_push(&tris, manifold_ivec3(idx[prev], idx[i], idx[next]));
        for (size_t k = i; k < remaining - 1; k++) {
          idx[k] = idx[k + 1];
          pts[k] = pts[k + 1];
        }
        remaining--;
        forced = true;
        break;
      }
      if (!forced) break;
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
// Helper: compute signed area of a 2D polygon (positive = CCW, negative = CW)
static inline double mtwh_signed_area(const ManifoldVec2 *v, int n) {
  double a = 0.0;
  for (int i = 0; i < n; i++) {
    int j = (i + 1) % n;
    a += v[i].x * v[j].y - v[j].x * v[i].y;
  }
  return a * 0.5;
}

// Helper: point-in-polygon ray casting test
static inline bool mtwh_pip(ManifoldVec2 p, const ManifoldVec2 *poly, int n) {
  int crossings = 0;
  for (int i = 0, j = n - 1; i < n; j = i++) {
    bool yi = poly[i].y > p.y;
    bool yj = poly[j].y > p.y;
    if (yi != yj) {
      double xint = (poly[j].x - poly[i].x) * (p.y - poly[i].y) /
                    (poly[j].y - poly[i].y) + poly[i].x;
      if (p.x < xint) crossings++;
    }
  }
  return (crossings & 1) == 1;
}

// Triangulate a single outer contour with its holes
static inline ManifoldVecIVec3 manifold_triangulate_one_group(
    const ManifoldVec2 *polyVerts, const int *polySizes, int nPolys,
    const int *polyOffsets, const int *groupIndices, int groupCount,
    int baseIndex) {
  if (groupCount <= 1) {
    int pi = groupIndices[0];
    int off = polyOffsets[pi];
    int n = polySizes[pi];
    // Need to create index mapping for baseIndex offset
    int *idxMap = (int *)malloc((size_t)n * sizeof(int));
    for (int i = 0; i < n; i++) idxMap[i] = baseIndex + off + i;
    ManifoldVecIVec3 tris = manifold_triangulate_polygon(
        polyVerts + off, idxMap, (size_t)n);
    free(idxMap);
    return tris;
  }

  // First polygon in group is the outer contour
  int outerPi = groupIndices[0];
  int outerOff = polyOffsets[outerPi];
  int outerN = polySizes[outerPi];

  size_t totalVerts = 0;
  for (int g = 0; g < groupCount; g++)
    totalVerts += (size_t)polySizes[groupIndices[g]];

  ManifoldVec2 *merged = (ManifoldVec2 *)malloc(
      (totalVerts + 2 * (size_t)(groupCount - 1)) * sizeof(ManifoldVec2));
  int *mergedIdx = (int *)malloc(
      (totalVerts + 2 * (size_t)(groupCount - 1)) * sizeof(int));
  size_t mergedLen = (size_t)outerN;

  for (int i = 0; i < outerN; i++) {
    merged[i] = polyVerts[outerOff + i];
    mergedIdx[i] = baseIndex + outerOff + i;
  }

  for (int h = 1; h < groupCount; h++) {
    int holePi = groupIndices[h];
    int holeOff = polyOffsets[holePi];
    int holeN = polySizes[holePi];
    const ManifoldVec2 *hole = polyVerts + holeOff;

    int rightmost = 0;
    for (int i = 1; i < holeN; i++) {
      if (hole[i].x > hole[rightmost].x ||
          (hole[i].x == hole[rightmost].x && hole[i].y > hole[rightmost].y))
        rightmost = i;
    }

    ManifoldVec2 holeVert = hole[rightmost];
    int bestOuter = 0;
    double bestDist = 1e308;
    for (size_t i = 0; i < mergedLen; i++) {
      double dx = merged[i].x - holeVert.x;
      double dy = merged[i].y - holeVert.y;
      double dist = dx * dx + dy * dy;
      if (dist < bestDist) { bestDist = dist; bestOuter = (int)i; }
    }

    size_t newLen = mergedLen + (size_t)holeN + 2;
    ManifoldVec2 *newMerged = (ManifoldVec2 *)malloc(newLen * sizeof(ManifoldVec2));
    int *newIdx = (int *)malloc(newLen * sizeof(int));
    size_t k = 0;

    for (int i = 0; i <= bestOuter; i++) {
      newMerged[k] = merged[i]; newIdx[k] = mergedIdx[i]; k++;
    }
    for (int i = 0; i < holeN; i++) {
      int hi = (rightmost + i) % holeN;
      newMerged[k] = hole[hi]; newIdx[k] = baseIndex + holeOff + hi; k++;
    }
    newMerged[k] = hole[rightmost]; newIdx[k] = baseIndex + holeOff + rightmost; k++;
    newMerged[k] = merged[bestOuter]; newIdx[k] = mergedIdx[bestOuter]; k++;
    for (size_t i = (size_t)(bestOuter + 1); i < mergedLen; i++) {
      newMerged[k] = merged[i]; newIdx[k] = mergedIdx[i]; k++;
    }

    free(merged); free(mergedIdx);
    merged = newMerged; mergedIdx = newIdx; mergedLen = k;
  }

  ManifoldVecIVec3 tris = manifold_triangulate_polygon(merged, mergedIdx, mergedLen);
  free(merged); free(mergedIdx);
  return tris;
}

// baseIndex: starting vertex index offset for the output triangles
// Returns triangles as ManifoldVecIVec3
// Handles multiple outer contours with their respective holes.
static inline ManifoldVecIVec3 manifold_triangulate_with_holes(
    const ManifoldVec2 *polyVerts, const int *polySizes, int nPolys,
    int baseIndex) {
  if (nPolys <= 0) {
    ManifoldVecIVec3 empty = {0};
    return empty;
  }
  if (nPolys == 1) {
    return manifold_triangulate_polygon(polyVerts, NULL, (size_t)polySizes[0]);
  }

  // Compute offsets and signed areas
  int *offsets = (int *)malloc((size_t)nPolys * sizeof(int));
  double *areas = (double *)malloc((size_t)nPolys * sizeof(double));
  offsets[0] = 0;
  for (int i = 0; i < nPolys; i++) {
    if (i > 0) offsets[i] = offsets[i-1] + polySizes[i-1];
    areas[i] = mtwh_signed_area(polyVerts + offsets[i], polySizes[i]);
  }

  // Identify outer contours (positive area = CCW)
  int nOuter = 0;
  int *outerIndices = (int *)malloc((size_t)nPolys * sizeof(int));
  for (int i = 0; i < nPolys; i++) {
    if (areas[i] > 0) outerIndices[nOuter++] = i;
  }

  if (nOuter == 0) {
    // Fallback: treat first as outer, rest as holes (original behavior)
    free(offsets); free(areas); free(outerIndices);
    // Use old single-group approach
    int *group = (int *)malloc((size_t)nPolys * sizeof(int));
    int *offs = (int *)malloc((size_t)nPolys * sizeof(int));
    offs[0] = 0;
    for (int i = 0; i < nPolys; i++) {
      group[i] = i;
      if (i > 0) offs[i] = offs[i-1] + polySizes[i-1];
    }
    ManifoldVecIVec3 tris = manifold_triangulate_one_group(
        polyVerts, polySizes, nPolys, offs, group, nPolys, baseIndex);
    free(group); free(offs);
    return tris;
  }

  if (nOuter == 1) {
    // Single outer with holes: use original approach
    int *group = (int *)malloc((size_t)nPolys * sizeof(int));
    for (int i = 0; i < nPolys; i++) group[i] = i;
    ManifoldVecIVec3 tris = manifold_triangulate_one_group(
        polyVerts, polySizes, nPolys, offsets, group, nPolys, baseIndex);
    free(group); free(offsets); free(areas); free(outerIndices);
    return tris;
  }

  // Multiple outers: assign each hole to its enclosing outer
  int *holeOwner = (int *)malloc((size_t)nPolys * sizeof(int));
  for (int i = 0; i < nPolys; i++) holeOwner[i] = -1;
  for (int i = 0; i < nPolys; i++) {
    if (areas[i] > 0) continue; // outer
    ManifoldVec2 pt = polyVerts[offsets[i]];
    double bestArea = 1e30;
    int bestOwner = 0;
    for (int j = 0; j < nOuter; j++) {
      int oi = outerIndices[j];
      if (mtwh_pip(pt, polyVerts + offsets[oi], polySizes[oi])) {
        if (areas[oi] < bestArea) {
          bestArea = areas[oi];
          bestOwner = j;
        }
      }
    }
    holeOwner[i] = bestOwner;
  }

  // Triangulate each group separately and merge results
  ManifoldVecIVec3 allTris = {0};
  int *groupBuf = (int *)malloc((size_t)nPolys * sizeof(int));

  for (int j = 0; j < nOuter; j++) {
    int gc = 0;
    groupBuf[gc++] = outerIndices[j]; // outer first
    for (int i = 0; i < nPolys; i++) {
      if (holeOwner[i] == j) groupBuf[gc++] = i;
    }

    ManifoldVecIVec3 tris = manifold_triangulate_one_group(
        polyVerts, polySizes, nPolys, offsets, groupBuf, gc, baseIndex);
    for (size_t t = 0; t < tris.len; t++) {
      vec_ivec3_push(&allTris, tris.data[t]);
    }
    vec_ivec3_free(&tris);
  }

  free(groupBuf); free(holeOwner); free(offsets); free(areas); free(outerIndices);
  return allTris;
}

#endif // MANIFOLD_POLYGON_H
