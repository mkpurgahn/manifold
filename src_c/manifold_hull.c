// Copyright 2024 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// Convex hull using incremental algorithm.
// For each point, find visible faces and replace them with new faces
// connecting the point to the horizon edges.

#include "manifold_hull.h"
#include <stdlib.h>
#include <string.h>
#include <float.h>

typedef struct {
  int v[3];     // vertex indices
  int adj[3];   // adjacent face indices (-1 = boundary)
  ManifoldVec3 normal;
  double d;     // plane offset: normal·p + d = 0
  bool removed;
} HullFace;

static void hull_face_init(HullFace *f, int v0, int v1, int v2,
                           const ManifoldVec3 *pts) {
  f->v[0] = v0; f->v[1] = v1; f->v[2] = v2;
  f->adj[0] = f->adj[1] = f->adj[2] = -1;
  f->removed = false;

  ManifoldVec3 ab = vec3_sub(pts[v1], pts[v0]);
  ManifoldVec3 ac = vec3_sub(pts[v2], pts[v0]);
  f->normal = vec3_cross(ab, ac);
  double len = vec3_length(f->normal);
  if (len > 1e-15) {
    f->normal = vec3_scale(f->normal, 1.0 / len);
  }
  f->d = -vec3_dot(f->normal, pts[v0]);
}

static double hull_face_dist(const HullFace *f, ManifoldVec3 p) {
  return vec3_dot(f->normal, p) + f->d;
}

// Find initial tetrahedron from extreme points
static bool find_initial_tet(const ManifoldVec3 *pts, size_t n,
                              int out[4]) {
  if (n < 4) return false;

  // Find extreme points along each axis
  int minX = 0, maxX = 0;
  for (size_t i = 1; i < n; i++) {
    if (pts[i].x < pts[minX].x) minX = (int)i;
    if (pts[i].x > pts[maxX].x) maxX = (int)i;
  }
  if (minX == maxX) return false;
  out[0] = minX;
  out[1] = maxX;

  // Find point farthest from line out[0]-out[1]
  ManifoldVec3 lineDir = vec3_sub(pts[out[1]], pts[out[0]]);
  double lineLen2 = vec3_dot(lineDir, lineDir);
  if (lineLen2 < 1e-30) return false;

  double bestDist2 = -1;
  out[2] = -1;
  for (size_t i = 0; i < n; i++) {
    if ((int)i == out[0] || (int)i == out[1]) continue;
    ManifoldVec3 v = vec3_sub(pts[i], pts[out[0]]);
    double t = vec3_dot(v, lineDir) / lineLen2;
    ManifoldVec3 proj = vec3_add(pts[out[0]], vec3_scale(lineDir, t));
    ManifoldVec3 diff = vec3_sub(pts[i], proj);
    double d2 = vec3_dot(diff, diff);
    if (d2 > bestDist2) { bestDist2 = d2; out[2] = (int)i; }
  }
  if (out[2] < 0 || bestDist2 < 1e-30) return false;

  // Find point farthest from plane of first 3 points
  ManifoldVec3 ab = vec3_sub(pts[out[1]], pts[out[0]]);
  ManifoldVec3 ac = vec3_sub(pts[out[2]], pts[out[0]]);
  ManifoldVec3 planeN = vec3_cross(ab, ac);
  double planeLen = vec3_length(planeN);
  if (planeLen < 1e-15) return false;
  planeN = vec3_scale(planeN, 1.0 / planeLen);
  double planeD = -vec3_dot(planeN, pts[out[0]]);

  double bestAbsDist = -1;
  out[3] = -1;
  for (size_t i = 0; i < n; i++) {
    if ((int)i == out[0] || (int)i == out[1] || (int)i == out[2]) continue;
    double dist = fabs(vec3_dot(planeN, pts[i]) + planeD);
    if (dist > bestAbsDist) { bestAbsDist = dist; out[3] = (int)i; }
  }
  if (out[3] < 0 || bestAbsDist < 1e-15) return false;

  // Ensure out[3] is above the plane (swap if needed for consistent orientation)
  double d = vec3_dot(planeN, pts[out[3]]) + planeD;
  if (d < 0) {
    int tmp = out[1]; out[1] = out[2]; out[2] = tmp;
  }

  return true;
}

// Dynamic array of HullFace
typedef struct {
  HullFace *data;
  size_t len;
  size_t cap;
} HullFaceVec;

static void hfv_push(HullFaceVec *v, HullFace f) {
  if (v->len >= v->cap) {
    v->cap = v->cap ? v->cap * 2 : 16;
    v->data = (HullFace *)realloc(v->data, v->cap * sizeof(HullFace));
  }
  v->data[v->len++] = f;
}

// Edge structure for horizon tracking
typedef struct {
  int v0, v1;
  int face;
  int edge; // which edge of the face (0,1,2)
} HorizonEdge;

void manifold_convex_hull(ManifoldImpl *impl, const ManifoldVec3 *points,
                           size_t numPoints) {
  manifold_impl_init(impl);

  if (numPoints < 4) {
    manifold_impl_make_empty(impl, MANIFOLD_ERROR_INVALID_CONSTRUCTION);
    return;
  }

  // Find initial tetrahedron
  int tetIdx[4];
  if (!find_initial_tet(points, numPoints, tetIdx)) {
    manifold_impl_make_empty(impl, MANIFOLD_ERROR_INVALID_CONSTRUCTION);
    return;
  }

  // Build initial tetrahedron (4 faces)
  HullFaceVec faces = {0};
  HullFace f;

  hull_face_init(&f, tetIdx[0], tetIdx[1], tetIdx[2], points);
  hfv_push(&faces, f);  // face 0
  hull_face_init(&f, tetIdx[0], tetIdx[2], tetIdx[3], points);
  hfv_push(&faces, f);  // face 1
  hull_face_init(&f, tetIdx[0], tetIdx[3], tetIdx[1], points);
  hfv_push(&faces, f);  // face 2
  hull_face_init(&f, tetIdx[1], tetIdx[3], tetIdx[2], points);
  hfv_push(&faces, f);  // face 3

  // Set adjacency for initial tetrahedron
  // Face 0: 012, Face 1: 023, Face 2: 031, Face 3: 132
  faces.data[0].adj[0] = 2; faces.data[0].adj[1] = 3; faces.data[0].adj[2] = 1;
  faces.data[1].adj[0] = 0; faces.data[1].adj[1] = 3; faces.data[1].adj[2] = 2;
  faces.data[2].adj[0] = 1; faces.data[2].adj[1] = 3; faces.data[2].adj[2] = 0;
  faces.data[3].adj[0] = 2; faces.data[3].adj[1] = 1; faces.data[3].adj[2] = 0;

  // For each remaining point, if it's outside any face, update the hull
  for (size_t pi = 0; pi < numPoints; pi++) {
    if ((int)pi == tetIdx[0] || (int)pi == tetIdx[1] ||
        (int)pi == tetIdx[2] || (int)pi == tetIdx[3])
      continue;

    // Find a visible face
    int visibleFace = -1;
    double maxDist = 1e-10;
    for (size_t fi = 0; fi < faces.len; fi++) {
      if (faces.data[fi].removed) continue;
      double d = hull_face_dist(&faces.data[fi], points[pi]);
      if (d > maxDist) { maxDist = d; visibleFace = (int)fi; }
    }
    if (visibleFace < 0) continue; // point inside hull

    // Mark all visible faces
    bool *visible = (bool *)calloc(faces.len, sizeof(bool));
    // BFS to find all visible faces starting from visibleFace
    int *stack = (int *)malloc(faces.len * sizeof(int));
    int stackLen = 0;
    stack[stackLen++] = visibleFace;
    visible[visibleFace] = true;

    while (stackLen > 0) {
      int fi = stack[--stackLen];
      for (int ei = 0; ei < 3; ei++) {
        int adj = faces.data[fi].adj[ei];
        if (adj >= 0 && !visible[adj] && !faces.data[adj].removed) {
          double d = hull_face_dist(&faces.data[adj], points[pi]);
          if (d > 1e-10) {
            visible[adj] = true;
            stack[stackLen++] = adj;
          }
        }
      }
    }

    // Find horizon edges (edges between visible and non-visible faces)
    HorizonEdge *horizon = (HorizonEdge *)malloc(faces.len * 3 * sizeof(HorizonEdge));
    int horizonLen = 0;

    for (size_t fi = 0; fi < faces.len; fi++) {
      if (!visible[fi]) continue;
      for (int ei = 0; ei < 3; ei++) {
        int adj = faces.data[fi].adj[ei];
        if (adj < 0 || !visible[adj]) {
          HorizonEdge he;
          he.v0 = faces.data[fi].v[ei];
          he.v1 = faces.data[fi].v[(ei + 1) % 3];
          he.face = adj;
          he.edge = ei;
          horizon[horizonLen++] = he;
        }
      }
    }

    if (horizonLen < 3) {
      free(visible);
      free(stack);
      free(horizon);
      continue;
    }

    // Remove visible faces
    for (size_t fi = 0; fi < faces.len; fi++) {
      if (visible[fi]) faces.data[fi].removed = true;
    }

    // Create new faces connecting point to horizon edges
    int *newFaceIdx = (int *)malloc(horizonLen * sizeof(int));
    for (int hi = 0; hi < horizonLen; hi++) {
      HullFace nf;
      hull_face_init(&nf, horizon[hi].v0, horizon[hi].v1, (int)pi, points);
      newFaceIdx[hi] = (int)faces.len;
      hfv_push(&faces, nf);

      // Set adjacency to the non-visible adjacent face
      if (horizon[hi].face >= 0) {
        HullFace *adjFace = &faces.data[horizon[hi].face];
        for (int aei = 0; aei < 3; aei++) {
          // Find the edge in the adjacent face that connects to this visible face
          int av0 = adjFace->v[aei];
          int av1 = adjFace->v[(aei + 1) % 3];
          if ((av0 == horizon[hi].v1 && av1 == horizon[hi].v0) ||
              (av0 == horizon[hi].v0 && av1 == horizon[hi].v1)) {
            adjFace->adj[aei] = newFaceIdx[hi];
            faces.data[newFaceIdx[hi]].adj[0] = horizon[hi].face;
            break;
          }
        }
      }
    }

    // Set adjacency between new faces
    for (int hi = 0; hi < horizonLen; hi++) {
      for (int hj = hi + 1; hj < horizonLen; hj++) {
        // Check if these two new faces share an edge
        int fi1 = newFaceIdx[hi];
        int fi2 = newFaceIdx[hj];
        for (int e1 = 0; e1 < 3; e1++) {
          for (int e2 = 0; e2 < 3; e2++) {
            int a0 = faces.data[fi1].v[e1];
            int a1 = faces.data[fi1].v[(e1 + 1) % 3];
            int b0 = faces.data[fi2].v[e2];
            int b1 = faces.data[fi2].v[(e2 + 1) % 3];
            if (a0 == b1 && a1 == b0) {
              faces.data[fi1].adj[e1] = fi2;
              faces.data[fi2].adj[e2] = fi1;
            }
          }
        }
      }
    }

    free(newFaceIdx);
    free(visible);
    free(stack);
    free(horizon);
  }

  // Collect surviving faces and build mesh
  // First, collect all used vertices and remap
  int *vertMap = (int *)calloc(numPoints, sizeof(int));
  for (size_t i = 0; i < numPoints; i++) vertMap[i] = -1;

  int vertCount = 0;
  ManifoldVecIVec3 triVerts = {0};
  for (size_t fi = 0; fi < faces.len; fi++) {
    if (faces.data[fi].removed) continue;
    for (int vi = 0; vi < 3; vi++) {
      if (vertMap[faces.data[fi].v[vi]] < 0) {
        vertMap[faces.data[fi].v[vi]] = vertCount++;
      }
    }
    vec_ivec3_push(&triVerts, manifold_ivec3(
        vertMap[faces.data[fi].v[0]],
        vertMap[faces.data[fi].v[1]],
        vertMap[faces.data[fi].v[2]]));
  }

  // Copy vertices
  impl->vertPos = vec_vec3_create_n((size_t)vertCount);
  for (size_t i = 0; i < numPoints; i++) {
    if (vertMap[i] >= 0) {
      impl->vertPos.data[vertMap[i]] = points[i];
    }
  }

  ManifoldVecIVec3 emptyTriVert = {0};
  manifold_impl_create_halfedges(impl, &triVerts, &emptyTriVert);
  manifold_impl_initialize_original(impl);
  manifold_impl_calculate_bbox(impl);
  manifold_impl_set_epsilon(impl, -1.0, false);
  manifold_impl_sort_geometry(impl);
  manifold_impl_set_normals_and_coplanar(impl);

  vec_ivec3_free(&triVerts);
  vec_ivec3_free(&emptyTriVert);
  free(vertMap);
  free(faces.data);
}
