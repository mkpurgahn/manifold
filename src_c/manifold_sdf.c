// Copyright 2023 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// SDF (Signed Distance Function) level set meshing for the C11 Manifold port.
// Implements marching tetrahedra on a body-centered cubic (BCC) grid,
// matching the C++ algorithm for correct manifold topology.

#include "manifold_impl.h"
#include <float.h>

// Constants matching C++
#define SDF_CROSSING (-2)
#define SDF_NONE     (-1)
#define SDF_HASH_OPEN UINT64_MAX
// Maximum fraction of spacing that a vert can move
#define SDF_KS 0.25
// Corresponding approximate distance ratio bound
#define SDF_KD (1.0 / SDF_KS - 1.0)
// Maximum number of opposed verts (of 7) to allow collapse
#define SDF_MAX_OPPOSED 3

// Tet triangulation tables (from C++ TetTri0/TetTri1)
static const int tetTri0[16][3] = {
  {-1,-1,-1}, {0,3,4}, {0,1,5}, {1,5,3},
  {1,4,2}, {1,0,3}, {2,5,0}, {5,3,2},
  {2,3,5}, {0,5,2}, {3,0,1}, {2,4,1},
  {3,5,1}, {5,1,0}, {4,3,0}, {-1,-1,-1}
};
static const int tetTri1[16][3] = {
  {-1,-1,-1}, {-1,-1,-1}, {-1,-1,-1}, {3,4,1},
  {-1,-1,-1}, {3,2,1}, {0,4,2}, {-1,-1,-1},
  {-1,-1,-1}, {2,4,0}, {1,2,3}, {-1,-1,-1},
  {1,4,3}, {-1,-1,-1}, {-1,-1,-1}, {-1,-1,-1}
};

// BCC neighbor offsets: 7 owned edges + 7 opposite edges
static const int sdf_neighbors[14][4] = {
  {0,0,0,1}, {1,0,0,0}, {0,1,0,0}, {0,0,1,0},
  {-1,0,0,1}, {0,-1,0,1}, {0,0,-1,1},
  {-1,-1,-1,1}, {-1,0,0,0}, {0,-1,0,0}, {0,0,-1,0},
  {0,-1,-1,1}, {-1,0,-1,1}, {-1,-1,0,1}
};

static inline int sdf_next3(int i) {
  static const int n[3] = {1, 2, 0};
  return n[i];
}

static inline int sdf_prev3(int i) {
  static const int p[3] = {2, 0, 1};
  return p[i];
}

static inline ManifoldIVec4 sdf_neighbor(ManifoldIVec4 base, int i) {
  ManifoldIVec4 n;
  n.x = base.x + sdf_neighbors[i][0];
  n.y = base.y + sdf_neighbors[i][1];
  n.z = base.z + sdf_neighbors[i][2];
  n.w = base.w + sdf_neighbors[i][3];
  if (n.w == 2) {
    n.x += 1; n.y += 1; n.z += 1;
    n.w = 0;
  }
  return n;
}

static inline int sdf_log2_ceil(int v) {
  int r = 0;
  while ((1 << r) < v) r++;
  return r;
}

static inline uint64_t sdf_encode(ManifoldIVec4 pos, ManifoldIVec3 gpow) {
  return (uint64_t)pos.w |
         ((uint64_t)pos.z << 1) |
         ((uint64_t)pos.y << (1 + gpow.z)) |
         ((uint64_t)pos.x << (1 + gpow.z + gpow.y));
}

static inline ManifoldIVec4 sdf_decode(uint64_t idx, ManifoldIVec3 gpow) {
  ManifoldIVec4 p;
  p.w = idx & 1; idx >>= 1;
  p.z = idx & ((1 << gpow.z) - 1); idx >>= gpow.z;
  p.y = idx & ((1 << gpow.y) - 1); idx >>= gpow.y;
  p.x = idx & ((1 << gpow.x) - 1);
  return p;
}

static inline ManifoldVec3 sdf_position(ManifoldIVec4 gi, ManifoldVec3 origin,
                                         ManifoldVec3 spacing) {
  double off = (gi.w == 1) ? 0.0 : -0.5;
  return manifold_vec3(origin.x + spacing.x * (gi.x + off),
                       origin.y + spacing.y * (gi.y + off),
                       origin.z + spacing.z * (gi.z + off));
}

static inline ManifoldVec3 sdf_bound_pos(ManifoldVec3 pos, ManifoldVec3 origin,
                                          ManifoldVec3 spacing,
                                          ManifoldIVec3 gridSize) {
  ManifoldVec3 maxPos = manifold_vec3(
    origin.x + spacing.x * (gridSize.x - 1),
    origin.y + spacing.y * (gridSize.y - 1),
    origin.z + spacing.z * (gridSize.z - 1));
  return manifold_vec3(
    fmin(fmax(pos.x, origin.x), maxPos.x),
    fmin(fmax(pos.y, origin.y), maxPos.y),
    fmin(fmax(pos.z, origin.z), maxPos.z));
}

static inline int ivec3_minelem(ManifoldIVec3 v) {
  int m = v.x < v.y ? v.x : v.y;
  return m < v.z ? m : v.z;
}

static double sdf_bounded_eval(ManifoldIVec4 gi, ManifoldVec3 origin,
                                ManifoldVec3 spacing, ManifoldIVec3 gridSize,
                                double level,
                                double (*sdf)(double, double, double, void*),
                                void *ctx) {
  ManifoldIVec3 xyz = {gi.x, gi.y, gi.z};
  int lowerBoundDist = ivec3_minelem(xyz);
  ManifoldIVec3 upper = {gridSize.x - xyz.x, gridSize.y - xyz.y,
                          gridSize.z - xyz.z};
  int upperBoundDist = ivec3_minelem(upper);
  int boundDist = lowerBoundDist < (upperBoundDist - gi.w) ?
                  lowerBoundDist : (upperBoundDist - gi.w);
  if (boundDist < 0) return 0.0;
  ManifoldVec3 pos = sdf_position(gi, origin, spacing);
  double d = sdf(pos.x, pos.y, pos.z, ctx) - level;
  return boundDist == 0 ? fmin(d, 0.0) : d;
}

// ITP root finding for surface location
static ManifoldVec3 sdf_find_surface(ManifoldVec3 pos0, double d0,
                                      ManifoldVec3 pos1, double d1,
                                      double tol, double level,
                                      double (*sdf)(double,double,double,void*),
                                      void *ctx) {
  if (d0 == 0) return pos0;
  if (d1 == 0) return pos1;

  const double k = 0.1;
  ManifoldVec3 diff = vec3_sub(pos0, pos1);
  double len = vec3_length(diff);
  double check = len > 0 ? 2 * tol / len : 1.0;
  double frac = 1;
  double biFrac = 1;
  while (frac > check) {
    double tInterp = d0 / (d0 - d1);
    double t = tInterp * (1.0 - k) + 0.5 * k;
    double r = biFrac / frac - 0.5;
    double x;
    if (fabs(t - 0.5) < r) x = t;
    else x = 0.5 - r * (t < 0.5 ? 1 : -1);

    ManifoldVec3 mid = vec3_add(vec3_scale(pos0, 1.0 - x),
                                 vec3_scale(pos1, x));
    double d = sdf(mid.x, mid.y, mid.z, ctx) - level;

    if ((d > 0) == (d0 > 0)) {
      d0 = d; pos0 = mid; frac *= 1 - x;
    } else {
      d1 = d; pos1 = mid; frac *= x;
    }
    biFrac /= 2;
  }
  double t = d0 / (d0 - d1);
  return vec3_add(vec3_scale(pos0, 1.0 - t), vec3_scale(pos1, t));
}

// GridVert: each vertex on the BCC grid
typedef struct {
  double distance;
  int movedVert;
  int edgeVerts[7];
} SdfGridVert;

static inline SdfGridVert sdf_gridvert_default(void) {
  SdfGridVert gv;
  gv.distance = NAN;
  gv.movedVert = SDF_NONE;
  for (int i = 0; i < 7; i++) gv.edgeVerts[i] = SDF_NONE;
  return gv;
}

static inline bool sdf_gv_has_moved(const SdfGridVert *gv) {
  return gv->movedVert >= 0;
}
static inline bool sdf_gv_same_side(const SdfGridVert *gv, double dist) {
  return (dist > 0) == (gv->distance > 0);
}
static inline int sdf_gv_inside(const SdfGridVert *gv) {
  return gv->distance > 0 ? 1 : -1;
}
static inline int sdf_gv_neighbor_inside(const SdfGridVert *gv, int i) {
  return sdf_gv_inside(gv) * (gv->edgeVerts[i] == SDF_NONE ? 1 : -1);
}

// Hash table for GridVert (open addressing)
typedef struct {
  uint64_t *keys;
  SdfGridVert *values;
  size_t size;
  size_t used;
} SdfGridHash;

static SdfGridHash sdf_hash_create(size_t size) {
  size_t s = 1;
  while (s < size) s <<= 1;
  SdfGridHash h;
  h.size = s;
  h.used = 0;
  h.keys = (uint64_t*)malloc(s * sizeof(uint64_t));
  h.values = (SdfGridVert*)malloc(s * sizeof(SdfGridVert));
  for (size_t i = 0; i < s; i++) {
    h.keys[i] = SDF_HASH_OPEN;
    h.values[i] = sdf_gridvert_default();
  }
  return h;
}

static void sdf_hash_free(SdfGridHash *h) {
  free(h->keys); free(h->values);
  h->keys = NULL; h->values = NULL;
  h->size = 0; h->used = 0;
}

static bool sdf_hash_full(const SdfGridHash *h) {
  return h->used * 2 > h->size;
}

static void sdf_hash_insert(SdfGridHash *h, uint64_t key, SdfGridVert val) {
  if (h->size == 0) return;
  uint32_t idx = (uint32_t)(manifold_hash64bit(key) & (h->size - 1));
  while (1) {
    if (sdf_hash_full(h)) return;
    if (h->keys[idx] == SDF_HASH_OPEN) {
      h->keys[idx] = key;
      h->values[idx] = val;
      h->used++;
      return;
    }
    if (h->keys[idx] == key) return;
    idx = (idx + 1) & (uint32_t)(h->size - 1);
  }
}

static const SdfGridVert* sdf_hash_find_const(const SdfGridHash *h,
                                               uint64_t key) {
  if (h->size == 0) return NULL;
  uint32_t idx = (uint32_t)(manifold_hash64bit(key) & (h->size - 1));
  while (1) {
    if (h->keys[idx] == key) return &h->values[idx];
    if (h->keys[idx] == SDF_HASH_OPEN) return &h->values[idx];
    idx = (idx + 1) & (uint32_t)(h->size - 1);
  }
}

// SDF level set meshing using marching tetrahedra on BCC grid
void manifold_impl_level_set(ManifoldImpl *impl,
                              double (*sdf)(double x, double y, double z,
                                            void *ctx),
                              void *ctx, ManifoldBox bounds,
                              double edgeLength, double level,
                              double tolerance_in) {
  manifold_impl_init(impl);

  double tolerance = tolerance_in > 0 ? tolerance_in : DBL_MAX;

  ManifoldVec3 dim = manifold_box_size(bounds);
  ManifoldIVec3 gridSize = {
    (int)(dim.x / edgeLength + 1.0),
    (int)(dim.y / edgeLength + 1.0),
    (int)(dim.z / edgeLength + 1.0)
  };
  if (gridSize.x < 2 || gridSize.y < 2 || gridSize.z < 2) {
    manifold_impl_make_empty(impl, MANIFOLD_ERROR_INVALID_CONSTRUCTION);
    return;
  }
  ManifoldVec3 spacing = manifold_vec3(
    dim.x / (gridSize.x - 1),
    dim.y / (gridSize.y - 1),
    dim.z / (gridSize.z - 1));

  ManifoldIVec3 gridPow = {
    sdf_log2_ceil(gridSize.x + 2) + 1,
    sdf_log2_ceil(gridSize.y + 2) + 1,
    sdf_log2_ceil(gridSize.z + 2) + 1
  };
  ManifoldIVec4 maxIV = {gridSize.x + 2, gridSize.y + 2, gridSize.z + 2, 1};
  uint64_t maxIndex = sdf_encode(maxIV, gridPow);
  ManifoldVec3 origin = bounds.min;

  // Evaluate bounded SDF on the full BCC grid
  ManifoldIVec4 voxelOffset = {1, 1, 1, 0};
  double *voxels = (double*)malloc(maxIndex * sizeof(double));
  for (uint64_t idx = 0; idx < maxIndex; idx++) {
    ManifoldIVec4 di = sdf_decode(idx, gridPow);
    ManifoldIVec4 gi = {di.x - voxelOffset.x, di.y - voxelOffset.y,
                         di.z - voxelOffset.z, di.w - voxelOffset.w};
    voxels[idx] = sdf_bounded_eval(gi, origin, spacing, gridSize, level,
                                    sdf, ctx);
  }

  // Phase 1: Identify near-surface grid vertices
  uint64_t encMax = sdf_encode(
    (ManifoldIVec4){gridSize.x, gridSize.y, gridSize.z, 1}, gridPow);
  double powApprox = pow((double)maxIndex, 0.667);
  size_t tableSize = (size_t)(2 * maxIndex);
  size_t altSize = (size_t)(10 * powApprox);
  if (altSize < tableSize) tableSize = altSize;

  SdfGridHash gridVerts = sdf_hash_create(tableSize);

  // Allocate vertex positions
  size_t maxVerts = gridVerts.size * 7;
  ManifoldVecVec3 vertPos = {0};
  vec_vec3_resize(&vertPos, maxVerts);
  int vertCount = 0;

  // NearSurface pass: find grid verts near the surface, optionally snap
  for (uint64_t index = 0; index < encMax; index++) {
    if (sdf_hash_full(&gridVerts)) break;
    ManifoldIVec4 gridIndex = sdf_decode(index, gridPow);
    if (gridIndex.x > gridSize.x || gridIndex.y > gridSize.y ||
        gridIndex.z > gridSize.z) continue;

    SdfGridVert gv = sdf_gridvert_default();
    ManifoldIVec4 gi_vo = {gridIndex.x + voxelOffset.x,
                            gridIndex.y + voxelOffset.y,
                            gridIndex.z + voxelOffset.z,
                            gridIndex.w + voxelOffset.w};
    gv.distance = voxels[sdf_encode(gi_vo, gridPow)];

    bool keep = false;
    double vMax = 0;
    int closestNeighbor = -1;
    int opposedVerts = 0;
    for (int i = 0; i < 7; i++) {
      ManifoldIVec4 ni = sdf_neighbor(gridIndex, i);
      ManifoldIVec4 ni_vo = {ni.x + voxelOffset.x, ni.y + voxelOffset.y,
                              ni.z + voxelOffset.z, ni.w + voxelOffset.w};
      double val = voxels[sdf_encode(ni_vo, gridPow)];

      ManifoldIVec4 oi = sdf_neighbor(gridIndex, i + 7);
      ManifoldIVec4 oi_vo = {oi.x + voxelOffset.x, oi.y + voxelOffset.y,
                              oi.z + voxelOffset.z, oi.w + voxelOffset.w};
      double valOp = voxels[sdf_encode(oi_vo, gridPow)];

      if (!sdf_gv_same_side(&gv, val)) {
        gv.edgeVerts[i] = SDF_CROSSING;
        keep = true;
        if (!sdf_gv_same_side(&gv, valOp)) ++opposedVerts;
        if (fabs(val) > SDF_KD * fabs(gv.distance) &&
            fabs(val) > fabs(vMax)) {
          vMax = val;
          closestNeighbor = i;
        }
      } else if (!sdf_gv_same_side(&gv, valOp) &&
                 fabs(valOp) > SDF_KD * fabs(gv.distance) &&
                 fabs(valOp) > fabs(vMax)) {
        vMax = valOp;
        closestNeighbor = i + 7;
      }
    }

    // Vertex snapping: collapse crossing edges into this grid vert
    if (closestNeighbor >= 0 && opposedVerts <= SDF_MAX_OPPOSED) {
      ManifoldVec3 gridPos = sdf_position(gridIndex, origin, spacing);
      ManifoldIVec4 neighborIndex = sdf_neighbor(gridIndex, closestNeighbor);
      ManifoldVec3 pos = sdf_find_surface(
        gridPos, gv.distance,
        sdf_position(neighborIndex, origin, spacing),
        vMax, tolerance, level, sdf, ctx);
      ManifoldVec3 delta = vec3_sub(pos, gridPos);
      if (fabs(delta.x) < SDF_KS * spacing.x &&
          fabs(delta.y) < SDF_KS * spacing.y &&
          fabs(delta.z) < SDF_KS * spacing.z) {
        int idx = vertCount++;
        if ((size_t)idx >= vertPos.len) vec_vec3_resize(&vertPos, vertPos.len*2);
        vertPos.data[idx] = sdf_bound_pos(pos, origin, spacing, gridSize);
        gv.movedVert = idx;
        for (int j = 0; j < 7; j++) {
          if (gv.edgeVerts[j] == SDF_CROSSING) gv.edgeVerts[j] = idx;
        }
        keep = true;
      }
    } else {
      for (int j = 0; j < 7; j++) gv.edgeVerts[j] = SDF_NONE;
    }

    if (keep) sdf_hash_insert(&gridVerts, index, gv);
  }

  // Phase 2: Compute edge vertices for crossing edges not yet assigned
  for (size_t hidx = 0; hidx < gridVerts.size; hidx++) {
    if (gridVerts.keys[hidx] == SDF_HASH_OPEN) continue;
    SdfGridVert *gv = &gridVerts.values[hidx];
    if (sdf_gv_has_moved(gv)) continue;

    ManifoldIVec4 gridIndex = sdf_decode(gridVerts.keys[hidx], gridPow);
    ManifoldVec3 position = sdf_position(gridIndex, origin, spacing);

    for (int i = 0; i < 7; i++) {
      ManifoldIVec4 neighborIndex = sdf_neighbor(gridIndex, i);
      uint64_t nkey = sdf_encode(neighborIndex, gridPow);
      const SdfGridVert *neighbor = sdf_hash_find_const(&gridVerts, nkey);

      double val;
      if (neighbor && isfinite(neighbor->distance)) {
        val = neighbor->distance;
      } else {
        ManifoldIVec4 ni_vo = {neighborIndex.x + voxelOffset.x,
                                neighborIndex.y + voxelOffset.y,
                                neighborIndex.z + voxelOffset.z,
                                neighborIndex.w + voxelOffset.w};
        val = voxels[sdf_encode(ni_vo, gridPow)];
      }
      if (sdf_gv_same_side(gv, val)) continue;

      if (neighbor && sdf_gv_has_moved(neighbor)) {
        gv->edgeVerts[i] = neighbor->movedVert;
        continue;
      }

      int idx = vertCount++;
      if ((size_t)idx >= vertPos.len) vec_vec3_resize(&vertPos, vertPos.len*2);
      ManifoldVec3 pos = sdf_find_surface(
        position, gv->distance,
        sdf_position(neighborIndex, origin, spacing),
        val, tolerance, level, sdf, ctx);
      vertPos.data[idx] = sdf_bound_pos(pos, origin, spacing, gridSize);
      gv->edgeVerts[i] = idx;
    }
  }

  // Phase 3: Build triangles from tetrahedra
  ManifoldVecIVec3 triVerts = {0};
  size_t estTris = gridVerts.used * 12;
  vec_ivec3_resize(&triVerts, estTris > 128 ? estTris : 128);
  int triCount = 0;

  for (size_t hidx = 0; hidx < gridVerts.size; hidx++) {
    if (gridVerts.keys[hidx] == SDF_HASH_OPEN) continue;
    const SdfGridVert *base = &gridVerts.values[hidx];
    ManifoldIVec4 baseIndex = sdf_decode(gridVerts.keys[hidx], gridPow);

    ManifoldIVec4 leadIndex = baseIndex;
    if (leadIndex.w == 0) {
      leadIndex.w = 1;
    } else {
      leadIndex.x += 1; leadIndex.y += 1; leadIndex.z += 1;
      leadIndex.w = 0;
    }

    // 6 tetrahedra around edge 0 in the (1,1,1) direction
    int tet[4];
    tet[0] = sdf_gv_neighbor_inside(base, 0);
    tet[1] = sdf_gv_inside(base);

    ManifoldIVec4 thisIndex = baseIndex;
    thisIndex.x += 1;
    SdfGridVert thisVert = *sdf_hash_find_const(&gridVerts,
                             sdf_encode(thisIndex, gridPow));
    tet[2] = sdf_gv_neighbor_inside(base, 1);

    for (int i = 0; i < 3; i++) {
      ManifoldIVec4 ti = leadIndex;
      int prev = sdf_prev3(i);
      int next = sdf_next3(i);
      // Decrement prev3 component
      if (prev == 0) ti.x -= 1;
      else if (prev == 1) ti.y -= 1;
      else ti.z -= 1;

      SdfGridVert nextVert;
      int prevCoord = (prev == 0) ? ti.x : ((prev == 1) ? ti.y : ti.z);
      if (prevCoord < 0) {
        nextVert = sdf_gridvert_default();
      } else {
        nextVert = *sdf_hash_find_const(&gridVerts,
                     sdf_encode(ti, gridPow));
      }
      tet[3] = sdf_gv_neighbor_inside(base, prev + 4);

      int edges1[6];
      edges1[0] = base->edgeVerts[0];
      edges1[1] = base->edgeVerts[i + 1];
      edges1[2] = nextVert.edgeVerts[next + 4];
      edges1[3] = nextVert.edgeVerts[prev + 1];
      edges1[4] = thisVert.edgeVerts[i + 4];
      edges1[5] = base->edgeVerts[prev + 4];
      thisVert = nextVert;

      // Create triangles from tet
      int ci = (tet[0] > 0 ? 1 : 0) + (tet[1] > 0 ? 2 : 0) +
               (tet[2] > 0 ? 4 : 0) + (tet[3] > 0 ? 8 : 0);
      // TetTri0
      if (tetTri0[ci][0] >= 0) {
        int v0 = edges1[tetTri0[ci][0]];
        int v1 = edges1[tetTri0[ci][1]];
        int v2 = edges1[tetTri0[ci][2]];
        if (v0 >= 0 && v1 >= 0 && v2 >= 0 &&
            v0 != v1 && v1 != v2 && v2 != v0) {
          if ((size_t)triCount >= triVerts.len)
            vec_ivec3_resize(&triVerts, triVerts.len * 2);
          triVerts.data[triCount++] = manifold_ivec3(v0, v1, v2);
        }
      }
      // TetTri1
      if (tetTri1[ci][0] >= 0) {
        int v0 = edges1[tetTri1[ci][0]];
        int v1 = edges1[tetTri1[ci][1]];
        int v2 = edges1[tetTri1[ci][2]];
        if (v0 >= 0 && v1 >= 0 && v2 >= 0 &&
            v0 != v1 && v1 != v2 && v2 != v0) {
          if ((size_t)triCount >= triVerts.len)
            vec_ivec3_resize(&triVerts, triVerts.len * 2);
          triVerts.data[triCount++] = manifold_ivec3(v0, v1, v2);
        }
      }

      // Second tetrahedron
      ManifoldIVec4 ti2 = baseIndex;
      if (next == 0) ti2.x += 1;
      else if (next == 1) ti2.y += 1;
      else ti2.z += 1;
      SdfGridVert nextVert2 = *sdf_hash_find_const(&gridVerts,
                                sdf_encode(ti2, gridPow));
      tet[2] = tet[3];
      tet[3] = sdf_gv_neighbor_inside(base, next + 1);

      int edges2[6];
      edges2[0] = base->edgeVerts[0];
      edges2[1] = edges1[5];
      edges2[2] = thisVert.edgeVerts[i + 4];
      edges2[3] = nextVert2.edgeVerts[next + 4];
      edges2[4] = edges1[3];
      edges2[5] = base->edgeVerts[next + 1];
      thisVert = nextVert2;

      ci = (tet[0] > 0 ? 1 : 0) + (tet[1] > 0 ? 2 : 0) +
           (tet[2] > 0 ? 4 : 0) + (tet[3] > 0 ? 8 : 0);
      if (tetTri0[ci][0] >= 0) {
        int v0 = edges2[tetTri0[ci][0]];
        int v1 = edges2[tetTri0[ci][1]];
        int v2 = edges2[tetTri0[ci][2]];
        if (v0 >= 0 && v1 >= 0 && v2 >= 0 &&
            v0 != v1 && v1 != v2 && v2 != v0) {
          if ((size_t)triCount >= triVerts.len)
            vec_ivec3_resize(&triVerts, triVerts.len * 2);
          triVerts.data[triCount++] = manifold_ivec3(v0, v1, v2);
        }
      }
      if (tetTri1[ci][0] >= 0) {
        int v0 = edges2[tetTri1[ci][0]];
        int v1 = edges2[tetTri1[ci][1]];
        int v2 = edges2[tetTri1[ci][2]];
        if (v0 >= 0 && v1 >= 0 && v2 >= 0 &&
            v0 != v1 && v1 != v2 && v2 != v0) {
          if ((size_t)triCount >= triVerts.len)
            vec_ivec3_resize(&triVerts, triVerts.len * 2);
          triVerts.data[triCount++] = manifold_ivec3(v0, v1, v2);
        }
      }

      tet[2] = tet[3];
    }
  }

  free(voxels);

  if (triCount == 0) {
    vec_vec3_free(&vertPos);
    vec_ivec3_free(&triVerts);
    sdf_hash_free(&gridVerts);
    manifold_impl_make_empty(impl, MANIFOLD_ERROR_NO_ERROR);
    return;
  }

  vec_vec3_resize(&vertPos, vertCount);
  vec_ivec3_resize(&triVerts, triCount);
  impl->vertPos = vertPos;

  ManifoldVecIVec3 emptyTriVert = {0};
  manifold_impl_create_halfedges(impl, &triVerts, &emptyTriVert);
  manifold_impl_cleanup_topology(impl);
  manifold_impl_remove_unreferenced_verts(impl);
  manifold_impl_initialize_original(impl);
  manifold_impl_calculate_bbox(impl);
  manifold_impl_set_epsilon(impl, -1.0, false);
  manifold_impl_sort_geometry(impl);
  manifold_impl_set_normals_and_coplanar(impl);

  vec_ivec3_free(&triVerts);
  vec_ivec3_free(&emptyTriVert);
  sdf_hash_free(&gridVerts);
}
