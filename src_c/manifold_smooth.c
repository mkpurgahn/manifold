// Copyright 2021 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// C11 port of smoothing.cpp and subdivision.cpp

#include "manifold_impl.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

// ===================== Helpers =====================

static bool bvec4_get(ManifoldBVec4 v, int i) {
  switch (i) {
    case 0: return v.x;
    case 1: return v.y;
    case 2: return v.z;
    case 3: return v.w;
    default: return false;
  }
}

static ManifoldBVec4 bvec4_not(ManifoldBVec4 v) {
  ManifoldBVec4 r = {!v.x, !v.y, !v.z, !v.w};
  return r;
}

static ManifoldVec3 orthogonal_to(ManifoldVec3 in, ManifoldVec3 altIn,
                                   ManifoldVec3 ref) {
  ManifoldVec3 out = vec3_sub(in, vec3_scale(ref, vec3_dot(in, ref)));
  if (vec3_dot(out, out) < MANIFOLD_PRECISION * vec3_dot(in, in)) {
    out = vec3_sub(altIn, vec3_scale(ref, vec3_dot(altIn, ref)));
  }
  return manifold_safe_normalize(out);
}

static double wrap_angle(double radians) {
  return radians < -MANIFOLD_PI  ? radians + MANIFOLD_TWO_PI
       : radians > MANIFOLD_PI  ? radians - MANIFOLD_TWO_PI
       : radians;
}

static double angle_between(ManifoldVec3 a, ManifoldVec3 b) {
  double d = vec3_dot(a, b);
  return d >= 1.0 ? 0.0 : (d <= -1.0 ? MANIFOLD_PI : acos(d));
}

static ManifoldVec4 circular_tangent(ManifoldVec3 tangent, ManifoldVec3 edgeVec) {
  ManifoldVec3 dir = manifold_safe_normalize(tangent);
  double weight = fmax(0.5, vec3_dot(dir, manifold_safe_normalize(edgeVec)));
  // Quadratic weighted bezier for circular interpolation
  ManifoldVec3 bz2_xyz = vec3_scale(dir, 0.5 * vec3_length(edgeVec));
  ManifoldVec4 bz2 = vec3_to_vec4(bz2_xyz, weight);
  // Equivalent cubic weighted bezier
  ManifoldVec4 bz3 = vec4_lerp(manifold_vec4(0, 0, 0, 1), bz2, 2.0 / 3.0);
  // Convert from homogeneous form to geometric form
  if (bz3.w != 0.0)
    return manifold_vec4(bz3.x / bz3.w, bz3.y / bz3.w, bz3.z / bz3.w, bz3.w);
  return manifold_vec4(bz3.x, bz3.y, bz3.z, bz3.w);
}

// ===================== Smoothing Queries =====================

bool manifold_impl_is_inside_quad(const ManifoldImpl *impl, int halfedge) {
  if (impl->halfedgeTangent.len > 0) {
    return impl->halfedgeTangent.data[halfedge].w < 0;
  }
  int tri = halfedge / 3;
  ManifoldTriRef ref = impl->meshRelation.triRef.data[tri];
  int pair = impl->halfedge.data[halfedge].pairedHalfedge;
  int pairTri = pair / 3;
  ManifoldTriRef pairRef = impl->meshRelation.triRef.data[pairTri];
  if (!manifold_triref_same_face(&ref, &pairRef)) return false;

  int neighbor = manifold_next_halfedge(halfedge);
  if (manifold_triref_same_face(&ref,
          &impl->meshRelation.triRef.data[impl->halfedge.data[neighbor].pairedHalfedge / 3]))
    return false;
  neighbor = manifold_next_halfedge(neighbor);
  if (manifold_triref_same_face(&ref,
          &impl->meshRelation.triRef.data[impl->halfedge.data[neighbor].pairedHalfedge / 3]))
    return false;
  neighbor = manifold_next_halfedge(pair);
  if (manifold_triref_same_face(&pairRef,
          &impl->meshRelation.triRef.data[impl->halfedge.data[neighbor].pairedHalfedge / 3]))
    return false;
  neighbor = manifold_next_halfedge(neighbor);
  if (manifold_triref_same_face(&pairRef,
          &impl->meshRelation.triRef.data[impl->halfedge.data[neighbor].pairedHalfedge / 3]))
    return false;
  return true;
}

bool manifold_impl_is_marked_inside_quad(const ManifoldImpl *impl, int halfedge) {
  return impl->halfedgeTangent.len > 0 && impl->halfedgeTangent.data[halfedge].w < 0;
}

// ===================== Subdivision =====================

// GetNeighbor: returns the tri side index connected to the quad partner, or -1
int manifold_impl_get_neighbor(const ManifoldImpl *impl, int tri) {
  int neighbor = -1;
  for (int i = 0; i < 3; i++) {
    if (manifold_impl_is_marked_inside_quad(impl, 3 * tri + i)) {
      neighbor = (neighbor == -1) ? i : -2;
    }
  }
  return neighbor;
}

// GetHalfedges: returns {h0,h1,h2,-1} for tri, or {h0,h1,h2,h3} for quad
ManifoldIVec4 manifold_impl_get_halfedges(const ManifoldImpl *impl, int tri) {
  ManifoldIVec4 halfedges = manifold_ivec4(-1, -1, -1, -1);
  for (int i = 0; i < 3; i++) halfedges = ivec4_set(halfedges, i, 3 * tri + i);

  int neighbor = manifold_impl_get_neighbor(impl, tri);
  if (neighbor >= 0) {
    int pair = impl->halfedge.data[3 * tri + neighbor].pairedHalfedge;
    if (pair / 3 < tri) return manifold_ivec4(-1, -1, -1, -1);
    int nh2 = manifold_next_halfedge(3 * tri + neighbor);
    int nh3 = manifold_next_halfedge(nh2);
    int np0 = manifold_next_halfedge(pair);
    int np1 = manifold_next_halfedge(np0);
    halfedges = manifold_ivec4(np0, np1, nh2, nh3);
  }
  return halfedges;
}

// GetIndices: returns BaryIndices for a halfedge
ManifoldBaryIndices manifold_impl_get_indices(const ManifoldImpl *impl, int halfedge) {
  int tri = halfedge / 3;
  int idx = halfedge % 3;
  int neighbor = manifold_impl_get_neighbor(impl, tri);
  if (idx == neighbor) {
    ManifoldBaryIndices r = {-1, -1, -1};
    return r;
  }
  if (neighbor < 0) {
    ManifoldBaryIndices r = {tri, idx, manifold_next3(idx)};
    return r;
  }
  int pair = impl->halfedge.data[3 * tri + neighbor].pairedHalfedge;
  if (pair / 3 < tri) {
    tri = pair / 3;
    idx = (manifold_next3(neighbor) == idx) ? 0 : 1;
  } else {
    idx = (manifold_next3(neighbor) == idx) ? 2 : 3;
  }
  ManifoldBaryIndices r = {tri, idx, (idx + 1) % 4};
  return r;
}

// ===================== Partition =====================

typedef struct {
  ManifoldIVec4 idx;
  ManifoldIVec4 sortedDivisions;
  ManifoldVecVec4 vertBary;
  ManifoldVecIVec3 triVert;
} Partition;

static void partition_init(Partition *p) {
  p->idx = manifold_ivec4(0, 0, 0, 0);
  p->sortedDivisions = manifold_ivec4(0, 0, 0, 0);
  p->vertBary = (ManifoldVecVec4)MANIFOLD_VEC_INIT;
  p->triVert = (ManifoldVecIVec3)MANIFOLD_VEC_INIT;
}

static void partition_free(Partition *p) {
  vec_vec4_free(&p->vertBary);
  vec_ivec3_free(&p->triVert);
}

static int partition_interior_offset(const Partition *p) {
  return ivec4_get(p->sortedDivisions, 0) + ivec4_get(p->sortedDivisions, 1) +
         ivec4_get(p->sortedDivisions, 2) + ivec4_get(p->sortedDivisions, 3);
}

static int partition_num_interior(const Partition *p) {
  return (int)p->vertBary.len - partition_interior_offset(p);
}

static void partition_fan(ManifoldVecIVec3 *triVert, ManifoldIVec3 cornerVerts,
                          int added, int edgeOffset) {
  int last = cornerVerts.x;
  for (int i = 0; i < added; i++) {
    int next = edgeOffset + i;
    vec_ivec3_push(triVert, manifold_ivec3(last, next, cornerVerts.z));
    last = next;
  }
  vec_ivec3_push(triVert, manifold_ivec3(last, cornerVerts.y, cornerVerts.z));
}

static int get_edge_vert(ManifoldIVec4 edgeOffsets, ManifoldBVec4 edgeFwd,
                         int edge, int idx) {
  bool fwd;
  int off;
  switch (edge) {
    case 0: off = edgeOffsets.x; fwd = edgeFwd.x; break;
    case 1: off = edgeOffsets.y; fwd = edgeFwd.y; break;
    case 2: off = edgeOffsets.z; fwd = edgeFwd.z; break;
    case 3: off = edgeOffsets.w; fwd = edgeFwd.w; break;
    default: off = 0; fwd = true; break;
  }
  return off + (fwd ? 1 : -1) * idx;
}

static void partition_quad(ManifoldVecIVec3 *triVert, ManifoldVecVec4 *vertBary,
                           ManifoldIVec4 cornerVerts, ManifoldIVec4 edgeOffsets,
                           ManifoldIVec4 edgeAdded, ManifoldBVec4 edgeFwd) {
  // Check that all edgeAdded >= 0
  for (int i = 0; i < 4; i++) {
    if (ivec4_get(edgeAdded, i) < 0) return;
  }

  int corner = -1;
  int last = 3;
  int maxEdge = -1;
  for (int i = 0; i < 4; i++) {
    if (corner == -1 && ivec4_get(edgeAdded, i) == 0 && ivec4_get(edgeAdded, last) == 0) {
      corner = i;
    }
    if (ivec4_get(edgeAdded, i) > 0) {
      maxEdge = (maxEdge == -1) ? i : -2;
    }
    last = i;
  }

  if (corner >= 0) {  // terminate
    if (maxEdge >= 0) {
      int edge[4];
      for (int i = 0; i < 4; i++) edge[i] = (maxEdge + i) % 4;
      int middle = ivec4_get(edgeAdded, maxEdge) / 2;
      vec_ivec3_push(triVert, manifold_ivec3(
          ivec4_get(cornerVerts, edge[2]), ivec4_get(cornerVerts, edge[3]),
          get_edge_vert(edgeOffsets, edgeFwd, maxEdge, middle)));
      int lst = ivec4_get(cornerVerts, edge[0]);
      for (int i = 0; i <= middle; i++) {
        int next = get_edge_vert(edgeOffsets, edgeFwd, maxEdge, i);
        vec_ivec3_push(triVert, manifold_ivec3(
            ivec4_get(cornerVerts, edge[3]), lst, next));
        lst = next;
      }
      lst = ivec4_get(cornerVerts, edge[1]);
      for (int i = ivec4_get(edgeAdded, maxEdge) - 1; i >= middle; i--) {
        int next = get_edge_vert(edgeOffsets, edgeFwd, maxEdge, i);
        vec_ivec3_push(triVert, manifold_ivec3(
            ivec4_get(cornerVerts, edge[2]), next, lst));
        lst = next;
      }
    } else {
      int sideVert = ivec4_get(cornerVerts, 0);
      for (int j = 1; j <= 2; j++) {
        int side = (corner + j) % 4;
        if (j == 2 && ivec4_get(edgeAdded, side) > 0) {
          vec_ivec3_push(triVert, manifold_ivec3(
              ivec4_get(cornerVerts, side),
              get_edge_vert(edgeOffsets, edgeFwd, side, 0), sideVert));
        } else {
          sideVert = ivec4_get(cornerVerts, side);
        }
        for (int i = 0; i < ivec4_get(edgeAdded, side); i++) {
          int nextVert = get_edge_vert(edgeOffsets, edgeFwd, side, i);
          vec_ivec3_push(triVert, manifold_ivec3(
              ivec4_get(cornerVerts, corner), sideVert, nextVert));
          sideVert = nextVert;
        }
        if (j == 2 || ivec4_get(edgeAdded, side) == 0) {
          vec_ivec3_push(triVert, manifold_ivec3(
              ivec4_get(cornerVerts, corner), sideVert,
              ivec4_get(cornerVerts, (corner + j + 1) % 4)));
        }
      }
    }
    return;
  }

  // recursively partition
  int ea1 = ivec4_get(edgeAdded, 1);
  int ea3 = ivec4_get(edgeAdded, 3);
  int partitions = 1 + (ea1 < ea3 ? ea1 : ea3);
  ManifoldIVec4 newCV = manifold_ivec4(ivec4_get(cornerVerts, 1), -1, -1,
                                        ivec4_get(cornerVerts, 0));
  ManifoldIVec4 newEO = manifold_ivec4(ivec4_get(edgeOffsets, 1), -1,
      get_edge_vert(edgeOffsets, edgeFwd, 3, ea3 + 1),
      ivec4_get(edgeOffsets, 0));
  ManifoldIVec4 newEA = manifold_ivec4(0, -1, 0, ivec4_get(edgeAdded, 0));
  ManifoldBVec4 newEF = {bvec4_get(edgeFwd, 1), true, bvec4_get(edgeFwd, 3),
                         bvec4_get(edgeFwd, 0)};

  for (int i = 1; i < partitions; i++) {
    int cornerOffset1 = (ea1 * i) / partitions;
    int cornerOffset3 = ea3 - 1 - (ea3 * i) / partitions;
    int nextOffset1 = get_edge_vert(edgeOffsets, edgeFwd, 1, cornerOffset1 + 1);
    int nextOffset3 = get_edge_vert(edgeOffsets, edgeFwd, 3, cornerOffset3 + 1);
    double ea0d = (double)ivec4_get(edgeAdded, 0);
    double ea2d = (double)ivec4_get(edgeAdded, 2);
    int added = (int)(0.5 + double_lerp(ea0d, ea2d, (double)i / partitions));

    newCV = ivec4_set(newCV, 1, get_edge_vert(edgeOffsets, edgeFwd, 1, cornerOffset1));
    newCV = ivec4_set(newCV, 2, get_edge_vert(edgeOffsets, edgeFwd, 3, cornerOffset3));
    newEA = ivec4_set(newEA, 0, abs(nextOffset1 - ivec4_get(newEO, 0)) - 1);
    newEA = ivec4_set(newEA, 1, added);
    newEA = ivec4_set(newEA, 2, abs(nextOffset3 - ivec4_get(newEO, 2)) - 1);
    newEO = ivec4_set(newEO, 1, (int)vertBary->len);
    newEO = ivec4_set(newEO, 2, nextOffset3);

    for (int j = 0; j < added; j++) {
      ManifoldVec4 b = vec4_lerp(vertBary->data[ivec4_get(newCV, 1)],
                                  vertBary->data[ivec4_get(newCV, 2)],
                                  (j + 1.0) / (added + 1.0));
      vec_vec4_push(vertBary, b);
    }

    partition_quad(triVert, vertBary, newCV, newEO, newEA, newEF);

    newCV = ivec4_set(newCV, 0, ivec4_get(newCV, 1));
    newCV = ivec4_set(newCV, 3, ivec4_get(newCV, 2));
    newEA = ivec4_set(newEA, 3, ivec4_get(newEA, 1));
    newEO = ivec4_set(newEO, 0, nextOffset1);
    newEO = ivec4_set(newEO, 3, ivec4_get(newEO, 1) + ivec4_get(newEA, 1) - 1);
    newEF.w = false;
  }

  newCV = ivec4_set(newCV, 1, ivec4_get(cornerVerts, 2));
  newCV = ivec4_set(newCV, 2, ivec4_get(cornerVerts, 3));
  newEO = ivec4_set(newEO, 1, ivec4_get(edgeOffsets, 2));
  newEA = ivec4_set(newEA, 0,
      ea1 - abs(ivec4_get(newEO, 0) - ivec4_get(edgeOffsets, 1)));
  newEA = ivec4_set(newEA, 1, ivec4_get(edgeAdded, 2));
  newEA = ivec4_set(newEA, 2,
      abs(ivec4_get(newEO, 2) - ivec4_get(edgeOffsets, 3)) - 1);
  newEO = ivec4_set(newEO, 2, ivec4_get(edgeOffsets, 3));
  newEF.y = bvec4_get(edgeFwd, 2);

  partition_quad(triVert, vertBary, newCV, newEO, newEA, newEF);
}

static Partition get_cached_partition(ManifoldIVec4 n) {
  Partition partition;
  partition_init(&partition);
  partition.sortedDivisions = n;

  if (n.w > 0) {  // quad
    vec_vec4_push(&partition.vertBary, manifold_vec4(1, 0, 0, 0));
    vec_vec4_push(&partition.vertBary, manifold_vec4(0, 1, 0, 0));
    vec_vec4_push(&partition.vertBary, manifold_vec4(0, 0, 1, 0));
    vec_vec4_push(&partition.vertBary, manifold_vec4(0, 0, 0, 1));
    ManifoldIVec4 edgeOffsets;
    edgeOffsets = ivec4_set(manifold_ivec4(0,0,0,0), 0, 4);
    for (int i = 0; i < 4; i++) {
      if (i > 0) {
        edgeOffsets = ivec4_set(edgeOffsets, i,
            ivec4_get(edgeOffsets, i - 1) + ivec4_get(n, i - 1) - 1);
      }
      ManifoldVec4 nextBary = partition.vertBary.data[(i + 1) % 4];
      for (int j = 1; j < ivec4_get(n, i); j++) {
        vec_vec4_push(&partition.vertBary,
            vec4_lerp(partition.vertBary.data[i], nextBary,
                      (double)j / ivec4_get(n, i)));
      }
    }
    ManifoldIVec4 nMinus1 = manifold_ivec4(n.x - 1, n.y - 1, n.z - 1, n.w - 1);
    ManifoldBVec4 allTrue = {true, true, true, true};
    partition_quad(&partition.triVert, &partition.vertBary,
                   manifold_ivec4(0, 1, 2, 3), edgeOffsets, nMinus1, allTrue);
  } else {  // tri
    vec_vec4_push(&partition.vertBary, manifold_vec4(1, 0, 0, 0));
    vec_vec4_push(&partition.vertBary, manifold_vec4(0, 1, 0, 0));
    vec_vec4_push(&partition.vertBary, manifold_vec4(0, 0, 1, 0));
    for (int i = 0; i < 3; i++) {
      ManifoldVec4 nextBary = partition.vertBary.data[(i + 1) % 3];
      for (int j = 1; j < ivec4_get(n, i); j++) {
        vec_vec4_push(&partition.vertBary,
            vec4_lerp(partition.vertBary.data[i], nextBary,
                      (double)j / ivec4_get(n, i)));
      }
    }
    ManifoldIVec3 edgeOffsets = manifold_ivec3(3, 3 + n.x - 1, 3 + n.x - 1 + n.y - 1);

    double f = (double)(n.z * n.z + n.x * n.x);
    if (n.y == 1) {
      if (n.x == 1) {
        vec_ivec3_push(&partition.triVert, manifold_ivec3(0, 1, 2));
      } else {
        partition_fan(&partition.triVert, manifold_ivec3(0, 1, 2),
                      n.x - 1, edgeOffsets.x);
      }
    } else if ((double)(n.y * n.y) > f - sqrt(2.0) * n.x * n.z) {  // acute-ish
      vec_ivec3_push(&partition.triVert,
          manifold_ivec3(edgeOffsets.y - 1, 1, edgeOffsets.y));
      ManifoldBVec4 allTrue = {true, true, true, true};
      partition_quad(&partition.triVert, &partition.vertBary,
          manifold_ivec4(edgeOffsets.y - 1, edgeOffsets.y, 2, 0),
          manifold_ivec4(-1, edgeOffsets.y + 1, edgeOffsets.z, edgeOffsets.x),
          manifold_ivec4(0, n.y - 2, n.z - 1, n.x - 2),
          allTrue);
    } else {  // obtuse -> split into two acute
      int ns = n.x - 2;
      int ns2 = (int)(0.5 + (f - (double)(n.y * n.y)) / (2.0 * n.x));
      if (ns2 < ns) ns = ns2;
      int nh = (int)(0.5 + sqrt(fmax(1.0, (double)(n.z * n.z) - (double)(ns * ns))));
      if (nh < 1) nh = 1;

      int hOffset = (int)partition.vertBary.len;
      ManifoldVec4 middleBary = partition.vertBary.data[edgeOffsets.x + ns - 1];
      for (int j = 1; j < nh; j++) {
        vec_vec4_push(&partition.vertBary,
            vec4_lerp(partition.vertBary.data[2], middleBary, (double)j / nh));
      }

      vec_ivec3_push(&partition.triVert,
          manifold_ivec3(edgeOffsets.y - 1, 1, edgeOffsets.y));
      ManifoldBVec4 allTrue = {true, true, true, true};
      partition_quad(&partition.triVert, &partition.vertBary,
          manifold_ivec4(edgeOffsets.y - 1, edgeOffsets.y, 2, edgeOffsets.x + ns - 1),
          manifold_ivec4(-1, edgeOffsets.y + 1, hOffset, edgeOffsets.x + ns),
          manifold_ivec4(0, n.y - 2, nh - 1, n.x - ns - 2),
          allTrue);

      if (n.z == 1) {
        partition_fan(&partition.triVert,
            manifold_ivec3(0, edgeOffsets.x + ns - 1, 2), ns - 1, edgeOffsets.x);
      } else {
        if (ns == 1) {
          vec_ivec3_push(&partition.triVert,
              manifold_ivec3(hOffset, 2, edgeOffsets.z));
          ManifoldBVec4 ef = {true, true, true, false};
          partition_quad(&partition.triVert, &partition.vertBary,
              manifold_ivec4(hOffset, edgeOffsets.z, 0, edgeOffsets.x),
              manifold_ivec4(-1, edgeOffsets.z + 1, -1, hOffset + nh - 2),
              manifold_ivec4(0, n.z - 2, ns - 1, nh - 2),
              ef);
        } else {
          vec_ivec3_push(&partition.triVert,
              manifold_ivec3(hOffset - 1, 0, edgeOffsets.x));
          ManifoldBVec4 ef = {true, true, false, true};
          partition_quad(&partition.triVert, &partition.vertBary,
              manifold_ivec4(hOffset - 1, edgeOffsets.x, edgeOffsets.x + ns - 1, 2),
              manifold_ivec4(-1, edgeOffsets.x + 1, hOffset + nh - 2, edgeOffsets.z),
              manifold_ivec4(0, ns - 2, nh - 1, n.z - 2),
              ef);
        }
      }
    }
  }
  return partition;
}

static Partition get_partition(ManifoldIVec4 divisions) {
  if (divisions.x == 0) {
    Partition p;
    partition_init(&p);
    return p;
  }

  ManifoldIVec4 sortedDiv = divisions;
  ManifoldIVec4 triIdx = manifold_ivec4(0, 1, 2, 3);
  if (divisions.w == 0) {  // triangle
    // Sort: sortedDiv[0] >= sortedDiv[1] >= sortedDiv[2]
    if (ivec4_get(sortedDiv, 2) > ivec4_get(sortedDiv, 1)) {
      int tmp = sortedDiv.z; sortedDiv.z = sortedDiv.y; sortedDiv.y = tmp;
      tmp = triIdx.z; triIdx.z = triIdx.y; triIdx.y = tmp;
    }
    if (ivec4_get(sortedDiv, 1) > ivec4_get(sortedDiv, 0)) {
      int tmp = sortedDiv.y; sortedDiv.y = sortedDiv.x; sortedDiv.x = tmp;
      tmp = triIdx.y; triIdx.y = triIdx.x; triIdx.x = tmp;
      if (ivec4_get(sortedDiv, 2) > ivec4_get(sortedDiv, 1)) {
        tmp = sortedDiv.z; sortedDiv.z = sortedDiv.y; sortedDiv.y = tmp;
        tmp = triIdx.z; triIdx.z = triIdx.y; triIdx.y = tmp;
      }
    }
  } else {  // quad
    int minIdx = 0;
    int min = ivec4_get(divisions, 0);
    int next = ivec4_get(divisions, 1);
    for (int i = 1; i < 4; i++) {
      int n = ivec4_get(divisions, (i + 1) % 4);
      if (ivec4_get(divisions, i) < min ||
          (ivec4_get(divisions, i) == min && n < next)) {
        minIdx = i;
        min = ivec4_get(divisions, i);
        next = n;
      }
    }
    ManifoldIVec4 tmp = sortedDiv;
    for (int i = 0; i < 4; i++) {
      triIdx = ivec4_set(triIdx, i, (i + minIdx) % 4);
      sortedDiv = ivec4_set(sortedDiv, i, ivec4_get(tmp, (i + minIdx) % 4));
    }
  }

  Partition partition = get_cached_partition(sortedDiv);
  partition.idx = triIdx;
  return partition;
}

static ManifoldVecIVec3 partition_reindex(const Partition *part,
    ManifoldIVec4 triVerts, ManifoldIVec4 edgeOffsets, ManifoldBVec4 edgeFwd,
    int interiorOffset) {

  // Allocate newVerts array
  int capacity = (int)part->vertBary.len;
  int *newVerts = (int *)malloc(capacity * sizeof(int));
  int nv = 0;

  ManifoldIVec4 triIdx = part->idx;
  ManifoldIVec4 outTri = manifold_ivec4(0, 1, 2, 3);
  if (ivec4_get(triVerts, 3) < 0 &&
      ivec4_get(triIdx, 1) != manifold_next3(ivec4_get(triIdx, 0))) {
    triIdx = manifold_ivec4(ivec4_get(part->idx, 2), ivec4_get(part->idx, 0),
                            ivec4_get(part->idx, 1), ivec4_get(part->idx, 3));
    edgeFwd = bvec4_not(edgeFwd);
    outTri = manifold_ivec4(1, 0, 2, 3);
  }
  for (int i = 0; i < 4; i++) {
    if (ivec4_get(triVerts, ivec4_get(triIdx, i)) >= 0)
      newVerts[nv++] = ivec4_get(triVerts, ivec4_get(triIdx, i));
  }
  for (int i = 0; i < 4; i++) {
    int sd = ivec4_get(part->sortedDivisions, i);
    int n = sd - 1;
    int idx_i = ivec4_get(part->idx, i);
    int off = ivec4_get(edgeOffsets, idx_i) +
              (bvec4_get(edgeFwd, idx_i) ? 0 : n - 1);
    for (int j = 0; j < n; j++) {
      newVerts[nv++] = off;
      off += bvec4_get(edgeFwd, idx_i) ? 1 : -1;
    }
  }
  int offset = interiorOffset - nv;
  for (int i = nv; i < capacity; i++) {
    newVerts[i] = i + offset;
  }

  int numTri = (int)part->triVert.len;
  ManifoldVecIVec3 result = MANIFOLD_VEC_INIT;
  vec_ivec3_resize(&result, numTri);
  for (int tri = 0; tri < numTri; tri++) {
    ManifoldIVec3 tv = part->triVert.data[tri];
    ManifoldIVec3 out = manifold_ivec3(0, 0, 0);
    out = ivec3_set(out, ivec4_get(outTri, 0), newVerts[tv.x]);
    out = ivec3_set(out, ivec4_get(outTri, 1), newVerts[tv.y]);
    out = ivec3_set(out, ivec4_get(outTri, 2), newVerts[tv.z]);
    result.data[tri] = out;
  }

  free(newVerts);
  return result;
}

// ===================== CreateTmpEdges =====================

static ManifoldVecTmpEdge create_tmp_edges(const ManifoldImpl *impl) {
  size_t numHE = impl->halfedge.len;
  ManifoldVecTmpEdge edges = MANIFOLD_VEC_INIT;
  for (size_t i = 0; i < numHE; i++) {
    ManifoldHalfedge h = impl->halfedge.data[i];
    if (manifold_halfedge_is_forward(&h)) {
      ManifoldTmpEdge e = manifold_tmpedge(h.startVert, h.endVert, (int)i);
      vec_tmpedge_push(&edges, e);
    }
  }
  return edges;
}

// ===================== FillRetainedVerts =====================

static void fill_retained_verts(const ManifoldImpl *impl,
                                 ManifoldVecBarycentric *vertBary) {
  int numTri = (int)(impl->halfedge.len / 3);
  for (int tri = 0; tri < numTri; tri++) {
    for (int i = 0; i < 3; i++) {
      ManifoldBaryIndices indices = manifold_impl_get_indices(impl, 3 * tri + i);
      if (indices.start4 < 0) continue;
      ManifoldVec4 uvw = manifold_vec4(0, 0, 0, 0);
      uvw = vec4_set(uvw, indices.start4, 1.0);
      int sv = impl->halfedge.data[3 * tri + i].startVert;
      if (sv >= 0 && sv < (int)vertBary->len) {
        vertBary->data[sv].tri = indices.tri;
        vertBary->data[sv].uvw = uvw;
      }
    }
  }
}

// ===================== Subdivide =====================

// Edge division callback type
typedef int (*EdgeDivisionsFn)(ManifoldVec3 edgeVec, ManifoldVec4 tangent0,
                                ManifoldVec4 tangent1, void *ctx);

ManifoldVecBarycentric manifold_impl_subdivide(
    ManifoldImpl *impl, EdgeDivisionsFn edgeDivisions, void *ctx,
    bool keepInterior) {

  ManifoldVecBarycentric emptyBary = MANIFOLD_VEC_INIT;
  if (manifold_impl_is_empty(impl)) return emptyBary;

  ManifoldVecTmpEdge edges = create_tmp_edges(impl);
  int numVert = (int)impl->vertPos.len;
  int numEdge = (int)edges.len;
  int numTri = (int)(impl->halfedge.len / 3);

  // half2Edge: map halfedge → edge index
  int *half2Edge = (int *)calloc(impl->halfedge.len, sizeof(int));
  for (int i = 0; i < numEdge; i++) {
    int idx = edges.data[i].halfedgeIdx;
    half2Edge[idx] = i;
    half2Edge[impl->halfedge.data[idx].pairedHalfedge] = i;
  }

  // faceHalfedges
  ManifoldIVec4 *faceHalfedges = (ManifoldIVec4 *)malloc(numTri * sizeof(ManifoldIVec4));
  for (int tri = 0; tri < numTri; tri++) {
    faceHalfedges[tri] = manifold_impl_get_halfedges(impl, tri);
  }

  // edgeAdded: how many verts to add on each edge
  int *edgeAdded = (int *)calloc(numEdge, sizeof(int));
  for (int i = 0; i < numEdge; i++) {
    ManifoldTmpEdge edge = edges.data[i];
    int hIdx = edge.halfedgeIdx;
    if (manifold_impl_is_marked_inside_quad(impl, hIdx)) {
      edgeAdded[i] = 0;
      continue;
    }
    ManifoldVec3 vec = vec3_sub(impl->vertPos.data[edge.first],
                                 impl->vertPos.data[edge.second]);
    ManifoldVec4 t0 = impl->halfedgeTangent.len > 0
                          ? impl->halfedgeTangent.data[hIdx]
                          : manifold_vec4(0, 0, 0, 0);
    ManifoldVec4 t1 = impl->halfedgeTangent.len > 0
                          ? impl->halfedgeTangent.data[impl->halfedge.data[hIdx].pairedHalfedge]
                          : manifold_vec4(0, 0, 0, 0);
    edgeAdded[i] = edgeDivisions(vec, t0, t1, ctx);
  }

  if (keepInterior) {
    int *tmp = (int *)malloc(numEdge * sizeof(int));
    memcpy(tmp, edgeAdded, numEdge * sizeof(int));
    for (int i = 0; i < numEdge; i++) {
      ManifoldTmpEdge edge = edges.data[i];
      int hIdx = edge.halfedgeIdx;
      if (manifold_impl_is_marked_inside_quad(impl, hIdx)) continue;

      int thisAdded = tmp[i];
      int maxExtra = 0;
      for (int side = 0; side < 2; side++) {
        int h = side == 0 ? hIdx : impl->halfedge.data[hIdx].pairedHalfedge;
        int longest = 0, total = 0;
        bool skip = false;
        for (int s = 0; s < 3; s++) {
          int added = edgeAdded[half2Edge[h]];
          if (added > longest) longest = added;
          total += added;
          h = manifold_next_halfedge(h);
          if (manifold_impl_is_marked_inside_quad(impl, h)) {
            skip = true;
            break;
          }
        }
        if (skip) continue;
        int minExtra = (int)(longest * 0.2) + 1;
        int extra = 2 * longest + minExtra - total;
        if (extra > 0 && longest > 0) {
          int e = (extra * (longest - thisAdded)) / longest;
          if (e > maxExtra) maxExtra = e;
        }
      }
      tmp[i] += maxExtra;
    }
    memcpy(edgeAdded, tmp, numEdge * sizeof(int));
    free(tmp);
  }

  // edgeOffset: exclusive scan of edgeAdded, starting from numVert
  int *edgeOffset = (int *)malloc(numEdge * sizeof(int));
  {
    int sum = numVert;
    for (int i = 0; i < numEdge; i++) {
      edgeOffset[i] = sum;
      sum += edgeAdded[i];
    }
  }

  int totalVerts = (numEdge > 0) ? edgeOffset[numEdge - 1] + edgeAdded[numEdge - 1] : numVert;

  ManifoldVecBarycentric vertBary = MANIFOLD_VEC_INIT;
  vec_bary_resize(&vertBary, totalVerts);
  int totalEdgeAdded = totalVerts - numVert;

  fill_retained_verts(impl, &vertBary);

  // Fill edge barycentrics
  for (int i = 0; i < numEdge; i++) {
    int n = edgeAdded[i];
    int offset = edgeOffset[i];
    ManifoldBaryIndices indices = manifold_impl_get_indices(impl, edges.data[i].halfedgeIdx);
    if (indices.tri < 0) continue;
    double frac = 1.0 / (n + 1);
    for (int j = 0; j < n; j++) {
      ManifoldVec4 uvw = manifold_vec4(0, 0, 0, 0);
      uvw = vec4_set(uvw, indices.end4, (j + 1) * frac);
      uvw = vec4_set(uvw, indices.start4, 1.0 - vec4_get(uvw, indices.end4));
      vertBary.data[offset + j].uvw = uvw;
      vertBary.data[offset + j].tri = indices.tri;
    }
  }

  // Build sub-triangulations for each tri
  Partition *subTris = (Partition *)malloc(numTri * sizeof(Partition));
  for (int tri = 0; tri < numTri; tri++) {
    ManifoldIVec4 he = faceHalfedges[tri];
    ManifoldIVec4 divisions = manifold_ivec4(0, 0, 0, 0);
    for (int i = 0; i < 4; i++) {
      if (ivec4_get(he, i) >= 0) {
        divisions = ivec4_set(divisions, i, edgeAdded[half2Edge[ivec4_get(he, i)]] + 1);
      }
    }
    subTris[tri] = get_partition(divisions);
  }

  // triOffset: exclusive scan of subTri sizes
  int *triOffset = (int *)malloc(numTri * sizeof(int));
  {
    int sum = 0;
    for (int tri = 0; tri < numTri; tri++) {
      triOffset[tri] = sum;
      sum += (int)subTris[tri].triVert.len;
    }
  }

  // interiorOffset: exclusive scan of interior vert counts
  int *interiorOffset = (int *)malloc(numTri * sizeof(int));
  {
    int sum = (int)vertBary.len;
    for (int tri = 0; tri < numTri; tri++) {
      interiorOffset[tri] = sum;
      sum += partition_num_interior(&subTris[tri]);
    }
  }

  int totalNewTris = triOffset[numTri - 1] + (int)subTris[numTri - 1].triVert.len;
  int totalNewVerts = interiorOffset[numTri - 1] + partition_num_interior(&subTris[numTri - 1]);

  ManifoldVecIVec3 triVerts = MANIFOLD_VEC_INIT;
  vec_ivec3_resize(&triVerts, totalNewTris);
  vec_bary_resize(&vertBary, totalNewVerts);
  ManifoldVecTriRef triRef = MANIFOLD_VEC_INIT;
  vec_triref_resize(&triRef, totalNewTris);
  ManifoldVecVec3 faceNormal = MANIFOLD_VEC_INIT;
  vec_vec3_resize(&faceNormal, totalNewTris);

  for (int tri = 0; tri < numTri; tri++) {
    ManifoldIVec4 he = faceHalfedges[tri];
    if (ivec4_get(he, 0) < 0) continue;

    ManifoldIVec4 tri3 = manifold_ivec4(0, 0, 0, 0);
    ManifoldIVec4 eo = manifold_ivec4(0, 0, 0, 0);
    ManifoldBVec4 ef = {false, false, false, false};
    for (int i = 0; i < 4; i++) {
      if (ivec4_get(he, i) < 0) {
        tri3 = ivec4_set(tri3, i, -1);
        continue;
      }
      ManifoldHalfedge *h = &impl->halfedge.data[ivec4_get(he, i)];
      tri3 = ivec4_set(tri3, i, h->startVert);
      eo = ivec4_set(eo, i, edgeOffset[half2Edge[ivec4_get(he, i)]]);
      bool fwd = manifold_halfedge_is_forward(h);
      switch (i) {
        case 0: ef.x = fwd; break;
        case 1: ef.y = fwd; break;
        case 2: ef.z = fwd; break;
        case 3: ef.w = fwd; break;
      }
    }

    ManifoldVecIVec3 newTris = partition_reindex(&subTris[tri], tri3, eo, ef,
                                                  interiorOffset[tri]);
    memcpy(&triVerts.data[triOffset[tri]], newTris.data,
           newTris.len * sizeof(ManifoldIVec3));
    for (size_t j = 0; j < newTris.len; j++) {
      triRef.data[triOffset[tri] + j] = impl->meshRelation.triRef.data[tri];
      faceNormal.data[triOffset[tri] + j] = impl->faceNormal.data[tri];
    }
    vec_ivec3_free(&newTris);

    // Fill interior barycentrics
    ManifoldIVec4 idx = subTris[tri].idx;
    ManifoldIVec4 vIdx;
    if (ivec4_get(he, 3) >= 0 || ivec4_get(idx, 1) == manifold_next3(ivec4_get(idx, 0))) {
      vIdx = idx;
    } else {
      vIdx = manifold_ivec4(ivec4_get(idx, 2), ivec4_get(idx, 0),
                            ivec4_get(idx, 1), ivec4_get(idx, 3));
    }
    ManifoldIVec4 rIdx = manifold_ivec4(0, 0, 0, 0);
    for (int i = 0; i < 4; i++) {
      rIdx = ivec4_set(rIdx, ivec4_get(vIdx, i), i);
    }

    int intOff = partition_interior_offset(&subTris[tri]);
    int numInt = partition_num_interior(&subTris[tri]);
    for (int j = 0; j < numInt; j++) {
      ManifoldVec4 bary = subTris[tri].vertBary.data[intOff + j];
      ManifoldVec4 rBary = manifold_vec4(
          vec4_get(bary, ivec4_get(rIdx, 0)),
          vec4_get(bary, ivec4_get(rIdx, 1)),
          vec4_get(bary, ivec4_get(rIdx, 2)),
          vec4_get(bary, ivec4_get(rIdx, 3)));
      vertBary.data[interiorOffset[tri] + j].tri = tri;
      vertBary.data[interiorOffset[tri] + j].uvw = rBary;
    }
  }

  // Replace mesh data
  vec_triref_free(&impl->meshRelation.triRef);
  impl->meshRelation.triRef = triRef;
  vec_vec3_free(&impl->faceNormal);
  impl->faceNormal = faceNormal;

  // Calculate new vertex positions
  ManifoldVecVec3 newVertPos = MANIFOLD_VEC_INIT;
  vec_vec3_resize(&newVertPos, vertBary.len);
  for (size_t vert = 0; vert < vertBary.len; vert++) {
    ManifoldBarycentric bary = vertBary.data[vert];
    ManifoldIVec4 he = faceHalfedges[bary.tri];
    if (ivec4_get(he, 3) < 0) {
      // triangle
      ManifoldVec3 p0 = impl->vertPos.data[impl->halfedge.data[ivec4_get(he, 0)].startVert];
      ManifoldVec3 p1 = impl->vertPos.data[impl->halfedge.data[ivec4_get(he, 1)].startVert];
      ManifoldVec3 p2 = impl->vertPos.data[impl->halfedge.data[ivec4_get(he, 2)].startVert];
      newVertPos.data[vert] = vec3_add(vec3_add(
          vec3_scale(p0, bary.uvw.x), vec3_scale(p1, bary.uvw.y)),
          vec3_scale(p2, bary.uvw.z));
    } else {
      // quad
      ManifoldVec3 p0 = impl->vertPos.data[impl->halfedge.data[ivec4_get(he, 0)].startVert];
      ManifoldVec3 p1 = impl->vertPos.data[impl->halfedge.data[ivec4_get(he, 1)].startVert];
      ManifoldVec3 p2 = impl->vertPos.data[impl->halfedge.data[ivec4_get(he, 2)].startVert];
      ManifoldVec3 p3 = impl->vertPos.data[impl->halfedge.data[ivec4_get(he, 3)].startVert];
      newVertPos.data[vert] = vec3_add(vec3_add(vec3_add(
          vec3_scale(p0, bary.uvw.x), vec3_scale(p1, bary.uvw.y)),
          vec3_scale(p2, bary.uvw.z)), vec3_scale(p3, bary.uvw.w));
    }
  }
  vec_vec3_free(&impl->vertPos);
  impl->vertPos = newVertPos;

  // Handle properties
  if (impl->numProp > 0) {
    int numPropVert = (int)manifold_impl_num_prop_vert(impl);
    int addedVerts = (int)impl->vertPos.len - numVert;
    int propOffset = numPropVert - numVert;
    int numProp = impl->numProp;

    ManifoldVecDouble prop = MANIFOLD_VEC_INIT;
    vec_double_resize(&prop, numProp * (numPropVert + addedVerts + totalEdgeAdded));

    // Copy retained prop verts
    memcpy(prop.data, impl->properties.data, impl->properties.len * sizeof(double));

    // Copy interior and forward edge prop verts
    for (int i = 0; i < addedVerts; i++) {
      int vert = numPropVert + i;
      ManifoldBarycentric bary = vertBary.data[numVert + i];
      ManifoldIVec4 he = faceHalfedges[bary.tri];

      for (int p2 = 0; p2 < numProp; p2++) {
        if (ivec4_get(he, 3) < 0) {
          double v0 = impl->properties.data[impl->halfedge.data[3 * bary.tri + 0].propVert * numProp + p2];
          double v1 = impl->properties.data[impl->halfedge.data[3 * bary.tri + 1].propVert * numProp + p2];
          double v2 = impl->properties.data[impl->halfedge.data[3 * bary.tri + 2].propVert * numProp + p2];
          prop.data[vert * numProp + p2] = v0 * bary.uvw.x + v1 * bary.uvw.y + v2 * bary.uvw.z;
        } else {
          double v0 = impl->properties.data[impl->halfedge.data[ivec4_get(he, 0)].propVert * numProp + p2];
          double v1 = impl->properties.data[impl->halfedge.data[ivec4_get(he, 1)].propVert * numProp + p2];
          double v2 = impl->properties.data[impl->halfedge.data[ivec4_get(he, 2)].propVert * numProp + p2];
          double v3 = impl->properties.data[impl->halfedge.data[ivec4_get(he, 3)].propVert * numProp + p2];
          prop.data[vert * numProp + p2] = v0 * bary.uvw.x + v1 * bary.uvw.y +
                                            v2 * bary.uvw.z + v3 * bary.uvw.w;
        }
      }
    }

    // Copy backward edge prop verts
    for (int i = 0; i < numEdge; i++) {
      int n = edgeAdded[i];
      int offset = edgeOffset[i] + propOffset + addedVerts;
      double frac = 1.0 / (n + 1);
      int halfedgeIdx = impl->halfedge.data[edges.data[i].halfedgeIdx].pairedHalfedge;
      int prop0 = impl->halfedge.data[halfedgeIdx].propVert;
      int prop1 = impl->halfedge.data[manifold_next_halfedge(halfedgeIdx)].propVert;
      for (int j = 0; j < n; j++) {
        for (int p2 = 0; p2 < numProp; p2++) {
          prop.data[(offset + j) * numProp + p2] = double_lerp(
              impl->properties.data[prop0 * numProp + p2],
              impl->properties.data[prop1 * numProp + p2], (j + 1) * frac);
        }
      }
    }

    // Build triProp
    ManifoldVecIVec3 triProp = MANIFOLD_VEC_INIT;
    vec_ivec3_resize(&triProp, totalNewTris);
    for (int tri = 0; tri < numTri; tri++) {
      ManifoldIVec4 he = faceHalfedges[tri];
      if (ivec4_get(he, 0) < 0) continue;

      ManifoldIVec4 tri3 = manifold_ivec4(0, 0, 0, 0);
      ManifoldIVec4 eo2 = manifold_ivec4(0, 0, 0, 0);
      ManifoldBVec4 ef2 = {true, true, true, true};
      for (int i = 0; i < 4; i++) {
        if (ivec4_get(he, i) < 0) {
          tri3 = ivec4_set(tri3, i, -1);
          continue;
        }
        ManifoldHalfedge *h = &impl->halfedge.data[ivec4_get(he, i)];
        tri3 = ivec4_set(tri3, i, h->propVert);
        eo2 = ivec4_set(eo2, i, edgeOffset[half2Edge[ivec4_get(he, i)]]);
        if (!manifold_halfedge_is_forward(h)) {
          int paired = h->pairedHalfedge;
          if (impl->halfedge.data[paired].propVert !=
                  impl->halfedge.data[manifold_next_halfedge(ivec4_get(he, i))].propVert ||
              impl->halfedge.data[manifold_next_halfedge(paired)].propVert !=
                  h->propVert) {
            eo2 = ivec4_set(eo2, i, ivec4_get(eo2, i) + addedVerts);
          } else {
            switch (i) {
              case 0: ef2.x = false; break;
              case 1: ef2.y = false; break;
              case 2: ef2.z = false; break;
              case 3: ef2.w = false; break;
            }
          }
        }
      }

      // Add propOffset to eo2
      ManifoldIVec4 eo2p = manifold_ivec4(
          ivec4_get(eo2, 0) + propOffset, ivec4_get(eo2, 1) + propOffset,
          ivec4_get(eo2, 2) + propOffset, ivec4_get(eo2, 3) + propOffset);

      ManifoldVecIVec3 newTris = partition_reindex(&subTris[tri], tri3, eo2p, ef2,
          interiorOffset[tri] + propOffset);
      memcpy(&triProp.data[triOffset[tri]], newTris.data,
             newTris.len * sizeof(ManifoldIVec3));
      vec_ivec3_free(&newTris);
    }

    vec_double_free(&impl->properties);
    impl->properties = prop;
    manifold_impl_create_halfedges(impl, &triProp, &triVerts);
    vec_ivec3_free(&triProp);
  } else {
    manifold_impl_create_halfedges(impl, &triVerts, &triVerts);
  }

  // Cleanup
  vec_ivec3_free(&triVerts);
  for (int i = 0; i < numTri; i++) partition_free(&subTris[i]);
  free(subTris);
  free(triOffset);
  free(interiorOffset);
  free(half2Edge);
  free(faceHalfedges);
  free(edgeAdded);
  free(edgeOffset);
  vec_tmpedge_free(&edges);

  return vertBary;
}

// ===================== InterpTri (Bezier interpolation) =====================

static ManifoldVec4 homogeneous_v4(ManifoldVec4 v) {
  return manifold_vec4(v.x * v.w, v.y * v.w, v.z * v.w, v.w);
}

static ManifoldVec4 homogeneous_v3(ManifoldVec3 v) {
  return vec3_to_vec4(v, 1.0);
}

static ManifoldVec3 hnormalize(ManifoldVec4 v) {
  if (v.w == 0) return vec4_to_vec3(v);
  return vec3_scale(vec4_to_vec3(v), 1.0 / v.w);
}

static ManifoldVec4 bezier_point_v4(ManifoldVec3 point, ManifoldVec4 tangent) {
  ManifoldVec4 sum = vec4_add(vec3_to_vec4(point, 0), tangent);
  return homogeneous_v4(sum);
}

// CubicBezier2Linear: returns two vec4 (points[0], points[1])
static void cubic_bezier_2linear(ManifoldVec4 p0, ManifoldVec4 p1,
                                  ManifoldVec4 p2, ManifoldVec4 p3,
                                  double x,
                                  ManifoldVec4 *out0, ManifoldVec4 *out1) {
  ManifoldVec4 p12 = vec4_lerp(p1, p2, x);
  *out0 = vec4_lerp(vec4_lerp(p0, p1, x), p12, x);
  *out1 = vec4_lerp(p12, vec4_lerp(p2, p3, x), x);
}

static ManifoldVec3 bezier_point(ManifoldVec4 pt0, ManifoldVec4 pt1, double x) {
  return hnormalize(vec4_lerp(pt0, pt1, x));
}

static ManifoldVec3 bezier_tangent(ManifoldVec4 pt0, ManifoldVec4 pt1) {
  return manifold_safe_normalize(vec3_sub(hnormalize(pt1), hnormalize(pt0)));
}

static ManifoldVec3 rotate_from_to(ManifoldVec3 v, ManifoldQuat start, ManifoldQuat end) {
  return quat_rotate(end, quat_rotate(quat_conj(start), v));
}

static ManifoldQuat smooth_slerp(ManifoldQuat x, ManifoldQuat y, double a, bool longWay) {
  ManifoldQuat z = y;
  double cosTheta = vec4_dot(x, y);
  if ((cosTheta < 0) != longWay) {
    z = vec4_neg(y);
    cosTheta = -cosTheta;
  }
  if (fabs(cosTheta) > 1.0 - DBL_EPSILON) {
    return vec4_lerp(x, z, a);
  }
  double angle = acos(cosTheta);
  ManifoldQuat r = vec4_scale(vec4_add(
      vec4_scale(x, sin((1.0 - a) * angle)),
      vec4_scale(z, sin(a * angle))), 1.0 / sin(angle));
  return r;
}

static void bezier2bezier(ManifoldVec3 c0, ManifoldVec3 c1,
                           ManifoldVec4 tx0, ManifoldVec4 tx1,
                           ManifoldVec4 ty0, ManifoldVec4 ty1,
                           double x, ManifoldVec3 anchor,
                           ManifoldVec4 *out0, ManifoldVec4 *out1) {
  ManifoldVec4 bz0, bz1;
  cubic_bezier_2linear(homogeneous_v3(c0), bezier_point_v4(c0, tx0),
                        bezier_point_v4(c1, tx1), homogeneous_v3(c1), x,
                        &bz0, &bz1);
  ManifoldVec3 end = bezier_point(bz0, bz1, x);
  ManifoldVec3 tangent = bezier_tangent(bz0, bz1);

  ManifoldVec3 nTX0 = manifold_safe_normalize(vec4_to_vec3(tx0));
  ManifoldVec3 nTX1 = vec3_neg(manifold_safe_normalize(vec4_to_vec3(tx1)));

  ManifoldVec3 biTan0 = orthogonal_to(vec4_to_vec3(ty0), vec3_sub(anchor, c0), nTX0);
  ManifoldVec3 biTan1 = orthogonal_to(vec4_to_vec3(ty1), vec3_sub(anchor, c1), nTX1);

  ManifoldMat3 m0;
  m0.cols[0] = nTX0;
  m0.cols[1] = biTan0;
  m0.cols[2] = vec3_cross(nTX0, biTan0);
  ManifoldQuat q0 = quat_from_mat3(m0);

  ManifoldMat3 m1;
  m1.cols[0] = nTX1;
  m1.cols[1] = biTan1;
  m1.cols[2] = vec3_cross(nTX1, biTan1);
  ManifoldQuat q1 = quat_from_mat3(m1);

  ManifoldVec3 edge = vec3_sub(c1, c0);
  bool longWay = vec3_dot(nTX0, edge) + vec3_dot(nTX1, edge) < 0;
  ManifoldQuat qTmp = smooth_slerp(q0, q1, x, longWay);
  ManifoldQuat q = quat_mul(quat_rotation(quat_xdir(qTmp), tangent), qTmp);

  ManifoldVec3 delta = vec3_lerp(
      rotate_from_to(vec4_to_vec3(ty0), q0, q),
      rotate_from_to(vec4_to_vec3(ty1), q1, q), x);
  double deltaW = double_lerp(ty0.w, ty1.w, x);

  *out0 = homogeneous_v3(end);
  *out1 = vec3_to_vec4(delta, deltaW);
}

static ManifoldVec3 bezier2d(ManifoldVec3 corners[4], ManifoldVec4 tangentsX[4],
                              ManifoldVec4 tangentsY[4], double x, double y,
                              ManifoldVec3 centroid) {
  ManifoldVec4 bz0_0, bz0_1, bz1_0, bz1_1;
  bezier2bezier(corners[0], corners[1], tangentsX[0], tangentsX[1],
                 tangentsY[0], tangentsY[1], x, centroid, &bz0_0, &bz0_1);
  bezier2bezier(corners[2], corners[3], tangentsX[2], tangentsX[3],
                 tangentsY[2], tangentsY[3], 1 - x, centroid, &bz1_0, &bz1_1);

  ManifoldVec4 final0, final1;
  cubic_bezier_2linear(bz0_0, bezier_point_v4(vec4_to_vec3(bz0_0), bz0_1),
                        bezier_point_v4(vec4_to_vec3(bz1_0), bz1_1), bz1_0, y,
                        &final0, &final1);
  return bezier_point(final0, final1, y);
}

static void interp_tri(ManifoldImpl *impl,
                        const ManifoldVecBarycentric *vertBary,
                        const ManifoldImpl *old, int vert) {
  ManifoldVec3 *pos = &impl->vertPos.data[vert];
  int tri = vertBary->data[vert].tri;
  ManifoldVec4 uvw = vertBary->data[vert].uvw;

  ManifoldIVec4 halfedges = manifold_impl_get_halfedges(old, tri);
  ManifoldVec3 corners[4];
  corners[0] = old->vertPos.data[old->halfedge.data[ivec4_get(halfedges, 0)].startVert];
  corners[1] = old->vertPos.data[old->halfedge.data[ivec4_get(halfedges, 1)].startVert];
  corners[2] = old->vertPos.data[old->halfedge.data[ivec4_get(halfedges, 2)].startVert];
  corners[3] = ivec4_get(halfedges, 3) < 0
      ? manifold_vec3(0, 0, 0)
      : old->vertPos.data[old->halfedge.data[ivec4_get(halfedges, 3)].startVert];

  for (int i = 0; i < 4; i++) {
    if (vec4_get(uvw, i) == 1.0) {
      *pos = corners[i];
      return;
    }
  }

  ManifoldVec4 posH = manifold_vec4(0, 0, 0, 0);

  if (ivec4_get(halfedges, 3) < 0) {  // tri
    ManifoldVec4 tangentR[3], tangentL[3];
    for (int i = 0; i < 3; i++)
      tangentR[i] = old->halfedgeTangent.data[ivec4_get(halfedges, i)];
    tangentL[0] = old->halfedgeTangent.data[old->halfedge.data[ivec4_get(halfedges, 2)].pairedHalfedge];
    tangentL[1] = old->halfedgeTangent.data[old->halfedge.data[ivec4_get(halfedges, 0)].pairedHalfedge];
    tangentL[2] = old->halfedgeTangent.data[old->halfedge.data[ivec4_get(halfedges, 1)].pairedHalfedge];
    ManifoldVec3 centroid = vec3_scale(
        vec3_add(vec3_add(corners[0], corners[1]), corners[2]), 1.0 / 3.0);

    for (int i = 0; i < 3; i++) {
      int j = manifold_next3(i);
      int k = manifold_prev3(i);
      double uvw_i = vec4_get(uvw, i);
      double uvw_j = vec4_get(uvw, j);
      double uvw_k = vec4_get(uvw, k);
      double x = uvw_k / (1.0 - uvw_i);

      ManifoldVec4 bz0, bz1;
      bezier2bezier(corners[j], corners[k], tangentR[j], tangentL[k],
                     tangentL[j], tangentR[k], x, centroid, &bz0, &bz1);

      ManifoldVec4 tanLerp = vec4_lerp(tangentR[i], tangentL[i], x);
      ManifoldVec4 bz2_0, bz2_1;
      cubic_bezier_2linear(bz0, bezier_point_v4(vec4_to_vec3(bz0), bz1),
                            bezier_point_v4(corners[i], tanLerp),
                            homogeneous_v3(corners[i]), uvw_i,
                            &bz2_0, &bz2_1);
      ManifoldVec3 p = bezier_point(bz2_0, bz2_1, uvw_i);
      ManifoldVec4 hp = homogeneous_v4(vec3_to_vec4(p, uvw_j * uvw_k));
      posH = vec4_add(posH, hp);
    }
  } else {  // quad
    ManifoldVec4 tangentsX[4], tangentsY[4];
    tangentsX[0] = old->halfedgeTangent.data[ivec4_get(halfedges, 0)];
    tangentsX[1] = old->halfedgeTangent.data[old->halfedge.data[ivec4_get(halfedges, 0)].pairedHalfedge];
    tangentsX[2] = old->halfedgeTangent.data[ivec4_get(halfedges, 2)];
    tangentsX[3] = old->halfedgeTangent.data[old->halfedge.data[ivec4_get(halfedges, 2)].pairedHalfedge];
    tangentsY[0] = old->halfedgeTangent.data[old->halfedge.data[ivec4_get(halfedges, 3)].pairedHalfedge];
    tangentsY[1] = old->halfedgeTangent.data[ivec4_get(halfedges, 1)];
    tangentsY[2] = old->halfedgeTangent.data[old->halfedge.data[ivec4_get(halfedges, 1)].pairedHalfedge];
    tangentsY[3] = old->halfedgeTangent.data[ivec4_get(halfedges, 3)];
    ManifoldVec3 centroid = vec3_scale(vec3_add(vec3_add(vec3_add(
        corners[0], corners[1]), corners[2]), corners[3]), 0.25);
    double x = uvw.y + uvw.z;
    double y = uvw.z + uvw.w;
    ManifoldVec3 pX = bezier2d(corners, tangentsX, tangentsY, x, y, centroid);
    ManifoldVec3 cornersY[4] = {corners[1], corners[2], corners[3], corners[0]};
    ManifoldVec4 txY[4] = {tangentsY[1], tangentsY[2], tangentsY[3], tangentsY[0]};
    ManifoldVec4 tyY[4] = {tangentsX[1], tangentsX[2], tangentsX[3], tangentsX[0]};
    ManifoldVec3 pY = bezier2d(cornersY, txY, tyY, y, 1 - x, centroid);
    posH = vec4_add(posH, homogeneous_v4(vec3_to_vec4(pX, x * (1 - x))));
    posH = vec4_add(posH, homogeneous_v4(vec3_to_vec4(pY, y * (1 - y))));
  }
  *pos = hnormalize(posH);
}

// ===================== Refine (using subdivision) =====================

typedef struct {
  int n;
} RefineNCtx;

MANIFOLD_UNUSED static int refine_n_divisions(ManifoldVec3 edge, ManifoldVec4 t0,
                               ManifoldVec4 t1, void *ctx) {
  (void)edge; (void)t0; (void)t1;
  return ((RefineNCtx *)ctx)->n - 1;
}

typedef struct {
  double length;
} RefineLengthCtx;

MANIFOLD_UNUSED static int refine_length_divisions(ManifoldVec3 edge, ManifoldVec4 t0,
                                    ManifoldVec4 t1, void *ctx) {
  (void)t0; (void)t1;
  double len = vec3_length(edge);
  return (int)(len / ((RefineLengthCtx *)ctx)->length);
}

typedef struct {
  double tolerance;
} RefineToleranceCtx;

MANIFOLD_UNUSED static int refine_tolerance_divisions(ManifoldVec3 edge, ManifoldVec4 t0,
                                       ManifoldVec4 t1, void *ctx) {
  double tolerance = ((RefineToleranceCtx *)ctx)->tolerance;
  ManifoldVec3 edgeNorm = manifold_safe_normalize(edge);
  ManifoldVec3 tStart = vec4_to_vec3(t0);
  ManifoldVec3 tEnd = vec4_to_vec3(t1);
  ManifoldVec3 start = vec3_sub(tStart, vec3_scale(edgeNorm, vec3_dot(edgeNorm, tStart)));
  ManifoldVec3 end = vec3_sub(tEnd, vec3_scale(edgeNorm, vec3_dot(edgeNorm, tEnd)));
  double d = 0.5 * (vec3_length(start) + vec3_length(end)) +
             vec3_length(vec3_sub(start, end));
  return (int)sqrt(3.0 * d / (4.0 * tolerance));
}

// Deep copy impl
static void impl_deep_copy(ManifoldImpl *dst, const ManifoldImpl *src) {
  manifold_impl_init(dst);
  dst->bBox = src->bBox;
  dst->epsilon = src->epsilon;
  dst->tolerance = src->tolerance;
  dst->numProp = src->numProp;
  dst->status = src->status;

  vec_vec3_resize(&dst->vertPos, src->vertPos.len);
  memcpy(dst->vertPos.data, src->vertPos.data, src->vertPos.len * sizeof(ManifoldVec3));

  vec_halfedge_resize(&dst->halfedge, src->halfedge.len);
  memcpy(dst->halfedge.data, src->halfedge.data, src->halfedge.len * sizeof(ManifoldHalfedge));

  if (src->properties.len > 0) {
    vec_double_resize(&dst->properties, src->properties.len);
    memcpy(dst->properties.data, src->properties.data, src->properties.len * sizeof(double));
  }

  if (src->vertNormal.len > 0) {
    vec_vec3_resize(&dst->vertNormal, src->vertNormal.len);
    memcpy(dst->vertNormal.data, src->vertNormal.data, src->vertNormal.len * sizeof(ManifoldVec3));
  }

  if (src->faceNormal.len > 0) {
    vec_vec3_resize(&dst->faceNormal, src->faceNormal.len);
    memcpy(dst->faceNormal.data, src->faceNormal.data, src->faceNormal.len * sizeof(ManifoldVec3));
  }

  if (src->halfedgeTangent.len > 0) {
    vec_vec4_resize(&dst->halfedgeTangent, src->halfedgeTangent.len);
    memcpy(dst->halfedgeTangent.data, src->halfedgeTangent.data, src->halfedgeTangent.len * sizeof(ManifoldVec4));
  }

  dst->meshRelation.originalID = src->meshRelation.originalID;
  vec_meshid_resize(&dst->meshRelation.meshIDtransform, src->meshRelation.meshIDtransform.len);
  memcpy(dst->meshRelation.meshIDtransform.data, src->meshRelation.meshIDtransform.data,
         src->meshRelation.meshIDtransform.len * sizeof(ManifoldMeshIDEntry));
  vec_triref_resize(&dst->meshRelation.triRef, src->meshRelation.triRef.len);
  memcpy(dst->meshRelation.triRef.data, src->meshRelation.triRef.data,
         src->meshRelation.triRef.len * sizeof(ManifoldTriRef));
}

void manifold_impl_refine(ManifoldImpl *impl, EdgeDivisionsFn edgeDivisions,
                           void *ctx, bool keepInterior) {
  if (manifold_impl_is_empty(impl)) return;

  ManifoldImpl old;
  impl_deep_copy(&old, impl);

  ManifoldVecBarycentric vertBary = manifold_impl_subdivide(
      impl, edgeDivisions, ctx, keepInterior);
  if (vertBary.len == 0) {
    vec_bary_free(&vertBary);
    manifold_impl_free(&old);
    return;
  }

  if (old.halfedgeTangent.len == old.halfedge.len) {
    for (size_t v = 0; v < impl->vertPos.len; v++) {
      interp_tri(impl, &vertBary, &old, (int)v);
    }
  }

  vec_vec4_free(&impl->halfedgeTangent);
  manifold_impl_sort_geometry(impl);
  if (old.halfedgeTangent.len == old.halfedge.len) {
    manifold_impl_set_normals_and_coplanar(impl);
  } else {
    manifold_impl_calculate_vert_normals(impl);
  }
  impl->meshRelation.originalID = -1;

  vec_bary_free(&vertBary);
  manifold_impl_free(&old);
}

// ===================== Smoothing Functions =====================

ManifoldVec3 manifold_impl_get_normal(const ManifoldImpl *impl, int halfedge,
                                       int normalIdx) {
  int prop = impl->halfedge.data[halfedge].propVert;
  ManifoldVec3 normal;
  normal.x = impl->properties.data[prop * impl->numProp + normalIdx];
  normal.y = impl->properties.data[prop * impl->numProp + normalIdx + 1];
  normal.z = impl->properties.data[prop * impl->numProp + normalIdx + 2];
  return normal;
}

ManifoldVec4 manifold_impl_tangent_from_normal(const ManifoldImpl *impl,
                                                ManifoldVec3 normal,
                                                int halfedge) {
  ManifoldHalfedge edge = impl->halfedge.data[halfedge];
  ManifoldVec3 edgeVec = vec3_sub(impl->vertPos.data[edge.endVert],
                                   impl->vertPos.data[edge.startVert]);
  ManifoldVec3 edgeNormal = vec3_add(impl->faceNormal.data[halfedge / 3],
                                      impl->faceNormal.data[edge.pairedHalfedge / 3]);
  ManifoldVec3 dir = vec3_cross(vec3_cross(edgeNormal, edgeVec), normal);
  return circular_tangent(dir, edgeVec);
}

// FlatFaces: find faces with at least 3 triangles
static ManifoldVecBool flat_faces(const ManifoldImpl *impl) {
  int numTri = (int)manifold_impl_num_tri(impl);
  ManifoldVecBool triIsFlatFace = MANIFOLD_VEC_INIT;
  vec_bool_resize(&triIsFlatFace, numTri);
  for (int i = 0; i < numTri; i++) triIsFlatFace.data[i] = false;

  for (int tri = 0; tri < numTri; tri++) {
    ManifoldTriRef ref = impl->meshRelation.triRef.data[tri];
    int faceNeighbors = 0;
    int faceTris[3] = {-1, -1, -1};
    for (int j = 0; j < 3; j++) {
      int neighborTri = impl->halfedge.data[3 * tri + j].pairedHalfedge / 3;
      ManifoldTriRef jRef = impl->meshRelation.triRef.data[neighborTri];
      if (manifold_triref_same_face(&jRef, &ref)) {
        faceNeighbors++;
        faceTris[j] = neighborTri;
      }
    }
    if (faceNeighbors > 1) {
      triIsFlatFace.data[tri] = true;
      for (int j = 0; j < 3; j++) {
        if (faceTris[j] >= 0) triIsFlatFace.data[faceTris[j]] = true;
      }
    }
  }
  return triIsFlatFace;
}

// VertFlatFace
static ManifoldVecInt vert_flat_face(const ManifoldImpl *impl,
                                      const ManifoldVecBool *flatFaces) {
  int numVert = (int)manifold_impl_num_vert(impl);
  ManifoldVecInt result = MANIFOLD_VEC_INIT;
  vec_int_resize(&result, numVert);
  for (int i = 0; i < numVert; i++) result.data[i] = -1;

  ManifoldVecTriRef vertRef = MANIFOLD_VEC_INIT;
  vec_triref_resize(&vertRef, numVert);
  ManifoldTriRef emptyRef = {-1, -1, -1, -1};
  for (int i = 0; i < numVert; i++) vertRef.data[i] = emptyRef;

  size_t numTri = manifold_impl_num_tri(impl);
  for (size_t tri = 0; tri < numTri; tri++) {
    if (!flatFaces->data[tri]) continue;
    for (int j = 0; j < 3; j++) {
      int vert = impl->halfedge.data[3 * tri + j].startVert;
      if (manifold_triref_same_face(&vertRef.data[vert],
                                     &impl->meshRelation.triRef.data[tri]))
        continue;
      vertRef.data[vert] = impl->meshRelation.triRef.data[tri];
      result.data[vert] = (result.data[vert] == -1) ? (int)tri : -2;
    }
  }
  vec_triref_free(&vertRef);
  return result;
}

// VertHalfedge: one halfedge per vertex
static ManifoldVecInt vert_halfedge(const ManifoldImpl *impl) {
  int numVert = (int)manifold_impl_num_vert(impl);
  ManifoldVecInt result = MANIFOLD_VEC_INIT;
  vec_int_resize(&result, numVert);
  // C++ overwrites, so last halfedge per vertex wins
  for (size_t i = 0; i < impl->halfedge.len; i++) {
    int v = impl->halfedge.data[i].startVert;
    if (v >= 0 && v < numVert) {
      result.data[v] = (int)i;
    }
  }
  return result;
}

// SharpenEdges
ManifoldVecSmoothness manifold_impl_sharpen_edges(const ManifoldImpl *impl,
                                                    double minSharpAngle,
                                                    double minSmoothness) {
  ManifoldVecSmoothness result = MANIFOLD_VEC_INIT;
  double minRadians = manifold_radians(minSharpAngle);
  for (size_t e = 0; e < impl->halfedge.len; e++) {
    if (!manifold_halfedge_is_forward(&impl->halfedge.data[e])) continue;
    size_t pair = impl->halfedge.data[e].pairedHalfedge;
    double d = vec3_dot(impl->faceNormal.data[e / 3], impl->faceNormal.data[pair / 3]);
    double dihedral = acos(fmax(-1.0, fmin(1.0, d)));
    if (dihedral > minRadians) {
      ManifoldSmoothness s1 = {e, minSmoothness};
      ManifoldSmoothness s2 = {pair, minSmoothness};
      vec_smooth_push(&result, s1);
      vec_smooth_push(&result, s2);
    }
  }
  return result;
}

// SharpenTangent
static void sharpen_tangent(ManifoldImpl *impl, int halfedge, double smoothness) {
  ManifoldVec4 *t = &impl->halfedgeTangent.data[halfedge];
  *t = manifold_vec4(smoothness * t->x, smoothness * t->y, smoothness * t->z,
                      smoothness == 0 ? 0 : t->w);
}

// LinearizeFlatTangents
static void linearize_flat_tangents(ManifoldImpl *impl) {
  int n = (int)impl->halfedgeTangent.len;
  for (int halfedge = 0; halfedge < n; halfedge++) {
    ManifoldVec4 *tangent = &impl->halfedgeTangent.data[halfedge];
    ManifoldVec4 *otherTangent =
        &impl->halfedgeTangent.data[impl->halfedge.data[halfedge].pairedHalfedge];

    bool flat0 = tangent->w == 0;
    bool flat1 = otherTangent->w == 0;
    if (!manifold_halfedge_is_forward(&impl->halfedge.data[halfedge]) || (!flat0 && !flat1))
      continue;

    ManifoldVec3 edgeVec = vec3_sub(
        impl->vertPos.data[impl->halfedge.data[halfedge].endVert],
        impl->vertPos.data[impl->halfedge.data[halfedge].startVert]);

    if (flat0 && flat1) {
      *tangent = vec3_to_vec4(vec3_scale(edgeVec, 1.0 / 3.0), 1.0);
      *otherTangent = vec3_to_vec4(vec3_scale(edgeVec, -1.0 / 3.0), 1.0);
    } else if (flat0) {
      *tangent = vec3_to_vec4(vec3_scale(vec3_add(edgeVec, vec4_to_vec3(*otherTangent)), 0.5), 1.0);
    } else {
      *otherTangent = vec3_to_vec4(vec3_scale(vec3_add(vec3_neg(edgeVec), vec4_to_vec3(*tangent)), 0.5), 1.0);
    }
  }
}

// DistributeTangents
static void distribute_tangents(ManifoldImpl *impl, const ManifoldVecBool *fixedHalfedges) {
  int numHalfedge = (int)fixedHalfedges->len;
  for (int halfedge = 0; halfedge < numHalfedge; halfedge++) {
    if (!fixedHalfedges->data[halfedge]) continue;

    int he = halfedge;
    if (manifold_impl_is_marked_inside_quad(impl, he)) {
      he = manifold_next_halfedge(impl->halfedge.data[he].pairedHalfedge);
    }

    ManifoldVec3 normal = manifold_vec3(0, 0, 0);
    ManifoldVec3 approxNormal = impl->vertNormal.data[impl->halfedge.data[he].startVert];
    ManifoldVec3 center = impl->vertPos.data[impl->halfedge.data[he].startVert];
    ManifoldVec3 lastEdgeVec = manifold_safe_normalize(
        vec3_sub(impl->vertPos.data[impl->halfedge.data[he].endVert], center));
    ManifoldVec3 firstTangent = manifold_safe_normalize(vec4_to_vec3(impl->halfedgeTangent.data[he]));
    ManifoldVec3 lastTangent = firstTangent;

    // Collect angles
    int maxDeg = (int)(impl->halfedge.len / 3) + 10;
    double *currentAngle = (double *)malloc(maxDeg * sizeof(double));
    double *desiredAngle = (double *)malloc(maxDeg * sizeof(double));
    int count = 0;

    int current = he;
    do {
      current = manifold_next_halfedge(impl->halfedge.data[current].pairedHalfedge);
      if (manifold_impl_is_marked_inside_quad(impl, current)) continue;
      ManifoldVec3 thisEdgeVec = manifold_safe_normalize(
          vec3_sub(impl->vertPos.data[impl->halfedge.data[current].endVert], center));
      ManifoldVec3 thisTangent = manifold_safe_normalize(vec4_to_vec3(impl->halfedgeTangent.data[current]));
      normal = vec3_add(normal, vec3_cross(thisTangent, lastTangent));

      desiredAngle[count] = angle_between(thisEdgeVec, lastEdgeVec) +
                            (count > 0 ? desiredAngle[count - 1] : 0.0);
      if (current == he) {
        currentAngle[count] = MANIFOLD_TWO_PI;
      } else {
        currentAngle[count] = angle_between(thisTangent, firstTangent);
        if (vec3_dot(approxNormal, vec3_cross(thisTangent, firstTangent)) < 0)
          currentAngle[count] = MANIFOLD_TWO_PI - currentAngle[count];
      }
      lastEdgeVec = thisEdgeVec;
      lastTangent = thisTangent;
      count++;
      if (count >= maxDeg) break;
    } while (!fixedHalfedges->data[current]);

    if (count <= 1 || vec3_dot(normal, normal) == 0) {
      free(currentAngle);
      free(desiredAngle);
      continue;
    }

    double scale = currentAngle[count - 1] / desiredAngle[count - 1];
    double offset = 0;
    if (current == he) {  // only one fixed
      for (int i = 0; i < count; i++)
        offset += wrap_angle(currentAngle[i] - scale * desiredAngle[i]);
      offset /= count;
    }

    current = he;
    int idx = 0;
    do {
      current = manifold_next_halfedge(impl->halfedge.data[current].pairedHalfedge);
      if (manifold_impl_is_marked_inside_quad(impl, current)) continue;
      desiredAngle[idx] *= scale;
      double lastAngle = idx > 0 ? desiredAngle[idx - 1] : 0;
      if (desiredAngle[idx] - lastAngle > MANIFOLD_PI) {
        desiredAngle[idx] = lastAngle + MANIFOLD_PI;
      } else if (idx + 1 < count &&
                 scale * desiredAngle[idx + 1] - desiredAngle[idx] > MANIFOLD_PI) {
        desiredAngle[idx] = scale * desiredAngle[idx + 1] - MANIFOLD_PI;
      }
      double angle = currentAngle[idx] - desiredAngle[idx] - offset;
      ManifoldVec3 tangent = vec4_to_vec3(impl->halfedgeTangent.data[current]);
      ManifoldQuat q = quat_from_axis_angle(vec3_normalize(normal), angle);
      impl->halfedgeTangent.data[current] =
          vec3_to_vec4(quat_rotate(q, tangent), impl->halfedgeTangent.data[current].w);
      idx++;
      if (idx >= count) break;
    } while (!fixedHalfedges->data[current]);

    free(currentAngle);
    free(desiredAngle);
  }
}

// SetNormals for SmoothOut with minSmoothness == 0
void manifold_impl_set_normals_smooth(ManifoldImpl *impl, int normalIdx,
                                       double minSharpAngle) {
  if (manifold_impl_is_empty(impl)) return;
  if (normalIdx < 0) return;

  int oldNumProp = impl->numProp;
  ManifoldVecBool triIsFlatFace = flat_faces(impl);
  ManifoldVecInt vertFlatFace_data = vert_flat_face(impl, &triIsFlatFace);
  int numVert = (int)manifold_impl_num_vert(impl);
  int *vertNumSharp = (int *)calloc(numVert, sizeof(int));

  for (size_t e = 0; e < impl->halfedge.len; e++) {
    if (!manifold_halfedge_is_forward(&impl->halfedge.data[e])) continue;
    int pair = impl->halfedge.data[e].pairedHalfedge;
    int tri1 = (int)(e / 3);
    int tri2 = pair / 3;
    double d = vec3_dot(impl->faceNormal.data[tri1], impl->faceNormal.data[tri2]);
    double dihedral = manifold_degrees(acos(fmax(-1.0, fmin(1.0, d))));
    if (dihedral > minSharpAngle) {
      vertNumSharp[impl->halfedge.data[e].startVert]++;
      vertNumSharp[impl->halfedge.data[e].endVert]++;
    } else {
      bool faceSplit =
          triIsFlatFace.data[tri1] != triIsFlatFace.data[tri2] ||
          (triIsFlatFace.data[tri1] && triIsFlatFace.data[tri2] &&
           !manifold_triref_same_face(&impl->meshRelation.triRef.data[tri1],
                                       &impl->meshRelation.triRef.data[tri2]));
      if (vertFlatFace_data.data[impl->halfedge.data[e].startVert] == -2 && faceSplit)
        vertNumSharp[impl->halfedge.data[e].startVert]++;
      if (vertFlatFace_data.data[impl->halfedge.data[e].endVert] == -2 && faceSplit)
        vertNumSharp[impl->halfedge.data[e].endVert]++;
    }
  }

  int numProp = oldNumProp > normalIdx + 3 ? oldNumProp : normalIdx + 3;
  int numPropVert = (int)manifold_impl_num_prop_vert(impl);
  ManifoldVecDouble oldProperties = MANIFOLD_VEC_INIT;
  vec_double_resize(&oldProperties, numProp * numPropVert);
  // Zero-fill
  memset(oldProperties.data, 0, oldProperties.len * sizeof(double));
  // Copy existing properties
  for (int v = 0; v < numPropVert; v++) {
    for (int p = 0; p < oldNumProp && p < numProp; p++) {
      oldProperties.data[v * numProp + p] = impl->properties.data[v * oldNumProp + p];
    }
  }
  vec_double_free(&impl->properties);
  vec_double_resize(&impl->properties, numProp * numPropVert);
  memset(impl->properties.data, 0, impl->properties.len * sizeof(double));
  impl->numProp = numProp;

  // Save old propVert and reset
  int *oldHalfedgeProp = (int *)malloc(impl->halfedge.len * sizeof(int));
  for (size_t i = 0; i < impl->halfedge.len; i++) {
    oldHalfedgeProp[i] = impl->halfedge.data[i].propVert;
    impl->halfedge.data[i].propVert = -1;
  }

  int numEdge = (int)impl->halfedge.len;
  for (int startEdge = 0; startEdge < numEdge; startEdge++) {
    if (impl->halfedge.data[startEdge].propVert >= 0) continue;
    int vert = impl->halfedge.data[startEdge].startVert;

    if (vertNumSharp[vert] < 2) {
      // single normal vertex
      ManifoldVec3 normal;
      if (vertFlatFace_data.data[vert] >= 0)
        normal = impl->faceNormal.data[vertFlatFace_data.data[vert]];
      else
        normal = impl->vertNormal.data[vert];

      int lastProp = -1;
      int current = startEdge;
      int iterations = 0;
      do {
        int prop = oldHalfedgeProp[current];
        impl->halfedge.data[current].propVert = prop;
        if (prop != lastProp) {
          lastProp = prop;
          for (int p = 0; p < oldNumProp && p < numProp; p++)
            impl->properties.data[prop * numProp + p] = oldProperties.data[prop * numProp + p];
          for (int i = 0; i < 3; i++)
            impl->properties.data[prop * numProp + normalIdx + i] =
                (i == 0 ? normal.x : (i == 1 ? normal.y : normal.z));
        }
        current = manifold_next_halfedge(impl->halfedge.data[current].pairedHalfedge);
        if (++iterations > numEdge) break;
      } while (current != startEdge);
    } else {
      // multi-normal vertex: compute per-group pseudo-normals
      ManifoldVec3 centerPos = impl->vertPos.data[vert];

      // Find a sharp edge to start on
      int current = startEdge;
      int prevFace = current / 3;
      int iterations2 = 0;
      do {
        int next = manifold_next_halfedge(impl->halfedge.data[current].pairedHalfedge);
        int face = next / 3;
        double d = vec3_dot(impl->faceNormal.data[face], impl->faceNormal.data[prevFace]);
        double dihedral = manifold_degrees(acos(fmax(-1.0, fmin(1.0, d))));
        bool isSharp = dihedral > minSharpAngle ||
            triIsFlatFace.data[face] != triIsFlatFace.data[prevFace] ||
            (triIsFlatFace.data[face] && triIsFlatFace.data[prevFace] &&
             !manifold_triref_same_face(&impl->meshRelation.triRef.data[face],
                                         &impl->meshRelation.triRef.data[prevFace]));
        if (isSharp) break;
        current = next;
        prevFace = face;
        if (++iterations2 > numEdge) break;
      } while (current != startEdge);
      int endEdge = current;

      // Phase 1: ForVert binary op to compute groups and pseudo-normals
      // Dynamic arrays for groups and normals
      int *group = (int *)malloc(numEdge * sizeof(int));
      ManifoldVec3 *normals = (ManifoldVec3 *)malloc(numEdge * sizeof(ManifoldVec3));
      int groupCount = 0;
      int numNormals = 0;

      // ForVert<FaceEdge> binary op
      // Transform: compute FaceEdge for each halfedge
      // First, compute the initial "here" value
      int curHE = endEdge;
      int hereFace;
      ManifoldVec3 hereEdgeVec;
      if (manifold_impl_is_inside_quad(impl, curHE)) {
        hereFace = curHE / 3;
        hereEdgeVec = manifold_vec3(NAN, NAN, NAN);
      } else {
        hereFace = curHE / 3;
        int endV = impl->halfedge.data[curHE].endVert;
        ManifoldVec3 pos = impl->vertPos.data[endV];
        if (vertNumSharp[endV] < 2) {
          ManifoldVec3 n2 = vertFlatFace_data.data[endV] >= 0
              ? impl->faceNormal.data[vertFlatFace_data.data[endV]]
              : impl->vertNormal.data[endV];
          ManifoldVec4 tFromN = manifold_impl_tangent_from_normal(impl, n2,
              impl->halfedge.data[curHE].pairedHalfedge);
          pos = vec3_add(pos, vec4_to_vec3(tFromN));
        }
        hereEdgeVec = manifold_safe_normalize(vec3_sub(pos, centerPos));
      }

      iterations2 = 0;
      current = endEdge;
      do {
        int nextHE = manifold_next_halfedge(impl->halfedge.data[current].pairedHalfedge);
        // Transform for next
        int nextFace;
        ManifoldVec3 nextEdgeVec;
        if (manifold_impl_is_inside_quad(impl, nextHE)) {
          nextFace = nextHE / 3;
          nextEdgeVec = manifold_vec3(NAN, NAN, NAN);
        } else {
          nextFace = nextHE / 3;
          int endV = impl->halfedge.data[nextHE].endVert;
          ManifoldVec3 pos = impl->vertPos.data[endV];
          if (vertNumSharp[endV] < 2) {
            ManifoldVec3 n2 = vertFlatFace_data.data[endV] >= 0
                ? impl->faceNormal.data[vertFlatFace_data.data[endV]]
                : impl->vertNormal.data[endV];
            ManifoldVec4 tFromN = manifold_impl_tangent_from_normal(impl, n2,
                impl->halfedge.data[nextHE].pairedHalfedge);
            pos = vec3_add(pos, vec4_to_vec3(tFromN));
          }
          nextEdgeVec = manifold_safe_normalize(vec3_sub(pos, centerPos));
        }

        // Binary op: check if edge is sharp between here and next
        double d2 = vec3_dot(impl->faceNormal.data[hereFace], impl->faceNormal.data[nextFace]);
        double dihedral2 = manifold_degrees(acos(fmax(-1.0, fmin(1.0, d2))));
        bool isSharp2 = dihedral2 > minSharpAngle ||
            triIsFlatFace.data[hereFace] != triIsFlatFace.data[nextFace] ||
            (triIsFlatFace.data[hereFace] && triIsFlatFace.data[nextFace] &&
             !manifold_triref_same_face(&impl->meshRelation.triRef.data[hereFace],
                                         &impl->meshRelation.triRef.data[nextFace]));
        if (isSharp2) {
          normals[numNormals] = manifold_vec3(0, 0, 0);
          numNormals++;
        }
        group[groupCount++] = numNormals - 1;
        if (isfinite(nextEdgeVec.x)) {
          ManifoldVec3 cr = manifold_safe_normalize(vec3_cross(nextEdgeVec, hereEdgeVec));
          double ang = angle_between(hereEdgeVec, nextEdgeVec);
          normals[numNormals - 1] = vec3_add(normals[numNormals - 1], vec3_scale(cr, ang));
        } else {
          nextEdgeVec = hereEdgeVec;
        }

        hereFace = nextFace;
        hereEdgeVec = nextEdgeVec;
        current = nextHE;
        if (++iterations2 > numEdge) break;
      } while (current != endEdge);

      // Normalize pseudo-normals
      for (int i = 0; i < numNormals; i++) {
        normals[i] = manifold_safe_normalize(normals[i]);
      }

      // Phase 2: assign normals to halfedges, splitting properties as needed
      // C++ ForVert advances before processing, so first halfedge processed
      // is next(endEdge), not endEdge itself
      int lastGroup2 = 0;
      int lastProp = -1;
      int newProp = -1;
      int idx = 0;
      current = endEdge;
      iterations2 = 0;
      do {
        current = manifold_next_halfedge(impl->halfedge.data[current].pairedHalfedge);
        int prop = oldHalfedgeProp[current];

        if (idx < groupCount && group[idx] != lastGroup2 && group[idx] != 0 && prop == lastProp) {
          // Split property vertex
          lastGroup2 = group[idx];
          newProp = (int)(impl->properties.len / numProp);
          vec_double_resize(&impl->properties, impl->properties.len + numProp);
          // Copy old property data
          for (int p = 0; p < oldNumProp && p < numProp; p++)
            impl->properties.data[newProp * numProp + p] = oldProperties.data[prop * numProp + p];
          for (int p = oldNumProp; p < numProp; p++)
            impl->properties.data[newProp * numProp + p] = 0;
          // Set normal
          ManifoldVec3 n = idx < groupCount ? normals[group[idx]] : manifold_vec3(0,0,0);
          impl->properties.data[newProp * numProp + normalIdx + 0] = n.x;
          impl->properties.data[newProp * numProp + normalIdx + 1] = n.y;
          impl->properties.data[newProp * numProp + normalIdx + 2] = n.z;
        } else if (prop != lastProp) {
          lastProp = prop;
          newProp = prop;
          for (int p = 0; p < oldNumProp && p < numProp; p++)
            impl->properties.data[prop * numProp + p] = oldProperties.data[prop * numProp + p];
          ManifoldVec3 n = idx < groupCount ? normals[group[idx]] : manifold_vec3(0,0,0);
          impl->properties.data[prop * numProp + normalIdx + 0] = n.x;
          impl->properties.data[prop * numProp + normalIdx + 1] = n.y;
          impl->properties.data[prop * numProp + normalIdx + 2] = n.z;
        }

        impl->halfedge.data[current].propVert = newProp;
        idx++;

        if (++iterations2 > numEdge) break;
      } while (current != endEdge);

      free(group);
      free(normals);
    }
  }

  free(oldHalfedgeProp);
  free(vertNumSharp);
  vec_double_free(&oldProperties);
  vec_bool_free(&triIsFlatFace);
  vec_int_free(&vertFlatFace_data);
}

// CreateTangents from normalIdx
void manifold_impl_create_tangents_normals(ManifoldImpl *impl, int normalIdx) {
  int numVert = (int)manifold_impl_num_vert(impl);
  int numHalfedge = (int)impl->halfedge.len;

  vec_vec4_free(&impl->halfedgeTangent);
  ManifoldVecVec4 tangent = MANIFOLD_VEC_INIT;
  vec_vec4_resize(&tangent, numHalfedge);
  ManifoldVecBool fixedHalfedge = MANIFOLD_VEC_INIT;
  vec_bool_resize(&fixedHalfedge, numHalfedge);
  for (int i = 0; i < numHalfedge; i++) {
    tangent.data[i] = manifold_vec4(0, 0, 0, 0);
    fixedHalfedge.data[i] = false;
  }

  ManifoldVecInt vertHE = vert_halfedge(impl);

  for (int v = 0; v < numVert; v++) {
    int e = vertHE.data[v];
    ManifoldIVec2 faceEdges = manifold_ivec2(-1, -1);

    // ForVert with binary op
    int current = e;
    int iterations = 0;
    ManifoldVec3 hereNormal = manifold_impl_get_normal(impl, e, normalIdx);
    ManifoldVec3 hereDiff = vec3_sub(impl->faceNormal.data[e / 3], hereNormal);
    bool hereIsFlat = vec3_dot(hereDiff, hereDiff) < MANIFOLD_PRECISION * MANIFOLD_PRECISION;

    do {
      int nextHalfedge = manifold_next_halfedge(impl->halfedge.data[current].pairedHalfedge);
      ManifoldVec3 nextNormal = manifold_impl_get_normal(impl, nextHalfedge, normalIdx);
      ManifoldVec3 nextDiff = vec3_sub(impl->faceNormal.data[nextHalfedge / 3], nextNormal);
      bool nextIsFlat = vec3_dot(nextDiff, nextDiff) < MANIFOLD_PRECISION * MANIFOLD_PRECISION;

      if (manifold_impl_is_inside_quad(impl, current)) {
        tangent.data[current] = manifold_vec4(0, 0, 0, -1);
      } else {
        ManifoldVec3 diff = vec3_sub(nextNormal, hereNormal);
        bool differentNormals = vec3_dot(diff, diff) > MANIFOLD_PRECISION * MANIFOLD_PRECISION;
        if (differentNormals || hereIsFlat != nextIsFlat) {
          fixedHalfedge.data[current] = true;
          if (faceEdges.x == -1) faceEdges.x = current;
          else if (faceEdges.y == -1) faceEdges.y = current;
          else faceEdges.x = -2;
        }
        if (differentNormals) {
          ManifoldVec3 edgeVec = vec3_sub(
              impl->vertPos.data[impl->halfedge.data[current].endVert],
              impl->vertPos.data[impl->halfedge.data[current].startVert]);
          ManifoldVec3 dir = vec3_cross(hereNormal, nextNormal);
          tangent.data[current] = circular_tangent(
              vec3_scale(dir, vec3_dot(dir, edgeVec) < 0 ? -1.0 : 1.0), edgeVec);
        } else {
          tangent.data[current] = manifold_impl_tangent_from_normal(impl, hereNormal, current);
        }
      }

      hereNormal = nextNormal;
      hereIsFlat = nextIsFlat;
      current = nextHalfedge;
      if (++iterations > numHalfedge) break;
    } while (current != e);

    if (faceEdges.x >= 0 && faceEdges.y >= 0) {
      ManifoldVec3 edge0 = vec3_sub(
          impl->vertPos.data[impl->halfedge.data[faceEdges.x].endVert],
          impl->vertPos.data[impl->halfedge.data[faceEdges.x].startVert]);
      ManifoldVec3 edge1 = vec3_sub(
          impl->vertPos.data[impl->halfedge.data[faceEdges.y].endVert],
          impl->vertPos.data[impl->halfedge.data[faceEdges.y].startVert]);
      ManifoldVec3 newTangent = vec3_sub(vec3_normalize(edge0), vec3_normalize(edge1));
      tangent.data[faceEdges.x] = circular_tangent(newTangent, edge0);
      tangent.data[faceEdges.y] = circular_tangent(vec3_neg(newTangent), edge1);
    } else if (faceEdges.x == -1 && faceEdges.y == -1) {
      fixedHalfedge.data[e] = true;
    }
  }

  impl->halfedgeTangent = tangent;
  distribute_tangents(impl, &fixedHalfedge);
  vec_bool_free(&fixedHalfedge);
  vec_int_free(&vertHE);
}

// CreateTangents from sharpenedEdges (for SmoothOut with minSmoothness > 0,
// and for Smooth constructor)
void manifold_impl_create_tangents_smooth(ManifoldImpl *impl,
    const ManifoldSmoothness *sharpenedEdges, int numSharpened) {

  int numHalfedge = (int)impl->halfedge.len;
  vec_vec4_free(&impl->halfedgeTangent);
  ManifoldVecVec4 tangent = MANIFOLD_VEC_INIT;
  vec_vec4_resize(&tangent, numHalfedge);
  ManifoldVecBool fixedHalfedge = MANIFOLD_VEC_INIT;
  vec_bool_resize(&fixedHalfedge, numHalfedge);
  for (int i = 0; i < numHalfedge; i++) fixedHalfedge.data[i] = false;

  ManifoldVecInt vertHE = vert_halfedge(impl);
  ManifoldVecBool triIsFlatFace = flat_faces(impl);
  ManifoldVecInt vertFlatFace_data = vert_flat_face(impl, &triIsFlatFace);
  int numVert = (int)manifold_impl_num_vert(impl);

  ManifoldVecVec3 vertNormal2 = MANIFOLD_VEC_INIT;
  vec_vec3_resize(&vertNormal2, numVert);
  for (int v = 0; v < numVert; v++) {
    if (vertFlatFace_data.data[v] >= 0)
      vertNormal2.data[v] = impl->faceNormal.data[vertFlatFace_data.data[v]];
    else
      vertNormal2.data[v] = impl->vertNormal.data[v];
  }

  for (int i = 0; i < numHalfedge; i++) {
    tangent.data[i] = manifold_impl_is_inside_quad(impl, i)
        ? manifold_vec4(0, 0, 0, -1)
        : manifold_impl_tangent_from_normal(impl,
              vertNormal2.data[impl->halfedge.data[i].startVert], i);
  }
  impl->halfedgeTangent = tangent;

  // Copy sharpened edges and add flat face edges
  ManifoldVecSmoothness allSharp = MANIFOLD_VEC_INIT;
  for (int i = 0; i < numSharpened; i++)
    vec_smooth_push(&allSharp, sharpenedEdges[i]);

  // Add flat face edges
  size_t numTri = manifold_impl_num_tri(impl);
  for (size_t tri = 0; tri < numTri; tri++) {
    if (!triIsFlatFace.data[tri]) continue;
    for (int j = 0; j < 3; j++) {
      int tri2 = impl->halfedge.data[3 * tri + j].pairedHalfedge / 3;
      if (!triIsFlatFace.data[tri2] ||
          !manifold_triref_same_face(&impl->meshRelation.triRef.data[tri],
                                      &impl->meshRelation.triRef.data[tri2])) {
        ManifoldSmoothness s = {3 * tri + (size_t)j, 0};
        vec_smooth_push(&allSharp, s);
      }
    }
  }

  // Build edge pairs map using simple arrays
  // For each forward halfedge, store both smoothness values
  typedef struct {
    int idx;  // forward halfedge index
    double smoothFwd;
    double smoothRev;
  } EdgePair;

  int numEdges2 = (int)(impl->halfedge.len / 2);
  EdgePair *edgePairs = (EdgePair *)calloc(numEdges2, sizeof(EdgePair));
  int numPairs = 0;

  // Use a simple hash map approach: array indexed by forward halfedge index
  int *pairMap = (int *)malloc(numHalfedge * sizeof(int));
  for (int i = 0; i < numHalfedge; i++) pairMap[i] = -1;

  for (size_t i = 0; i < allSharp.len; i++) {
    ManifoldSmoothness s = allSharp.data[i];
    if (s.smoothness >= 1.0) continue;
    int he_idx = (int)s.halfedge;
    bool forward = manifold_halfedge_is_forward(&impl->halfedge.data[he_idx]);
    int pair_he = impl->halfedge.data[he_idx].pairedHalfedge;
    int fwd_idx = forward ? he_idx : pair_he;

    if (pairMap[fwd_idx] < 0) {
      pairMap[fwd_idx] = numPairs;
      edgePairs[numPairs].idx = fwd_idx;
      edgePairs[numPairs].smoothFwd = forward ? s.smoothness : 1.0;
      edgePairs[numPairs].smoothRev = forward ? 1.0 : s.smoothness;
      numPairs++;
    } else {
      int pi = pairMap[fwd_idx];
      if (forward) {
        if (s.smoothness < edgePairs[pi].smoothFwd)
          edgePairs[pi].smoothFwd = s.smoothness;
      } else {
        if (s.smoothness < edgePairs[pi].smoothRev)
          edgePairs[pi].smoothRev = s.smoothness;
      }
    }
  }

  // Build per-vertex sharp edge list
  typedef struct {
    int firstHE;
    double firstSmooth;
    int secondHE;
    double secondSmooth;
  } VertEdgePair;

  // vertTangents: list of VertEdgePair per vertex
  int *vertPairCount = (int *)calloc(numVert, sizeof(int));
  VertEdgePair *vertPairBuf = (VertEdgePair *)malloc(numPairs * 2 * sizeof(VertEdgePair));
  int *vertPairStart = (int *)malloc(numVert * sizeof(int));

  // First pass: count
  for (int i = 0; i < numPairs; i++) {
    int fwdHE = edgePairs[i].idx;
    int revHE = impl->halfedge.data[fwdHE].pairedHalfedge;
    int v0 = impl->halfedge.data[fwdHE].startVert;
    int v1 = impl->halfedge.data[revHE].startVert;
    vertPairCount[v0]++;
    vertPairCount[v1]++;
  }
  // Prefix sum
  {
    int sum = 0;
    for (int v = 0; v < numVert; v++) {
      vertPairStart[v] = sum;
      sum += vertPairCount[v];
      vertPairCount[v] = 0;
    }
  }
  // Second pass: fill
  for (int i = 0; i < numPairs; i++) {
    int fwdHE = edgePairs[i].idx;
    int revHE = impl->halfedge.data[fwdHE].pairedHalfedge;
    int v0 = impl->halfedge.data[fwdHE].startVert;
    int v1 = impl->halfedge.data[revHE].startVert;
    // v0 sees (fwd, rev)
    int off0 = vertPairStart[v0] + vertPairCount[v0]++;
    vertPairBuf[off0] = (VertEdgePair){fwdHE, edgePairs[i].smoothFwd,
                                        revHE, edgePairs[i].smoothRev};
    // v1 sees (rev, fwd)
    int off1 = vertPairStart[v1] + vertPairCount[v1]++;
    vertPairBuf[off1] = (VertEdgePair){revHE, edgePairs[i].smoothRev,
                                        fwdHE, edgePairs[i].smoothFwd};
  }

  for (int v = 0; v < numVert; v++) {
    int cnt = vertPairCount[v];
    if (cnt == 0) {
      fixedHalfedge.data[vertHE.data[v]] = true;
      continue;
    }
    if (cnt == 1) continue;  // smooth terminal
    VertEdgePair *pairs = &vertPairBuf[vertPairStart[v]];
    if (cnt == 2) {
      // Make continuous edge
      int first = pairs[0].firstHE;
      int second = pairs[1].firstHE;
      fixedHalfedge.data[first] = true;
      fixedHalfedge.data[second] = true;
      ManifoldVec3 newTangent = vec3_normalize(
          vec3_sub(vec4_to_vec3(impl->halfedgeTangent.data[first]),
                   vec4_to_vec3(impl->halfedgeTangent.data[second])));
      ManifoldVec3 pos = impl->vertPos.data[impl->halfedge.data[first].startVert];
      impl->halfedgeTangent.data[first] = circular_tangent(
          newTangent, vec3_sub(impl->vertPos.data[impl->halfedge.data[first].endVert], pos));
      impl->halfedgeTangent.data[second] = circular_tangent(
          vec3_neg(newTangent),
          vec3_sub(impl->vertPos.data[impl->halfedge.data[second].endVert], pos));

      double smoothness = (pairs[0].secondSmooth + pairs[1].firstSmooth) / 2.0;
      int cur = first;
      int iter = 0;
      do {
        cur = manifold_next_halfedge(impl->halfedge.data[cur].pairedHalfedge);
        if (cur == second) {
          smoothness = (pairs[1].secondSmooth + pairs[0].firstSmooth) / 2.0;
        } else if (cur != first && !manifold_impl_is_marked_inside_quad(impl, cur)) {
          sharpen_tangent(impl, cur, smoothness);
        }
        if (++iter > numHalfedge) break;
      } while (cur != first);
    } else {
      // Sharpen vertex uniformly
      double smoothness = 0;
      double denom = 0;
      for (int i = 0; i < cnt; i++) {
        smoothness += pairs[i].firstSmooth + pairs[i].secondSmooth;
        denom += (pairs[i].firstSmooth == 0 ? 0 : 1) +
                 (pairs[i].secondSmooth == 0 ? 0 : 1);
      }
      if (denom > 0) smoothness /= denom;

      int cur = pairs[0].firstHE;
      int iter = 0;
      do {
        cur = manifold_next_halfedge(impl->halfedge.data[cur].pairedHalfedge);
        if (!manifold_impl_is_marked_inside_quad(impl, cur)) {
          int pair_he = impl->halfedge.data[cur].pairedHalfedge;
          sharpen_tangent(impl, cur,
              triIsFlatFace.data[cur / 3] || triIsFlatFace.data[pair_he / 3]
                  ? 0 : smoothness);
        }
        if (++iter > numHalfedge) break;
      } while (cur != pairs[0].firstHE);
    }
  }

  linearize_flat_tangents(impl);
  distribute_tangents(impl, &fixedHalfedge);

  free(pairMap);
  free(edgePairs);
  free(vertPairCount);
  free(vertPairBuf);
  free(vertPairStart);
  vec_smooth_free(&allSharp);
  vec_bool_free(&fixedHalfedge);
  vec_int_free(&vertHE);
  vec_bool_free(&triIsFlatFace);
  vec_int_free(&vertFlatFace_data);
  vec_vec3_free(&vertNormal2);
}

// UpdateSharpenedEdges: remap halfedge indices after sorting
ManifoldVecSmoothness manifold_impl_update_sharpened_edges(
    const ManifoldImpl *impl, const ManifoldSmoothness *edges, int numEdges) {
  // Build old->new halfedge mapping using faceID
  int *oldHE2New = (int *)calloc(impl->halfedge.len, sizeof(int));
  size_t numTri = manifold_impl_num_tri(impl);
  for (size_t tri = 0; tri < numTri; tri++) {
    int oldTri = impl->meshRelation.triRef.data[tri].faceID;
    for (int i = 0; i < 3; i++) {
      if (oldTri >= 0 && (size_t)(3 * oldTri + i) < impl->halfedge.len)
        oldHE2New[3 * oldTri + i] = (int)(3 * tri + i);
    }
  }

  ManifoldVecSmoothness result = MANIFOLD_VEC_INIT;
  for (int i = 0; i < numEdges; i++) {
    ManifoldSmoothness s = edges[i];
    if (s.halfedge < impl->halfedge.len)
      s.halfedge = oldHE2New[s.halfedge];
    vec_smooth_push(&result, s);
  }
  free(oldHE2New);
  return result;
}
