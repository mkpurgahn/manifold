// Copyright 2021 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// Full boolean algorithm for the C11 Manifold port.
// Ported from boolean3.cpp and boolean_result.cpp.

#include "manifold_boolean.h"
#include "manifold_polygon.h"
#include <string.h>
#include <limits.h>
#include <stdio.h>

// ============== Intersection helpers ==============

static inline double bool_with_sign(bool pos, double v) {
  return pos ? v : -v;
}

// Interpolate along an edge at a given x-coordinate
static inline ManifoldVec2 bool_interpolate(ManifoldVec3 aL, ManifoldVec3 aR,
                                             double x) {
  double dxL = x - aL.x;
  double dxR = x - aR.x;
  bool useL = fabs(dxL) < fabs(dxR);
  ManifoldVec3 dLR = vec3_sub(aR, aL);
  double lambda = (useL ? dxL : dxR) / dLR.x;
  if (!isfinite(lambda) || !isfinite(dLR.y) || !isfinite(dLR.z))
    return manifold_vec2(aL.y, aL.z);
  ManifoldVec2 yz;
  yz.x = lambda * dLR.y + (useL ? aL.y : aR.y);
  yz.y = lambda * dLR.z + (useL ? aL.z : aR.z);
  return yz;
}

// Intersect two interpolated lines
static inline ManifoldVec4 bool_intersect(ManifoldVec3 aL, ManifoldVec3 aR,
                                           ManifoldVec3 bL, ManifoldVec3 bR) {
  double dyL = bL.y - aL.y;
  double dyR = bR.y - aR.y;
  bool useL = fabs(dyL) < fabs(dyR);
  double dx = aR.x - aL.x;
  double lambda = (useL ? dyL : dyR) / (dyL - dyR);
  if (!isfinite(lambda)) lambda = 0.0;
  ManifoldVec4 xyzz;
  xyzz.x = lambda * dx + (useL ? aL.x : aR.x);
  double aDy = aR.y - aL.y;
  double bDy = bR.y - bL.y;
  bool useA = fabs(aDy) < fabs(bDy);
  xyzz.y = lambda * (useA ? aDy : bDy) +
            (useL ? (useA ? aL.y : bL.y) : (useA ? aR.y : bR.y));
  xyzz.z = lambda * (aR.z - aL.z) + (useL ? aL.z : aR.z);
  xyzz.w = lambda * (bR.z - bL.z) + (useL ? bL.z : bR.z);
  return xyzz;
}

static inline bool bool_shadows(double p, double q, double dir) {
  return p == q ? dir < 0 : p < q;
}

// Shadow01: test vertex a0 against edge b1
// expandP: symbolic perturbation direction
// forward: whether computing p→q or q→p
typedef struct {
  int s01;
  ManifoldVec2 yz01;
} Shadow01Result;

static Shadow01Result shadow01(int a0, int b1,
                                const ManifoldImpl *inA,
                                const ManifoldImpl *inB,
                                bool expandP, bool forward) {
  Shadow01Result r = {0, manifold_vec2(NAN, NAN)};
  if (b1 < 0 || (size_t)b1 >= inB->halfedge.len) return r;
  int b1s = inB->halfedge.data[b1].startVert;
  int b1e = inB->halfedge.data[b1].endVert;
  if (b1s < 0 || (size_t)b1s >= inB->vertPos.len ||
      b1e < 0 || (size_t)b1e >= inB->vertPos.len ||
      a0 < 0 || (size_t)a0 >= inA->vertPos.len) return r;
  if ((size_t)a0 >= inA->vertNormal.len ||
      (size_t)b1s >= inB->vertNormal.len ||
      (size_t)b1e >= inB->vertNormal.len) return r;
  double a0x = inA->vertPos.data[a0].x;
  double b1sx = inB->vertPos.data[b1s].x;
  double b1ex = inB->vertPos.data[b1e].x;
  double a0xp = inA->vertNormal.data[a0].x;
  double b1sxp = inB->vertNormal.data[b1s].x;
  double b1exp = inB->vertNormal.data[b1e].x;

  if (forward) {
    r.s01 = (int)bool_shadows(a0x, b1ex, bool_with_sign(expandP, a0xp) - b1exp) -
            (int)bool_shadows(a0x, b1sx, bool_with_sign(expandP, a0xp) - b1sxp);
  } else {
    r.s01 = (int)bool_shadows(b1sx, a0x, bool_with_sign(expandP, b1sxp) - a0xp) -
            (int)bool_shadows(b1ex, a0x, bool_with_sign(expandP, b1exp) - a0xp);
  }

  if (r.s01 != 0) {
    r.yz01 = bool_interpolate(inB->vertPos.data[b1s],
                               inB->vertPos.data[b1e],
                               inA->vertPos.data[a0].x);
    int b1pair = inB->halfedge.data[b1].pairedHalfedge;
    if (b1pair < 0 || (size_t)(b1pair / 3) >= inB->faceNormal.len) {
      r.s01 = 0;
      return r;
    }
    double dir = inB->faceNormal.data[b1 / 3].y +
                 inB->faceNormal.data[b1pair / 3].y;
    if (forward) {
      if (!bool_shadows(inA->vertPos.data[a0].y, r.yz01.x, -dir))
        r.s01 = 0;
    } else {
      if (!bool_shadows(r.yz01.x, inA->vertPos.data[a0].y,
                         bool_with_sign(expandP, dir)))
        r.s01 = 0;
    }
  }
  return r;
}

// Kernel11: edge-edge intersection
typedef struct {
  int s11;
  ManifoldVec4 xyzz11;
} Kernel11Result;

static Kernel11Result kernel11(int p1, int q1,
                                const ManifoldImpl *inP,
                                const ManifoldImpl *inQ,
                                bool expandP) {
  Kernel11Result res = {0, manifold_vec4(NAN, NAN, NAN, NAN)};
  int k = 0;
  ManifoldVec3 pRL[2], qRL[2];
  bool shadows = false;

  int p0[2] = {inP->halfedge.data[p1].startVert,
               inP->halfedge.data[p1].endVert};
  for (int i = 0; i < 2; i++) {
    Shadow01Result sr = shadow01(p0[i], q1, inP, inQ, expandP, true);
    if (isfinite(sr.yz01.x)) {
      res.s11 += sr.s01 * (i == 0 ? -1 : 1);
      if (k < 2 && (k == 0 || (sr.s01 != 0) != shadows)) {
        shadows = sr.s01 != 0;
        pRL[k] = inP->vertPos.data[p0[i]];
        qRL[k] = manifold_vec3(pRL[k].x, sr.yz01.x, sr.yz01.y);
        k++;
      }
    }
  }

  int q0[2] = {inQ->halfedge.data[q1].startVert,
               inQ->halfedge.data[q1].endVert};
  for (int i = 0; i < 2; i++) {
    Shadow01Result sr = shadow01(q0[i], p1, inQ, inP, expandP, false);
    if (isfinite(sr.yz01.x)) {
      res.s11 += sr.s01 * (i == 0 ? -1 : 1);
      if (k < 2 && (k == 0 || (sr.s01 != 0) != shadows)) {
        shadows = sr.s01 != 0;
        qRL[k] = inQ->vertPos.data[q0[i]];
        pRL[k] = manifold_vec3(qRL[k].x, sr.yz01.x, sr.yz01.y);
        k++;
      }
    }
  }

  if (res.s11 == 0) {
    res.xyzz11 = manifold_vec4(NAN, NAN, NAN, NAN);
  } else {
    if (k != 2) { res.s11 = 0; return res; }
    res.xyzz11 = bool_intersect(pRL[0], pRL[1], qRL[0], qRL[1]);
    int p1pair = inP->halfedge.data[p1].pairedHalfedge;
    int q1pair = inQ->halfedge.data[q1].pairedHalfedge;
    if (p1pair < 0 || (size_t)(p1pair / 3) >= inP->faceNormal.len ||
        q1pair < 0 || (size_t)(q1pair / 3) >= inQ->faceNormal.len) {
      res.s11 = 0; return res;
    }
    double dirP = inP->faceNormal.data[p1 / 3].z +
                  inP->faceNormal.data[p1pair / 3].z;
    double dirQ = inQ->faceNormal.data[q1 / 3].z +
                  inQ->faceNormal.data[q1pair / 3].z;
    if (!bool_shadows(res.xyzz11.z, res.xyzz11.w,
                       bool_with_sign(expandP, dirP) - dirQ))
      res.s11 = 0;
  }
  return res;
}

// Kernel02: vertex-face intersection
typedef struct {
  int s02;
  double z02;
} Kernel02Result;

static Kernel02Result kernel02(int a0, int b2,
                                const ManifoldImpl *inA,
                                const ManifoldImpl *inB,
                                bool expandP, bool forward) {
  Kernel02Result res = {0, 0.0};
  int k = 0;
  ManifoldVec3 yzzRL[2];
  bool shadows_val = false;

  for (int i = 0; i < 3; i++) {
    int b1 = 3 * b2 + i;
    ManifoldHalfedge edgeB = inB->halfedge.data[b1];
    int b1F = manifold_halfedge_is_forward(&edgeB) ? b1 : edgeB.pairedHalfedge;
    if (b1F < 0 || (size_t)b1F >= inB->halfedge.len) continue;

    Shadow01Result sr = shadow01(a0, b1F, inA, inB, expandP, forward);
    if (isfinite(sr.yz01.x)) {
      res.s02 += sr.s01 * ((forward == manifold_halfedge_is_forward(&edgeB)) ? -1 : 1);
      if (k < 2 && (k == 0 || (sr.s01 != 0) != shadows_val)) {
        shadows_val = sr.s01 != 0;
        yzzRL[k] = manifold_vec3(sr.yz01.x, sr.yz01.y, sr.yz01.y);
        k++;
      }
    }
  }

  if (res.s02 == 0) {
    res.z02 = NAN;
  } else {
    if (k != 2) { res.s02 = 0; res.z02 = NAN; return res; }
    ManifoldVec3 vertPosA = inA->vertPos.data[a0];
    res.z02 = bool_interpolate(yzzRL[0], yzzRL[1], vertPosA.y).y;
    if (forward) {
      if (!bool_shadows(vertPosA.z, res.z02, -inB->faceNormal.data[b2].z))
        res.s02 = 0;
    } else {
      if (!bool_shadows(res.z02, vertPosA.z,
                         bool_with_sign(expandP, inB->faceNormal.data[b2].z)))
        res.s02 = 0;
    }
  }
  return res;
}

// Kernel12: edge-face intersection
typedef struct {
  int x12;
  ManifoldVec3 v12;
} Kernel12Result;

static Kernel12Result kernel12(int a1, int b2,
                                const ManifoldImpl *inA,
                                const ManifoldImpl *inB,
                                const ManifoldImpl *inP,
                                const ManifoldImpl *inQ,
                                bool expandP, bool forward) {
  Kernel12Result res = {0, manifold_vec3(NAN, NAN, NAN)};
  int k = 0;
  ManifoldVec3 xzyLR0[2], xzyLR1[2];
  bool shadows_val = false;

  ManifoldHalfedge edgeA = inA->halfedge.data[a1];
  int vertAList[2] = {edgeA.startVert, edgeA.endVert};

  for (int vi = 0; vi < 2; vi++) {
    int vertA = vertAList[vi];
    Kernel02Result kr = kernel02(vertA, b2, inA, inB, expandP, forward);
    if (isfinite(kr.z02)) {
      res.x12 += kr.s02 * ((vertA == edgeA.startVert) == forward ? 1 : -1);
      if (k < 2 && (k == 0 || (kr.s02 != 0) != shadows_val)) {
        shadows_val = kr.s02 != 0;
        xzyLR0[k] = inA->vertPos.data[vertA];
        double tmpY = xzyLR0[k].y;
        xzyLR0[k].y = xzyLR0[k].z;
        xzyLR0[k].z = tmpY;
        xzyLR1[k] = xzyLR0[k];
        xzyLR1[k].y = kr.z02;
        k++;
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    int b1 = 3 * b2 + i;
    ManifoldHalfedge edgeB = inB->halfedge.data[b1];
    int b1F = manifold_halfedge_is_forward(&edgeB) ? b1 : edgeB.pairedHalfedge;
    if (b1F < 0 || (size_t)b1F >= inB->halfedge.len) continue;
    Kernel11Result kr;
    if (forward)
      kr = kernel11(a1, b1F, inP, inQ, expandP);
    else
      kr = kernel11(b1F, a1, inP, inQ, expandP);
    if (isfinite(kr.xyzz11.x)) {
      res.x12 -= kr.s11 * (manifold_halfedge_is_forward(&edgeB) ? 1 : -1);
      if (k < 2 && (k == 0 || (kr.s11 != 0) != shadows_val)) {
        shadows_val = kr.s11 != 0;
        xzyLR0[k].x = kr.xyzz11.x;
        xzyLR0[k].y = kr.xyzz11.z;
        xzyLR0[k].z = kr.xyzz11.y;
        xzyLR1[k] = xzyLR0[k];
        xzyLR1[k].y = kr.xyzz11.w;
        if (!forward) {
          double tmp = xzyLR0[k].y;
          xzyLR0[k].y = xzyLR1[k].y;
          xzyLR1[k].y = tmp;
        }
        k++;
      }
    }
  }

  if (res.x12 == 0) {
    res.v12 = manifold_vec3(NAN, NAN, NAN);
  } else {
    if (k != 2) { res.x12 = 0; return res; }
    ManifoldVec4 xzyy = bool_intersect(xzyLR0[0], xzyLR0[1],
                                        xzyLR1[0], xzyLR1[1]);
    res.v12.x = xzyy.x;
    res.v12.y = xzyy.z;
    res.v12.z = xzyy.y;
  }
  return res;
}

// ============== Collision callback for edge-face intersections ==============

typedef struct {
  const ManifoldImpl *inP;
  const ManifoldImpl *inQ;
  const ManifoldImpl *inA;
  const ManifoldImpl *inB;
  bool expandP;
  bool forward;
  ManifoldIntersections *result;
} Intersect12Ctx;

static void intersect12_callback(int queryIdx, int leafIdx, void *ctx) {
  Intersect12Ctx *c = (Intersect12Ctx *)ctx;
  int a1 = queryIdx;
  int b2 = leafIdx;

  Kernel12Result kr = kernel12(a1, b2, c->inA, c->inB,
                                c->inP, c->inQ, c->expandP, c->forward);
  if (isfinite(kr.v12.x)) {
    ManifoldIntArr2 pair;
    if (c->forward) {
      pair.v[0] = a1;
      pair.v[1] = b2;
    } else {
      pair.v[0] = b2;
      pair.v[1] = a1;
    }
    vec_intarr2_push(&c->result->p1q2, pair);
    vec_int_push(&c->result->x12, kr.x12);
    vec_vec3_push(&c->result->v12, kr.v12);
  }
}

// Build edge bounding boxes for one mesh
static ManifoldVecBox build_edge_boxes(const ManifoldImpl *impl) {
  size_t nHe = impl->halfedge.len;
  ManifoldVecBox boxes = vec_box_create_n(nHe);
  for (size_t i = 0; i < nHe; i++) {
    ManifoldHalfedge he = impl->halfedge.data[i];
    if (manifold_halfedge_is_forward(&he)) {
      ManifoldVec3 s = impl->vertPos.data[he.startVert];
      ManifoldVec3 e = impl->vertPos.data[he.endVert];
      boxes.data[i].min = vec3_min(s, e);
      boxes.data[i].max = vec3_max(s, e);
    } else {
      boxes.data[i] = manifold_box_empty();
    }
  }
  return boxes;
}

// Compute edge-face intersections
static ManifoldIntersections intersect12(const ManifoldImpl *inP,
                                          const ManifoldImpl *inQ,
                                          bool expandP, bool forward) {
  ManifoldIntersections result = {{0}, {0}, {0}};

  const ManifoldImpl *a = forward ? inP : inQ;
  const ManifoldImpl *b = forward ? inQ : inP;

  if (a->halfedge.len == 0 || !b->colliderBuilt) return result;

  // Build query boxes from edges of A
  ManifoldVecBox edgeBoxes = build_edge_boxes(a);

  Intersect12Ctx ctx = {inP, inQ, a, b, expandP, forward, &result};

  // Query each edge of A against B's face collider
  for (size_t qi = 0; qi < edgeBoxes.len; qi++) {
    if (!manifold_halfedge_is_forward(&a->halfedge.data[qi])) continue;
    ManifoldBox qbox = edgeBoxes.data[qi];
    if (qbox.min.x > qbox.max.x) continue;  // empty

    // Walk the BVH
    if (b->collider.internalChildren.len == 0) continue;
    int stack[64];
    int top = -1;
    int node = MANIFOLD_COLLIDER_ROOT;

    while (1) {
      int internal = collider_node2internal(node);
      int child1 = b->collider.internalChildren.data[internal].first;
      int child2 = b->collider.internalChildren.data[internal].second;

      bool overlap1 = manifold_box_overlaps(b->collider.nodeBBox.data[child1], qbox);
      bool overlap2 = manifold_box_overlaps(b->collider.nodeBBox.data[child2], qbox);

      bool traverse1 = false, traverse2 = false;
      if (overlap1 && collider_is_leaf(child1)) {
        int leafIdx = collider_node2leaf(child1);
        intersect12_callback((int)qi, leafIdx, &ctx);
      } else {
        traverse1 = overlap1 && collider_is_internal(child1);
      }
      if (overlap2 && collider_is_leaf(child2)) {
        int leafIdx = collider_node2leaf(child2);
        intersect12_callback((int)qi, leafIdx, &ctx);
      } else {
        traverse2 = overlap2 && collider_is_internal(child2);
      }

      if (!traverse1 && !traverse2) {
        if (top < 0) break;
        node = stack[top--];
      } else {
        node = traverse1 ? child1 : child2;
        if (traverse1 && traverse2) {
          if (top < 63) stack[++top] = child2;
        }
      }
    }
  }

  vec_box_free(&edgeBoxes);

  // Sort by edge index
  if (result.p1q2.len > 0) {
    int index = forward ? 0 : 1;
    size_t n = result.p1q2.len;
    ManifoldVecSize order = vec_size_create_n(n);
    vec_size_sequence(&order);

    // Simple insertion sort (good enough for typical intersection counts)
    for (size_t i = 1; i < n; i++) {
      size_t key = order.data[i];
      int keyVal0 = result.p1q2.data[key].v[index];
      int keyVal1 = result.p1q2.data[key].v[1 - index];
      size_t j = i;
      while (j > 0) {
        size_t prev = order.data[j - 1];
        int pv0 = result.p1q2.data[prev].v[index];
        int pv1 = result.p1q2.data[prev].v[1 - index];
        if (pv0 < keyVal0 || (pv0 == keyVal0 && pv1 <= keyVal1)) break;
        order.data[j] = order.data[j - 1];
        j--;
      }
      order.data[j] = key;
    }

    // Apply permutation
    ManifoldVecIntArr2 newP1q2 = vec_intarr2_create_n(n);
    ManifoldVecInt newX12 = vec_int_create_n(n);
    ManifoldVecVec3 newV12 = vec_vec3_create_n(n);
    for (size_t i = 0; i < n; i++) {
      size_t idx = order.data[i];
      newP1q2.data[i] = result.p1q2.data[idx];
      newX12.data[i] = result.x12.data[idx];
      newV12.data[i] = result.v12.data[idx];
    }
    vec_intarr2_free(&result.p1q2);
    vec_int_free(&result.x12);
    vec_vec3_free(&result.v12);
    result.p1q2 = newP1q2;
    result.x12 = newX12;
    result.v12 = newV12;
    vec_size_free(&order);
  }

  return result;
}

// ============== Winding number computation ==============

typedef struct {
  const ManifoldImpl *inA;
  const ManifoldImpl *inB;
  bool expandP;
  bool forward;
  int *verts;
  int *w03;
} Winding03Ctx;

static void winding03_callback(int queryIdx, int leafIdx, void *ctx) {
  Winding03Ctx *c = (Winding03Ctx *)ctx;
  Kernel02Result kr = kernel02(c->verts[queryIdx], leafIdx,
                                c->inA, c->inB, c->expandP, c->forward);
  if (isfinite(kr.z02)) {
    c->w03[c->verts[queryIdx]] += kr.s02 * (c->forward ? 1 : -1);
  }
}

static ManifoldVecInt winding03(const ManifoldImpl *inP,
                                 const ManifoldImpl *inQ,
                                 const ManifoldIntersections *xv,
                                 bool expandP, bool forward) {
  const ManifoldImpl *a = forward ? inP : inQ;
  const ManifoldImpl *b = forward ? inQ : inP;
  int index = forward ? 0 : 1;

  size_t numVertA = a->vertPos.len;
  ManifoldVecInt w03 = vec_int_create_fill(numVertA, 0);

  // Use disjoint sets to find connected components
  ManifoldDisjointSets uf = manifold_disjoint_sets_create((uint32_t)numVertA);

  // Unite vertices along unbroken edges
  for (size_t edge = 0; edge < a->halfedge.len; edge++) {
    ManifoldHalfedge he = a->halfedge.data[edge];
    if (!manifold_halfedge_is_forward(&he)) continue;

    // Check if edge is broken (has intersection)
    bool broken = false;
    for (size_t j = 0; j < xv->p1q2.len; j++) {
      if (xv->p1q2.data[j].v[index] == (int)edge) {
        broken = true;
        break;
      }
      if (xv->p1q2.data[j].v[index] > (int)edge) break;
    }
    if (!broken) {
      manifold_disjoint_sets_unite(&uf, he.startVert, he.endVert);
    }
  }

  // Find component representatives
  ManifoldVecInt verts = {0};
  bool *seen = (bool *)calloc(numVertA, sizeof(bool));
  for (size_t v = 0; v < numVertA; v++) {
    int root = manifold_disjoint_sets_find(&uf, (int)v);
    if (!seen[root]) {
      seen[root] = true;
      vec_int_push(&verts, root);
    }
  }

  // Compute winding number for each component rep
  if (b->colliderBuilt && b->collider.internalChildren.len > 0) {
    Winding03Ctx ctx = {a, b, expandP, forward, verts.data, w03.data};

    for (size_t qi = 0; qi < verts.len; qi++) {
      int vert = verts.data[qi];
      ManifoldVec3 pos = a->vertPos.data[vert];

      // Walk the BVH to find faces whose x-y bounding box covers the vertex.
      // Use z-projected overlap (x-y only) matching C++ DoesOverlap(vec3).
      int stack[64];
      int top2 = -1;
      int node = MANIFOLD_COLLIDER_ROOT;

      while (1) {
        int internal = collider_node2internal(node);
        int child1 = b->collider.internalChildren.data[internal].first;
        int child2 = b->collider.internalChildren.data[internal].second;

        bool overlap1 = manifold_box_overlaps_point(b->collider.nodeBBox.data[child1], pos);
        bool overlap2 = manifold_box_overlaps_point(b->collider.nodeBBox.data[child2], pos);

        bool traverse1 = false, traverse2 = false;
        if (overlap1 && collider_is_leaf(child1)) {
          int leafIdx = collider_node2leaf(child1);
          winding03_callback((int)qi, leafIdx, &ctx);
        } else {
          traverse1 = overlap1 && collider_is_internal(child1);
        }
        if (overlap2 && collider_is_leaf(child2)) {
          int leafIdx = collider_node2leaf(child2);
          winding03_callback((int)qi, leafIdx, &ctx);
        } else {
          traverse2 = overlap2 && collider_is_internal(child2);
        }

        if (!traverse1 && !traverse2) {
          if (top2 < 0) break;
          node = stack[top2--];
        } else {
          node = traverse1 ? child1 : child2;
          if (traverse1 && traverse2) {
            if (top2 < 63) stack[++top2] = child2;
          }
        }
      }
    }
  }

  // Flood fill: propagate winding numbers from roots to all verts
  for (size_t v = 0; v < numVertA; v++) {
    int root = manifold_disjoint_sets_find(&uf, (int)v);
    if ((int)v != root) {
      w03.data[v] = w03.data[root];
    }
  }

  free(seen);
  vec_int_free(&verts);
  manifold_disjoint_sets_free(&uf);
  return w03;
}

// ============== Edge position for result construction ==============

typedef struct {
  double edgePos;
  int vert;
  int collisionId;
  bool isStart;
} EdgePos;

static int edgepos_cmp(const void *a, const void *b) {
  const EdgePos *ea = (const EdgePos *)a;
  const EdgePos *eb = (const EdgePos *)b;
  if (ea->edgePos < eb->edgePos) return -1;
  if (ea->edgePos > eb->edgePos) return 1;
  if (ea->collisionId < eb->collisionId) return -1;
  if (ea->collisionId > eb->collisionId) return 1;
  return 0;
}

MANIFOLD_VEC_STRUCT(EdgePos, VecEdgePos)
MANIFOLD_VEC_FUNCS(EdgePos, VecEdgePos, vec_edgepos)

// Simple key-value map for edges: int → VecEdgePos
typedef struct {
  int key;
  VecEdgePos value;
} EdgeMapEntry;

MANIFOLD_VEC_STRUCT(EdgeMapEntry, EdgeMap)
MANIFOLD_VEC_FUNCS(EdgeMapEntry, EdgeMap, vec_edgemap)

static VecEdgePos *edgemap_get_or_create(EdgeMap *map, int key) {
  for (size_t i = 0; i < map->len; i++) {
    if (map->data[i].key == key) return &map->data[i].value;
  }
  EdgeMapEntry entry;
  entry.key = key;
  entry.value = (VecEdgePos){0};
  vec_edgemap_push(map, entry);
  return &map->data[map->len - 1].value;
}

// Pair key for new edges: (faceP, faceQ) → VecEdgePos
typedef struct {
  int faceP, faceQ;
  VecEdgePos value;
} NewEdgeEntry;

MANIFOLD_VEC_STRUCT(NewEdgeEntry, NewEdgeMap)
MANIFOLD_VEC_FUNCS(NewEdgeEntry, NewEdgeMap, vec_newedgemap)

static VecEdgePos *newedgemap_get_or_create(NewEdgeMap *map,
                                             int faceP, int faceQ) {
  for (size_t i = 0; i < map->len; i++) {
    if (map->data[i].faceP == faceP && map->data[i].faceQ == faceQ)
      return &map->data[i].value;
  }
  NewEdgeEntry entry;
  entry.faceP = faceP;
  entry.faceQ = faceQ;
  entry.value = (VecEdgePos){0};
  vec_newedgemap_push(map, entry);
  return &map->data[map->len - 1].value;
}

// ============== Result construction ==============

// Add new edge verts from intersections
static void add_new_edge_verts(EdgeMap *edgesP,
                                NewEdgeMap *edgesNew,
                                const ManifoldIntersections *xv,
                                const ManifoldVecInt *i12,
                                const ManifoldVecInt *v12R,
                                const ManifoldVecHalfedge *halfedgeP,
                                bool forward, size_t offset) {
  for (size_t i = 0; i < xv->p1q2.len; i++) {
    int edgeP = xv->p1q2.data[i].v[forward ? 0 : 1];
    int faceQ = xv->p1q2.data[i].v[forward ? 1 : 0];
    int vert = v12R->data[i];
    int inclusion = i12->data[i];

    ManifoldHalfedge halfedge = halfedgeP->data[edgeP];
    int keyRightFace = (halfedge.pairedHalfedge >= 0) ? halfedge.pairedHalfedge / 3 : edgeP / 3;
    int keyLeftFace = edgeP / 3;

    int newKeyRightP, newKeyRightQ, newKeyLeftP, newKeyLeftQ;
    if (forward) {
      newKeyRightP = keyRightFace;
      newKeyRightQ = faceQ;
      newKeyLeftP = keyLeftFace;
      newKeyLeftQ = faceQ;
    } else {
      newKeyRightP = faceQ;
      newKeyRightQ = keyRightFace;
      newKeyLeftP = faceQ;
      newKeyLeftQ = keyLeftFace;
    }

    bool direction = inclusion < 0;

    // Add to edgesP
    VecEdgePos *edgePList = edgemap_get_or_create(edgesP, edgeP);
    for (int j = 0; j < abs(inclusion); j++) {
      EdgePos ep = {0.0, vert + j, (int)(i + offset), direction};
      vec_edgepos_push(edgePList, ep);
    }

    // Add to right new edge
    VecEdgePos *rightList = newedgemap_get_or_create(edgesNew,
                                                      newKeyRightP,
                                                      newKeyRightQ);
    bool dirRight = direction ^ !forward;
    for (int j = 0; j < abs(inclusion); j++) {
      EdgePos ep = {0.0, vert + j, (int)(i + offset), dirRight};
      vec_edgepos_push(rightList, ep);
    }

    // Add to left new edge
    VecEdgePos *leftList = newedgemap_get_or_create(edgesNew,
                                                     newKeyLeftP,
                                                     newKeyLeftQ);
    bool dirLeft = direction ^ forward;
    for (int j = 0; j < abs(inclusion); j++) {
      EdgePos ep = {0.0, vert + j, (int)(i + offset), dirLeft};
      vec_edgepos_push(leftList, ep);
    }
  }
}

// ============== Face2Tri: triangulate general faces ==============

// Assemble halfedges into polygon loops (matching C++ AssembleHalfedges).
// Returns number of loops. loops[i] contains edge local indices for loop i,
// loopLens[i] is the length. Caller must free loops[i] and loops/loopLens.
static int assemble_halfedge_loops(const ManifoldHalfedge *edges, int numEdge,
                                    int ***loops_out, int **loopLens_out) {
  // Build multimap: startVert → local edge index
  // Use a simple array of {key, value, next} linked list per vertex
  int *nextEntry = (int *)malloc((size_t)numEdge * sizeof(int));
  int *entryEdge = (int *)malloc((size_t)numEdge * sizeof(int));
  int *entryKey = (int *)malloc((size_t)numEdge * sizeof(int));
  // Hash table: vert → first entry index (-1 = empty)
  int hashSize = numEdge * 2 + 1;
  int *hashHead = (int *)malloc((size_t)hashSize * sizeof(int));
  for (int i = 0; i < hashSize; i++) hashHead[i] = -1;

  for (int i = 0; i < numEdge; i++) {
    int sv = edges[i].startVert;
    int bucket = ((unsigned)sv) % (unsigned)hashSize;
    entryKey[i] = sv;
    entryEdge[i] = i;
    nextEntry[i] = hashHead[bucket];
    hashHead[bucket] = i;
  }

  int numLoops = 0;
  int loopCap = 4;
  int **loops = (int **)malloc((size_t)loopCap * sizeof(int *));
  int *loopLens = (int *)malloc((size_t)loopCap * sizeof(int));

  while (1) {
    // Find first remaining edge
    int startEdge = -1;
    for (int i = 0; i < numEdge; i++) {
      if (entryKey[i] >= 0) { startEdge = i; break; }
    }
    if (startEdge < 0) break;

    // Start a new loop
    int *loop = (int *)malloc((size_t)numEdge * sizeof(int));
    int loopLen = 0;
    int thisEdge = startEdge;

    // Follow the C++ approach: don't erase startEdge until we circle back
    // Add startEdge to the loop and mark it used
    loop[loopLen++] = startEdge;
    entryKey[startEdge] = -1;

    int endV = edges[startEdge].endVert;
    while (endV != edges[startEdge].startVert) {
      int bucket = ((unsigned)endV) % (unsigned)hashSize;
      int prev = -1;
      int cur = hashHead[bucket];
      int found = -1;
      while (cur >= 0) {
        if (entryKey[cur] == endV) {
          found = entryEdge[cur];
          // Remove from hash chain
          if (prev < 0) hashHead[bucket] = nextEntry[cur];
          else nextEntry[prev] = nextEntry[cur];
          break;
        }
        prev = cur;
        cur = nextEntry[cur];
      }
      if (found < 0) break;  // broken loop
      loop[loopLen++] = found;
      entryKey[found] = -1;
      endV = edges[found].endVert;
      if (loopLen > numEdge) break;  // safety
    }

    if (numLoops >= loopCap) {
      loopCap *= 2;
      loops = (int **)realloc(loops, (size_t)loopCap * sizeof(int *));
      loopLens = (int *)realloc(loopLens, (size_t)loopCap * sizeof(int));
    }
    loops[numLoops] = loop;
    loopLens[numLoops] = loopLen;
    numLoops++;
  }

  free(nextEntry); free(entryEdge); free(entryKey); free(hashHead);
  *loops_out = loops;
  *loopLens_out = loopLens;
  return numLoops;
}

static void face2tri(ManifoldImpl *impl, const ManifoldVecInt *faceEdge,
                      const ManifoldVecTriRef *halfedgeRef, bool allowConvex) {
  (void)allowConvex;
  ManifoldVecIVec3 triVerts = {0};
  ManifoldVecVec3 triNormal = {0};
  ManifoldVecIVec3 triProp = {0};
  ManifoldVecTriRef triRef = {0};

  size_t numFaces = faceEdge->len - 1;
  for (size_t face = 0; face < numFaces; face++) {
    int firstEdge = faceEdge->data[face];
    int lastEdge = faceEdge->data[face + 1];
    int numEdge = lastEdge - firstEdge;
    if (numEdge < 3) continue;

    ManifoldVec3 normal = impl->faceNormal.data[face];
    ManifoldTriRef ref = halfedgeRef->data[firstEdge];

    if (numEdge == 3) {
      // Already a triangle - check orientation
      ManifoldHalfedge *he = impl->halfedge.data + firstEdge;
      ManifoldIVec3 tri = manifold_ivec3(he[0].startVert, he[1].startVert,
                                          he[2].startVert);
      ManifoldIVec3 prop = manifold_ivec3(he[0].propVert, he[1].propVert,
                                           he[2].propVert);

      // Verify edges form a cycle, swap if needed
      if (he[0].endVert == tri.z) {
        int tmp;
        tmp = tri.y; tri.y = tri.z; tri.z = tmp;
        tmp = prop.y; prop.y = prop.z; prop.z = tmp;
      }

      vec_ivec3_push(&triVerts, tri);
      vec_ivec3_push(&triProp, prop);
      vec_vec3_push(&triNormal, normal);
      vec_triref_push(&triRef, ref);
    } else {
      // General case: assemble into loops, project, and triangulate
      ManifoldMat2x3 projection = manifold_get_axis_aligned_projection(normal);

      // Assemble halfedges into polygon loops
      int **loops = NULL;
      int *loopLens = NULL;
      int numLoops = assemble_halfedge_loops(
          impl->halfedge.data + firstEdge, numEdge, &loops, &loopLens);

      if (numLoops == 0) {
        free(loops); free(loopLens);
        continue;
      }

      // Project all loops into 2D and build per-loop point arrays
      // Also build a global index array mapping flat position → halfedge index
      int totalPts = 0;
      for (int li = 0; li < numLoops; li++) totalPts += loopLens[li];

      ManifoldVec2 *allPts = (ManifoldVec2 *)malloc((size_t)totalPts * sizeof(ManifoldVec2));
      int *heMap = (int *)malloc((size_t)totalPts * sizeof(int));  // flat idx → local edge idx
      int *polySizes = (int *)malloc((size_t)numLoops * sizeof(int));

      int flatIdx = 0;
      for (int li = 0; li < numLoops; li++) {
        polySizes[li] = loopLens[li];
        for (int pi = 0; pi < loopLens[li]; pi++) {
          int localEdge = loops[li][pi];
          int sv = impl->halfedge.data[firstEdge + localEdge].startVert;
          if (sv < 0 || (size_t)sv >= impl->vertPos.len) {
            // Invalid vertex reference - skip this face
            goto skip_face;
          }
          allPts[flatIdx] = mat2x3_mul_vec3(projection, impl->vertPos.data[sv]);
          heMap[flatIdx] = localEdge;
          flatIdx++;
        }
      }

      // Triangulate (handles multiple loops = polygon with holes)
      ManifoldVecIVec3 tris = {0};
      if (numLoops == 1) {
        tris = manifold_triangulate_polygon(allPts, NULL, (size_t)polySizes[0]);
      } else {
        // Find the outer polygon (largest absolute signed area) and put it first
        double *areas = (double *)malloc((size_t)numLoops * sizeof(double));
        int offset = 0;
        int maxIdx = 0;
        double maxArea = 0;
        for (int li = 0; li < numLoops; li++) {
          double area = 0;
          for (int pi = 0; pi < polySizes[li]; pi++) {
            int ni = (pi + 1) % polySizes[li];
            area += allPts[offset + pi].x * allPts[offset + ni].y
                  - allPts[offset + ni].x * allPts[offset + pi].y;
          }
          areas[li] = area * 0.5;
          if (fabs(areas[li]) > maxArea) {
            maxArea = fabs(areas[li]);
            maxIdx = li;
          }
          offset += polySizes[li];
        }
        // If largest area loop is not first, swap it to position 0
        if (maxIdx != 0) {
          // Swap polySizes
          int tmpS = polySizes[0]; polySizes[0] = polySizes[maxIdx]; polySizes[maxIdx] = tmpS;
          // Swap areas
          double tmpA = areas[0]; areas[0] = areas[maxIdx]; areas[maxIdx] = tmpA;
          // Rebuild allPts/heMap with outer polygon first
          int **newLoops = (int **)malloc((size_t)numLoops * sizeof(int *));
          int *newLoopLens = (int *)malloc((size_t)numLoops * sizeof(int));
          // Put maxIdx first, then others
          newLoops[0] = loops[maxIdx]; newLoopLens[0] = loopLens[maxIdx];
          int ni = 1;
          for (int li = 0; li < numLoops; li++) {
            if (li != maxIdx) {
              newLoops[ni] = loops[li]; newLoopLens[ni] = loopLens[li];
              ni++;
            }
          }
          // Rebuild allPts and heMap
          flatIdx = 0;
          bool badVert2 = false;
          for (int li = 0; li < numLoops; li++) {
            polySizes[li] = newLoopLens[li];
            for (int pi = 0; pi < newLoopLens[li]; pi++) {
              int localEdge = newLoops[li][pi];
              int sv = impl->halfedge.data[firstEdge + localEdge].startVert;
              if (sv < 0 || (size_t)sv >= impl->vertPos.len) {
                badVert2 = true; break;
              }
              allPts[flatIdx] = mat2x3_mul_vec3(projection, impl->vertPos.data[sv]);
              heMap[flatIdx] = localEdge;
              flatIdx++;
            }
            if (badVert2) break;
          }
          // Use newLoops for cleanup later
          for (int li = 0; li < numLoops; li++) loops[li] = newLoops[li];
          free(newLoops); free(newLoopLens);
          if (badVert2) { free(areas); goto skip_face; }
        }
        // Ensure outer polygon is CCW (positive area) — flip if needed
        {
          double outerArea = 0;
          for (int pi = 0; pi < polySizes[0]; pi++) {
            int ni2 = (pi + 1) % polySizes[0];
            outerArea += allPts[pi].x * allPts[ni2].y - allPts[ni2].x * allPts[pi].y;
          }
          if (outerArea < 0) {
            // Reverse the outer loop
            for (int i = 0; i < polySizes[0] / 2; i++) {
              int j = polySizes[0] - 1 - i;
              ManifoldVec2 tmpP = allPts[i]; allPts[i] = allPts[j]; allPts[j] = tmpP;
              int tmpH = heMap[i]; heMap[i] = heMap[j]; heMap[j] = tmpH;
            }
          }
        }
        free(areas);
        tris = manifold_triangulate_with_holes(allPts, polySizes, numLoops, 0);
      }

      for (size_t ti = 0; ti < tris.len; ti++) {
        ManifoldIVec3 t = tris.data[ti];
        int he0 = firstEdge + heMap[t.x];
        int he1 = firstEdge + heMap[t.y];
        int he2 = firstEdge + heMap[t.z];
        ManifoldIVec3 triV = manifold_ivec3(
            impl->halfedge.data[he0].startVert,
            impl->halfedge.data[he1].startVert,
            impl->halfedge.data[he2].startVert);
        ManifoldIVec3 triP = manifold_ivec3(
            impl->halfedge.data[he0].propVert,
            impl->halfedge.data[he1].propVert,
            impl->halfedge.data[he2].propVert);
        vec_ivec3_push(&triVerts, triV);
        vec_ivec3_push(&triProp, triP);
        vec_vec3_push(&triNormal, normal);
        vec_triref_push(&triRef, ref);
      }

      vec_ivec3_free(&tris);
skip_face:
      for (int li = 0; li < numLoops; li++) free(loops[li]);
      free(loops); free(loopLens);
      free(allPts); free(heMap); free(polySizes);
    }
  }

  // Remove degenerate face pairs: two triangles with same 3 verts in reversed
  // winding order AND same face normal (opposite sign). These create
  // zero-volume flaps that corrupt topology.
  {
    size_t nTri = triVerts.len;
    bool *remove = (bool *)calloc(nTri, sizeof(bool));
    if (remove) {
      for (size_t i = 0; i < nTri; i++) {
        if (remove[i]) continue;
        int a0 = triVerts.data[i].x, a1 = triVerts.data[i].y, a2 = triVerts.data[i].z;
        for (size_t j = i + 1; j < nTri; j++) {
          if (remove[j]) continue;
          int b0 = triVerts.data[j].x, b1 = triVerts.data[j].y, b2 = triVerts.data[j].z;
          // Check if same 3 verts in reversed winding
          if ((a0 == b0 && a1 == b2 && a2 == b1) ||
              (a0 == b1 && a1 == b0 && a2 == b2) ||
              (a0 == b2 && a1 == b1 && a2 == b0)) {
            // Verify they come from the same original face (coplanar)
            // Same normal → degenerate pair from split face
            ManifoldVec3 ni = triNormal.data[i];
            ManifoldVec3 nj = triNormal.data[j];
            double dot = ni.x * nj.x + ni.y * nj.y + ni.z * nj.z;
            if (dot > 0.9) {
              remove[i] = remove[j] = true;
              break;
            }
          }
        }
      }
      // Compact arrays
      size_t dst = 0;
      for (size_t src = 0; src < nTri; src++) {
        if (!remove[src]) {
          triVerts.data[dst] = triVerts.data[src];
          triProp.data[dst] = triProp.data[src];
          triNormal.data[dst] = triNormal.data[src];
          triRef.data[dst] = triRef.data[src];
          dst++;
        }
      }
      triVerts.len = triProp.len = triNormal.len = triRef.len = dst;
      free(remove);
    }
  }

  // Replace faceNormal
  vec_vec3_free(&impl->faceNormal);
  impl->faceNormal = triNormal;

  // Replace meshRelation triRef
  vec_triref_free(&impl->meshRelation.triRef);
  impl->meshRelation.triRef = triRef;

  // Create halfedges from triangles
  manifold_impl_create_halfedges(impl, &triProp, &triVerts);

  vec_ivec3_free(&triVerts);
  vec_ivec3_free(&triProp);
}

// ============== Boolean3 Result ==============

static void intersections_free(ManifoldIntersections *x) {
  vec_intarr2_free(&x->p1q2);
  vec_int_free(&x->x12);
  vec_vec3_free(&x->v12);
}

// Copy manifold with triRef, mesh ID transforms
static void copy_impl_full(ManifoldImpl *dst, const ManifoldImpl *src) {
  manifold_impl_init(dst);
  dst->vertPos = vec_vec3_copy(&src->vertPos);
  dst->halfedge = vec_halfedge_copy(&src->halfedge);
  dst->faceNormal = vec_vec3_copy(&src->faceNormal);
  dst->vertNormal = vec_vec3_copy(&src->vertNormal);
  dst->halfedgeTangent = vec_vec4_copy(&src->halfedgeTangent);
  dst->bBox = src->bBox;
  dst->epsilon = src->epsilon;
  dst->tolerance = src->tolerance;
  dst->numProp = src->numProp;
  dst->properties = vec_double_copy(&src->properties);
  dst->meshRelation.originalID = src->meshRelation.originalID;
  dst->meshRelation.triRef = vec_triref_copy(&src->meshRelation.triRef);
  for (size_t i = 0; i < src->meshRelation.meshIDtransform.len; i++) {
    vec_meshid_push(&dst->meshRelation.meshIDtransform,
                    src->meshRelation.meshIDtransform.data[i]);
  }
}

static ManifoldError boolean3_result(const ManifoldBoolean3 *b3,
                                      ManifoldImpl *outR, ManifoldOpType op) {
  (void)b3->expandP;
  const ManifoldImpl *inP = b3->inP;
  const ManifoldImpl *inQ = b3->inQ;
  const ManifoldIntersections *xv12 = &b3->xv12;
  const ManifoldIntersections *xv21 = &b3->xv21;

  int c1 = (op == MANIFOLD_OP_INTERSECT) ? 0 : 1;
  int c2 = (op == MANIFOLD_OP_ADD) ? 1 : 0;
  int c3 = (op == MANIFOLD_OP_INTERSECT) ? 1 : -1;

  if (inP->status != MANIFOLD_ERROR_NO_ERROR) {
    manifold_impl_init(outR);
    outR->status = inP->status;
    return inP->status;
  }
  if (inQ->status != MANIFOLD_ERROR_NO_ERROR) {
    manifold_impl_init(outR);
    outR->status = inQ->status;
    return inQ->status;
  }

  if (manifold_impl_is_empty(inP)) {
    if (!manifold_impl_is_empty(inQ) && op == MANIFOLD_OP_ADD) {
      copy_impl_full(outR, inQ);
      return MANIFOLD_ERROR_NO_ERROR;
    }
    manifold_impl_init(outR);
    return MANIFOLD_ERROR_NO_ERROR;
  }
  if (manifold_impl_is_empty(inQ)) {
    if (op == MANIFOLD_OP_INTERSECT) {
      manifold_impl_init(outR);
      return MANIFOLD_ERROR_NO_ERROR;
    }
    copy_impl_full(outR, inP);
    return MANIFOLD_ERROR_NO_ERROR;
  }

  if (!b3->valid) {
    manifold_impl_init(outR);
    outR->status = MANIFOLD_ERROR_RESULT_TOO_LARGE;
    return MANIFOLD_ERROR_RESULT_TOO_LARGE;
  }

  bool invertQ = (op == MANIFOLD_OP_SUBTRACT);

  // Convert winding numbers to inclusion values
  ManifoldVecInt i12 = vec_int_create_n(xv12->x12.len);
  ManifoldVecInt i21 = vec_int_create_n(xv21->x12.len);
  ManifoldVecInt i03 = vec_int_create_n(b3->w03.len);
  ManifoldVecInt i30 = vec_int_create_n(b3->w30.len);

  for (size_t i = 0; i < xv12->x12.len; i++)
    i12.data[i] = c3 * xv12->x12.data[i];
  for (size_t i = 0; i < xv21->x12.len; i++)
    i21.data[i] = c3 * xv21->x12.data[i];
  for (size_t i = 0; i < b3->w03.len; i++)
    i03.data[i] = c1 + c3 * b3->w03.data[i];
  for (size_t i = 0; i < b3->w30.len; i++)
    i30.data[i] = c2 + c3 * b3->w30.data[i];


  // Calculate vertex mappings (exclusive scan with abs sum)
  ManifoldVecInt vP2R = vec_int_create_n(inP->vertPos.len);
  int numVertR = 0;
  for (size_t i = 0; i < inP->vertPos.len; i++) {
    vP2R.data[i] = numVertR;
    numVertR += abs(i03.data[i]);
  }
  int nPv = numVertR;

  ManifoldVecInt vQ2R = vec_int_create_n(inQ->vertPos.len);
  for (size_t i = 0; i < inQ->vertPos.len; i++) {
    vQ2R.data[i] = numVertR;
    numVertR += abs(i30.data[i]);
  }
  int nQv = numVertR - nPv;

  ManifoldVecInt v12R = vec_int_create_n(xv12->v12.len);
  if (xv12->v12.len > 0) {
    for (size_t i = 0; i < xv12->v12.len; i++) {
      v12R.data[i] = numVertR;
      numVertR += abs(i12.data[i]);
    }
  }

  ManifoldVecInt v21R = vec_int_create_n(xv21->v12.len);
  if (xv21->v12.len > 0) {
    for (size_t i = 0; i < xv21->v12.len; i++) {
      v21R.data[i] = numVertR;
      numVertR += abs(i21.data[i]);
    }
  }

  manifold_impl_init(outR);
  if (numVertR == 0) {
    vec_int_free(&i12); vec_int_free(&i21);
    vec_int_free(&i03); vec_int_free(&i30);
    vec_int_free(&vP2R); vec_int_free(&vQ2R);
    vec_int_free(&v12R); vec_int_free(&v21R);
    return MANIFOLD_ERROR_NO_ERROR;
  }

  outR->epsilon = fmax(inP->epsilon, inQ->epsilon);
  outR->tolerance = fmax(inP->tolerance, inQ->tolerance);

  // Create output vertices
  outR->vertPos = vec_vec3_create_n((size_t)numVertR);

  // Duplicate verts from P
  for (size_t v = 0; v < inP->vertPos.len; v++) {
    int n = abs(i03.data[v]);
    for (int i = 0; i < n; i++)
      outR->vertPos.data[vP2R.data[v] + i] = inP->vertPos.data[v];
  }
  // Duplicate verts from Q
  for (size_t v = 0; v < inQ->vertPos.len; v++) {
    int n = abs(i30.data[v]);
    for (int i = 0; i < n; i++)
      outR->vertPos.data[vQ2R.data[v] + i] = inQ->vertPos.data[v];
  }
  // New verts from P edges × Q faces
  for (size_t v = 0; v < i12.len; v++) {
    int n = abs(i12.data[v]);
    for (int i = 0; i < n; i++)
      outR->vertPos.data[v12R.data[v] + i] = xv12->v12.data[v];
  }
  // New verts from Q edges × P faces
  for (size_t v = 0; v < i21.len; v++) {
    int n = abs(i21.data[v]);
    for (int i = 0; i < n; i++)
      outR->vertPos.data[v21R.data[v] + i] = xv21->v12.data[v];
  }

  // Build edge maps for partial and new edges
  EdgeMap edgesP = {0}, edgesQ = {0};
  NewEdgeMap edgesNew = {0};

  add_new_edge_verts(&edgesP, &edgesNew, xv12, &i12, &v12R,
                      &inP->halfedge, true, 0);
  add_new_edge_verts(&edgesQ, &edgesNew, xv21, &i21, &v21R,
                      &inQ->halfedge, false, xv12->p1q2.len);

  // Size output: count sides per face
  size_t totalFaces = manifold_impl_num_tri(inP) + manifold_impl_num_tri(inQ);
  ManifoldVecInt sidesPerFace = vec_int_create_fill(totalFaces, 0);

  // Count sides from P's retained vertices
  for (size_t i = 0; i < inP->halfedge.len; i++) {
    int sv = inP->halfedge.data[i].startVert;
    sidesPerFace.data[i / 3] += abs(i03.data[sv]);
  }
  // Count sides from Q's retained vertices
  for (size_t i = 0; i < inQ->halfedge.len; i++) {
    int sv = inQ->halfedge.data[i].startVert;
    sidesPerFace.data[manifold_impl_num_tri(inP) + i / 3] += abs(i30.data[sv]);
  }
  // Count new vertices on each face from intersections
  for (size_t idx = 0; idx < i12.len; idx++) {
    int edgeP = xv12->p1q2.data[idx].v[0];
    int faceQ = xv12->p1q2.data[idx].v[1];
    int inc = abs(i12.data[idx]);
    if ((size_t)(faceQ + manifold_impl_num_tri(inP)) < sidesPerFace.len)
      sidesPerFace.data[faceQ + manifold_impl_num_tri(inP)] += inc;
    if ((size_t)edgeP < inP->halfedge.len) {
      ManifoldHalfedge he = inP->halfedge.data[edgeP];
      if ((size_t)(edgeP / 3) < sidesPerFace.len)
        sidesPerFace.data[edgeP / 3] += inc;
      if (he.pairedHalfedge >= 0 && (size_t)(he.pairedHalfedge / 3) < sidesPerFace.len)
        sidesPerFace.data[he.pairedHalfedge / 3] += inc;
    }
  }
  for (size_t idx = 0; idx < i21.len; idx++) {
    int edgeQ = xv21->p1q2.data[idx].v[1];
    int faceP = xv21->p1q2.data[idx].v[0];
    int inc = abs(i21.data[idx]);
    if ((size_t)faceP < sidesPerFace.len)
      sidesPerFace.data[faceP] += inc;
    if ((size_t)edgeQ < inQ->halfedge.len) {
      ManifoldHalfedge he = inQ->halfedge.data[edgeQ];
      if ((size_t)(manifold_impl_num_tri(inP) + edgeQ / 3) < sidesPerFace.len)
        sidesPerFace.data[manifold_impl_num_tri(inP) + edgeQ / 3] += inc;
      if (he.pairedHalfedge >= 0 && (size_t)(manifold_impl_num_tri(inP) + he.pairedHalfedge / 3) < sidesPerFace.len)
        sidesPerFace.data[manifold_impl_num_tri(inP) + he.pairedHalfedge / 3] += inc;
    }
  }

  // Build face mapping PQ → R (exclusive scan of keepFace)
  ManifoldVecInt facePQ2R = vec_int_create_n(totalFaces);
  int numFaceR = 0;
  for (size_t i = 0; i < totalFaces; i++) {
    facePQ2R.data[i] = numFaceR;
    if (sidesPerFace.data[i] > 0) numFaceR++;
  }

  // Build face normals for output
  outR->faceNormal = vec_vec3_create_n((size_t)numFaceR);
  int fIdx = 0;
  for (size_t i = 0; i < manifold_impl_num_tri(inP); i++) {
    if (sidesPerFace.data[i] > 0) {
      outR->faceNormal.data[fIdx++] = inP->faceNormal.data[i];
    }
  }
  for (size_t i = 0; i < manifold_impl_num_tri(inQ); i++) {
    if (sidesPerFace.data[manifold_impl_num_tri(inP) + i] > 0) {
      ManifoldVec3 n = inQ->faceNormal.data[i];
      if (invertQ) n = vec3_neg(n);
      outR->faceNormal.data[fIdx++] = n;
    }
  }

  // Build faceEdge (prefix sum of sides per kept face)
  ManifoldVecInt faceEdge = vec_int_create_n((size_t)(numFaceR + 1));
  faceEdge.data[0] = 0;
  int feIdx = 0;
  for (size_t i = 0; i < totalFaces; i++) {
    if (sidesPerFace.data[i] > 0) {
      faceEdge.data[feIdx + 1] = faceEdge.data[feIdx] + sidesPerFace.data[i];
      feIdx++;
    }
  }

  // Allocate halfedges
  int totalHe = faceEdge.data[numFaceR];
  outR->halfedge = vec_halfedge_create_n((size_t)totalHe);
  ManifoldVecTriRef halfedgeRef = vec_triref_create_n((size_t)totalHe);

  // Working face pointer (incremented as edges are added)
  ManifoldVecInt facePtrR = vec_int_copy(&faceEdge);

  // Reset face pointers to start of each face
  for (int i = 0; i < numFaceR; i++)
    facePtrR.data[i] = faceEdge.data[i];

  // Mark which halfedges are "whole" (not split by intersections)
  ManifoldVecChar wholeP = vec_char_create_fill(inP->halfedge.len, 1);
  ManifoldVecChar wholeQ = vec_char_create_fill(inQ->halfedge.len, 1);

  // ---- Append partial edges from P ----
  for (size_t em = 0; em < edgesP.len; em++) {
    int edgeP = edgesP.data[em].key;
    VecEdgePos *edgePosP = &edgesP.data[em].value;

    ManifoldHalfedge he = inP->halfedge.data[edgeP];
    wholeP.data[edgeP] = 0;
    if (he.pairedHalfedge >= 0 && (size_t)he.pairedHalfedge < wholeP.len)
      wholeP.data[he.pairedHalfedge] = 0;

    int vStart = he.startVert;
    int vEnd = he.endVert;
    ManifoldVec3 edgeVec = vec3_sub(inP->vertPos.data[vEnd],
                                     inP->vertPos.data[vStart]);

    // Fill in edge positions
    for (size_t j = 0; j < edgePosP->len; j++) {
      edgePosP->data[j].edgePos = vec3_dot(outR->vertPos.data[edgePosP->data[j].vert], edgeVec);
    }

    // Add original start/end verts
    int inclusion = i03.data[vStart];
    EdgePos startEP = {vec3_dot(outR->vertPos.data[vP2R.data[vStart]], edgeVec),
                        vP2R.data[vStart], INT_MAX, inclusion > 0};
    for (int j = 0; j < abs(inclusion); j++) {
      vec_edgepos_push(edgePosP, startEP);
      startEP.vert++;
    }

    inclusion = i03.data[vEnd];
    EdgePos endEP = {vec3_dot(outR->vertPos.data[vP2R.data[vEnd]], edgeVec),
                      vP2R.data[vEnd], INT_MAX, inclusion < 0};
    for (int j = 0; j < abs(inclusion); j++) {
      vec_edgepos_push(edgePosP, endEP);
      endEP.vert++;
    }

    // Pair up and add to faces
    int faceLeftP = edgeP / 3;
    if ((size_t)faceLeftP >= facePQ2R.len) continue;
    int faceLeft = facePQ2R.data[faceLeftP];
    int faceRightP = (he.pairedHalfedge >= 0) ? he.pairedHalfedge / 3 : faceLeftP;
    if ((size_t)faceRightP >= facePQ2R.len) faceRightP = faceLeftP;
    int faceRight = facePQ2R.data[faceRightP];
    ManifoldTriRef forwardRef = {0, -1, faceLeftP, -1};
    ManifoldTriRef backwardRef = {0, -1, faceRightP, -1};

    if (edgePosP->len >= 2 && edgePosP->len % 2 == 0) {
      size_t nEdges = edgePosP->len / 2;
      // Partition
      EdgePos *temp = (EdgePos *)malloc(edgePosP->len * sizeof(EdgePos));
      size_t si = 0, ei = nEdges;
      for (size_t j = 0; j < edgePosP->len; j++) {
        if (edgePosP->data[j].isStart && si < nEdges) temp[si++] = edgePosP->data[j];
        else if (ei < edgePosP->len) temp[ei++] = edgePosP->data[j];
      }
      if (si == nEdges) {
        qsort(temp, nEdges, sizeof(EdgePos), edgepos_cmp);
        qsort(temp + nEdges, nEdges, sizeof(EdgePos), edgepos_cmp);
        for (size_t j = 0; j < nEdges; j++) {
          int forwardEdge = facePtrR.data[faceLeft]++;
          int backwardEdge = facePtrR.data[faceRight]++;

          ManifoldHalfedge fwd = {temp[j].vert, temp[j + nEdges].vert,
                                   backwardEdge, temp[j].vert};
          outR->halfedge.data[forwardEdge] = fwd;
          halfedgeRef.data[forwardEdge] = forwardRef;

          ManifoldHalfedge bwd = {temp[j + nEdges].vert, temp[j].vert,
                                   forwardEdge, temp[j + nEdges].vert};
          outR->halfedge.data[backwardEdge] = bwd;
          halfedgeRef.data[backwardEdge] = backwardRef;
        }
      }
      free(temp);
    }
  }

  // ---- Append partial edges from Q ----
  for (size_t em = 0; em < edgesQ.len; em++) {
    int edgeQ = edgesQ.data[em].key;
    VecEdgePos *edgePosQ = &edgesQ.data[em].value;

    ManifoldHalfedge he = inQ->halfedge.data[edgeQ];
    wholeQ.data[edgeQ] = 0;
    if (he.pairedHalfedge >= 0 && (size_t)he.pairedHalfedge < wholeQ.len)
      wholeQ.data[he.pairedHalfedge] = 0;

    ManifoldVec3 edgeVec = vec3_sub(inQ->vertPos.data[he.endVert],
                                     inQ->vertPos.data[he.startVert]);

    for (size_t j = 0; j < edgePosQ->len; j++) {
      edgePosQ->data[j].edgePos = vec3_dot(outR->vertPos.data[edgePosQ->data[j].vert], edgeVec);
    }

    int inclusion = i30.data[he.startVert];
    EdgePos startEP = {vec3_dot(outR->vertPos.data[vQ2R.data[he.startVert]], edgeVec),
                        vQ2R.data[he.startVert], INT_MAX, inclusion > 0};
    for (int j = 0; j < abs(inclusion); j++) {
      vec_edgepos_push(edgePosQ, startEP);
      startEP.vert++;
    }

    inclusion = i30.data[he.endVert];
    EdgePos endEP = {vec3_dot(outR->vertPos.data[vQ2R.data[he.endVert]], edgeVec),
                      vQ2R.data[he.endVert], INT_MAX, inclusion < 0};
    for (int j = 0; j < abs(inclusion); j++) {
      vec_edgepos_push(edgePosQ, endEP);
      endEP.vert++;
    }

    int faceLeftQ = edgeQ / 3;
    size_t faceLeftQIdx = manifold_impl_num_tri(inP) + faceLeftQ;
    if (faceLeftQIdx >= facePQ2R.len) continue;
    int faceLeft = facePQ2R.data[faceLeftQIdx];
    int faceRightQ = (he.pairedHalfedge >= 0) ? he.pairedHalfedge / 3 : faceLeftQ;
    size_t faceRightQIdx = manifold_impl_num_tri(inP) + faceRightQ;
    if (faceRightQIdx >= facePQ2R.len) faceRightQ = faceLeftQ;
    faceRightQIdx = manifold_impl_num_tri(inP) + faceRightQ;
    int faceRight = facePQ2R.data[faceRightQIdx];
    ManifoldTriRef forwardRef = {1, -1, faceLeftQ, -1};
    ManifoldTriRef backwardRef = {1, -1, faceRightQ, -1};

    if (edgePosQ->len >= 2 && edgePosQ->len % 2 == 0) {
      size_t nEdges = edgePosQ->len / 2;
      EdgePos *temp = (EdgePos *)malloc(edgePosQ->len * sizeof(EdgePos));
      size_t si = 0, ei = nEdges;
      for (size_t j = 0; j < edgePosQ->len; j++) {
        if (edgePosQ->data[j].isStart && si < nEdges) temp[si++] = edgePosQ->data[j];
        else if (ei < edgePosQ->len) temp[ei++] = edgePosQ->data[j];
      }
      if (si == nEdges) {
        qsort(temp, nEdges, sizeof(EdgePos), edgepos_cmp);
        qsort(temp + nEdges, nEdges, sizeof(EdgePos), edgepos_cmp);
        for (size_t j = 0; j < nEdges; j++) {
          int forwardEdge = facePtrR.data[faceLeft]++;
          int backwardEdge = facePtrR.data[faceRight]++;

          ManifoldHalfedge fwd = {temp[j].vert, temp[j + nEdges].vert,
                                   backwardEdge, temp[j].vert};
          outR->halfedge.data[forwardEdge] = fwd;
          halfedgeRef.data[forwardEdge] = forwardRef;

          ManifoldHalfedge bwd = {temp[j + nEdges].vert, temp[j].vert,
                                   forwardEdge, temp[j + nEdges].vert};
          outR->halfedge.data[backwardEdge] = bwd;
          halfedgeRef.data[backwardEdge] = backwardRef;
        }
      }
      free(temp);
    }
  }

  // ---- Append new edges ----
  for (size_t em = 0; em < edgesNew.len; em++) {
    int faceP = edgesNew.data[em].faceP;
    int faceQ = edgesNew.data[em].faceQ;
    VecEdgePos *edgePos = &edgesNew.data[em].value;

    // Sort by longest dimension
    ManifoldBox bbox = manifold_box_empty();
    for (size_t j = 0; j < edgePos->len; j++) {
      ManifoldVec3 vp = outR->vertPos.data[edgePos->data[j].vert];
      bbox.min = vec3_min(bbox.min, vp);
      bbox.max = vec3_max(bbox.max, vp);
    }
    ManifoldVec3 size = vec3_sub(bbox.max, bbox.min);
    int dim = (size.x > size.y && size.x > size.z) ? 0
              : (size.y > size.z ? 1 : 2);
    for (size_t j = 0; j < edgePos->len; j++) {
      ManifoldVec3 vp = outR->vertPos.data[edgePos->data[j].vert];
      edgePos->data[j].edgePos = (dim == 0) ? vp.x : (dim == 1) ? vp.y : vp.z;
    }

    if ((size_t)faceP >= facePQ2R.len) continue;
    size_t faceQIdx = manifold_impl_num_tri(inP) + faceQ;
    if (faceQIdx >= facePQ2R.len) continue;
    int faceLeft = facePQ2R.data[faceP];
    int faceRight = facePQ2R.data[faceQIdx];
    ManifoldTriRef forwardRef = {0, -1, faceP, -1};
    ManifoldTriRef backwardRef = {1, -1, faceQ, -1};

    if (edgePos->len >= 2 && edgePos->len % 2 == 0) {
      size_t nEdges = edgePos->len / 2;
      EdgePos *temp = (EdgePos *)malloc(edgePos->len * sizeof(EdgePos));
      size_t si = 0, ei = nEdges;
      for (size_t j = 0; j < edgePos->len; j++) {
        if (edgePos->data[j].isStart && si < nEdges) temp[si++] = edgePos->data[j];
        else if (ei < edgePos->len) temp[ei++] = edgePos->data[j];
      }
      if (si == nEdges) {
        qsort(temp, nEdges, sizeof(EdgePos), edgepos_cmp);
        qsort(temp + nEdges, nEdges, sizeof(EdgePos), edgepos_cmp);
        for (size_t j = 0; j < nEdges; j++) {
          int forwardEdge = facePtrR.data[faceLeft]++;
          int backwardEdge = facePtrR.data[faceRight]++;

          ManifoldHalfedge fwd = {temp[j].vert, temp[j + nEdges].vert,
                                   backwardEdge, temp[j].vert};
          outR->halfedge.data[forwardEdge] = fwd;
          halfedgeRef.data[forwardEdge] = forwardRef;

          ManifoldHalfedge bwd = {temp[j + nEdges].vert, temp[j].vert,
                                   forwardEdge, temp[j + nEdges].vert};
          outR->halfedge.data[backwardEdge] = bwd;
          halfedgeRef.data[backwardEdge] = backwardRef;
        }
      }
      free(temp);
    }
  }

  // ---- Append whole edges from P ----
  for (size_t idx = 0; idx < inP->halfedge.len; idx++) {
    if (!wholeP.data[idx]) continue;
    ManifoldHalfedge he = inP->halfedge.data[idx];
    if (!manifold_halfedge_is_forward(&he)) continue;

    int inclusion = i03.data[he.startVert];
    if (inclusion == 0) continue;

    ManifoldHalfedge outHe = he;
    if (inclusion < 0) {
      int tmp = outHe.startVert;
      outHe.startVert = outHe.endVert;
      outHe.endVert = tmp;
    }
    outHe.startVert = vP2R.data[outHe.startVert];
    outHe.endVert = vP2R.data[outHe.endVert];

    int faceLeftP = (int)(idx / 3);
    if ((size_t)faceLeftP >= facePQ2R.len) continue;
    int newFace = facePQ2R.data[faceLeftP];
    int faceRightP = (he.pairedHalfedge >= 0) ? he.pairedHalfedge / 3 : faceLeftP;
    if ((size_t)faceRightP >= facePQ2R.len) faceRightP = faceLeftP;
    int faceRight = facePQ2R.data[faceRightP];
    ManifoldTriRef forwardRef = {0, -1, faceLeftP, -1};
    ManifoldTriRef backwardRef = {0, -1, faceRightP, -1};

    for (int i = 0; i < abs(inclusion); i++) {
      int forwardEdge = facePtrR.data[newFace]++;
      int backwardEdge = facePtrR.data[faceRight]++;

      outHe.pairedHalfedge = backwardEdge;
      outHe.propVert = outHe.startVert;
      outR->halfedge.data[forwardEdge] = outHe;
      halfedgeRef.data[forwardEdge] = forwardRef;

      ManifoldHalfedge bwd = {outHe.endVert, outHe.startVert,
                               forwardEdge, outHe.endVert};
      outR->halfedge.data[backwardEdge] = bwd;
      halfedgeRef.data[backwardEdge] = backwardRef;

      outHe.startVert++;
      outHe.endVert++;
    }
  }

  // ---- Append whole edges from Q ----
  for (size_t idx = 0; idx < inQ->halfedge.len; idx++) {
    if (!wholeQ.data[idx]) continue;
    ManifoldHalfedge he = inQ->halfedge.data[idx];
    if (!manifold_halfedge_is_forward(&he)) continue;

    int inclusion = i30.data[he.startVert];
    if (inclusion == 0) continue;

    ManifoldHalfedge outHe = he;
    if (inclusion < 0) {
      int tmp = outHe.startVert;
      outHe.startVert = outHe.endVert;
      outHe.endVert = tmp;
    }
    outHe.startVert = vQ2R.data[outHe.startVert];
    outHe.endVert = vQ2R.data[outHe.endVert];

    int faceLeftQ = (int)(idx / 3);
    size_t fLQIdx = manifold_impl_num_tri(inP) + faceLeftQ;
    if (fLQIdx >= facePQ2R.len) continue;
    int newFace = facePQ2R.data[fLQIdx];
    int faceRightQ = (he.pairedHalfedge >= 0) ? he.pairedHalfedge / 3 : faceLeftQ;
    size_t fRQIdx = manifold_impl_num_tri(inP) + faceRightQ;
    if (fRQIdx >= facePQ2R.len) { faceRightQ = faceLeftQ; fRQIdx = fLQIdx; }
    int faceRight = facePQ2R.data[fRQIdx];
    ManifoldTriRef forwardRef = {1, -1, faceLeftQ, -1};
    ManifoldTriRef backwardRef = {1, -1, faceRightQ, -1};

    for (int i = 0; i < abs(inclusion); i++) {
      int forwardEdge = facePtrR.data[newFace]++;
      int backwardEdge = facePtrR.data[faceRight]++;

      outHe.pairedHalfedge = backwardEdge;
      outHe.propVert = outHe.startVert;
      outR->halfedge.data[forwardEdge] = outHe;
      halfedgeRef.data[forwardEdge] = forwardRef;

      ManifoldHalfedge bwd = {outHe.endVert, outHe.startVert,
                               forwardEdge, outHe.endVert};
      outR->halfedge.data[backwardEdge] = bwd;
      halfedgeRef.data[backwardEdge] = backwardRef;

      outHe.startVert++;
      outHe.endVert++;
    }
  }


  // Merge coincident vertices before face2tri to avoid degenerate polygons.
  // The boolean creates many duplicate vertices at the same position (multiple
  // intersection verts, duplicated P/Q verts). Merging them first simplifies
  // the polygon faces so the ear-clipper can triangulate correctly.
  {
    int nv = (int)outR->vertPos.len;
    int *vertMap = (int *)malloc((size_t)nv * sizeof(int));
    for (int i = 0; i < nv; i++) vertMap[i] = i;
    double eps = fmax(outR->epsilon, 1e-12);
    for (int i = 0; i < nv; i++) {
      if (vertMap[i] != i) continue;
      for (int j = i + 1; j < nv; j++) {
        if (vertMap[j] != j) continue;
        ManifoldVec3 d = vec3_sub(outR->vertPos.data[i], outR->vertPos.data[j]);
        if (fabs(d.x) <= eps && fabs(d.y) <= eps && fabs(d.z) <= eps) {
          vertMap[j] = i;
        }
      }
    }
    // Update halfedge vertex references
    for (size_t i = 0; i < outR->halfedge.len; i++) {
      int sv = outR->halfedge.data[i].startVert;
      int ev = outR->halfedge.data[i].endVert;
      if (sv >= 0 && sv < nv) outR->halfedge.data[i].startVert = vertMap[sv];
      if (ev >= 0 && ev < nv) outR->halfedge.data[i].endVert = vertMap[ev];
    }
    free(vertMap);

    // Remove degenerate edges (start == end) and rebuild faceEdge
    int writeIdx = 0;
    ManifoldVecInt newFaceEdge = vec_int_create_n((size_t)(numFaceR + 1));
    ManifoldVecTriRef newHalfedgeRef = {0};
    for (int f = 0; f < numFaceR; f++) {
      newFaceEdge.data[f] = writeIdx;
      int start = faceEdge.data[f];
      int end = faceEdge.data[f + 1];
      for (int e = start; e < end; e++) {
        if (outR->halfedge.data[e].startVert != outR->halfedge.data[e].endVert) {
          if (writeIdx != e) {
            outR->halfedge.data[writeIdx] = outR->halfedge.data[e];
            // halfedgeRef is indexed the same way
          }
          vec_triref_push(&newHalfedgeRef, halfedgeRef.data[e]);
          writeIdx++;
        }
      }
    }
    newFaceEdge.data[numFaceR] = writeIdx;
    outR->halfedge.len = (size_t)writeIdx;
    vec_int_free(&faceEdge);
    faceEdge = newFaceEdge;
    vec_triref_free(&halfedgeRef);
    halfedgeRef = newHalfedgeRef;

    // Remove empty faces
    int newNumFaceR = 0;
    ManifoldVecInt cleanFaceEdge = vec_int_create_n((size_t)(numFaceR + 1));
    ManifoldVecVec3 cleanNormals = {0};
    cleanFaceEdge.data[0] = 0;
    for (int f = 0; f < numFaceR; f++) {
      int nEdges = faceEdge.data[f + 1] - faceEdge.data[f];
      if (nEdges >= 3) {
        // Compact halfedges
        int src = faceEdge.data[f];
        int dst = cleanFaceEdge.data[newNumFaceR];
        if (src != dst) {
          for (int e = 0; e < nEdges; e++) {
            outR->halfedge.data[dst + e] = outR->halfedge.data[src + e];
          }
          // Also shift halfedgeRef
          for (int e = 0; e < nEdges; e++) {
            halfedgeRef.data[dst + e] = halfedgeRef.data[src + e];
          }
        }
        vec_vec3_push(&cleanNormals, outR->faceNormal.data[f]);
        cleanFaceEdge.data[newNumFaceR + 1] = dst + nEdges;
        newNumFaceR++;
      }
    }
    outR->halfedge.len = (size_t)cleanFaceEdge.data[newNumFaceR];
    vec_int_free(&faceEdge);
    faceEdge = cleanFaceEdge;
    numFaceR = newNumFaceR;
    vec_vec3_free(&outR->faceNormal);
    outR->faceNormal = cleanNormals;
  }

  // Triangulate the faces
  faceEdge.len = (size_t)(numFaceR + 1);
  face2tri(outR, &faceEdge, &halfedgeRef, true);


  // Reorder halfedges for determinism
  manifold_impl_reorder_halfedges(outR);


  // Update references
  size_t offsetQ = manifold_mesh_id_counter;
  for (size_t i = 0; i < outR->meshRelation.triRef.len; i++) {
    ManifoldTriRef *ref = &outR->meshRelation.triRef.data[i];
    int tri = ref->faceID;
    bool isPQ = (ref->meshID == 0);
    if (isPQ && tri >= 0 && tri < (int)inP->meshRelation.triRef.len) {
      *ref = inP->meshRelation.triRef.data[tri];
    } else if (!isPQ && tri >= 0 && tri < (int)inQ->meshRelation.triRef.len) {
      *ref = inQ->meshRelation.triRef.data[tri];
      ref->meshID += (int)offsetQ;
    }
  }

  // Copy mesh ID transforms
  for (size_t i = 0; i < inP->meshRelation.meshIDtransform.len; i++) {
    vec_meshid_push(&outR->meshRelation.meshIDtransform,
                    inP->meshRelation.meshIDtransform.data[i]);
  }
  for (size_t i = 0; i < inQ->meshRelation.meshIDtransform.len; i++) {
    ManifoldMeshIDEntry entry = inQ->meshRelation.meshIDtransform.data[i];
    entry.key += (int)offsetQ;
    entry.value.backSide ^= invertQ;
    vec_meshid_push(&outR->meshRelation.meshIDtransform, entry);
  }

  // Simplify topology
  manifold_impl_simplify_topology(outR, nPv + nQv);


  manifold_impl_remove_unreferenced_verts(outR);


  // Finalize
  manifold_impl_calculate_bbox(outR);
  manifold_impl_sort_geometry(outR);
  manifold_reserve_ids(1);

  // Cleanup
  vec_int_free(&i12); vec_int_free(&i21);
  vec_int_free(&i03); vec_int_free(&i30);
  vec_int_free(&vP2R); vec_int_free(&vQ2R);
  vec_int_free(&v12R); vec_int_free(&v21R);
  vec_int_free(&sidesPerFace);
  vec_int_free(&facePQ2R);
  vec_int_free(&faceEdge);
  vec_int_free(&facePtrR);
  vec_char_free(&wholeP);
  vec_char_free(&wholeQ);
  vec_triref_free(&halfedgeRef);

  for (size_t i = 0; i < edgesP.len; i++) vec_edgepos_free(&edgesP.data[i].value);
  vec_edgemap_free(&edgesP);
  for (size_t i = 0; i < edgesQ.len; i++) vec_edgepos_free(&edgesQ.data[i].value);
  vec_edgemap_free(&edgesQ);
  for (size_t i = 0; i < edgesNew.len; i++) vec_edgepos_free(&edgesNew.data[i].value);
  vec_newedgemap_free(&edgesNew);

  return MANIFOLD_ERROR_NO_ERROR;
}

// ============== Main entry point ==============

ManifoldError manifold_boolean_op(ManifoldImpl *result,
                                   const ManifoldImpl *p,
                                   const ManifoldImpl *q,
                                   ManifoldOpType op) {
  // Handle empty cases
  if (manifold_impl_is_empty(p) && manifold_impl_is_empty(q)) {
    manifold_impl_init(result);
    return MANIFOLD_ERROR_NO_ERROR;
  }

  if (manifold_impl_is_empty(p)) {
    if (op == MANIFOLD_OP_ADD) {
      copy_impl_full(result, q);
      return MANIFOLD_ERROR_NO_ERROR;
    }
    manifold_impl_init(result);
    return MANIFOLD_ERROR_NO_ERROR;
  }

  if (manifold_impl_is_empty(q)) {
    if (op == MANIFOLD_OP_INTERSECT) {
      manifold_impl_init(result);
      return MANIFOLD_ERROR_NO_ERROR;
    }
    copy_impl_full(result, p);
    return MANIFOLD_ERROR_NO_ERROR;
  }

  // Identical geometry shortcut: if both meshes have the same vertices and
  // triangles, handle the boolean result directly without the full algorithm.
  // This avoids issues with perfectly coplanar faces.
  {
    size_t pNv = manifold_impl_num_vert(p);
    size_t qNv = manifold_impl_num_vert(q);
    size_t pNt = manifold_impl_num_tri(p);
    size_t qNt = manifold_impl_num_tri(q);
    if (pNv == qNv && pNt == qNt && pNv > 0) {
      bool same = true;
      double eps = fmax(p->epsilon, q->epsilon);
      if (eps < 1e-12) eps = 1e-12;
      for (size_t i = 0; i < pNv && same; i++) {
        ManifoldVec3 a = p->vertPos.data[i];
        ManifoldVec3 b = q->vertPos.data[i];
        double d = (a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) +
                   (a.z-b.z)*(a.z-b.z);
        if (d > eps * eps) same = false;
      }
      if (same) {
        if (op == MANIFOLD_OP_ADD || op == MANIFOLD_OP_INTERSECT) {
          copy_impl_full(result, p);
          return MANIFOLD_ERROR_NO_ERROR;
        }
        if (op == MANIFOLD_OP_SUBTRACT) {
          manifold_impl_init(result);
          return MANIFOLD_ERROR_NO_ERROR;
        }
      }
    }
  }

  // Non-overlapping shortcut
  if (!manifold_box_overlaps(p->bBox, q->bBox)) {
    if (op == MANIFOLD_OP_ADD) {
      // Simple merge for non-overlapping union
      manifold_impl_init(result);
      size_t pNv = manifold_impl_num_vert(p);
      size_t qNv = manifold_impl_num_vert(q);
      size_t pNhe = p->halfedge.len;
      size_t qNhe = q->halfedge.len;
      result->vertPos = vec_vec3_create_n(pNv + qNv);
      memcpy(result->vertPos.data, p->vertPos.data, pNv * sizeof(ManifoldVec3));
      memcpy(result->vertPos.data + pNv, q->vertPos.data, qNv * sizeof(ManifoldVec3));
      result->halfedge = vec_halfedge_create_n(pNhe + qNhe);
      memcpy(result->halfedge.data, p->halfedge.data, pNhe * sizeof(ManifoldHalfedge));
      for (size_t i = 0; i < qNhe; i++) {
        ManifoldHalfedge he = q->halfedge.data[i];
        he.startVert += (int)pNv;
        he.endVert += (int)pNv;
        he.pairedHalfedge += (int)pNhe;
        he.propVert = he.startVert;
        result->halfedge.data[pNhe + i] = he;
      }
      size_t pNt = manifold_impl_num_tri(p);
      size_t qNt = manifold_impl_num_tri(q);
      result->meshRelation.triRef = vec_triref_create_n(pNt + qNt);
      for (size_t i = 0; i < pNt && i < p->meshRelation.triRef.len; i++)
        result->meshRelation.triRef.data[i] = p->meshRelation.triRef.data[i];
      for (size_t i = 0; i < qNt && i < q->meshRelation.triRef.len; i++)
        result->meshRelation.triRef.data[pNt + i] = q->meshRelation.triRef.data[i];
      for (size_t i = 0; i < p->meshRelation.meshIDtransform.len; i++)
        vec_meshid_push(&result->meshRelation.meshIDtransform,
                        p->meshRelation.meshIDtransform.data[i]);
      for (size_t i = 0; i < q->meshRelation.meshIDtransform.len; i++)
        vec_meshid_push(&result->meshRelation.meshIDtransform,
                        q->meshRelation.meshIDtransform.data[i]);
      manifold_impl_calculate_bbox(result);
      manifold_impl_set_epsilon(result, fmax(p->epsilon, q->epsilon), false);
      result->tolerance = fmax(p->tolerance, q->tolerance);
      manifold_impl_set_normals_and_coplanar(result);
      return MANIFOLD_ERROR_NO_ERROR;
    }
    if (op == MANIFOLD_OP_SUBTRACT) {
      copy_impl_full(result, p);
      return MANIFOLD_ERROR_NO_ERROR;
    }
    // Intersect with no overlap → empty
    manifold_impl_init(result);
    return MANIFOLD_ERROR_NO_ERROR;
  }

  // Full boolean algorithm for overlapping meshes
  // Ensure colliders are built - make local copies if needed to avoid
  // modifying const inputs (sort_geometry reorders vertices/faces in place)
  ManifoldImpl pLocal, qLocal;
  bool pCopied = false, qCopied = false;
  if (!p->colliderBuilt) {
    copy_impl_full(&pLocal, p);
    manifold_impl_sort_geometry(&pLocal);
    p = &pLocal;
    pCopied = true;
  }
  if (!q->colliderBuilt) {
    copy_impl_full(&qLocal, q);
    manifold_impl_sort_geometry(&qLocal);
    q = &qLocal;
    qCopied = true;
  }

  // Build Boolean3 state
  ManifoldBoolean3 b3;
  b3.inP = p;
  b3.inQ = q;
  b3.expandP = (op == MANIFOLD_OP_ADD);
  b3.valid = true;

  // If no bounding box overlap, early out with zero winding
  if (!manifold_box_overlaps(p->bBox, q->bBox)) {
    b3.w03 = vec_int_create_fill(p->vertPos.len, 0);
    b3.w30 = vec_int_create_fill(q->vertPos.len, 0);
    b3.xv12 = (ManifoldIntersections){{0}, {0}, {0}};
    b3.xv21 = (ManifoldIntersections){{0}, {0}, {0}};
  } else {
    // Compute intersections
    b3.xv12 = intersect12(p, q, b3.expandP, true);
    b3.xv21 = intersect12(p, q, b3.expandP, false);

    if (b3.xv12.x12.len > (size_t)INT_MAX || b3.xv21.x12.len > (size_t)INT_MAX) {
      b3.valid = false;
    }

    // Compute winding numbers
    b3.w03 = winding03(p, q, &b3.xv12, b3.expandP, true);
    b3.w30 = winding03(p, q, &b3.xv21, b3.expandP, false);
  }

  ManifoldError err = boolean3_result(&b3, result, op);

  // Cleanup Boolean3
  intersections_free(&b3.xv12);
  intersections_free(&b3.xv21);
  vec_int_free(&b3.w03);
  vec_int_free(&b3.w30);

  // Free local copies if we made them
  if (pCopied) manifold_impl_free(&pLocal);
  if (qCopied) manifold_impl_free(&qLocal);

  return err;
}
