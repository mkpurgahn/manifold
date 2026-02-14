// Copyright 2021 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// Ear-clipping polygon triangulation - C11 port of src/polygon.cpp.
// Faithfully ports the C++ EarClip class with cost-based ear priority,
// keyhole cutting for holes, and convex polygon fast path.

#ifndef MANIFOLD_EARCLIP_H
#define MANIFOLD_EARCLIP_H

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <stdbool.h>

#include "manifold_vec_math.h"
#include "manifold_vec.h"

// ==================== Data Structures ====================

typedef struct {
  double min_x, min_y, max_x, max_y;
} ECRect;

static inline ECRect ec_rect_empty(void) {
  return (ECRect){INFINITY, INFINITY, -INFINITY, -INFINITY};
}

static inline void ec_rect_union_pt(ECRect *r, ManifoldVec2 p) {
  if (p.x < r->min_x) r->min_x = p.x;
  if (p.y < r->min_y) r->min_y = p.y;
  if (p.x > r->max_x) r->max_x = p.x;
  if (p.y > r->max_y) r->max_y = p.y;
}

static inline double ec_rect_scale(ECRect r) {
  double dx = r.max_x - r.min_x;
  double dy = r.max_y - r.min_y;
  return fmax(dx, dy);
}

#define EC_NIL (-1)
static const double EC_KBEST = -INFINITY;

typedef struct {
  int mesh_idx;
  double cost;
  ManifoldVec2 pos;
  ManifoldVec2 rightDir;
  int left;   // index in verts array
  int right;  // index in verts array
  int ear_heap_idx;  // index in ear heap, or EC_NIL
} ECVert;

// Simple min-heap for ear priority queue (sorted by cost)
typedef struct {
  int *data;    // array of vert indices
  int len;
  int cap;
  ECVert *verts;  // pointer to verts array (for cost lookup)
} ECEarHeap;

static inline void ec_heap_init(ECEarHeap *h, ECVert *verts, int cap) {
  h->data = (int *)malloc((size_t)cap * sizeof(int));
  h->len = 0;
  h->cap = cap;
  h->verts = verts;
}

static inline void ec_heap_free(ECEarHeap *h) {
  free(h->data);
  h->data = NULL;
  h->len = 0;
}

static inline void ec_heap_swap(ECEarHeap *h, int i, int j) {
  int a = h->data[i], b = h->data[j];
  h->data[i] = b;
  h->data[j] = a;
  h->verts[a].ear_heap_idx = j;
  h->verts[b].ear_heap_idx = i;
}

static inline void ec_heap_sift_up(ECEarHeap *h, int i) {
  while (i > 0) {
    int parent = (i - 1) / 2;
    if (h->verts[h->data[i]].cost < h->verts[h->data[parent]].cost) {
      ec_heap_swap(h, i, parent);
      i = parent;
    } else break;
  }
}

static inline void ec_heap_sift_down(ECEarHeap *h, int i) {
  while (1) {
    int smallest = i;
    int l = 2 * i + 1, r = 2 * i + 2;
    if (l < h->len && h->verts[h->data[l]].cost < h->verts[h->data[smallest]].cost)
      smallest = l;
    if (r < h->len && h->verts[h->data[r]].cost < h->verts[h->data[smallest]].cost)
      smallest = r;
    if (smallest != i) {
      ec_heap_swap(h, i, smallest);
      i = smallest;
    } else break;
  }
}

static inline void ec_heap_insert(ECEarHeap *h, int vi) {
  if (h->len >= h->cap) {
    h->cap = h->cap * 2 + 1;
    h->data = (int *)realloc(h->data, (size_t)h->cap * sizeof(int));
  }
  int i = h->len++;
  h->data[i] = vi;
  h->verts[vi].ear_heap_idx = i;
  ec_heap_sift_up(h, i);
}

static inline int ec_heap_pop(ECEarHeap *h) {
  if (h->len == 0) return EC_NIL;
  int top = h->data[0];
  h->verts[top].ear_heap_idx = EC_NIL;
  h->len--;
  if (h->len > 0) {
    h->data[0] = h->data[h->len];
    h->verts[h->data[0]].ear_heap_idx = 0;
    ec_heap_sift_down(h, 0);
  }
  return top;
}

static inline void ec_heap_remove(ECEarHeap *h, int vi) {
  int i = h->verts[vi].ear_heap_idx;
  if (i == EC_NIL || i >= h->len) return;
  h->verts[vi].ear_heap_idx = EC_NIL;
  h->len--;
  if (i == h->len) return;
  h->data[i] = h->data[h->len];
  h->verts[h->data[i]].ear_heap_idx = i;
  ec_heap_sift_up(h, i);
  ec_heap_sift_down(h, i);
}

static inline void ec_heap_update(ECEarHeap *h, int vi) {
  int i = h->verts[vi].ear_heap_idx;
  if (i == EC_NIL) return;
  ec_heap_sift_up(h, i);
  ec_heap_sift_down(h, i);
}

// ==================== Vert Operations ====================

static inline ManifoldVec2 ec_safe_normalize(ManifoldVec2 v) {
  ManifoldVec2 n = vec2_normalize(v);
  return isfinite(n.x) ? n : manifold_vec2(0, 0);
}

static inline bool ec_clipped(const ECVert *verts, int v) {
  return verts[verts[v].right].left != v;
}

static inline void ec_link(ECVert *verts, int left, int right) {
  verts[left].right = right;
  verts[right].left = left;
  verts[left].rightDir = ec_safe_normalize(
      vec2_sub(verts[right].pos, verts[left].pos));
}

static inline bool ec_is_short(const ECVert *verts, int v, double epsilon) {
  ManifoldVec2 edge = vec2_sub(verts[verts[v].right].pos, verts[v].pos);
  return vec2_dot(edge, edge) * 4 < epsilon * epsilon;
}

// Like CCW, returns 1 if pt is on the inside of the angle at v
static inline int ec_interior(const ECVert *verts, int v, ManifoldVec2 pt, double epsilon) {
  ManifoldVec2 diff = vec2_sub(pt, verts[v].pos);
  if (vec2_dot(diff, diff) < epsilon * epsilon) return 0;
  return manifold_ccw(verts[v].pos, verts[verts[v].left].pos, verts[verts[v].right].pos, epsilon) +
         manifold_ccw(verts[v].pos, verts[verts[v].right].pos, pt, epsilon) +
         manifold_ccw(verts[v].pos, pt, verts[verts[v].left].pos, epsilon);
}

// InsideEdge: returns true if vert vi is on the inside of the edge from tail to tail->right
static inline bool ec_inside_edge(const ECVert *verts, int vi, int tail, double epsilon, bool toLeft) {
  double p2 = epsilon * epsilon;
  int nextL = verts[verts[vi].left].right;
  int nextR = verts[tail].right;
  int center = tail;
  int last = center;

  int maxIter = 10000; // safety
  while (nextL != nextR && tail != nextR &&
         nextL != (toLeft ? verts[vi].right : verts[vi].left) && --maxIter > 0) {
    ManifoldVec2 edgeL = vec2_sub(verts[nextL].pos, verts[center].pos);
    double l2 = vec2_dot(edgeL, edgeL);
    if (l2 <= p2) {
      nextL = toLeft ? verts[nextL].left : verts[nextL].right;
      continue;
    }

    ManifoldVec2 edgeR = vec2_sub(verts[nextR].pos, verts[center].pos);
    double r2 = vec2_dot(edgeR, edgeR);
    if (r2 <= p2) {
      nextR = verts[nextR].right;
      continue;
    }

    ManifoldVec2 vecLR = vec2_sub(verts[nextR].pos, verts[nextL].pos);
    double lr2 = vec2_dot(vecLR, vecLR);
    if (lr2 <= p2) {
      last = center;
      center = nextL;
      nextL = toLeft ? verts[nextL].left : verts[nextL].right;
      if (nextL == nextR) break;
      nextR = verts[nextR].right;
      continue;
    }

    int convexity = manifold_ccw(verts[nextL].pos, verts[center].pos, verts[nextR].pos, epsilon);
    if (center != last) {
      convexity += manifold_ccw(verts[last].pos, verts[center].pos, verts[nextL].pos, epsilon) +
                   manifold_ccw(verts[nextR].pos, verts[center].pos, verts[last].pos, epsilon);
    }
    if (convexity != 0) return convexity > 0;

    if (l2 < r2) {
      center = nextL;
      nextL = toLeft ? verts[nextL].left : verts[nextL].right;
    } else {
      center = nextR;
      nextR = verts[nextR].right;
    }
    last = center;
  }
  return true;
}

static inline bool ec_is_convex(const ECVert *verts, int v, double epsilon) {
  return manifold_ccw(verts[verts[v].left].pos, verts[v].pos,
                      verts[verts[v].right].pos, epsilon) >= 0;
}

static inline bool ec_is_reflex(const ECVert *verts, int v, double epsilon) {
  int tail = verts[v].left;
  return !ec_inside_edge(verts, tail, verts[tail].right, epsilon, true);
}

// InterpY2X: returns x-value on edge at v corresponding to start.y
static inline double ec_interp_y2x(const ECVert *verts, int v, ManifoldVec2 start,
                                    int onTop, double epsilon) {
  int r = verts[v].right;
  if (fabs(verts[v].pos.y - start.y) <= epsilon) {
    if (verts[r].pos.y <= start.y + epsilon || onTop == 1) return NAN;
    else return verts[v].pos.x;
  } else if (verts[v].pos.y < start.y - epsilon) {
    if (verts[r].pos.y > start.y + epsilon) {
      return verts[v].pos.x + (start.y - verts[v].pos.y) *
             (verts[r].pos.x - verts[v].pos.x) /
             (verts[r].pos.y - verts[v].pos.y);
    } else if (verts[r].pos.y < start.y - epsilon || onTop == -1) return NAN;
    else return verts[r].pos.x;
  }
  return NAN;
}

// SignedDist for ear cost
static inline double ec_signed_dist(const ECVert *verts, int from, int test,
                                     ManifoldVec2 unit, double epsilon) {
  double d = vec2_cross(unit, vec2_sub(verts[test].pos, verts[from].pos));
  if (fabs(d) < epsilon) {
    double dR = vec2_cross(unit, vec2_sub(verts[verts[test].right].pos, verts[from].pos));
    if (fabs(dR) > epsilon) return dR;
    double dL = vec2_cross(unit, vec2_sub(verts[verts[test].left].pos, verts[from].pos));
    if (fabs(dL) > epsilon) return dL;
  }
  return d;
}

// Cost of test vert within ear at v
static inline double ec_cost(const ECVert *verts, int v, int test,
                              ManifoldVec2 openSide, double epsilon) {
  double cost = fmin(ec_signed_dist(verts, v, test, verts[v].rightDir, epsilon),
                     ec_signed_dist(verts, verts[v].left, test,
                                    verts[verts[v].left].rightDir, epsilon));
  double openCost = vec2_cross(openSide,
      vec2_sub(verts[test].pos, verts[verts[v].right].pos));
  return fmin(cost, openCost);
}

static inline double ec_delaunay_cost(ManifoldVec2 diff, double scale, double epsilon) {
  return -epsilon - scale * vec2_dot(diff, diff);
}

// EarCost: compute cost for ear at v by checking all verts in polygon
// Uses brute-force O(n) check instead of BVH for simplicity
static inline double ec_ear_cost(ECVert *verts, int v, int numVerts, double epsilon) {
  ManifoldVec2 openSide = vec2_sub(verts[verts[v].left].pos, verts[verts[v].right].pos);
  ManifoldVec2 center = vec2_scale(
      vec2_add(verts[verts[v].left].pos, verts[verts[v].right].pos), 0.5);
  double scale = 4.0 / vec2_dot(openSide, openSide);
  openSide = ec_safe_normalize(openSide);

  double totalCost = vec2_dot(verts[verts[v].left].rightDir, verts[v].rightDir) - 1 - epsilon;
  if (manifold_ccw(verts[v].pos, verts[verts[v].left].pos,
                   verts[verts[v].right].pos, epsilon) == 0) {
    return totalCost;
  }

  int lid = verts[verts[v].left].mesh_idx;
  int rid = verts[verts[v].right].mesh_idx;

  for (int i = 0; i < numVerts; i++) {
    if (ec_clipped(verts, i)) continue;
    if (verts[i].mesh_idx == verts[v].mesh_idx ||
        verts[i].mesh_idx == lid || verts[i].mesh_idx == rid) continue;
    double cost = ec_cost(verts, v, i, openSide, epsilon);
    if (cost < -epsilon) {
      cost = ec_delaunay_cost(vec2_sub(verts[i].pos, center), scale, epsilon);
    }
    if (cost > totalCost) totalCost = cost;
  }
  return totalCost;
}

// ==================== Main EarClip Algorithm ====================

typedef struct {
  ECVert *verts;
  int numVerts;
  int cap;
  ECEarHeap ears;
  ManifoldVecIVec3 triangles;
  double epsilon;
  ECRect bBox;
  // Holes (sorted by maxX, descending)
  int *holes;   // start vert indices
  int numHoles;
  ECRect *holeBBox;
  // Outers (positive area)
  int *outers;
  int numOuters;
  // Simples (to triangulate)
  int *simples;
  int numSimples;
} EarClipState;

static void ec_clip_ear(EarClipState *st, int ear) {
  ec_link(st->verts, st->verts[ear].left, st->verts[ear].right);
  int li = st->verts[ear].left;
  int ri = st->verts[ear].right;
  if (st->verts[li].mesh_idx != st->verts[ear].mesh_idx &&
      st->verts[ear].mesh_idx != st->verts[ri].mesh_idx &&
      st->verts[ri].mesh_idx != st->verts[li].mesh_idx) {
    vec_ivec3_push(&st->triangles,
        manifold_ivec3(st->verts[li].mesh_idx,
                       st->verts[ear].mesh_idx,
                       st->verts[ri].mesh_idx));
  }
}

static void ec_clip_if_degenerate(EarClipState *st, int ear) {
  if (ec_clipped(st->verts, ear)) return;
  if (st->verts[ear].left == st->verts[ear].right) return;
  if (ec_is_short(st->verts, ear, st->epsilon) ||
      (manifold_ccw(st->verts[st->verts[ear].left].pos, st->verts[ear].pos,
                    st->verts[st->verts[ear].right].pos, st->epsilon) == 0 &&
       vec2_dot(vec2_sub(st->verts[st->verts[ear].left].pos, st->verts[ear].pos),
                vec2_sub(st->verts[st->verts[ear].right].pos, st->verts[ear].pos)) > 0)) {
    ec_clip_ear(st, ear);
    ec_clip_if_degenerate(st, st->verts[ear].left);
    ec_clip_if_degenerate(st, st->verts[ear].right);
  }
}

// Process ear: compute cost and update heap
static void ec_process_ear(EarClipState *st, int v) {
  if (st->verts[v].ear_heap_idx != EC_NIL) {
    ec_heap_remove(&st->ears, v);
  }
  if (ec_is_short(st->verts, v, st->epsilon)) {
    st->verts[v].cost = EC_KBEST;
    ec_heap_insert(&st->ears, v);
  } else if (ec_is_convex(st->verts, v, 2 * st->epsilon)) {
    st->verts[v].cost = ec_ear_cost(st->verts, v, st->numVerts, st->epsilon);
    ec_heap_insert(&st->ears, v);
  } else {
    st->verts[v].cost = 1; // reflex
  }
}

// Loop over unclipped verts starting from first, calling func on each
// Returns an unclipped vert, or EC_NIL if polygon degenerates
typedef void (*ec_loop_fn)(EarClipState *, int);

static int ec_loop(EarClipState *st, int first, ec_loop_fn func) {
  int v = first;
  int maxIter = st->numVerts + 10;
  do {
    if (ec_clipped(st->verts, v)) {
      first = st->verts[st->verts[v].right].left;
      if (!ec_clipped(st->verts, first)) {
        v = first;
        if (st->verts[v].right == st->verts[v].left) return EC_NIL;
        func(st, v);
      }
    } else {
      if (st->verts[v].right == st->verts[v].left) return EC_NIL;
      func(st, v);
    }
    v = st->verts[v].right;
  } while (v != first && --maxIter > 0);
  return maxIter > 0 ? v : EC_NIL;
}

static int ec_loop_nofunc(EarClipState *st, int first, void (*func)(EarClipState*, int, void*), void *ctx) {
  int v = first;
  int maxIter = st->numVerts + 10;
  do {
    if (ec_clipped(st->verts, v)) {
      first = st->verts[st->verts[v].right].left;
      if (!ec_clipped(st->verts, first)) {
        v = first;
        if (st->verts[v].right == st->verts[v].left) return EC_NIL;
        func(st, v, ctx);
      }
    } else {
      if (st->verts[v].right == st->verts[v].left) return EC_NIL;
      func(st, v, ctx);
    }
    v = st->verts[v].right;
  } while (v != first && --maxIter > 0);
  return maxIter > 0 ? v : EC_NIL;
}

// FindStart: categorize polygon as hole, outer, or simple
typedef struct {
  double maxX;
  int start;
  ECRect bBox;
  double area;
  double areaComp;
  ManifoldVec2 origin;
} FindStartCtx;

static void ec_find_start_fn(EarClipState *st, int v, void *ctx) {
  FindStartCtx *c = (FindStartCtx *)ctx;
  ec_rect_union_pt(&c->bBox, st->verts[v].pos);
  double area1 = vec2_cross(vec2_sub(st->verts[v].pos, c->origin),
                            vec2_sub(st->verts[st->verts[v].right].pos, c->origin));
  double t1 = c->area + area1;
  c->areaComp += (c->area - t1) + area1;
  c->area = t1;
  if (st->verts[v].pos.x > c->maxX) {
    c->maxX = st->verts[v].pos.x;
    c->start = v;
  }
}

static void ec_find_start(EarClipState *st, int first) {
  FindStartCtx ctx;
  ctx.maxX = -INFINITY;
  ctx.start = first;
  ctx.bBox = ec_rect_empty();
  ctx.area = 0;
  ctx.areaComp = 0;
  ctx.origin = st->verts[first].pos;

  if (ec_loop_nofunc(st, first, ec_find_start_fn, &ctx) == EC_NIL) return;

  ctx.area += ctx.areaComp;
  double sx = ctx.bBox.max_x - ctx.bBox.min_x;
  double sy = ctx.bBox.max_y - ctx.bBox.min_y;
  double minArea = st->epsilon * fmax(sx, sy);

  if (isfinite(ctx.maxX) && ctx.area < -minArea) {
    // Hole - insert sorted by maxX descending
    int insertIdx = st->numHoles;
    for (int i = 0; i < st->numHoles; i++) {
      if (st->verts[ctx.start].pos.x > st->verts[st->holes[i]].pos.x) {
        insertIdx = i;
        break;
      }
    }
    // Shift
    for (int i = st->numHoles; i > insertIdx; i--) {
      st->holes[i] = st->holes[i-1];
      st->holeBBox[i] = st->holeBBox[i-1];
    }
    st->holes[insertIdx] = ctx.start;
    st->holeBBox[insertIdx] = ctx.bBox;
    st->numHoles++;
  } else {
    st->simples[st->numSimples++] = ctx.start;
    if (ctx.area > minArea) {
      st->outers[st->numOuters++] = ctx.start;
    }
  }
}

// FindCloserBridge
static int ec_find_closer_bridge(EarClipState *st, int start, int edge) {
  int connector;
  int eright = st->verts[edge].right;
  if (st->verts[edge].pos.x < st->verts[start].pos.x)
    connector = eright;
  else if (st->verts[eright].pos.x < st->verts[start].pos.x)
    connector = edge;
  else if (st->verts[eright].pos.y - st->verts[start].pos.y >
           st->verts[start].pos.y - st->verts[edge].pos.y)
    connector = edge;
  else
    connector = eright;

  if (fabs(st->verts[connector].pos.y - st->verts[start].pos.y) <= st->epsilon)
    return connector;

  double above = st->verts[connector].pos.y > st->verts[start].pos.y ? 1 : -1;

  for (int oi = 0; oi < st->numOuters; oi++) {
    int first = st->outers[oi];
    int v = first;
    int maxIter = st->numVerts + 10;
    do {
      if (!ec_clipped(st->verts, v) && st->verts[v].right != st->verts[v].left) {
        double inside = above * manifold_ccw(st->verts[start].pos, st->verts[v].pos,
                                              st->verts[connector].pos, st->epsilon);
        if (st->verts[v].pos.x > st->verts[start].pos.x - st->epsilon &&
            st->verts[v].pos.y * above > st->verts[start].pos.y * above - st->epsilon &&
            (inside > 0 || (inside == 0 && st->verts[v].pos.x < st->verts[connector].pos.x &&
                            st->verts[v].pos.y * above < st->verts[connector].pos.y * above)) &&
            ec_inside_edge(st->verts, v, edge, st->epsilon, true) &&
            ec_is_reflex(st->verts, v, st->epsilon)) {
          connector = v;
        }
      }
      v = st->verts[v].right;
    } while (v != first && --maxIter > 0);
  }
  return connector;
}

// CutKeyhole: connect a hole to an outer polygon
static void ec_cut_keyhole(EarClipState *st, int startIdx) {
  int start = st->holes[startIdx];
  ECRect bBox = st->holeBBox[startIdx];
  int onTop = st->verts[start].pos.y >= bBox.max_y - st->epsilon ? 1 :
              st->verts[start].pos.y <= bBox.min_y + st->epsilon ? -1 : 0;
  int connector = EC_NIL;

  for (int oi = 0; oi < st->numOuters; oi++) {
    int first = st->outers[oi];
    int v = first;
    int maxIter = st->numVerts + 10;
    do {
      if (!ec_clipped(st->verts, v) && st->verts[v].right != st->verts[v].left) {
        double x = ec_interp_y2x(st->verts, v, st->verts[start].pos, onTop, st->epsilon);
        if (isfinite(x) && ec_inside_edge(st->verts, start, v, st->epsilon, true) &&
            (connector == EC_NIL ||
             manifold_ccw(manifold_vec2(x, st->verts[start].pos.y),
                          st->verts[connector].pos,
                          st->verts[st->verts[connector].right].pos, st->epsilon) == 1 ||
             (st->verts[connector].pos.y < st->verts[v].pos.y
                  ? ec_inside_edge(st->verts, v, connector, st->epsilon, false)
                  : !ec_inside_edge(st->verts, connector, v, st->epsilon, false)))) {
          connector = v;
        }
      }
      v = st->verts[v].right;
    } while (v != first && --maxIter > 0);
  }

  if (connector == EC_NIL) {
    st->simples[st->numSimples++] = start;
    return;
  }

  connector = ec_find_closer_bridge(st, start, connector);

  // JoinPolygons: duplicate start and connector
  // Ensure capacity
  if (st->numVerts + 2 > st->cap) {
    st->cap = st->cap * 2 + 2;
    st->verts = (ECVert *)realloc(st->verts, (size_t)st->cap * sizeof(ECVert));
    st->ears.verts = st->verts;
  }

  int newStart = st->numVerts++;
  int newConnector = st->numVerts++;
  st->verts[newStart] = st->verts[start];
  st->verts[newConnector] = st->verts[connector];
  st->verts[newStart].ear_heap_idx = EC_NIL;
  st->verts[newConnector].ear_heap_idx = EC_NIL;

  st->verts[st->verts[start].right].left = newStart;
  st->verts[st->verts[connector].left].right = newConnector;
  ec_link(st->verts, start, connector);
  ec_link(st->verts, newConnector, newStart);

  ec_clip_if_degenerate(st, start);
  ec_clip_if_degenerate(st, newStart);
  ec_clip_if_degenerate(st, connector);
  ec_clip_if_degenerate(st, newConnector);
}

// QueueVert callback
static void ec_queue_vert(EarClipState *st, int v) {
  ec_process_ear(st, v);
}

// TriangulatePoly: main ear-clipping for a simple polygon
static void ec_triangulate_poly(EarClipState *st, int start) {
  st->ears.len = 0; // clear heap

  int numTri = -2;
  // Count unclipped verts and queue them
  int v = start;
  int maxIter = st->numVerts + 10;
  int first2 = start;
  do {
    if (ec_clipped(st->verts, v)) {
      first2 = st->verts[st->verts[v].right].left;
      if (!ec_clipped(st->verts, first2)) {
        v = first2;
        if (st->verts[v].right == st->verts[v].left) return;
        ec_process_ear(st, v);
        numTri++;
      }
    } else {
      if (st->verts[v].right == st->verts[v].left) return;
      ec_process_ear(st, v);
      numTri++;
    }
    v = st->verts[v].right;
  } while (v != first2 && --maxIter > 0);
  if (maxIter <= 0) return;

  int backup = v;
  while (numTri > 0) {
    int ear = ec_heap_pop(&st->ears);
    if (ear != EC_NIL) {
      backup = st->verts[ear].right;
    } else {
      ear = backup;
    }
    ec_clip_ear(st, ear);
    numTri--;
    ec_process_ear(st, st->verts[ear].left);
    ec_process_ear(st, st->verts[ear].right);
    backup = st->verts[ear].right;
  }
}

// ==================== IsConvex ====================

static inline bool ec_poly_is_convex(const ManifoldVec2 *pts, const int *indices,
                                      int n, double epsilon) {
  ManifoldVec2 firstEdge = vec2_sub(pts[0], pts[n - 1]);
  ManifoldVec2 lastEdge = ec_safe_normalize(firstEdge);
  for (int v = 0; v < n; v++) {
    ManifoldVec2 edge = v + 1 < n ? vec2_sub(pts[v + 1], pts[v]) : firstEdge;
    double det = vec2_cross(lastEdge, edge);
    if (det <= 0 || (fabs(det) < epsilon && vec2_dot(lastEdge, edge) < 0))
      return false;
    lastEdge = ec_safe_normalize(edge);
  }
  return true;
}

// ==================== TriangulateConvex (zigzag) ====================

static inline ManifoldVecIVec3 ec_triangulate_convex(const int *indices, int n) {
  ManifoldVecIVec3 tris = {0};
  int i = 0, k = n - 1;
  bool right = true;
  while (i + 1 < k) {
    int j = right ? i + 1 : k - 1;
    vec_ivec3_push(&tris, manifold_ivec3(indices[i], indices[j], indices[k]));
    if (right) i = j;
    else k = j;
    right = !right;
  }
  return tris;
}

// ==================== Main Entry Point ====================

// Triangulate a set of polygon loops (first is outer, rest are holes).
// polyVerts: flat array of all vertices (outer vertices first, then holes)
// polySizes: array of vertex counts [outerN, hole1N, hole2N, ...]
// nPolys: number of polygons (1 outer + n-1 holes)
// indices: mapping from flat vertex index to halfedge index (or NULL for 0..n-1)
// epsilon: tolerance
// allowConvex: if true, use fast path for convex polygons
//
// Returns triangle indices referencing the flat vertex positions in polyVerts.
static inline ManifoldVecIVec3 ec_triangulate(
    const ManifoldVec2 *polyVerts, const int *polySizes, int nPolys,
    const int *indices, double epsilon, bool allowConvex) {
  // Total verts
  int totalVerts = 0;
  for (int i = 0; i < nPolys; i++) totalVerts += polySizes[i];
  if (totalVerts < 3) return (ManifoldVecIVec3){0};

  // Build index array
  int *idx = (int *)malloc((size_t)totalVerts * sizeof(int));
  for (int i = 0; i < totalVerts; i++) {
    idx[i] = indices ? indices[i] : i;
  }

  // Check convex fast path (single polygon only)
  if (allowConvex && nPolys == 1) {
    if (ec_poly_is_convex(polyVerts, idx, totalVerts, epsilon)) {
      ManifoldVecIVec3 result = ec_triangulate_convex(idx, totalVerts);
      free(idx);
      return result;
    }
  }

  // Full EarClip
  int cap = totalVerts + 2 * nPolys + 4;
  EarClipState st;
  st.cap = cap;
  st.verts = (ECVert *)malloc((size_t)cap * sizeof(ECVert));
  st.numVerts = totalVerts;
  st.triangles = (ManifoldVecIVec3){0};
  st.bBox = ec_rect_empty();

  // Initialize circular linked list
  int offset = 0;
  int *starts = (int *)malloc((size_t)nPolys * sizeof(int));
  for (int p = 0; p < nPolys; p++) {
    int polyN = polySizes[p];
    int first = offset;
    for (int i = 0; i < polyN; i++) {
      int vi = offset + i;
      st.verts[vi].mesh_idx = idx[vi];
      st.verts[vi].cost = 0;
      st.verts[vi].ear_heap_idx = EC_NIL;
      st.verts[vi].pos = polyVerts[vi];
      st.verts[vi].rightDir = manifold_vec2(0, 0);
      st.verts[vi].left = EC_NIL;
      st.verts[vi].right = EC_NIL;
      ec_rect_union_pt(&st.bBox, polyVerts[vi]);
    }
    // Link into circular list
    for (int i = 0; i < polyN; i++) {
      int vi = offset + i;
      int next = offset + (i + 1) % polyN;
      ec_link(st.verts, vi, next);
    }
    starts[p] = first;
    offset += polyN;
  }

  if (epsilon < 0) epsilon = ec_rect_scale(st.bBox) * 1.1920928955078125e-7; // kPrecision

  st.epsilon = epsilon;
  ec_heap_init(&st.ears, st.verts, totalVerts);

  // Clip degenerate ears
  for (int i = 0; i < totalVerts; i++) {
    ec_clip_if_degenerate(&st, i);
  }

  // Allocate space for holes, outers, simples
  st.holes = (int *)malloc((size_t)nPolys * sizeof(int));
  st.holeBBox = (ECRect *)malloc((size_t)nPolys * sizeof(ECRect));
  st.numHoles = 0;
  st.outers = (int *)malloc((size_t)nPolys * sizeof(int));
  st.numOuters = 0;
  st.simples = (int *)malloc((size_t)nPolys * sizeof(int));
  st.numSimples = 0;

  // FindStart for each polygon
  for (int p = 0; p < nPolys; p++) {
    ec_find_start(&st, starts[p]);
  }

  // CutKeyholes (process holes)
  for (int h = 0; h < st.numHoles; h++) {
    ec_cut_keyhole(&st, h);
  }

  // TriangulatePoly for each simple polygon
  for (int s = 0; s < st.numSimples; s++) {
    ec_triangulate_poly(&st, st.simples[s]);
  }

  // Cleanup
  ManifoldVecIVec3 result = st.triangles;
  ec_heap_free(&st.ears);
  free(st.verts);
  free(st.holes);
  free(st.holeBBox);
  free(st.outers);
  free(st.simples);
  free(starts);
  free(idx);

  return result;
}

#endif // MANIFOLD_EARCLIP_H
