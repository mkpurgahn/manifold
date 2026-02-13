// Copyright 2025 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// 2D KD-tree for polygon triangulation - C11 port of src/tree2d.cpp.
// Simple alternating x/y partition tree for spatial queries on polygon vertices.

#ifndef MANIFOLD_TREE2D_H_
#define MANIFOLD_TREE2D_H_

#include "manifold_vec_math.h"

// PolyVert: a 2D polygon vertex with position and metadata.
typedef struct {
  ManifoldVec2 pos;
  int idx;      // original vertex index
  int mesh_idx; // mesh/polygon index
} ManifoldPolyVert;

// ManifoldRect and manifold_rect_overlaps are already defined in
// manifold_types.h and manifold_vec_math.h respectively.
// manifold_rect_contains_point is also in manifold_vec_math.h.

// Compare functions for sorting
static inline int polyvert_cmp_x(const void *a, const void *b) {
  const ManifoldPolyVert *pa = (const ManifoldPolyVert *)a;
  const ManifoldPolyVert *pb = (const ManifoldPolyVert *)b;
  if (pa->pos.x < pb->pos.x) return -1;
  if (pa->pos.x > pb->pos.x) return 1;
  return 0;
}

static inline int polyvert_cmp_y(const void *a, const void *b) {
  const ManifoldPolyVert *pa = (const ManifoldPolyVert *)a;
  const ManifoldPolyVert *pb = (const ManifoldPolyVert *)b;
  if (pa->pos.y < pb->pos.y) return -1;
  if (pa->pos.y > pb->pos.y) return 1;
  return 0;
}

// Build 2D tree recursively by sorting alternating on x/y.
static inline void manifold_build_2d_tree_impl(ManifoldPolyVert *points,
                                                size_t count, bool sortX) {
  if (count < 2) return;
  if (sortX)
    qsort(points, count, sizeof(ManifoldPolyVert), polyvert_cmp_x);
  else
    qsort(points, count, sizeof(ManifoldPolyVert), polyvert_cmp_y);
  size_t mid = count / 2;
  manifold_build_2d_tree_impl(points, mid, !sortX);
  if (mid + 1 < count)
    manifold_build_2d_tree_impl(points + mid + 1, count - mid - 1, !sortX);
}

// Build 2D KD-tree. Points array is sorted in-place.
static inline void manifold_build_2d_tree(ManifoldPolyVert *points,
                                           size_t count) {
  if (count <= 8) return;
  manifold_build_2d_tree_impl(points, count, true);
}

// Query 2D tree for all points within rectangle r, calling callback for each.
typedef void (*ManifoldPolyVertCallback)(const ManifoldPolyVert *p, void *ctx);

static inline void manifold_query_2d_tree(ManifoldPolyVert *points,
                                           size_t count, ManifoldRect r,
                                           ManifoldPolyVertCallback callback,
                                           void *ctx) {
  if (count <= 8) {
    for (size_t i = 0; i < count; i++) {
      if (manifold_rect_contains_point(r, points[i].pos))
        callback(&points[i], ctx);
    }
    return;
  }

  ManifoldRect current;
  current.min.x = current.min.y = -1e308;
  current.max.x = current.max.y = 1e308;

  int level = 0;

  // Explicit stack for traversal (matching C++ max depth of 64)
  ManifoldRect rectStack[64];
  ManifoldPolyVert *viewStart[64];
  size_t viewCount[64];
  int levelStack[64];
  int sp = 0;

  ManifoldPolyVert *curStart = points;
  size_t curCount = count;

  for (;;) {
    if (curCount <= 8) {
      for (size_t i = 0; i < curCount; i++) {
        if (manifold_rect_contains_point(r, curStart[i].pos))
          callback(&curStart[i], ctx);
      }
      if (--sp < 0) break;
      level = levelStack[sp];
      curStart = viewStart[sp];
      curCount = viewCount[sp];
      current = rectStack[sp];
      continue;
    }

    size_t mid = curCount / 2;
    ManifoldPolyVert middle = curStart[mid];

    ManifoldRect left = current;
    ManifoldRect right = current;
    if (level % 2 == 0) {
      left.max.x = middle.pos.x;
      right.min.x = middle.pos.x;
    } else {
      left.max.y = middle.pos.y;
      right.min.y = middle.pos.y;
    }

    if (manifold_rect_contains_point(r, middle.pos))
      callback(&middle, ctx);

    bool goLeft = manifold_rect_overlaps(left, r);
    bool goRight = manifold_rect_overlaps(right, r);

    if (goLeft) {
      if (goRight && sp < 64) {
        rectStack[sp] = right;
        viewStart[sp] = curStart + mid + 1;
        viewCount[sp] = curCount - mid - 1;
        levelStack[sp] = level + 1;
        sp++;
      }
      current = left;
      curCount = mid;
      level++;
    } else if (goRight) {
      current = right;
      curStart = curStart + mid + 1;
      curCount = curCount - mid - 1;
      level++;
    } else {
      if (--sp < 0) break;
      level = levelStack[sp];
      curStart = viewStart[sp];
      curCount = viewCount[sp];
      current = rectStack[sp];
    }
  }
}

#endif  // MANIFOLD_TREE2D_H_
