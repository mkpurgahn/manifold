// Copyright 2021 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// BVH collider for the C11 Manifold port.

#ifndef MANIFOLD_COLLIDER_H
#define MANIFOLD_COLLIDER_H

#include "manifold_vec.h"
#include "manifold_vec_math.h"

// Internal node layout: even = leaf, odd = internal, root = 1
#define MANIFOLD_COLLIDER_ROOT 1

static inline bool collider_is_leaf(int node) { return node % 2 == 0; }
static inline bool collider_is_internal(int node) { return node % 2 == 1; }
static inline int collider_node2internal(int node) { return (node - 1) / 2; }
static inline int collider_internal2node(int internal) { return internal * 2 + 1; }
static inline int collider_node2leaf(int node) { return node / 2; }
static inline int collider_leaf2node(int leaf) { return leaf * 2; }

static inline uint32_t collider_spread_bits3(uint32_t v) {
  v = 0xFF0000FFu & (v * 0x00010001u);
  v = 0x0F00F00Fu & (v * 0x00000101u);
  v = 0xC30C30C3u & (v * 0x00000011u);
  v = 0x49249249u & (v * 0x00000005u);
  return v;
}

static inline uint32_t manifold_morton_code(ManifoldVec3 position,
                                             ManifoldBox bBox) {
  ManifoldVec3 range = vec3_sub(bBox.max, bBox.min);
  ManifoldVec3 xyz = vec3_div(vec3_sub(position, bBox.min), range);
  xyz = vec3_min(manifold_vec3_splat(1023.0),
                 vec3_max(manifold_vec3_splat(0.0),
                          vec3_scale(xyz, 1024.0)));
  uint32_t x = collider_spread_bits3((uint32_t)xyz.x);
  uint32_t y = collider_spread_bits3((uint32_t)xyz.y);
  uint32_t z = collider_spread_bits3((uint32_t)xyz.z);
  return x * 4 + y * 2 + z;
}

typedef struct { int first, second; } ManifoldIntPair;

MANIFOLD_VEC_STRUCT(ManifoldIntPair, ManifoldVecIntPair)
MANIFOLD_VEC_FUNCS(ManifoldIntPair, ManifoldVecIntPair, vec_intpair)

typedef struct {
  ManifoldVecBox nodeBBox;
  ManifoldVecInt nodeParent;
  ManifoldVecIntPair internalChildren;
} ManifoldCollider;

static inline size_t collider_num_internal(const ManifoldCollider *c) {
  return c->internalChildren.len;
}

static inline size_t collider_num_leaves(const ManifoldCollider *c) {
  return c->internalChildren.len == 0 ? 0 : c->internalChildren.len + 1;
}

static inline int collider_clz(uint32_t v) {
  if (v == 0) return 32;
  return __builtin_clz(v);
}

static inline int collider_prefix_length_raw(uint32_t a, uint32_t b) {
  return collider_clz(a ^ b);
}

static inline int collider_prefix_length(const uint32_t *morton, int size,
                                          int i, int j) {
  if (j < 0 || j >= size) return -1;
  if (morton[i] == morton[j])
    return 32 + collider_prefix_length_raw((uint32_t)i, (uint32_t)j);
  return collider_prefix_length_raw(morton[i], morton[j]);
}

static inline int collider_range_end(const uint32_t *morton, int size, int i) {
  int dir = collider_prefix_length(morton, size, i, i + 1) -
            collider_prefix_length(morton, size, i, i - 1);
  dir = (dir > 0) - (dir < 0);
  int commonPrefix = collider_prefix_length(morton, size, i, i - dir);
  int max_length = 128;
  while (collider_prefix_length(morton, size, i, i + dir * max_length) > commonPrefix)
    max_length *= 4;
  int length = 0;
  for (int step = max_length / 2; step > 0; step /= 2) {
    if (collider_prefix_length(morton, size, i, i + dir * (length + step)) > commonPrefix)
      length += step;
  }
  return i + dir * length;
}

static inline int collider_find_split(const uint32_t *morton, int size,
                                       int first, int last) {
  int commonPrefix = collider_prefix_length(morton, size, first, last);
  int split = first;
  int step = last - first;
  do {
    step = (step + 1) >> 1;
    int newSplit = split + step;
    if (newSplit < last) {
      if (collider_prefix_length(morton, size, first, newSplit) > commonPrefix)
        split = newSplit;
    }
  } while (step > 1);
  return split;
}

static inline void collider_build_radix_tree(ManifoldCollider *c,
                                              const uint32_t *morton,
                                              int numLeaves) {
  for (int internal = 0; internal < numLeaves - 1; internal++) {
    int first = internal;
    int last = collider_range_end(morton, numLeaves, first);
    if (first > last) { int t = first; first = last; last = t; }
    int split = collider_find_split(morton, numLeaves, first, last);
    int child1 = (split == first) ? collider_leaf2node(split)
                                  : collider_internal2node(split);
    int s2 = split + 1;
    int child2 = (s2 == last) ? collider_leaf2node(s2)
                              : collider_internal2node(s2);
    c->internalChildren.data[internal].first = child1;
    c->internalChildren.data[internal].second = child2;
    int node = collider_internal2node(internal);
    c->nodeParent.data[child1] = node;
    c->nodeParent.data[child2] = node;
  }
}

static inline void collider_build_internal_boxes(ManifoldCollider *c,
                                                   int numLeaves) {
  ManifoldVecInt counter = vec_int_create_fill((size_t)(numLeaves - 1), 0);
  for (int leaf = 0; leaf < numLeaves; leaf++) {
    int node = collider_leaf2node(leaf);
    while (node != MANIFOLD_COLLIDER_ROOT) {
      node = c->nodeParent.data[node];
      int internal = collider_node2internal(node);
      if (counter.data[internal] == 0) {
        counter.data[internal] = 1;
        break;  // first thread to reach this node, stop
      }
      int c1 = c->internalChildren.data[internal].first;
      int c2 = c->internalChildren.data[internal].second;
      c->nodeBBox.data[node] = manifold_box_union(
          c->nodeBBox.data[c1], c->nodeBBox.data[c2]);
    }
  }
  vec_int_free(&counter);
}

static inline ManifoldCollider manifold_collider_create(
    const ManifoldBox *leafBB, const uint32_t *leafMorton, size_t numLeaves) {
  ManifoldCollider c;
  if (numLeaves == 0) {
    c.nodeBBox = (ManifoldVecBox)MANIFOLD_VEC_INIT;
    c.nodeParent = (ManifoldVecInt)MANIFOLD_VEC_INIT;
    c.internalChildren = (ManifoldVecIntPair)MANIFOLD_VEC_INIT;
    return c;
  }
  size_t numNodes = 2 * numLeaves - 1;
  c.nodeBBox = vec_box_create_n(numNodes);
  c.nodeParent = vec_int_create_fill(numNodes, -1);
  ManifoldIntPair neg1 = {-1, -1};
  c.internalChildren = vec_intpair_create_fill(numLeaves - 1, neg1);

  // Copy leaf boxes into even-indexed positions
  for (size_t i = 0; i < numLeaves; i++) {
    c.nodeBBox.data[i * 2] = leafBB[i];
  }

  collider_build_radix_tree(&c, leafMorton, (int)numLeaves);
  collider_build_internal_boxes(&c, (int)numLeaves);
  return c;
}

static inline void manifold_collider_free(ManifoldCollider *c) {
  vec_box_free(&c->nodeBBox);
  vec_int_free(&c->nodeParent);
  vec_intpair_free(&c->internalChildren);
}

// Deep copy a collider
static inline ManifoldCollider manifold_collider_copy(
    const ManifoldCollider *src) {
  ManifoldCollider c;
  c.nodeBBox = vec_box_copy(&src->nodeBBox);
  c.nodeParent = vec_int_copy(&src->nodeParent);
  c.internalChildren = vec_intpair_copy(&src->internalChildren);
  return c;
}

// Collision query callback type: called with (queryIdx, leafIdx, userdata)
typedef void (*ManifoldCollisionCallback)(int queryIdx, int leafIdx, void *ctx);

static inline void manifold_collider_collisions(
    const ManifoldCollider *c,
    const ManifoldBox *queries, size_t numQueries,
    ManifoldCollisionCallback callback, void *ctx,
    bool selfCollision) {
  if (c->internalChildren.len == 0) return;

  for (size_t qi = 0; qi < numQueries; qi++) {
    int stack[64];
    int top = -1;
    int node = MANIFOLD_COLLIDER_ROOT;
    ManifoldBox queryBox = queries[qi];

    while (1) {
      int internal = collider_node2internal(node);
      int child1 = c->internalChildren.data[internal].first;
      int child2 = c->internalChildren.data[internal].second;

      bool overlap1 = manifold_box_overlaps(c->nodeBBox.data[child1], queryBox);
      bool overlap2 = manifold_box_overlaps(c->nodeBBox.data[child2], queryBox);

      bool traverse1 = false, traverse2 = false;
      if (overlap1 && collider_is_leaf(child1)) {
        int leafIdx = collider_node2leaf(child1);
        if (!selfCollision || leafIdx != (int)qi)
          callback((int)qi, leafIdx, ctx);
      } else {
        traverse1 = overlap1 && collider_is_internal(child1);
      }
      if (overlap2 && collider_is_leaf(child2)) {
        int leafIdx = collider_node2leaf(child2);
        if (!selfCollision || leafIdx != (int)qi)
          callback((int)qi, leafIdx, ctx);
      } else {
        traverse2 = overlap2 && collider_is_internal(child2);
      }

      if (!traverse1 && !traverse2) {
        if (top < 0) break;
        node = stack[top--];
      } else {
        node = traverse1 ? child1 : child2;
        if (traverse1 && traverse2)
          stack[++top] = child2;
      }
    }
  }
}

#endif // MANIFOLD_COLLIDER_H
