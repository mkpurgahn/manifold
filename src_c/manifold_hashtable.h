// Copyright 2022 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// Hash table for the C11 Manifold port.

#ifndef MANIFOLD_HASHTABLE_H
#define MANIFOLD_HASHTABLE_H

#include "manifold_vec.h"
#include "manifold_vec_math.h"

#define MANIFOLD_HASH_OPEN UINT64_MAX

// Hash table storing uint64_t keys → int values (most common usage)
typedef struct {
  ManifoldVecU64 keys;
  ManifoldVecInt values;
  size_t used;
  uint32_t step;
} ManifoldHashTable;

static inline size_t manifold_hash_next_pow2(size_t n) {
  if (n == 0) return 0;
  size_t p = 1;
  while (p < n) p <<= 1;
  return p;
}

static inline ManifoldHashTable manifold_hashtable_create(size_t size,
                                                          uint32_t step) {
  ManifoldHashTable ht;
  size_t sz = size == 0 ? 0 : manifold_hash_next_pow2(size);
  ht.keys = vec_u64_create_fill(sz, MANIFOLD_HASH_OPEN);
  int zero = 0;
  ht.values = vec_int_create_fill(sz, zero);
  ht.used = 0;
  ht.step = step;
  return ht;
}

static inline void manifold_hashtable_free(ManifoldHashTable *ht) {
  vec_u64_free(&ht->keys);
  vec_int_free(&ht->values);
  ht->used = 0;
}

static inline bool manifold_hashtable_full(const ManifoldHashTable *ht) {
  return ht->used * 2 > ht->keys.len;
}

static inline void manifold_hashtable_insert(ManifoldHashTable *ht,
                                              uint64_t key, int val) {
  if (ht->keys.len == 0) return;
  uint32_t idx = (uint32_t)(manifold_hash64bit(key) & (ht->keys.len - 1));
  while (1) {
    if (manifold_hashtable_full(ht)) return;
    uint64_t k = ht->keys.data[idx];
    if (k == MANIFOLD_HASH_OPEN) {
      ht->keys.data[idx] = key;
      ht->used++;
      ht->values.data[idx] = val;
      return;
    }
    if (k == key) return;
    idx = (idx + ht->step) & (uint32_t)(ht->keys.len - 1);
  }
}

static inline int manifold_hashtable_lookup(const ManifoldHashTable *ht,
                                             uint64_t key) {
  if (ht->keys.len == 0) return -1;
  uint32_t idx = (uint32_t)(manifold_hash64bit(key) & (ht->keys.len - 1));
  while (1) {
    uint64_t k = ht->keys.data[idx];
    if (k == key) return ht->values.data[idx];
    if (k == MANIFOLD_HASH_OPEN) return -1;
    idx = (idx + ht->step) & (uint32_t)(ht->keys.len - 1);
  }
}

#endif // MANIFOLD_HASHTABLE_H
