// Disjoint sets (Union-Find) for C11 Manifold port.
// Based on dset by Wenzel Jakob, with path compression and union by rank.
// Simplified: single-threaded (no atomics needed for serial).

#ifndef MANIFOLD_DISJOINT_SETS_H
#define MANIFOLD_DISJOINT_SETS_H

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  uint32_t *parent;
  uint32_t *rank_;
  uint32_t size;
} ManifoldDisjointSets;

static inline ManifoldDisjointSets manifold_disjoint_sets_create(uint32_t size) {
  ManifoldDisjointSets ds;
  ds.size = size;
  ds.parent = (uint32_t *)malloc(size * sizeof(uint32_t));
  ds.rank_ = (uint32_t *)calloc(size, sizeof(uint32_t));
  for (uint32_t i = 0; i < size; i++) ds.parent[i] = i;
  return ds;
}

static inline void manifold_disjoint_sets_free(ManifoldDisjointSets *ds) {
  free(ds->parent);
  free(ds->rank_);
  ds->parent = NULL;
  ds->rank_ = NULL;
  ds->size = 0;
}

static inline uint32_t manifold_disjoint_sets_find(ManifoldDisjointSets *ds,
                                                    uint32_t id) {
  while (ds->parent[id] != id) {
    ds->parent[id] = ds->parent[ds->parent[id]]; // path compression
    id = ds->parent[id];
  }
  return id;
}

static inline uint32_t manifold_disjoint_sets_unite(ManifoldDisjointSets *ds,
                                                     uint32_t id1,
                                                     uint32_t id2) {
  id1 = manifold_disjoint_sets_find(ds, id1);
  id2 = manifold_disjoint_sets_find(ds, id2);
  if (id1 == id2) return id1;

  if (ds->rank_[id1] < ds->rank_[id2] ||
      (ds->rank_[id1] == ds->rank_[id2] && id1 > id2)) {
    uint32_t tmp = id1; id1 = id2; id2 = tmp;
  }
  ds->parent[id2] = id1;
  if (ds->rank_[id1] == ds->rank_[id2]) ds->rank_[id1]++;
  return id1;
}

static inline bool manifold_disjoint_sets_same(ManifoldDisjointSets *ds,
                                                uint32_t id1, uint32_t id2) {
  return manifold_disjoint_sets_find(ds, id1) ==
         manifold_disjoint_sets_find(ds, id2);
}

static inline int manifold_disjoint_sets_connected_components(
    ManifoldDisjointSets *ds, int *components) {
  // Map each root to a component label
  int *labels = (int *)malloc(ds->size * sizeof(int));
  memset(labels, -1, ds->size * sizeof(int));
  int numComponents = 0;
  for (uint32_t i = 0; i < ds->size; i++) {
    uint32_t root = manifold_disjoint_sets_find(ds, i);
    if (labels[root] == -1) {
      labels[root] = numComponents++;
    }
    components[i] = labels[root];
  }
  free(labels);
  return numComponents;
}

#endif // MANIFOLD_DISJOINT_SETS_H
