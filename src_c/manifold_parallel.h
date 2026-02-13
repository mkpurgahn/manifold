// Copyright 2022 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// Parallel execution framework for the C11 Manifold port.
// C equivalent of src/parallel.h - uses pthreads for parallelism.

#ifndef MANIFOLD_PARALLEL_H_
#define MANIFOLD_PARALLEL_H_

#include <stdlib.h>
#include <string.h>
#include <pthread.h>

// Execution policy: determines whether to run in parallel or sequentially.
typedef enum {
  MANIFOLD_SEQ = 0,
  MANIFOLD_PAR = 1
} ManifoldExecPolicy;

// Threshold for auto-selecting parallel vs sequential execution.
#define MANIFOLD_SEQ_THRESHOLD 10000

static inline ManifoldExecPolicy manifold_auto_policy(size_t size) {
  return (size > MANIFOLD_SEQ_THRESHOLD) ? MANIFOLD_PAR : MANIFOLD_SEQ;
}

// ===== Parallel for_each =====

typedef struct {
  void (*func)(size_t idx, void *ctx);
  void *ctx;
  size_t start;
  size_t end;
} ManifoldParallelRange;

static inline void *manifold_parallel_worker(void *arg) {
  ManifoldParallelRange *range = (ManifoldParallelRange *)arg;
  for (size_t i = range->start; i < range->end; i++) {
    range->func(i, range->ctx);
  }
  return NULL;
}

// Parallel for_each over a range [start, end).
// Uses pthreads when policy is MANIFOLD_PAR, otherwise runs sequentially.
static inline void manifold_parallel_for(ManifoldExecPolicy policy,
                                         size_t start, size_t end,
                                         void (*func)(size_t, void *),
                                         void *ctx) {
  if (start >= end) return;
  size_t count = end - start;

  if (policy == MANIFOLD_SEQ || count < MANIFOLD_SEQ_THRESHOLD) {
    for (size_t i = start; i < end; i++) {
      func(i, ctx);
    }
    return;
  }

  // Determine number of threads (up to 8)
  int numThreads = 8;
  if (count < (size_t)numThreads * 2) numThreads = (int)(count / 2);
  if (numThreads < 1) numThreads = 1;

  ManifoldParallelRange *ranges =
      (ManifoldParallelRange *)malloc((size_t)numThreads * sizeof(ManifoldParallelRange));
  pthread_t *threads = (pthread_t *)malloc((size_t)numThreads * sizeof(pthread_t));

  size_t chunkSize = count / (size_t)numThreads;
  size_t remainder = count % (size_t)numThreads;
  size_t offset = start;

  for (int t = 0; t < numThreads; t++) {
    ranges[t].func = func;
    ranges[t].ctx = ctx;
    ranges[t].start = offset;
    ranges[t].end = offset + chunkSize + (t < (int)remainder ? 1 : 0);
    offset = ranges[t].end;
  }

  // Launch threads (first range runs on current thread)
  for (int t = 1; t < numThreads; t++) {
    pthread_create(&threads[t], NULL, manifold_parallel_worker, &ranges[t]);
  }
  manifold_parallel_worker(&ranges[0]);

  for (int t = 1; t < numThreads; t++) {
    pthread_join(threads[t], NULL);
  }

  free(ranges);
  free(threads);
}

// ===== Parallel reduce =====

typedef struct {
  void *(*reduce_func)(const void *a, const void *b, void *ctx);
  void *(*map_func)(size_t idx, void *ctx);
  void *ctx;
  size_t start;
  size_t end;
  void *result;
  size_t elem_size;
} ManifoldReduceRange;

// Sequential reduce with map: result = reduce(map(start), map(start+1), ...)
static inline void manifold_seq_map_reduce(
    size_t start, size_t end, void *init, size_t elem_size,
    void *(*map_func)(size_t idx, void *ctx),
    void *(*reduce_func)(const void *a, const void *b, void *ctx),
    void *ctx, void *out) {
  memcpy(out, init, elem_size);
  for (size_t i = start; i < end; i++) {
    void *mapped = map_func(i, ctx);
    void *reduced = reduce_func(out, mapped, ctx);
    memcpy(out, reduced, elem_size);
  }
}

// ===== Parallel exclusive scan =====

// Sequential exclusive scan: out[i] = init op in[0] op ... op in[i-1]
static inline void manifold_seq_exclusive_scan(
    const void *in, void *out, size_t count, size_t elem_size,
    const void *init,
    void (*op)(const void *a, const void *b, void *result)) {
  if (count == 0) return;
  char *cin = (char *)in;
  char *cout = (char *)out;
  char *acc = (char *)malloc(elem_size);
  char *tmp = (char *)malloc(elem_size);
  memcpy(acc, init, elem_size);
  for (size_t i = 0; i < count; i++) {
    memcpy(cout + i * elem_size, acc, elem_size);
    op(acc, cin + i * elem_size, tmp);
    memcpy(acc, tmp, elem_size);
  }
  free(acc);
  free(tmp);
}

#endif  // MANIFOLD_PARALLEL_H_
