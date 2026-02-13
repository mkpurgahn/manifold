// Copyright 2021 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// Generic dynamic array for the C11 Manifold port.
// Uses X-macro pattern for type-specific instantiation.

#ifndef MANIFOLD_VEC_H
#define MANIFOLD_VEC_H

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "manifold_types.h"

// Generic dynamic array structure and operations via macros.
// Usage: MANIFOLD_VEC_DECLARE(int) creates ManifoldVecInt type.

#ifdef __GNUC__
#define MANIFOLD_UNUSED __attribute__((unused))
#else
#define MANIFOLD_UNUSED
#endif

#define MANIFOLD_VEC_STRUCT(T, NAME) \
  typedef struct { \
    T *data; \
    size_t len; \
    size_t cap; \
  } NAME;

#define MANIFOLD_VEC_INIT {NULL, 0, 0}

#define MANIFOLD_VEC_FUNCS(T, NAME, PREFIX) \
\
MANIFOLD_UNUSED static inline NAME PREFIX##_create(size_t cap) { \
  NAME v; \
  v.len = 0; \
  v.cap = cap; \
  v.data = (cap > 0) ? (T*)malloc(cap * sizeof(T)) : NULL; \
  return v; \
} \
\
MANIFOLD_UNUSED static inline NAME PREFIX##_create_n(size_t n) { \
  NAME v; \
  v.len = n; \
  v.cap = n; \
  v.data = (n > 0) ? (T*)malloc(n * sizeof(T)) : NULL; \
  return v; \
} \
\
MANIFOLD_UNUSED static inline NAME PREFIX##_create_fill(size_t n, T val) { \
  NAME v = PREFIX##_create_n(n); \
  for (size_t i = 0; i < n; i++) v.data[i] = val; \
  return v; \
} \
\
MANIFOLD_UNUSED static inline void PREFIX##_free(NAME *v) { \
  free(v->data); \
  v->data = NULL; \
  v->len = 0; \
  v->cap = 0; \
} \
\
MANIFOLD_UNUSED static inline void PREFIX##_reserve(NAME *v, size_t cap) { \
  if (cap > v->cap) { \
    T *newbuf = (T*)realloc(v->data, cap * sizeof(T)); \
    assert(newbuf != NULL); \
    v->data = newbuf; \
    v->cap = cap; \
  } \
} \
\
MANIFOLD_UNUSED static inline void PREFIX##_push(NAME *v, T val) { \
  if (v->len >= v->cap) { \
    size_t newcap = v->cap == 0 ? 128 : v->cap * 2; \
    PREFIX##_reserve(v, newcap); \
  } \
  v->data[v->len++] = val; \
} \
\
MANIFOLD_UNUSED static inline void PREFIX##_resize(NAME *v, size_t n) { \
  PREFIX##_reserve(v, n); \
  v->len = n; \
} \
\
MANIFOLD_UNUSED static inline void PREFIX##_resize_fill(NAME *v, size_t n, T val) { \
  size_t old_len = v->len; \
  PREFIX##_resize(v, n); \
  for (size_t i = old_len; i < n; i++) v->data[i] = val; \
} \
\
MANIFOLD_UNUSED static inline void PREFIX##_clear(NAME *v) { \
  v->len = 0; \
} \
\
MANIFOLD_UNUSED static inline void PREFIX##_pop(NAME *v) { \
  assert(v->len > 0); \
  v->len--; \
} \
\
MANIFOLD_UNUSED static inline T PREFIX##_back(const NAME *v) { \
  assert(v->len > 0); \
  return v->data[v->len - 1]; \
} \
\
MANIFOLD_UNUSED static inline void PREFIX##_shrink_to_fit(NAME *v) { \
  if (v->len == 0) { \
    free(v->data); \
    v->data = NULL; \
    v->cap = 0; \
  } else if (v->cap > v->len) { \
    T *newbuf = (T*)realloc(v->data, v->len * sizeof(T)); \
    if (newbuf) { v->data = newbuf; v->cap = v->len; } \
  } \
} \
\
MANIFOLD_UNUSED static inline NAME PREFIX##_copy(const NAME *src) { \
  NAME v; \
  v.len = src->len; \
  v.cap = src->len; \
  if (v.len > 0) { \
    v.data = (T*)malloc(v.len * sizeof(T)); \
    assert(v.data != NULL); \
    memcpy(v.data, src->data, v.len * sizeof(T)); \
  } else { \
    v.data = NULL; \
  } \
  return v; \
} \
\
MANIFOLD_UNUSED static inline void PREFIX##_swap(NAME *a, NAME *b) { \
  NAME tmp = *a; \
  *a = *b; \
  *b = tmp; \
}

// Declare common Vec types

MANIFOLD_VEC_STRUCT(int, ManifoldVecInt)
MANIFOLD_VEC_FUNCS(int, ManifoldVecInt, vec_int)

MANIFOLD_VEC_STRUCT(uint32_t, ManifoldVecU32)
MANIFOLD_VEC_FUNCS(uint32_t, ManifoldVecU32, vec_u32)

MANIFOLD_VEC_STRUCT(uint64_t, ManifoldVecU64)
MANIFOLD_VEC_FUNCS(uint64_t, ManifoldVecU64, vec_u64)

MANIFOLD_VEC_STRUCT(size_t, ManifoldVecSize)
MANIFOLD_VEC_FUNCS(size_t, ManifoldVecSize, vec_size)

MANIFOLD_VEC_STRUCT(double, ManifoldVecDouble)
MANIFOLD_VEC_FUNCS(double, ManifoldVecDouble, vec_double)

MANIFOLD_VEC_STRUCT(ManifoldVec2, ManifoldVecVec2)
MANIFOLD_VEC_FUNCS(ManifoldVec2, ManifoldVecVec2, vec_vec2)

MANIFOLD_VEC_STRUCT(ManifoldVec3, ManifoldVecVec3)
MANIFOLD_VEC_FUNCS(ManifoldVec3, ManifoldVecVec3, vec_vec3)

MANIFOLD_VEC_STRUCT(ManifoldVec4, ManifoldVecVec4)
MANIFOLD_VEC_FUNCS(ManifoldVec4, ManifoldVecVec4, vec_vec4)

MANIFOLD_VEC_STRUCT(ManifoldIVec3, ManifoldVecIVec3)
MANIFOLD_VEC_FUNCS(ManifoldIVec3, ManifoldVecIVec3, vec_ivec3)

MANIFOLD_VEC_STRUCT(ManifoldIVec4, ManifoldVecIVec4)
MANIFOLD_VEC_FUNCS(ManifoldIVec4, ManifoldVecIVec4, vec_ivec4)

MANIFOLD_VEC_STRUCT(ManifoldHalfedge, ManifoldVecHalfedge)
MANIFOLD_VEC_FUNCS(ManifoldHalfedge, ManifoldVecHalfedge, vec_halfedge)

MANIFOLD_VEC_STRUCT(ManifoldTriRef, ManifoldVecTriRef)
MANIFOLD_VEC_FUNCS(ManifoldTriRef, ManifoldVecTriRef, vec_triref)

MANIFOLD_VEC_STRUCT(ManifoldBarycentric, ManifoldVecBarycentric)
MANIFOLD_VEC_FUNCS(ManifoldBarycentric, ManifoldVecBarycentric, vec_bary)

MANIFOLD_VEC_STRUCT(ManifoldTmpEdge, ManifoldVecTmpEdge)
MANIFOLD_VEC_FUNCS(ManifoldTmpEdge, ManifoldVecTmpEdge, vec_tmpedge)

MANIFOLD_VEC_STRUCT(ManifoldBox, ManifoldVecBox)
MANIFOLD_VEC_FUNCS(ManifoldBox, ManifoldVecBox, vec_box)

MANIFOLD_VEC_STRUCT(ManifoldSmoothness, ManifoldVecSmoothness)
MANIFOLD_VEC_FUNCS(ManifoldSmoothness, ManifoldVecSmoothness, vec_smooth)

MANIFOLD_VEC_STRUCT(bool, ManifoldVecBool)
MANIFOLD_VEC_FUNCS(bool, ManifoldVecBool, vec_bool)

MANIFOLD_VEC_STRUCT(char, ManifoldVecChar)
MANIFOLD_VEC_FUNCS(char, ManifoldVecChar, vec_char)

// Int pair type for boolean edges  
typedef struct { int v[2]; } ManifoldIntArr2;
MANIFOLD_VEC_STRUCT(ManifoldIntArr2, ManifoldVecIntArr2)
MANIFOLD_VEC_FUNCS(ManifoldIntArr2, ManifoldVecIntArr2, vec_intarr2)

// Sequence fill: fills v->data[0..n-1] with 0, 1, 2, ...
static inline void vec_int_sequence(ManifoldVecInt *v) {
  for (size_t i = 0; i < v->len; i++) v->data[i] = (int)i;
}

static inline void vec_size_sequence(ManifoldVecSize *v) {
  for (size_t i = 0; i < v->len; i++) v->data[i] = i;
}

// Permute: rearrange src into dst according to new2old mapping
#define MANIFOLD_VEC_PERMUTE(T, NAME, PREFIX) \
static inline void PREFIX##_permute(NAME *inOut, const ManifoldVecInt *new2Old) { \
  NAME tmp = PREFIX##_copy(inOut); \
  PREFIX##_resize(inOut, new2Old->len); \
  for (size_t i = 0; i < new2Old->len; i++) { \
    inOut->data[i] = tmp.data[new2Old->data[i]]; \
  } \
  PREFIX##_free(&tmp); \
}

MANIFOLD_VEC_PERMUTE(ManifoldVec3, ManifoldVecVec3, vec_vec3)
MANIFOLD_VEC_PERMUTE(ManifoldVec4, ManifoldVecVec4, vec_vec4)
MANIFOLD_VEC_PERMUTE(ManifoldBox, ManifoldVecBox, vec_box)
MANIFOLD_VEC_PERMUTE(uint32_t, ManifoldVecU32, vec_u32)
MANIFOLD_VEC_PERMUTE(ManifoldTriRef, ManifoldVecTriRef, vec_triref)
MANIFOLD_VEC_PERMUTE(ManifoldHalfedge, ManifoldVecHalfedge, vec_halfedge)
MANIFOLD_VEC_PERMUTE(int, ManifoldVecInt, vec_int)
MANIFOLD_VEC_PERMUTE(ManifoldIntArr2, ManifoldVecIntArr2, vec_intarr2)
MANIFOLD_VEC_PERMUTE(size_t, ManifoldVecSize, vec_size)

#endif // MANIFOLD_VEC_H
