// Copyright 2021 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// Boolean operations for the C11 Manifold port.
// This implements the core mesh boolean algorithm.
//
// TODO: Full implementation of the Manifold boolean algorithm requires:
// 1. Edge-face intersection detection (using collider BVH)
// 2. Winding number computation
// 3. Result mesh construction from intersection curves
//
// For now, this provides a basic union via mesh combination,
// which works for non-overlapping meshes.

#include "manifold_boolean.h"
#include <string.h>

// Simple union: combine vertices and triangles from both meshes
// This only works correctly for non-overlapping meshes
static ManifoldError manifold_simple_union(ManifoldImpl *result,
                                            const ManifoldImpl *p,
                                            const ManifoldImpl *q) {
  manifold_impl_init(result);

  size_t pNumVert = manifold_impl_num_vert(p);
  size_t qNumVert = manifold_impl_num_vert(q);
  size_t pNumHe = p->halfedge.len;
  size_t qNumHe = q->halfedge.len;

  // Copy vertices
  result->vertPos = vec_vec3_create_n(pNumVert + qNumVert);
  for (size_t i = 0; i < pNumVert; i++)
    result->vertPos.data[i] = p->vertPos.data[i];
  for (size_t i = 0; i < qNumVert; i++)
    result->vertPos.data[pNumVert + i] = q->vertPos.data[i];

  // Copy halfedges, offsetting Q's vertex indices
  result->halfedge = vec_halfedge_create_n(pNumHe + qNumHe);
  for (size_t i = 0; i < pNumHe; i++)
    result->halfedge.data[i] = p->halfedge.data[i];
  for (size_t i = 0; i < qNumHe; i++) {
    ManifoldHalfedge he = q->halfedge.data[i];
    he.startVert += (int)pNumVert;
    he.endVert += (int)pNumVert;
    he.pairedHalfedge += (int)pNumHe;
    if (p->numProp == 0) he.propVert = he.startVert;
    result->halfedge.data[pNumHe + i] = he;
  }

  // Copy triRefs
  size_t pNumTri = manifold_impl_num_tri(p);
  size_t qNumTri = manifold_impl_num_tri(q);
  result->meshRelation.triRef = vec_triref_create_n(pNumTri + qNumTri);
  for (size_t i = 0; i < pNumTri && i < p->meshRelation.triRef.len; i++)
    result->meshRelation.triRef.data[i] = p->meshRelation.triRef.data[i];
  for (size_t i = 0; i < qNumTri && i < q->meshRelation.triRef.len; i++)
    result->meshRelation.triRef.data[pNumTri + i] = q->meshRelation.triRef.data[i];

  // Copy mesh ID transforms
  for (size_t i = 0; i < p->meshRelation.meshIDtransform.len; i++) {
    vec_meshid_push(&result->meshRelation.meshIDtransform,
                    p->meshRelation.meshIDtransform.data[i]);
  }
  for (size_t i = 0; i < q->meshRelation.meshIDtransform.len; i++) {
    vec_meshid_push(&result->meshRelation.meshIDtransform,
                    q->meshRelation.meshIDtransform.data[i]);
  }

  result->numProp = p->numProp > q->numProp ? p->numProp : q->numProp;

  manifold_impl_calculate_bbox(result);
  manifold_impl_set_epsilon(result, fmax(p->epsilon, q->epsilon), false);
  result->tolerance = fmax(p->tolerance, q->tolerance);

  manifold_impl_set_normals_and_coplanar(result);

  return MANIFOLD_ERROR_NO_ERROR;
}

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
    switch (op) {
      case MANIFOLD_OP_ADD:
      case MANIFOLD_OP_INTERSECT: {
        // Return copy of q for union, or empty for intersect
        if (op == MANIFOLD_OP_ADD) {
          manifold_impl_init(result);
          result->vertPos = vec_vec3_copy(&q->vertPos);
          result->halfedge = vec_halfedge_copy(&q->halfedge);
          result->faceNormal = vec_vec3_copy(&q->faceNormal);
          result->vertNormal = vec_vec3_copy(&q->vertNormal);
          result->bBox = q->bBox;
          result->epsilon = q->epsilon;
          result->tolerance = q->tolerance;
          return MANIFOLD_ERROR_NO_ERROR;
        }
        manifold_impl_init(result);
        return MANIFOLD_ERROR_NO_ERROR;
      }
      case MANIFOLD_OP_SUBTRACT:
        manifold_impl_init(result);
        return MANIFOLD_ERROR_NO_ERROR;
    }
  }

  if (manifold_impl_is_empty(q)) {
    switch (op) {
      case MANIFOLD_OP_ADD:
      case MANIFOLD_OP_SUBTRACT: {
        manifold_impl_init(result);
        result->vertPos = vec_vec3_copy(&p->vertPos);
        result->halfedge = vec_halfedge_copy(&p->halfedge);
        result->faceNormal = vec_vec3_copy(&p->faceNormal);
        result->vertNormal = vec_vec3_copy(&p->vertNormal);
        result->bBox = p->bBox;
        result->epsilon = p->epsilon;
        result->tolerance = p->tolerance;
        return MANIFOLD_ERROR_NO_ERROR;
      }
      case MANIFOLD_OP_INTERSECT:
        manifold_impl_init(result);
        return MANIFOLD_ERROR_NO_ERROR;
    }
  }

  // Check if bounding boxes don't overlap - shortcut for non-overlapping
  if (!manifold_box_overlaps(p->bBox, q->bBox)) {
    switch (op) {
      case MANIFOLD_OP_ADD:
        return manifold_simple_union(result, p, q);
      case MANIFOLD_OP_SUBTRACT: {
        // Q doesn't touch P, so result is just P
        manifold_impl_init(result);
        result->vertPos = vec_vec3_copy(&p->vertPos);
        result->halfedge = vec_halfedge_copy(&p->halfedge);
        result->faceNormal = vec_vec3_copy(&p->faceNormal);
        result->vertNormal = vec_vec3_copy(&p->vertNormal);
        result->bBox = p->bBox;
        result->epsilon = p->epsilon;
        result->tolerance = p->tolerance;
        return MANIFOLD_ERROR_NO_ERROR;
      }
      case MANIFOLD_OP_INTERSECT:
        manifold_impl_init(result);
        return MANIFOLD_ERROR_NO_ERROR;
    }
  }

  // For overlapping meshes, use simple union for now
  // TODO: Implement full boolean algorithm with edge-face intersections
  if (op == MANIFOLD_OP_ADD) {
    return manifold_simple_union(result, p, q);
  }

  // Fallback: return p for subtract, empty for intersect
  manifold_impl_init(result);
  if (op == MANIFOLD_OP_SUBTRACT) {
    result->vertPos = vec_vec3_copy(&p->vertPos);
    result->halfedge = vec_halfedge_copy(&p->halfedge);
    result->faceNormal = vec_vec3_copy(&p->faceNormal);
    result->vertNormal = vec_vec3_copy(&p->vertNormal);
    result->bBox = p->bBox;
    result->epsilon = p->epsilon;
    result->tolerance = p->tolerance;
  }
  return MANIFOLD_ERROR_NO_ERROR;
}
