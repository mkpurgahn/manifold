// Copyright 2022 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// CSG tree operations for the C11 Manifold port.
// C equivalent of src/csg_tree.cpp.
//
// In the C++ version, the CSG tree uses lazy evaluation with virtual dispatch,
// shared_ptr, and caching. In this C11 port, boolean operations are evaluated
// eagerly (directly). This file implements the Compose function which merges
// pairwise disjoint meshes without boolean operations.
//
// Key C++ → C translations:
// - CsgLeafNode::Compose → manifold_compose_impls
// - BatchBoolean → manifold_batch_boolean (in manifold_api.c)
// - BatchUnion → uses compose + batch_boolean
// - SimpleBoolean → manifold_boolean3 (in manifold_boolean.c)
// - CsgOpNode::ToLeafNode → replaced by eager evaluation

#include "manifold_csg_tree.h"
#include <stdlib.h>
#include <string.h>

// Compose: merge pairwise disjoint meshes into one without boolean operations.
// This is the C equivalent of CsgLeafNode::Compose in csg_tree.cpp.
void manifold_compose_impls(const ManifoldImpl **impls, int count,
                            ManifoldImpl *out) {
  if (count <= 0) {
    manifold_impl_init(out);
    return;
  }
  if (count == 1) {
    // Simple copy of the single impl
    manifold_impl_init(out);
    size_t nv = manifold_impl_num_vert(impls[0]);
    size_t nt = manifold_impl_num_tri(impls[0]);
    size_t nhe = impls[0]->halfedge.len;
    vec_vec3_resize(&out->vertPos, nv);
    vec_vec3_resize(&out->vertNormal, nv);
    vec_halfedge_resize(&out->halfedge, nhe);
    vec_vec3_resize(&out->faceNormal, nt);
    vec_triref_resize(&out->meshRelation.triRef, nt);
    if (nv > 0) {
      memcpy(out->vertPos.data, impls[0]->vertPos.data, nv * sizeof(ManifoldVec3));
      if (impls[0]->vertNormal.len >= nv)
        memcpy(out->vertNormal.data, impls[0]->vertNormal.data, nv * sizeof(ManifoldVec3));
    }
    if (nhe > 0)
      memcpy(out->halfedge.data, impls[0]->halfedge.data, nhe * sizeof(ManifoldHalfedge));
    if (nt > 0) {
      if (impls[0]->faceNormal.len >= nt)
        memcpy(out->faceNormal.data, impls[0]->faceNormal.data, nt * sizeof(ManifoldVec3));
      if (impls[0]->meshRelation.triRef.len >= nt)
        memcpy(out->meshRelation.triRef.data, impls[0]->meshRelation.triRef.data, nt * sizeof(ManifoldTriRef));
    }
    out->epsilon = impls[0]->epsilon;
    out->tolerance = impls[0]->tolerance;
    if (impls[0]->numProp > 0) {
      out->numProp = impls[0]->numProp;
      out->properties = vec_double_copy(&impls[0]->properties);
    }
    manifold_impl_calculate_bbox(out);
    return;
  }

  // Calculate total sizes
  size_t totalVert = 0;
  size_t totalEdge = 0;
  size_t totalTri = 0;
  size_t totalPropVert = 0;
  int numPropOut = 0;
  double epsilon = -1;
  double tolerance = -1;

  for (int i = 0; i < count; i++) {
    totalVert += manifold_impl_num_vert(impls[i]);
    totalEdge += manifold_impl_num_edge(impls[i]);
    totalTri += manifold_impl_num_tri(impls[i]);
    int np = impls[i]->numProp;
    if (np > numPropOut) numPropOut = np;
    totalPropVert += (np == 0) ? 1 : (impls[i]->properties.len / (size_t)np);
    if (impls[i]->epsilon > epsilon) epsilon = impls[i]->epsilon;
    if (impls[i]->tolerance > tolerance) tolerance = impls[i]->tolerance;
  }

  manifold_impl_init(out);
  out->epsilon = epsilon;
  out->tolerance = tolerance;

  // Allocate combined arrays
  vec_vec3_resize(&out->vertPos, totalVert);
  vec_vec3_resize(&out->vertNormal, totalVert);
  vec_halfedge_resize(&out->halfedge, 2 * totalEdge);
  vec_vec3_resize(&out->faceNormal, totalTri);
  vec_triref_resize(&out->meshRelation.triRef, totalTri);
  if (numPropOut > 0) {
    out->numProp = numPropOut;
    vec_double_resize(&out->properties, (size_t)numPropOut * totalPropVert);
    memset(out->properties.data, 0, out->properties.len * sizeof(double));
  }

  // Copy data from each impl
  size_t vertOff = 0;
  size_t edgeOff = 0;
  size_t triOff = 0;
  size_t propOff = 0;

  for (int i = 0; i < count; i++) {
    const ManifoldImpl *src = impls[i];
    size_t nv = manifold_impl_num_vert(src);
    size_t ne = manifold_impl_num_edge(src);
    size_t nt = manifold_impl_num_tri(src);
    size_t nhe = 2 * ne;

    // Copy vertex positions and normals
    if (nv > 0) {
      memcpy(&out->vertPos.data[vertOff], src->vertPos.data,
             nv * sizeof(ManifoldVec3));
      if (src->vertNormal.len >= nv) {
        memcpy(&out->vertNormal.data[vertOff], src->vertNormal.data,
               nv * sizeof(ManifoldVec3));
      }
    }

    // Copy face normals
    if (nt > 0 && src->faceNormal.len >= nt) {
      memcpy(&out->faceNormal.data[triOff], src->faceNormal.data,
             nt * sizeof(ManifoldVec3));
    }

    // Copy halfedges with offset
    for (size_t j = 0; j < nhe; j++) {
      ManifoldHalfedge h = src->halfedge.data[j];
      h.startVert += (int)vertOff;
      h.endVert += (int)vertOff;
      h.pairedHalfedge += (int)edgeOff;
      if (numPropOut > 0) {
        if (src->numProp == 0) h.propVert = 0;
        h.propVert += (int)propOff;
      }
      out->halfedge.data[edgeOff + j] = h;
    }

    // Copy properties with stride conversion
    if (src->numProp > 0 && numPropOut > 0) {
      int srcNP = src->numProp;
      size_t srcPropVert = src->properties.len / (size_t)srcNP;
      for (size_t v = 0; v < srcPropVert; v++) {
        for (int p = 0; p < srcNP; p++) {
          out->properties.data[(propOff + v) * (size_t)numPropOut + p] =
              src->properties.data[v * (size_t)srcNP + p];
        }
      }
    }

    // Copy triRef with mesh ID offset
    for (size_t j = 0; j < nt; j++) {
      ManifoldTriRef ref = src->meshRelation.triRef.data[j];
      out->meshRelation.triRef.data[triOff + j] = ref;
    }

    // Copy meshIDtransform entries
    for (size_t j = 0; j < src->meshRelation.meshIDtransform.len; j++) {
      ManifoldMeshIDEntry entry = src->meshRelation.meshIDtransform.data[j];
      // Check if already present
      bool found = false;
      for (size_t k = 0; k < out->meshRelation.meshIDtransform.len; k++) {
        if (out->meshRelation.meshIDtransform.data[k].key == entry.key) {
          found = true;
          break;
        }
      }
      if (!found) {
        vec_meshid_push(&out->meshRelation.meshIDtransform, entry);
      }
    }

    vertOff += nv;
    edgeOff += nhe;
    triOff += nt;
    propOff += (src->numProp == 0) ? 1 : (src->properties.len / (size_t)src->numProp);
  }

  // Update bounding box
  manifold_impl_calculate_bbox(out);

  // Sort geometry for proper collider construction
  manifold_impl_sort_geometry(out);
}
