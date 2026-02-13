// Triangle-triangle distance for the C11 Manifold port.
// Ported from src/tri_dist.h (NVIDIA PhysX derived).

#include "manifold_tri_dist.h"
#include <math.h>
#include <float.h>

static inline double v3dot(ManifoldVec3 a, ManifoldVec3 b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline ManifoldVec3 v3sub(ManifoldVec3 a, ManifoldVec3 b) {
  return (ManifoldVec3){a.x - b.x, a.y - b.y, a.z - b.z};
}

static inline ManifoldVec3 v3add(ManifoldVec3 a, ManifoldVec3 b) {
  return (ManifoldVec3){a.x + b.x, a.y + b.y, a.z + b.z};
}

static inline ManifoldVec3 v3scale(ManifoldVec3 a, double s) {
  return (ManifoldVec3){a.x * s, a.y * s, a.z * s};
}

static inline ManifoldVec3 v3cross(ManifoldVec3 a, ManifoldVec3 b) {
  return (ManifoldVec3){
    a.y * b.z - a.z * b.y,
    a.z * b.x - a.x * b.z,
    a.x * b.y - a.y * b.x
  };
}

static inline double clamp01(double x) {
  return x < 0.0 ? 0.0 : (x > 1.0 ? 1.0 : x);
}

void manifold_edge_edge_dist(ManifoldVec3 *x, ManifoldVec3 *y,
                             ManifoldVec3 p, ManifoldVec3 a,
                             ManifoldVec3 q, ManifoldVec3 b) {
  ManifoldVec3 T = v3sub(q, p);
  double ADotA = v3dot(a, a);
  double BDotB = v3dot(b, b);
  double ADotB = v3dot(a, b);
  double ADotT = v3dot(a, T);
  double BDotT = v3dot(b, T);

  double Denom = ADotA * BDotB - ADotB * ADotB;
  double t = (Denom != 0.0)
    ? clamp01((ADotT * BDotB - BDotT * ADotB) / Denom)
    : 0.0;

  double u;
  if (BDotB != 0.0) {
    u = (t * ADotB - BDotT) / BDotB;
    if (u < 0.0) {
      u = 0.0;
      t = (ADotA != 0.0) ? clamp01(ADotT / ADotA) : 0.0;
    } else if (u > 1.0) {
      u = 1.0;
      t = (ADotA != 0.0) ? clamp01((ADotB + ADotT) / ADotA) : 0.0;
    }
  } else {
    u = 0.0;
    t = (ADotA != 0.0) ? clamp01(ADotT / ADotA) : 0.0;
  }

  *x = v3add(p, v3scale(a, t));
  *y = v3add(q, v3scale(b, u));
}

double manifold_distance_tri_tri_squared(const ManifoldVec3 p[3],
                                         const ManifoldVec3 q[3]) {
  ManifoldVec3 Sv[3], Tv[3];
  Sv[0] = v3sub(p[1], p[0]);
  Sv[1] = v3sub(p[2], p[1]);
  Sv[2] = v3sub(p[0], p[2]);
  Tv[0] = v3sub(q[1], q[0]);
  Tv[1] = v3sub(q[2], q[1]);
  Tv[2] = v3sub(q[0], q[2]);

  int shown_disjoint = 0;
  double mindd = DBL_MAX;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ManifoldVec3 cp, cq;
      manifold_edge_edge_dist(&cp, &cq, p[i], Sv[i], q[j], Tv[j]);
      ManifoldVec3 V = v3sub(cq, cp);
      double dd = v3dot(V, V);

      if (dd <= mindd) {
        mindd = dd;

        int id = i + 2;
        if (id >= 3) id -= 3;
        ManifoldVec3 Z = v3sub(p[id], cp);
        double a = v3dot(Z, V);
        id = j + 2;
        if (id >= 3) id -= 3;
        Z = v3sub(q[id], cq);
        double b = v3dot(Z, V);

        if ((a <= 0.0) && (b >= 0.0)) {
          return v3dot(V, V);
        }

        if (a <= 0.0)
          a = 0.0;
        else if (b > 0.0)
          b = 0.0;

        if ((mindd - a + b) > 0.0) shown_disjoint = 1;
      }
    }
  }

  // Check if q vertex projects inside p triangle
  ManifoldVec3 Sn = v3cross(Sv[0], Sv[1]);
  double Snl = v3dot(Sn, Sn);

  if (Snl > 1e-15) {
    double Tp[3];
    Tp[0] = v3dot(v3sub(p[0], q[0]), Sn);
    Tp[1] = v3dot(v3sub(p[0], q[1]), Sn);
    Tp[2] = v3dot(v3sub(p[0], q[2]), Sn);

    int index = -1;
    if ((Tp[0] > 0.0) && (Tp[1] > 0.0) && (Tp[2] > 0.0)) {
      index = Tp[0] < Tp[1] ? 0 : 1;
      if (Tp[2] < Tp[index]) index = 2;
    } else if ((Tp[0] < 0.0) && (Tp[1] < 0.0) && (Tp[2] < 0.0)) {
      index = Tp[0] > Tp[1] ? 0 : 1;
      if (Tp[2] > Tp[index]) index = 2;
    }

    if (index >= 0) {
      shown_disjoint = 1;
      ManifoldVec3 qIdx = q[index];

      ManifoldVec3 V = v3sub(qIdx, p[0]);
      ManifoldVec3 Z = v3cross(Sn, Sv[0]);
      if (v3dot(V, Z) > 0.0) {
        V = v3sub(qIdx, p[1]);
        Z = v3cross(Sn, Sv[1]);
        if (v3dot(V, Z) > 0.0) {
          V = v3sub(qIdx, p[2]);
          Z = v3cross(Sn, Sv[2]);
          if (v3dot(V, Z) > 0.0) {
            ManifoldVec3 cp = v3add(qIdx, v3scale(Sn, Tp[index] / Snl));
            ManifoldVec3 cq = qIdx;
            ManifoldVec3 diff = v3sub(cp, cq);
            return v3dot(diff, diff);
          }
        }
      }
    }
  }

  // Check if p vertex projects inside q triangle
  ManifoldVec3 Tn = v3cross(Tv[0], Tv[1]);
  double Tnl = v3dot(Tn, Tn);

  if (Tnl > 1e-15) {
    double Sp[3];
    Sp[0] = v3dot(v3sub(q[0], p[0]), Tn);
    Sp[1] = v3dot(v3sub(q[0], p[1]), Tn);
    Sp[2] = v3dot(v3sub(q[0], p[2]), Tn);

    int index = -1;
    if ((Sp[0] > 0.0) && (Sp[1] > 0.0) && (Sp[2] > 0.0)) {
      index = Sp[0] < Sp[1] ? 0 : 1;
      if (Sp[2] < Sp[index]) index = 2;
    } else if ((Sp[0] < 0.0) && (Sp[1] < 0.0) && (Sp[2] < 0.0)) {
      index = Sp[0] > Sp[1] ? 0 : 1;
      if (Sp[2] > Sp[index]) index = 2;
    }

    if (index >= 0) {
      shown_disjoint = 1;
      ManifoldVec3 pIdx = p[index];

      ManifoldVec3 V = v3sub(pIdx, q[0]);
      ManifoldVec3 Z = v3cross(Tn, Tv[0]);
      if (v3dot(V, Z) > 0.0) {
        V = v3sub(pIdx, q[1]);
        Z = v3cross(Tn, Tv[1]);
        if (v3dot(V, Z) > 0.0) {
          V = v3sub(pIdx, q[2]);
          Z = v3cross(Tn, Tv[2]);
          if (v3dot(V, Z) > 0.0) {
            ManifoldVec3 cp = pIdx;
            ManifoldVec3 cq = v3add(pIdx, v3scale(Tn, Sp[index] / Tnl));
            ManifoldVec3 diff = v3sub(cp, cq);
            return v3dot(diff, diff);
          }
        }
      }
    }
  }

  return shown_disjoint ? mindd : 0.0;
}
