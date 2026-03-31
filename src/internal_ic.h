#ifndef INTERNAL_IC_H
#define INTERNAL_IC_H

#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <limits>
#include <sstream>

#include "zmat/covalent_radii.h"

// Lightweight interpolation backends used by neb_interpolator / xyzgeom.
//
// LIIC implementation notes:
// - The interpolated target lives in primitive internal coordinates (bond / angle / torsion).
// - The primitive B matrix is computed via finite differences (fd_step).
// - A DLC-style basis is obtained by Jacobi diagonalization of G = Bp * Bp^T.
// - Back-transform uses: dx = Bp^T * (sum_i v_i * ((v_i·dq) / (lambda_i + damp))),
//   where (lambda_i, v_i) are selected eigenpairs of G.
// - Torsion differences and FD differences are wrapped to [-pi, pi] to avoid long-way interpolation.
//
// DM implementation notes:
// - The target distance matrix is linearly interpolated between endpoints.
// - Cartesian coordinates are then optimized against that target distance matrix.

namespace ICInterp {

// Constants and simple utilities
static inline double pi() { return 3.141592653589793238462643383279502884; }

// Configuration structures
struct LIICOptions {
    // Geometry / graph
    double bond_factor = 1.25;     // bond cutoff multiplier of (r_cov(i)+r_cov(j))
    // Finite difference
    double fd_step = 1.0e-4;       // Angstrom
    // DLC selection
    double ev_thresh = 1.0e-3;     // eigenvalue threshold for selecting DLC coordinates
    double damp = 1.0e-8;          // damping added to lambda in back-transform
    // Back-transform iteration
    int max_iter = 50;
    double tol = 1.0e-4;           // RMS of primitive residual (Ang/rad)
    double max_cart_step = 0.20;   // max atom displacement per iter (Ang)
    int line_search_steps = 6;     // reduce step if residual increases
    // Alignment per iteration
    bool kabsch_each_iter = true;
    bool kabsch_each_image = true;
    // Diagnostics
    int verbose = 0;
};


struct DMOptions {
    int max_iter = 800;
    double step = 5.0e-3;          // gradient descent step size
    double tol = 1.0e-3;           // RMS distance error (Ang)
    double max_cart_step = 0.20;   // max atom displacement per iter (Ang)
    int line_search_steps = 6;
    bool kabsch_each_iter = false; // usually not necessary for DM
    bool kabsch_each_image = true;
    int verbose = 0;
};

// Vector3 structure and simple operations
struct Vec3 {
    double x, y, z;
};

static inline Vec3 v_add(const Vec3& a, const Vec3& b) { return {a.x+b.x, a.y+b.y, a.z+b.z}; }
static inline Vec3 v_sub(const Vec3& a, const Vec3& b) { return {a.x-b.x, a.y-b.y, a.z-b.z}; }
static inline Vec3 v_mul(const Vec3& a, double s) { return {a.x*s, a.y*s, a.z*s}; }
static inline double v_dot(const Vec3& a, const Vec3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static inline Vec3 v_cross(const Vec3& a, const Vec3& b) {
    return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
}
static inline double v_norm(const Vec3& a) { return std::sqrt(v_dot(a,a)); }

// Simple coordinate access and manipulation
static inline Vec3 get_atom(const std::vector<double>& xyz, int i) {
    return {xyz[3*i+0], xyz[3*i+1], xyz[3*i+2]};
}
static inline void set_atom(std::vector<double>& xyz, int i, const Vec3& v) {
    xyz[3*i+0] = v.x; xyz[3*i+1] = v.y; xyz[3*i+2] = v.z;
}

static inline Vec3 centroid(const std::vector<double>& xyz) {
    const int n = static_cast<int>(xyz.size()/3);
    Vec3 c{0.0,0.0,0.0};
    if (n <= 0) return c;
    for (int i=0;i<n;++i) {
        c.x += xyz[3*i+0];
        c.y += xyz[3*i+1];
        c.z += xyz[3*i+2];
    }
    c.x /= n; c.y /= n; c.z /= n;
    return c;
}

static inline void translate(std::vector<double>& xyz, const Vec3& t) {
    const int n = static_cast<int>(xyz.size()/3);
    for (int i=0;i<n;++i) {
        xyz[3*i+0] += t.x;
        xyz[3*i+1] += t.y;
        xyz[3*i+2] += t.z;
    }
}

// Utility functions
static inline double wrap_to_pi(double a) {
    const double twopi = 2.0*pi();
    while (a >  pi()) a -= twopi;
    while (a <= -pi()) a += twopi;
    return a;
}

static inline double clamp01(double x) {
    if (x < -1.0) return -1.0;
    if (x >  1.0) return  1.0;
    return x;
}

// Primitive internal coordinate types
enum class PrimType { BOND=0, ANGLE=1, TORSION=2 };

struct PrimIC {
    PrimType type;
    int a, b, c, d; // use only needed; for bond: a,b; angle: a,b,c; torsion: a,b,c,d
};

// Primitive comparison and creation
static inline bool prim_less(const PrimIC& p1, const PrimIC& p2) {
    if (p1.type != p2.type) return static_cast<int>(p1.type) < static_cast<int>(p2.type);
    if (p1.a != p2.a) return p1.a < p2.a;
    if (p1.b != p2.b) return p1.b < p2.b;
    if (p1.c != p2.c) return p1.c < p2.c;
    return p1.d < p2.d;
}

static inline bool prim_equal(const PrimIC& p1, const PrimIC& p2) {
    return p1.type==p2.type && p1.a==p2.a && p1.b==p2.b && p1.c==p2.c && p1.d==p2.d;
}

static inline PrimIC make_bond(int i, int j) {
    if (i>j) std::swap(i,j);
    return {PrimType::BOND, i,j,-1,-1};
}

static inline PrimIC make_angle(int i, int j, int k) {
    if (i>k) std::swap(i,k);
    return {PrimType::ANGLE, i,j,k,-1};
}

static inline PrimIC make_torsion(int i, int j, int k, int l) {
    // canonicalize by comparing tuple with reversed tuple (l,k,j,i)
    int t1[4] = {i,j,k,l};
    int t2[4] = {l,k,j,i};
    bool use_rev = false;
    for (int m=0;m<4;++m) {
        if (t2[m] < t1[m]) { use_rev=true; break; }
        if (t2[m] > t1[m]) { use_rev=false; break; }
    }
    if (use_rev) { i=t2[0]; j=t2[1]; k=t2[2]; l=t2[3]; }
    return {PrimType::TORSION, i,j,k,l};
}

// Simple distance and internal coordinate computation
static inline double dist(const Vec3& a, const Vec3& b) {
    return v_norm(v_sub(a,b));
}

static inline double compute_bond(const std::vector<double>& xyz, int i, int j) {
    return dist(get_atom(xyz,i), get_atom(xyz,j));
}

static inline double compute_angle(const std::vector<double>& xyz, int i, int j, int k) {
    Vec3 vi = v_sub(get_atom(xyz,i), get_atom(xyz,j));
    Vec3 vk = v_sub(get_atom(xyz,k), get_atom(xyz,j));
    double ni = v_norm(vi);
    double nk = v_norm(vk);
    if (ni < 1e-15 || nk < 1e-15) return 0.0;
    double c = v_dot(vi,vk)/(ni*nk);
    return std::acos(clamp01(c));
}

static inline double compute_torsion(const std::vector<double>& xyz, int i, int j, int k, int l) {
    Vec3 r1 = get_atom(xyz,i);
    Vec3 r2 = get_atom(xyz,j);
    Vec3 r3 = get_atom(xyz,k);
    Vec3 r4 = get_atom(xyz,l);

    Vec3 b1 = v_sub(r2,r1);
    Vec3 b2 = v_sub(r3,r2);
    Vec3 b3 = v_sub(r4,r3);

    Vec3 n1 = v_cross(b1,b2);
    Vec3 n2 = v_cross(b2,b3);
    double n1n = v_norm(n1);
    double n2n = v_norm(n2);
    double b2n = v_norm(b2);
    if (n1n < 1e-15 || n2n < 1e-15 || b2n < 1e-15) return 0.0;

    Vec3 b2u = v_mul(b2, 1.0/b2n);
    Vec3 m1 = v_cross(n1, b2u);

    double x = v_dot(n1,n2);
    double y = v_dot(m1,n2);
    return std::atan2(y, x); // already in (-pi,pi]
}

static inline void compute_primitives(const std::vector<double>& xyz,
                                      const std::vector<PrimIC>& prims,
                                      std::vector<double>& q_out)
{
    q_out.assign(prims.size(), 0.0);
    for (size_t p=0;p<prims.size();++p) {
        const PrimIC& ic = prims[p];
        if (ic.type == PrimType::BOND) {
            q_out[p] = compute_bond(xyz, ic.a, ic.b);
        } else if (ic.type == PrimType::ANGLE) {
            q_out[p] = compute_angle(xyz, ic.a, ic.b, ic.c);
        } else {
            q_out[p] = compute_torsion(xyz, ic.a, ic.b, ic.c, ic.d);
        }
    }
}

// Forward declarations for complex functions (implemented in internal_ic.cpp)
bool kabsch_align(const std::vector<double>& reference, std::vector<double>& mobile);
void build_union_bonds(const std::vector<std::string>& symbols,
                      const std::vector<double>& x0,
                      const std::vector<double>& x1,
                      double bond_factor,
                      std::vector<PrimIC>& bonds_out,
                      std::vector<std::vector<int>>& nbrs_out);
void build_angles_torsions(const std::vector<std::vector<int>>& nbrs,
                           std::vector<PrimIC>& angles_out,
                           std::vector<PrimIC>& torsions_out);
void build_union_primitives(const std::vector<std::string>& symbols,
                           const std::vector<double>& x0,
                           const std::vector<double>& x1,
                           const LIICOptions& opt,
                           std::vector<PrimIC>& prims_out);
bool jacobi_eigen(std::vector<double> A, int n,
                  std::vector<double>& d,
                  std::vector<double>& V);
void compute_Bp_fd(const std::vector<double>& xyz,
                   const std::vector<PrimIC>& prims,
                   double h,
                   std::vector<double>& Bp_out);
double rms_residual(const std::vector<double>& q_target,
                    const std::vector<double>& q_cur,
                    const std::vector<PrimIC>& prims);
bool interpolate_liic(const std::vector<std::string>& symbols,
                     const std::vector<double>& x0,
                     const std::vector<double>& x1,
                     int nimages,
                     std::vector<std::vector<double>>& out_images,
                     const LIICOptions& opt,
                     std::string* err_msg);

bool interpolate_dm(const std::vector<std::string>& symbols,
                    const std::vector<double>& x0,
                    const std::vector<double>& x1,
                    int nimages,
                    std::vector<std::vector<double>>& out_images,
                    const DMOptions& opt,
                    std::string* err_msg);

} // namespace ICInterp

#endif // INTERNAL_IC_H
