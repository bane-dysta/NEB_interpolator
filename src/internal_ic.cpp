#include "internal_ic.h"
#include <iostream>
#include <algorithm>
#include <numeric>
#include <iomanip>

namespace ICInterp {

// Quaternion-based Kabsch alignment: rotate 'mobile' onto 'reference'.
// Both vectors are xyz length 3N. Returns true if succeeded.
bool kabsch_align(const std::vector<double>& reference,
                  std::vector<double>& mobile)
{
    const int n = static_cast<int>(reference.size()/3);
    if (n <= 0 || mobile.size() != reference.size()) return false;

    Vec3 c_ref = centroid(reference);
    Vec3 c_mob = centroid(mobile);

    // Compute covariance S = sum (y_i)(x_i)^T, where x=ref centered, y=mob centered
    double Sxx=0,Sxy=0,Sxz=0,Syx=0,Syy=0,Syz=0,Szx=0,Szy=0,Szz=0;
    for (int i=0;i<n;++i) {
        Vec3 x = v_sub(get_atom(reference,i), c_ref);
        Vec3 y = v_sub(get_atom(mobile,i), c_mob);
        Sxx += y.x*x.x; Sxy += y.x*x.y; Sxz += y.x*x.z;
        Syx += y.y*x.x; Syy += y.y*x.y; Syz += y.y*x.z;
        Szx += y.z*x.x; Szy += y.z*x.y; Szz += y.z*x.z;
    }

    // Build 4x4 symmetric K matrix (Horn 1987)
    double K[16];
    K[0]  = Sxx + Syy + Szz;
    K[1]  = Syz - Szy;
    K[2]  = Szx - Sxz;
    K[3]  = Sxy - Syx;

    K[4]  = Syz - Szy;
    K[5]  = Sxx - Syy - Szz;
    K[6]  = Sxy + Syx;
    K[7]  = Szx + Sxz;

    K[8]  = Szx - Sxz;
    K[9]  = Sxy + Syx;
    K[10] = -Sxx + Syy - Szz;
    K[11] = Syz + Szy;

    K[12] = Sxy - Syx;
    K[13] = Szx + Sxz;
    K[14] = Syz + Szy;
    K[15] = -Sxx - Syy + Szz;

    // Power iteration to get the largest eigenvector (quaternion)
    double q[4] = {1,0,0,0};
    for (int it=0; it<60; ++it) {
        double nq[4] = {
            K[0]*q[0] + K[1]*q[1] + K[2]*q[2] + K[3]*q[3],
            K[4]*q[0] + K[5]*q[1] + K[6]*q[2] + K[7]*q[3],
            K[8]*q[0] + K[9]*q[1] + K[10]*q[2] + K[11]*q[3],
            K[12]*q[0] + K[13]*q[1] + K[14]*q[2] + K[15]*q[3]
        };
        double norm = std::sqrt(nq[0]*nq[0] + nq[1]*nq[1] + nq[2]*nq[2] + nq[3]*nq[3]);
        if (norm < 1e-18) break;
        for (int k=0;k<4;++k) nq[k] /= norm;

        double diff = std::fabs(nq[0]-q[0]) + std::fabs(nq[1]-q[1]) + std::fabs(nq[2]-q[2]) + std::fabs(nq[3]-q[3]);
        for (int k=0;k<4;++k) q[k] = nq[k];
        if (diff < 1e-12) break;
    }

    // Convert quaternion to rotation matrix
    const double w=q[0], x=q[1], y=q[2], z=q[3];
    const double R00 = 1.0 - 2.0*(y*y + z*z);
    const double R01 = 2.0*(x*y - w*z);
    const double R02 = 2.0*(x*z + w*y);

    const double R10 = 2.0*(x*y + w*z);
    const double R11 = 1.0 - 2.0*(x*x + z*z);
    const double R12 = 2.0*(y*z - w*x);

    const double R20 = 2.0*(x*z - w*y);
    const double R21 = 2.0*(y*z + w*x);
    const double R22 = 1.0 - 2.0*(x*x + y*y);

    // Apply: y_rot = R * (y - c_mob) + c_ref
    for (int i=0;i<n;++i) {
        Vec3 yy = v_sub(get_atom(mobile,i), c_mob);
        Vec3 yr{R00*yy.x + R01*yy.y + R02*yy.z,
                R10*yy.x + R11*yy.y + R12*yy.z,
                R20*yy.x + R21*yy.y + R22*yy.z};
        yr = v_add(yr, c_ref);
        set_atom(mobile, i, yr);
    }
    return true;
}

// Build a union bond graph using covalent radii and bond_factor, based on both endpoints.
void build_union_bonds(const std::vector<std::string>& symbols,
                       const std::vector<double>& x0,
                       const std::vector<double>& x1,
                       double bond_factor,
                       std::vector<PrimIC>& bonds_out,
                       std::vector<std::vector<int>>& nbrs_out)
{
    const int n = static_cast<int>(symbols.size());
    bonds_out.clear();
    nbrs_out.assign(n, std::vector<int>());

    std::vector<double> rad(n, 0.0);
    for (int i=0;i<n;++i) {
        rad[i] = ChemData::CovalentRadii::getRadius(symbols[i]);
        if (rad[i] <= 1e-8) rad[i] = 0.77; // fallback: ~C
    }

    auto add_bond = [&](int i, int j) {
        PrimIC b = make_bond(i,j);
        bonds_out.push_back(b);
    };

    for (int i=0;i<n;++i) {
        for (int j=i+1;j<n;++j) {
            const double cutoff = bond_factor * (rad[i] + rad[j]);
            const double d0 = compute_bond(x0, i, j);
            const double d1 = compute_bond(x1, i, j);
            if (d0 <= cutoff || d1 <= cutoff) {
                add_bond(i,j);
            }
        }
    }

    // Dedup bonds (should be unique already)
    std::sort(bonds_out.begin(), bonds_out.end(), prim_less);
    bonds_out.erase(std::unique(bonds_out.begin(), bonds_out.end(), prim_equal), bonds_out.end());

    // Build adjacency
    for (const auto& b : bonds_out) {
        int i=b.a, j=b.b;
        nbrs_out[i].push_back(j);
        nbrs_out[j].push_back(i);
    }
    for (int i=0;i<n;++i) {
        std::sort(nbrs_out[i].begin(), nbrs_out[i].end());
        nbrs_out[i].erase(std::unique(nbrs_out[i].begin(), nbrs_out[i].end()), nbrs_out[i].end());
    }
}

void build_angles_torsions(const std::vector<std::vector<int>>& nbrs,
                           std::vector<PrimIC>& angles_out,
                           std::vector<PrimIC>& torsions_out)
{
    const int n = static_cast<int>(nbrs.size());
    angles_out.clear();
    torsions_out.clear();

    // Angles: for each center j, choose i<k in neighbors
    for (int j=0;j<n;++j) {
        const auto& nn = nbrs[j];
        for (size_t a=0;a<nn.size();++a) {
            for (size_t b=a+1;b<nn.size();++b) {
                int i = nn[a];
                int k = nn[b];
                angles_out.push_back(make_angle(i,j,k));
            }
        }
    }

    // Torsions: for each bond j-k, for i in nbrs[j]\{k}, l in nbrs[k]\{j}
    // To avoid duplicates, iterate j<k for central bond.
    for (int j=0;j<n;++j) {
        for (int k : nbrs[j]) {
            if (j >= k) continue;
            for (int i : nbrs[j]) {
                if (i == k) continue;
                for (int l : nbrs[k]) {
                    if (l == j) continue;
                    torsions_out.push_back(make_torsion(i,j,k,l));
                }
            }
        }
    }

    std::sort(angles_out.begin(), angles_out.end(), prim_less);
    angles_out.erase(std::unique(angles_out.begin(), angles_out.end(), prim_equal), angles_out.end());

    std::sort(torsions_out.begin(), torsions_out.end(), prim_less);
    torsions_out.erase(std::unique(torsions_out.begin(), torsions_out.end(), prim_equal), torsions_out.end());
}

void build_union_primitives(const std::vector<std::string>& symbols,
                            const std::vector<double>& x0,
                            const std::vector<double>& x1,
                            const LIICOptions& opt,
                            std::vector<PrimIC>& prims_out)
{
    std::vector<PrimIC> bonds, angles, torsions;
    std::vector<std::vector<int>> nbrs;
    build_union_bonds(symbols, x0, x1, opt.bond_factor, bonds, nbrs);
    build_angles_torsions(nbrs, angles, torsions);

    prims_out.clear();
    prims_out.reserve(bonds.size()+angles.size()+torsions.size());
    prims_out.insert(prims_out.end(), bonds.begin(), bonds.end());
    prims_out.insert(prims_out.end(), angles.begin(), angles.end());
    prims_out.insert(prims_out.end(), torsions.begin(), torsions.end());

    std::sort(prims_out.begin(), prims_out.end(), prim_less);
    prims_out.erase(std::unique(prims_out.begin(), prims_out.end(), prim_equal), prims_out.end());
}

// Jacobi diagonalization for symmetric matrix A (n x n, row-major).
// Returns eigenvalues in d (size n) and eigenvectors in V (n x n, columns are eigenvectors).
bool jacobi_eigen(std::vector<double> A, int n,
                  std::vector<double>& d,
                  std::vector<double>& V)
{
    V.assign(n*n, 0.0);
    for (int i=0;i<n;++i) V[i*n+i] = 1.0;

    d.assign(n, 0.0);
    if (n <= 0) return false;

    const int max_sweeps = std::max(30, 6*n);
    for (int sweep=0; sweep<max_sweeps; ++sweep) {
        double sm = 0.0;
        for (int p=0;p<n-1;++p) {
            for (int q=p+1;q<n;++q) sm += std::fabs(A[p*n+q]);
        }
        if (sm < 1e-14) break;

        for (int p=0;p<n-1;++p) {
            for (int q=p+1;q<n;++q) {
                double apq = A[p*n+q];
                if (std::fabs(apq) < 1e-14) continue;

                double app = A[p*n+p];
                double aqq = A[q*n+q];

                double tau = (aqq - app) / (2.0*apq);
                double t;
                if (tau >= 0.0) t = 1.0 / (tau + std::sqrt(1.0 + tau*tau));
                else            t = -1.0 / (-tau + std::sqrt(1.0 + tau*tau));
                double c = 1.0 / std::sqrt(1.0 + t*t);
                double s = t * c;

                // Rotate A: A' = J^T A J
                // Update diagonal elements
                double app_new = app - t*apq;
                double aqq_new = aqq + t*apq;
                A[p*n+p] = app_new;
                A[q*n+q] = aqq_new;
                A[p*n+q] = 0.0;
                A[q*n+p] = 0.0;

                for (int r=0;r<n;++r) {
                    if (r==p || r==q) continue;
                    double arp = A[r*n+p];
                    double arq = A[r*n+q];
                    double new_rp = c*arp - s*arq;
                    double new_rq = s*arp + c*arq;
                    A[r*n+p] = new_rp;
                    A[p*n+r] = new_rp;
                    A[r*n+q] = new_rq;
                    A[q*n+r] = new_rq;
                }

                // Update eigenvectors
                for (int r=0;r<n;++r) {
                    double vrp = V[r*n+p];
                    double vrq = V[r*n+q];
                    V[r*n+p] = c*vrp - s*vrq;
                    V[r*n+q] = s*vrp + c*vrq;
                }
            }
        }
    }

    for (int i=0;i<n;++i) d[i] = A[i*n+i];

    // sanity: eigenvalues could be tiny negative due to numeric; clamp
    for (int i=0;i<n;++i) if (d[i] < 0.0 && std::fabs(d[i]) < 1e-12) d[i] = 0.0;

    return true;
}

// Compute primitive B matrix by central finite difference: Bp (m x 3N).
void compute_Bp_fd(const std::vector<double>& xyz,
                    const std::vector<PrimIC>& prims,
                    double h,
                    std::vector<double>& Bp_out)
{
    const int ncart = static_cast<int>(xyz.size());
    const int m = static_cast<int>(prims.size());
    Bp_out.assign(m*ncart, 0.0);

    std::vector<double> q_plus, q_minus;
    std::vector<double> x_plus = xyz;
    std::vector<double> x_minus = xyz;

    for (int k=0;k<ncart;++k) {
        x_plus[k]  = xyz[k] + h;
        x_minus[k] = xyz[k] - h;

        compute_primitives(x_plus, prims, q_plus);
        compute_primitives(x_minus, prims, q_minus);

        for (int p=0;p<m;++p) {
            double diff = q_plus[p] - q_minus[p];
            if (prims[p].type == PrimType::TORSION) diff = wrap_to_pi(diff);
            Bp_out[p*ncart + k] = diff / (2.0*h);
        }

        x_plus[k]  = xyz[k];
        x_minus[k] = xyz[k];
    }
}

// Helper: compute RMS of primitive residual, with torsion residual wrapped.
double rms_residual(const std::vector<double>& q_target,
                    const std::vector<double>& q_cur,
                    const std::vector<PrimIC>& prims)
{
    const int m = static_cast<int>(prims.size());
    if ((int)q_target.size()!=m || (int)q_cur.size()!=m || m==0) return 0.0;
    double ss = 0.0;
    for (int p=0;p<m;++p) {
        double r = q_target[p] - q_cur[p];
        if (prims[p].type == PrimType::TORSION) r = wrap_to_pi(r);
        ss += r*r;
    }
    return std::sqrt(ss / m);
}

// Main LIIC interpolation: linear targets in primitive internal coordinates plus DLC-style back-transform.
// Returns false if any image fails; callers may apply fallback policy outside this function.
bool interpolate_liic(const std::vector<std::string>& symbols,
                     const std::vector<double>& x0,
                     const std::vector<double>& x1,
                     int nimages,
                     std::vector<std::vector<double>>& out_images,
                     const LIICOptions& opt,
                     std::string* err_msg)
{
    out_images.clear();
    if (nimages <= 0) return true;
    const int n = static_cast<int>(symbols.size());
    if ((int)x0.size() != 3*n || (int)x1.size() != 3*n) {
        if (err_msg) *err_msg = "LIIC: xyz size mismatch.";
        return false;
    }

    // Build union primitive list
    std::vector<PrimIC> prims;
    build_union_primitives(symbols, x0, x1, opt, prims);
    const int m = static_cast<int>(prims.size());
    if (m <= 0) {
        if (err_msg) *err_msg = "LIIC: primitive list empty (bond_factor too small?).";
        return false;
    }

    std::vector<double> q0, q1;
    compute_primitives(x0, prims, q0);
    compute_primitives(x1, prims, q1);

    // Build wrapped difference dq (torsions shortest path)
    std::vector<double> dq(m,0.0);
    for (int p=0;p<m;++p) {
        double d = q1[p] - q0[p];
        if (prims[p].type == PrimType::TORSION) d = wrap_to_pi(d);
        dq[p] = d;
    }

    // Sequential growth: each image starts from previous converged geometry
    std::vector<double> x_prev = x0;
    out_images.resize(nimages);

    for (int img=0; img<nimages; ++img) {
        const double f = static_cast<double>(img+1) / static_cast<double>(nimages+1);

        std::vector<double> q_target(m,0.0);
        for (int p=0;p<m;++p) q_target[p] = q0[p] + f*dq[p];

        std::vector<double> x = x_prev; // start from previous
        if (img==0) {
            // a small nudge towards product can help if starting residual is large
            for (int k=0;k<3*n;++k) x[k] = x0[k] + f*(x1[k]-x0[k]);
        }

        // Optional: align guess to previous to keep orientation consistent
        if (opt.kabsch_each_image && img>0) {
            std::vector<double> ref = x_prev;
            kabsch_align(ref, x);
        }

        std::vector<double> q_cur;

        double last_rms = std::numeric_limits<double>::infinity();
        int it_converged = -1;

        if (opt.verbose >= 1) {
            std::cerr << "    [LIIC] image " << (img + 1) << "/" << nimages
                      << "  f=" << f
                      << "  nprim=" << m
                      << "  tol=" << opt.tol
                      << "  max_iter=" << opt.max_iter
                      << std::endl;
        }

        bool converged = false;
        for (int it=0; it<opt.max_iter; ++it) {
            compute_primitives(x, prims, q_cur);
            double rms = rms_residual(q_target, q_cur, prims);
            last_rms = rms;

            // Verbose levels:
            //   0: silent
            //   1: print every ~5 iterations (plus first/last/converged)
            //   2: print every iteration
            if (opt.verbose >= 2 ||
                (opt.verbose >= 1 && (it == 0 || ((it + 1) % 5 == 0) || rms < opt.tol))) {
                std::cerr << "      iter " << it
                          << "  rms=" << std::scientific << rms << std::defaultfloat
                          << std::endl;
            }

            if (rms < opt.tol) { converged = true; it_converged = it; break; }

            // Finite difference Bp
            std::vector<double> Bp;
            compute_Bp_fd(x, prims, opt.fd_step, Bp);

            // G = Bp * Bp^T (m x m)
            std::vector<double> G(m*m, 0.0);
            const int ncart = 3*n;
            for (int i=0;i<m;++i) {
                for (int j=i;j<m;++j) {
                    double s = 0.0;
                    const double* Bi = &Bp[i*ncart];
                    const double* Bj = &Bp[j*ncart];
                    for (int k=0;k<ncart;++k) s += Bi[k]*Bj[k];
                    G[i*m+j] = s;
                    G[j*m+i] = s;
                }
            }

            // Diagonalize G
            std::vector<double> evals, evecs;
            if (!jacobi_eigen(G, m, evals, evecs)) {
                if (err_msg) *err_msg = "LIIC: Jacobi eigen failed.";
                return false;
            }

            // Sort eigenpairs by eval descending
            std::vector<int> idx(m);
            for (int i=0;i<m;++i) idx[i]=i;
            std::sort(idx.begin(), idx.end(), [&](int a, int b){ return evals[a] > evals[b]; });

            // Build residual r_p (with torsion wrap)
            std::vector<double> r_p(m,0.0);
            for (int p=0;p<m;++p) {
                double r = q_target[p] - q_cur[p];
                if (prims[p].type == PrimType::TORSION) r = wrap_to_pi(r);
                r_p[p] = r;
            }

            // Compute primitive-space weight vector t = sum_sel v_i * ( (v_i·r_p) / (lambda_i + damp) )
            std::vector<double> t(m, 0.0);
            int nsel = 0;
            for (int kk=0; kk<m; ++kk) {
                const int ev_i = idx[kk];
                const double lam = evals[ev_i];
                if (lam <= opt.ev_thresh) continue;

                // eigenvector v is column ev_i of evecs (size m)
                double dotvr = 0.0;
                for (int j=0;j<m;++j) dotvr += evecs[j*m + ev_i] * r_p[j];
                const double scale = dotvr / (lam + opt.damp);

                for (int j=0;j<m;++j) t[j] += evecs[j*m + ev_i] * scale;
                ++nsel;
            }
            if (nsel == 0) {
                if (err_msg) {
                    double max_eval = 0.0;
                    for (double v : evals) max_eval = std::max(max_eval, v);
                    std::ostringstream oss;
                    oss << "LIIC: no DLC eigenvectors selected (ev_thresh=" << opt.ev_thresh
                        << ", max_eval=" << max_eval << ").";
                    *err_msg = oss.str();
                }
                return false;
            }

            // dx = Bp^T * t
            std::vector<double> dx(ncart, 0.0);
            for (int k=0;k<ncart;++k) {
                double s = 0.0;
                for (int p=0;p<m;++p) s += Bp[p*ncart + k] * t[p];
                dx[k] = s;
            }

            // Limit step size
            double max_disp = 0.0;
            for (int iatom=0;iatom<n;++iatom) {
                Vec3 d{dx[3*iatom+0], dx[3*iatom+1], dx[3*iatom+2]};
                max_disp = std::max(max_disp, v_norm(d));
            }
            double scale_step = 1.0;
            if (max_disp > opt.max_cart_step && max_disp > 1e-18) {
                scale_step = opt.max_cart_step / max_disp;
            }

            if (opt.verbose >= 2) {
                std::cerr << "        nsel=" << nsel
                          << "  max_disp=" << max_disp
                          << "  step_scale=" << scale_step
                          << std::endl;
            }

            // Simple line search: ensure residual decreases
            bool accepted = false;
            std::vector<double> x_trial = x;
            for (int ls=0; ls<=opt.line_search_steps; ++ls) {
                for (int k=0;k<ncart;++k) x_trial[k] = x[k] + scale_step*dx[k];

                if (opt.kabsch_each_iter) {
                    // Keep translation/rotation stable relative to x_prev
                    std::vector<double> ref = (img==0 ? x0 : x_prev);
                    kabsch_align(ref, x_trial);
                }

                std::vector<double> q_trial;
                compute_primitives(x_trial, prims, q_trial);
                double rms_trial = rms_residual(q_target, q_trial, prims);

                // Accept only non-worsening steps (within a tiny tolerance).
                // The previous implementation also compared against a global
                // rms0 initialized to 1e100, which could accept a worse first
                // step unconditionally.
                if (rms_trial <= rms*(1.0 + 1e-6)) {
                    x.swap(x_trial);
                    accepted = true;
                    break;
                }
                scale_step *= 0.5;
            }

            if (!accepted) {
                if (err_msg) {
                    std::ostringstream oss;
                    oss << "LIIC: back-transform diverged (image " << (img+1)
                        << ", iter " << it
                        << ", last_rms=" << last_rms
                        << ").";
                    *err_msg = oss.str();
                }
                return false;
            }
        }

        if (!converged) {
            if (err_msg) {
                std::ostringstream oss;
                oss << "LIIC: did not converge within max_iter (image " << (img+1)
                    << ", last_rms=" << last_rms
                    << ", tol=" << opt.tol
                    << ", max_iter=" << opt.max_iter
                    << ").";
                *err_msg = oss.str();
            }
            return false;
        }

        if (opt.verbose >= 1) {
            std::cerr << "    [LIIC] image " << (img + 1)
                      << " converged in " << (it_converged + 1) << " iterations"
                      << "  final_rms=" << std::scientific << last_rms << std::defaultfloat
                      << std::endl;
        }

        out_images[img] = x;
        x_prev = x;
    }

    return true;
}

// Distance-matrix (DM) interpolation: optimize xyz to match interpolated pairwise distances.
// Useful as a fallback when LIIC back-transform fails.
bool interpolate_dm(const std::vector<std::string>& symbols,
                    const std::vector<double>& x0,
                    const std::vector<double>& x1,
                    int nimages,
                    std::vector<std::vector<double>>& out_images,
                    const DMOptions& opt,
                    std::string* err_msg)
{
    (void)symbols; // currently unused (could be used for weighting)
    out_images.clear();
    if (nimages <= 0) return true;

    const int n = static_cast<int>(x0.size()/3);
    if ((int)x1.size() != 3*n) {
        if (err_msg) *err_msg = "DM: xyz size mismatch.";
        return false;
    }

    // Precompute distance matrices
    std::vector<double> D0(n*n, 0.0), D1(n*n, 0.0);
    auto fill_dist = [&](const std::vector<double>& x, std::vector<double>& D) {
        for (int i=0;i<n;++i) {
            Vec3 ri = get_atom(x,i);
            for (int j=0;j<n;++j) {
                if (i==j) { D[i*n+j] = 0.0; continue; }
                Vec3 rj = get_atom(x,j);
                D[i*n+j] = dist(ri,rj);
            }
        }
    };
    fill_dist(x0, D0);
    fill_dist(x1, D1);

    std::vector<double> x_prev = x0;
    out_images.resize(nimages);

    for (int img=0; img<nimages; ++img) {
        const double f = static_cast<double>(img+1) / static_cast<double>(nimages+1);

        // Target distances
        std::vector<double> Dt(n*n, 0.0);
        for (int i=0;i<n*n;++i) Dt[i] = D0[i] + f*(D1[i] - D0[i]);

        // Initial guess: LIC
        std::vector<double> x(3*n, 0.0);
        for (int k=0;k<3*n;++k) x[k] = x0[k] + f*(x1[k] - x0[k]);

        if (opt.kabsch_each_image && img>0) {
            std::vector<double> ref = x_prev;
            kabsch_align(ref, x);
        }

        auto rms_dist_err = [&](const std::vector<double>& xcur) -> double {
            double ss = 0.0;
            long long cnt = 0;
            for (int i=0;i<n;++i) {
                Vec3 ri = get_atom(xcur,i);
                for (int j=i+1;j<n;++j) {
                    Vec3 rj = get_atom(xcur,j);
                    double dij = dist(ri,rj);
                    double r = dij - Dt[i*n+j];
                    ss += r*r;
                    ++cnt;
                }
            }
            return (cnt>0) ? std::sqrt(ss / static_cast<double>(cnt)) : 0.0;
        };

        double rms0 = rms_dist_err(x);
        std::vector<double> ref_geom = (img==0 ? x0 : x_prev);
        Vec3 c_ref = centroid(ref_geom);
        bool converged = (rms0 < opt.tol);

        for (int it=0; it<opt.max_iter && !converged; ++it) {
            std::vector<double> grad(3*n, 0.0);

            // Gradient of sum_{i<j} (d_ij - t_ij)^2
            for (int i=0;i<n;++i) {
                Vec3 ri = get_atom(x,i);
                for (int j=i+1;j<n;++j) {
                    Vec3 rj = get_atom(x,j);
                    Vec3 rij = v_sub(ri,rj);
                    double dij = v_norm(rij);
                    if (dij < 1e-12) continue;
                    double diff = dij - Dt[i*n+j];
                    Vec3 u = v_mul(rij, 1.0/dij);
                    // Derivative: 2*diff*u on i, -2*diff*u on j
                    Vec3 gi = v_mul(u, 2.0*diff);
                    grad[3*i+0] += gi.x; grad[3*i+1] += gi.y; grad[3*i+2] += gi.z;
                    grad[3*j+0] -= gi.x; grad[3*j+1] -= gi.y; grad[3*j+2] -= gi.z;
                }
            }

            // Step limiting
            double max_disp = 0.0;
            for (int iatom=0;iatom<n;++iatom) {
                Vec3 g{grad[3*iatom+0], grad[3*iatom+1], grad[3*iatom+2]};
                max_disp = std::max(max_disp, v_norm(g));
            }
            double step = opt.step;
            if (max_disp > 1e-18) {
                double trial_disp = step * max_disp;
                if (trial_disp > opt.max_cart_step) step *= (opt.max_cart_step / trial_disp);
            }

            bool accepted = false;
            std::vector<double> x_trial = x;
            for (int ls=0; ls<=opt.line_search_steps; ++ls) {
                for (int k=0;k<3*n;++k) x_trial[k] = x[k] - step * grad[k];

                // Keep COM stable: match reference centroid (translation is a null mode for distances)
                Vec3 c = centroid(x_trial);
                translate(x_trial, v_sub(c_ref, c));

                if (opt.kabsch_each_iter) {
                    kabsch_align(ref_geom, x_trial);
                }

                double rms1 = rms_dist_err(x_trial);
                if (rms1 <= rms0*(1.0 + 1e-6) || rms1 < rms0) {
                    x.swap(x_trial);
                    rms0 = rms1;
                    accepted = true;
                    break;
                }
                step *= 0.5;
            }

            if (!accepted) break;

            if (rms0 < opt.tol) converged = true;
        }

        if (!converged) {
            if (err_msg) {
                std::ostringstream oss;
                oss << "DM: did not converge (image " << (img+1) << "), final RMS(dist) = " << rms0;
                *err_msg = oss.str();
            }
            return false;
        }

        if (opt.kabsch_each_image) {
            kabsch_align(ref_geom, x);
        } else {
            Vec3 ccur = centroid(x);
            translate(x, v_sub(c_ref, ccur));
        }

        out_images[img] = x;
        x_prev = x;
    }

    return true;
}

} // namespace ICInterp

