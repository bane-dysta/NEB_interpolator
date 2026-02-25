#include "neb_core.h"

#include <algorithm>
#include <cmath>

namespace neb {

bool assembleNEBForces(const geom::Structure& initial,
                       const geom::Structure& final,
                       const std::vector<geom::Structure>& images,
                       const std::vector<std::vector<double>>& real_forces,
                       double spring_constant,
                       std::vector<std::vector<double>>& neb_forces,
                       double* max_force,
                       std::string* err) {
    const int num_images = static_cast<int>(images.size());
    if (num_images == 0) {
        neb_forces.clear();
        if (max_force) *max_force = 0.0;
        return true;
    }

    const int natoms = static_cast<int>(initial.size());
    if (natoms <= 0) {
        if (err) *err = "NEBCore: natoms is zero";
        return false;
    }
    if (static_cast<int>(final.size()) != natoms) {
        if (err) *err = "NEBCore: initial/final natoms mismatch";
        return false;
    }
    if (static_cast<int>(real_forces.size()) != num_images) {
        if (err) *err = "NEBCore: real_forces size mismatch";
        return false;
    }

    neb_forces.assign(num_images, std::vector<double>(3 * static_cast<size_t>(natoms), 0.0));
    double local_max = 0.0;

    for (int img = 0; img < num_images; ++img) {
        const geom::Structure& current = images[img];
        const geom::Structure* prev = (img == 0) ? &initial : &images[img - 1];
        const geom::Structure* next = (img == num_images - 1) ? &final : &images[img + 1];

        const size_t n = current.size();
        if (static_cast<int>(n) != natoms) {
            if (err) *err = "NEBCore: natoms mismatch in intermediate image";
            return false;
        }
        if (real_forces[img].size() != 3 * n) {
            if (err) *err = "NEBCore: real force length mismatch";
            return false;
        }

        // Tangent vector (global): t = R_{i+1} - R_{i-1}
        std::vector<double> t(3 * n, 0.0);
        for (size_t i = 0; i < n; ++i) {
            t[3 * i + 0] = next->atoms[i].x - prev->atoms[i].x;
            t[3 * i + 1] = next->atoms[i].y - prev->atoms[i].y;
            t[3 * i + 2] = next->atoms[i].z - prev->atoms[i].z;
        }

        double tnorm2 = 0.0;
        for (double v : t) tnorm2 += v * v;
        const double tnorm = std::sqrt(tnorm2);

        std::vector<double> t_unit(3 * n, 0.0);
        if (tnorm > 1e-16) {
            for (size_t k = 0; k < t_unit.size(); ++k) t_unit[k] = t[k] / tnorm;
        }

        // Spring force along tangent (standard NEB spacing term)
        double d_next2 = 0.0, d_prev2 = 0.0;
        for (size_t i = 0; i < n; ++i) {
            const double dxn = next->atoms[i].x - current.atoms[i].x;
            const double dyn = next->atoms[i].y - current.atoms[i].y;
            const double dzn = next->atoms[i].z - current.atoms[i].z;
            d_next2 += dxn * dxn + dyn * dyn + dzn * dzn;

            const double dxp = current.atoms[i].x - prev->atoms[i].x;
            const double dyp = current.atoms[i].y - prev->atoms[i].y;
            const double dzp = current.atoms[i].z - prev->atoms[i].z;
            d_prev2 += dxp * dxp + dyp * dyp + dzp * dzp;
        }
        const double d_next = std::sqrt(d_next2);
        const double d_prev = std::sqrt(d_prev2);
        const double spring_scalar = (tnorm > 1e-16) ? (spring_constant * (d_next - d_prev)) : 0.0;

        // Project real force perpendicular to tangent
        const std::vector<double>& freal = real_forces[img];
        double dot_ft = 0.0;
        for (size_t k = 0; k < freal.size(); ++k) {
            dot_ft += freal[k] * t_unit[k];
        }

        for (size_t k = 0; k < freal.size(); ++k) {
            const double f_perp = freal[k] - dot_ft * t_unit[k];
            neb_forces[img][k] = f_perp + spring_scalar * t_unit[k];
        }

        for (size_t i = 0; i < n; ++i) {
            const double fx = neb_forces[img][3 * i + 0];
            const double fy = neb_forces[img][3 * i + 1];
            const double fz = neb_forces[img][3 * i + 2];
            const double force_mag = std::sqrt(fx * fx + fy * fy + fz * fz);
            local_max = std::max(local_max, force_mag);
        }
    }

    if (max_force) *max_force = local_max;
    return true;
}

} // namespace neb
