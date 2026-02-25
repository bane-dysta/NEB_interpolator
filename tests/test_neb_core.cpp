#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "geometry.h"
#include "neb_core.h"

static bool almost_equal(double a, double b, double tol = 1e-10) {
    return std::fabs(a - b) <= tol;
}

static void require(bool cond, const std::string& msg) {
    if (!cond) {
        std::cerr << "[FAIL] " << msg << std::endl;
        std::exit(1);
    }
}

static geom::Structure make_structure_1atom(double x, double y, double z) {
    geom::Structure s;
    s.atoms.push_back(geom::Atom("H", x, y, z));
    return s;
}

static geom::Structure make_structure_2atoms(double x1, double y1, double z1,
                                             double x2, double y2, double z2) {
    geom::Structure s;
    s.atoms.push_back(geom::Atom("H", x1, y1, z1));
    s.atoms.push_back(geom::Atom("H", x2, y2, z2));
    return s;
}

static void test_spring_only_parallel() {
    // 1 atom, 1 image on x axis.
    // initial: x=0, image: x=1, final: x=3
    // tangent unit = +x
    // d_next=2, d_prev=1 => spring_scalar = k*(2-1) = k
    const double k = 1.5;
    geom::Structure initial = make_structure_1atom(0.0, 0.0, 0.0);
    geom::Structure final = make_structure_1atom(3.0, 0.0, 0.0);
    std::vector<geom::Structure> images = {make_structure_1atom(1.0, 0.0, 0.0)};

    // Real force is parallel to tangent => removed by perpendicular projection.
    std::vector<std::vector<double>> real_forces = {{5.0, 0.0, 0.0}};

    std::vector<std::vector<double>> neb_forces;
    double max_force = -1.0;
    std::string err;

    bool ok = neb::assembleNEBForces(initial, final, images, real_forces, k, neb_forces, &max_force, &err);
    require(ok, "assembleNEBForces returned false: " + err);
    require(neb_forces.size() == 1 && neb_forces[0].size() == 3, "unexpected output shape");

    require(almost_equal(neb_forces[0][0], k), "spring-only x component mismatch");
    require(almost_equal(neb_forces[0][1], 0.0), "spring-only y should be 0");
    require(almost_equal(neb_forces[0][2], 0.0), "spring-only z should be 0");
    require(almost_equal(max_force, std::fabs(k)), "max_force mismatch for spring-only case");
}

static void test_projection_perp() {
    // Same geometry as above, but with a perpendicular (y) real force.
    // Expected: y component preserved, plus spring along x.
    const double k = 1.5;
    geom::Structure initial = make_structure_1atom(0.0, 0.0, 0.0);
    geom::Structure final = make_structure_1atom(3.0, 0.0, 0.0);
    std::vector<geom::Structure> images = {make_structure_1atom(1.0, 0.0, 0.0)};

    std::vector<std::vector<double>> real_forces = {{0.0, 7.0, 0.0}};

    std::vector<std::vector<double>> neb_forces;
    double max_force = -1.0;
    std::string err;

    bool ok = neb::assembleNEBForces(initial, final, images, real_forces, k, neb_forces, &max_force, &err);
    require(ok, "assembleNEBForces returned false: " + err);
    require(neb_forces.size() == 1 && neb_forces[0].size() == 3, "unexpected output shape");

    require(almost_equal(neb_forces[0][0], k), "x component mismatch (spring)");
    require(almost_equal(neb_forces[0][1], 7.0), "y component mismatch (perp projection)");
    require(almost_equal(neb_forces[0][2], 0.0), "z component mismatch");

    const double expected_mag = std::sqrt(k * k + 7.0 * 7.0);
    require(almost_equal(max_force, expected_mag), "max_force mismatch for perp case");
}

static void test_global_projection_multiatom() {
    // 2 atoms, 1 image. Tangent in 3N space mixes atoms.
    // Make tangent direction: atom1 moves along x, atom2 moves along y.
    // Choose real force exactly parallel to tangent => perpendicular projection yields ~0.
    geom::Structure initial = make_structure_2atoms(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    geom::Structure final = make_structure_2atoms(2.0, 0.0, 0.0, 0.0, 2.0, 0.0);
    std::vector<geom::Structure> images = {make_structure_2atoms(1.0, 0.0, 0.0, 0.0, 1.0, 0.0)};

    // Parallel to tangent (scaled)
    std::vector<std::vector<double>> real_forces = {{10.0, 0.0, 0.0, 0.0, 10.0, 0.0}};

    std::vector<std::vector<double>> neb_forces;
    double max_force = -1.0;
    std::string err;

    // Here d_next == d_prev, so spring term is also zero.
    const double k = 1.0;
    bool ok = neb::assembleNEBForces(initial, final, images, real_forces, k, neb_forces, &max_force, &err);
    require(ok, "assembleNEBForces returned false: " + err);
    require(neb_forces.size() == 1 && neb_forces[0].size() == 6, "unexpected output shape");

    for (size_t i = 0; i < 6; ++i) {
        require(std::fabs(neb_forces[0][i]) <= 1e-9, "multi-atom projection should yield ~0");
    }
    require(std::fabs(max_force) <= 1e-9, "max_force should be ~0 for multi-atom projection test");
}

int main() {
    std::cout << "Running NEB core unit tests..." << std::endl;

    test_spring_only_parallel();
    std::cout << "  [OK] test_spring_only_parallel" << std::endl;

    test_projection_perp();
    std::cout << "  [OK] test_projection_perp" << std::endl;

    test_global_projection_multiatom();
    std::cout << "  [OK] test_global_projection_multiatom" << std::endl;

    std::cout << "All NEB core unit tests passed." << std::endl;
    return 0;
}
