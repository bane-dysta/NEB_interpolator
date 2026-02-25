#ifndef NEB_ENGINE_H
#define NEB_ENGINE_H

#include <string>
#include <vector>

#include "geometry.h"

namespace neb {

// External engine configuration.
//
// The external engine is expected to read an input file (XYZ blocks) and write an
// output file (gradient/force blocks). The exact text format is documented in
// neb_engine.cpp.
struct ExternalEngineConfig {
    bool enabled = false;
    std::string cmd;                         // command line
    std::string in_file = "neb_engine_in.dat";
    std::string out_file = "neb_engine_out.dat";
    std::string units = "Angstrom";          // label written to input/output headers
    bool keep_cycle_files = false;           // keep _cycle#### files
    bool output_is_force = false;            // default: output is gradient; if true, interpret as force
    double spring_constant = 1.0;            // NEB spring constant (only used by caller/NEB core)
    int run_every = 1;                       // run engine every N cycles
};

// Abstract engine interface: returns *real* forces (or virtual "physical" forces)
// for each intermediate image.
class NEBEngine {
public:
    virtual ~NEBEngine() = default;

    // Given the current intermediate images, return the *real* forces for each image.
    // forces[img] has length 3*natoms in the order [x1,y1,z1,x2,y2,z2,...].
    virtual bool compute_forces(const std::vector<geom::Structure>& images,
                                int cycle,
                                std::vector<std::vector<double>>& forces,
                                std::string* err) = 0;
};

// ExternalEngine: calls an external program each cycle (or every N cycles) and
// reads gradients/forces back from a text file.
class ExternalEngine : public NEBEngine {
public:
    explicit ExternalEngine(const ExternalEngineConfig& cfg);

    bool compute_forces(const std::vector<geom::Structure>& images,
                        int cycle,
                        std::vector<std::vector<double>>& forces,
                        std::string* err) override;

private:
    ExternalEngineConfig cfg_;
    bool have_cache_ = false;
    std::vector<std::vector<double>> last_forces_; // [img][3N]
};

// DistancePenaltyEngine: a dependency-free "virtual physics" engine that penalizes
// deviations from a linearly-interpolated distance matrix. Produces *real* forces
// only (no NEB spring, no projection).
class DistancePenaltyEngine : public NEBEngine {
public:
    DistancePenaltyEngine(const geom::Structure& initial, const geom::Structure& final);

    bool compute_forces(const std::vector<geom::Structure>& images,
                        int cycle,
                        std::vector<std::vector<double>>& forces,
                        std::string* err) override;

private:
    int natoms_ = 0;
    std::vector<double> d0_; // initial distance matrix (natoms*natoms)
    std::vector<double> d1_; // final distance matrix
};

} // namespace neb

#endif // NEB_ENGINE_H
