#ifndef NEB_DRIVER_H
#define NEB_DRIVER_H

#include <string>
#include <vector>

#include "geometry.h"
#include "internal_ic.h"
#include "neb_core.h"
#include "neb_engine.h"
#include "rmsd_align.h"

namespace neb {

// Public method names (canonical terminology)
enum class Method {
    LIC,       // Linear interpolation in Cartesian coordinates
    LIIC,      // Linear interpolation in internal coordinates (DLC-style)
    DM,        // Distance-matrix interpolation
    NEB,       // Nudged elastic band optimization (LIC initialization)
    NEB_LIIC   // NEB optimization (LIIC initialization)
};

class NEBDriver {
public:
    NEBDriver(int num_img = 5,
              double step = 0.0001,
              double conv_thresh = 0.01,
              int max_iter = 10000,
              bool align = true,
              std::string rmsd_exec = "./calc_rmsd_xyz");

    // Optional external-engine mode (used by NEB / NEB-LIIC)
    void setExternalEngine(const ExternalEngineConfig& cfg);

    // LIIC/DM configuration (dependency-free)
    void setLIICOptions(const ICInterp::LIICOptions& opt);
    void setDMOptions(const ICInterp::DMOptions& opt);
    void setDMFallback(bool enable);

    // Alignment (Fortran calc_rmsd_xyz)
    void setAlignment(bool enable);
    void setRMSDExecutable(const std::string& path);

    // Set endpoint structures
    bool setStructuresFromFiles(const std::string& initial_file,
                                const std::string& final_file,
                                std::string* err = nullptr);

    // For interactive tools (e.g., xyzgeom): set initial from memory then final from file.
    void setInitialFromMemory(const std::vector<geom::Atom>& atoms_in,
                              const std::string& comment = "Initial structure from memory");
    bool setFinalFromFile(const std::string& final_file, std::string* err = nullptr);

    // Run one of the supported methods.
    // Returns false only for *fatal* errors that prevent generating any path.
    bool run(Method method, std::string* err = nullptr);

    // Write path (initial + intermediates + final).
    //  - separate: <prefix>00.xyz, <prefix>01.xyz, ..., <prefix>NN.xyz
    //  - multiframe: <prefix>trajectory.xyz
    bool writeResults(const std::string& prefix = "", bool multiframe = false) const;

    // Accessors
    const geom::Structure& initial() const { return initial_; }
    const geom::Structure& final() const { return final_; }
    const std::vector<geom::Structure>& images() const { return images_; }

private:
    geom::Structure initial_;
    geom::Structure final_;
    std::vector<geom::Structure> images_; // intermediate images only

    int num_images_;
    double step_init_;
    double convergence_threshold_;
    int max_iterations_;

    bool use_alignment_;
    rmsd::FortranRMSDAligner aligner_;

    ExternalEngineConfig engine_cfg_;
    bool use_external_engine_ = false;

    ICInterp::LIICOptions liic_options_;
    ICInterp::DMOptions dm_options_;
    bool enable_dm_fallback_ = true;

    // Internal helpers
    void performLIC_();

    // Returns true if LIIC succeeded (or DM fallback succeeded).
    // Returns false if it had to fall back to Cartesian LIC.
    // Throws/returns fatal error via err only for invalid input (e.g. unknown element).
    bool performLIIC_(std::string* err);

    // Returns true if DM succeeded; false if it had to fall back to Cartesian LIC.
    bool performDM_(std::string* err);

    bool performNEB_(bool init_with_liic, std::string* err);

    bool writeSeparateXYZ_(const std::string& prefix) const;
    bool writeMultiframeXYZ_(const std::string& prefix) const;
};

} // namespace neb

#endif // NEB_DRIVER_H
