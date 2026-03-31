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

// User-facing run modes.
enum class Method {
    LIC,       // Linear interpolation in Cartesian coordinates
    LIIC,      // Linear interpolation in primitive internal coordinates (DLC-style back-transform)
    DM,        // Distance-matrix interpolation
    NEB,       // Nudged elastic band optimization (LIC initialization)
    NEB_LIIC   // Nudged elastic band optimization (LIIC initialization)
};

// Actual interpolation backend used to generate the path initializer.
enum class PathMethod {
    NONE,
    LIC,
    LIIC,
    DM
};

enum class RunStatus {
    NOT_RUN,
    GENERATED_EXACT,
    GENERATED_WITH_FALLBACK,
    FATAL_ERROR
};

std::string toString(Method method);
std::string toString(PathMethod method);
std::string toString(RunStatus status);
PathMethod requestedPathMethod(Method method);
bool modeUsesNEB(Method method);
std::string formatPathMethodChain(const std::vector<PathMethod>& methods);

struct RunReport {
    Method requested_mode = Method::LIC;
    PathMethod requested_path_method = PathMethod::LIC;
    PathMethod actual_path_method = PathMethod::NONE;
    RunStatus status = RunStatus::NOT_RUN;
    bool path_generated = false;
    bool used_fallback = false;
    bool optimized_with_neb = false;
    std::vector<PathMethod> attempted_path_methods;
    std::string detail;
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

    // LIIC / DM configuration (dependency-free)
    void setLIICOptions(const ICInterp::LIICOptions& opt);
    void setDMOptions(const ICInterp::DMOptions& opt);
    void setLIICToDMFallback(bool enable);

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

    // Run one of the supported modes.
    // Return value means: did the requested mode produce a usable path / optimized path?
    // Whether the requested interpolation backend was actually used is reported via lastRunReport().
    bool run(Method method, std::string* err = nullptr);

    // Write path (initial + intermediates + final).
    //  - separate: <prefix>00.xyz, <prefix>01.xyz, ..., <prefix>NN.xyz
    //  - multiframe: <prefix>trajectory.xyz
    bool writeResults(const std::string& prefix = "", bool multiframe = false) const;

    // Accessors
    const geom::Structure& initial() const { return initial_; }
    const geom::Structure& final() const { return final_; }
    const std::vector<geom::Structure>& images() const { return images_; }
    const RunReport& lastRunReport() const { return last_run_report_; }
    bool lastRunUsedRequestedPathMethod() const {
        return last_run_report_.path_generated && !last_run_report_.used_fallback;
    }

private:
    struct PathBuildResult {
        bool path_generated = false;
        bool fatal = false;
        PathMethod actual_method = PathMethod::NONE;
        std::vector<PathMethod> attempted_methods;
        std::string detail;
    };

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
    bool enable_liic_to_dm_fallback_ = true;
    RunReport last_run_report_;

    // Internal helpers
    void collectEndpointArrays_(std::vector<std::string>& symbols,
                                std::vector<double>& x0,
                                std::vector<double>& x1) const;
    void loadImagesFromCartesian_(const std::vector<std::string>& symbols,
                                  const std::vector<std::vector<double>>& imgs_xyz,
                                  PathMethod actual_method,
                                  PathMethod requested_method);

    PathBuildResult generateLICPath_(PathMethod requested_method = PathMethod::LIC,
                                     const std::string& intro_message = "");
    PathBuildResult generateLIICPath_(std::string* err);
    PathBuildResult generateDMPath_(std::string* err);
    bool optimizeNEB_(std::string* err);
    void updateNEBImageComments_(Method requested_mode, const PathBuildResult& init_result);
    void resetRunReport_(Method method);
    void finalizeRunReport_(Method method,
                            const PathBuildResult& result,
                            bool optimized_with_neb,
                            const std::string& detail_override = "");

    bool writeSeparateXYZ_(const std::string& prefix) const;
    bool writeMultiframeXYZ_(const std::string& prefix) const;
};

} // namespace neb

#endif // NEB_DRIVER_H
