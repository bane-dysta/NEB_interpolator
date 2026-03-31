#include "neb_driver.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>

#include "util.h"
#include "zmat/covalent_radii.h"

namespace neb {
namespace {

std::string makeDirectImageComment(PathMethod actual_method,
                                   PathMethod requested_method,
                                   int image_index)
{
    std::ostringstream oss;
    oss << toString(actual_method) << " intermediate image " << image_index;
    if (requested_method != PathMethod::NONE && requested_method != actual_method) {
        oss << " [fallback from " << toString(requested_method) << "]";
    }
    return oss.str();
}

std::string makeNEBImageComment(Method requested_mode,
                                PathMethod actual_method,
                                PathMethod requested_path_method,
                                int image_index)
{
    std::ostringstream oss;
    oss << toString(requested_mode) << " optimized image " << image_index
        << " [initialized with " << toString(actual_method);
    if (actual_method != PathMethod::NONE && actual_method != requested_path_method) {
        oss << ", fallback from " << toString(requested_path_method);
    }
    oss << "]";
    return oss.str();
}

std::string decorateFailureReason(const std::string& prefix, const std::string& reason) {
    if (reason.empty()) return prefix;
    return prefix + " (" + reason + ")";
}

} // namespace

std::string toString(Method method) {
    switch (method) {
        case Method::LIC:      return "LIC";
        case Method::LIIC:     return "LIIC";
        case Method::DM:       return "DM";
        case Method::NEB:      return "NEB";
        case Method::NEB_LIIC: return "NEB-LIIC";
        default:               return "UNKNOWN";
    }
}

std::string toString(PathMethod method) {
    switch (method) {
        case PathMethod::NONE: return "NONE";
        case PathMethod::LIC:  return "LIC";
        case PathMethod::LIIC: return "LIIC";
        case PathMethod::DM:   return "DM";
        default:               return "UNKNOWN";
    }
}

std::string toString(RunStatus status) {
    switch (status) {
        case RunStatus::NOT_RUN:                 return "NOT_RUN";
        case RunStatus::GENERATED_EXACT:         return "GENERATED_EXACT";
        case RunStatus::GENERATED_WITH_FALLBACK: return "GENERATED_WITH_FALLBACK";
        case RunStatus::FATAL_ERROR:             return "FATAL_ERROR";
        default:                                 return "UNKNOWN";
    }
}

PathMethod requestedPathMethod(Method method) {
    switch (method) {
        case Method::LIC:
        case Method::NEB:
            return PathMethod::LIC;
        case Method::LIIC:
        case Method::NEB_LIIC:
            return PathMethod::LIIC;
        case Method::DM:
            return PathMethod::DM;
        default:
            return PathMethod::NONE;
    }
}

bool modeUsesNEB(Method method) {
    return method == Method::NEB || method == Method::NEB_LIIC;
}

std::string formatPathMethodChain(const std::vector<PathMethod>& methods) {
    if (methods.empty()) return "(none)";
    std::ostringstream oss;
    for (size_t i = 0; i < methods.size(); ++i) {
        if (i > 0) oss << " -> ";
        oss << toString(methods[i]);
    }
    return oss.str();
}

NEBDriver::NEBDriver(int num_img,
                     double step,
                     double conv_thresh,
                     int max_iter,
                     bool align,
                     std::string rmsd_exec)
    : num_images_(num_img),
      step_init_(step),
      convergence_threshold_(conv_thresh),
      max_iterations_(max_iter),
      use_alignment_(align),
      aligner_(std::move(rmsd_exec))
{
}

void NEBDriver::setExternalEngine(const ExternalEngineConfig& cfg) {
    engine_cfg_ = cfg;
    use_external_engine_ = cfg.enabled;
}

void NEBDriver::setLIICOptions(const ICInterp::LIICOptions& opt) {
    liic_options_ = opt;
}

void NEBDriver::setDMOptions(const ICInterp::DMOptions& opt) {
    dm_options_ = opt;
}

void NEBDriver::setLIICToDMFallback(bool enable) {
    enable_liic_to_dm_fallback_ = enable;
}

void NEBDriver::setAlignment(bool enable) {
    use_alignment_ = enable;
}

void NEBDriver::setRMSDExecutable(const std::string& path) {
    aligner_.setExecutablePath(path);
}

bool NEBDriver::setStructuresFromFiles(const std::string& initial_file,
                                       const std::string& final_file,
                                       std::string* err)
{
    if (!initial_.readXYZ(initial_file)) {
        if (err) *err = "Error: failed to read initial structure: " + initial_file;
        return false;
    }
    if (!final_.readXYZ(final_file)) {
        if (err) *err = "Error: failed to read final structure: " + final_file;
        return false;
    }

    if (!initial_.isCompatible(final_)) {
        if (err) *err = "Error: initial and final structures are not compatible (atom symbols/order mismatch).";
        return false;
    }

    if (use_alignment_) {
        std::cout << "Aligning structures using Fortran RMSD..." << std::endl;
        const double rmsd_before = initial_.calculateRMSD(final_);
        std::cout << "RMSD before alignment: " << std::fixed << std::setprecision(6) << rmsd_before << std::endl;

        // Create temporary file for alignment
        const std::string temp_final = util::makeTempFilePath("neb_final_", ".xyz").string();
        if (!final_.writeXYZ(temp_final)) {
            if (err) *err = "Error: cannot create temporary file for alignment.";
            return false;
        }

        if (aligner_.alignAndReplace(initial_file, temp_final)) {
            geom::Structure final_aligned;
            if (final_aligned.readXYZ(temp_final)) {
                const double rmsd_after = initial_.calculateRMSD(final_aligned);
                std::cout << "RMSD after alignment: " << std::fixed << std::setprecision(6) << rmsd_after << std::endl;

                if (rmsd_after < rmsd_before) {
                    final_ = final_aligned;
                    std::cout << "Alignment improved RMSD, using aligned final structure." << std::endl;
                } else {
                    std::cout << "Alignment did not improve RMSD, using original final structure." << std::endl;
                }
            } else {
                std::cerr << "Warning: cannot read aligned structure, using original." << std::endl;
            }
        } else {
            std::cerr << "Warning: alignment failed, using original structures." << std::endl;
        }

        // Cleanup temp file (non-fatal)
        std::error_code ec;
        std::filesystem::remove(temp_final, ec);
        (void)ec;
    }

    return true;
}

void NEBDriver::setInitialFromMemory(const std::vector<geom::Atom>& atoms_in,
                                     const std::string& comment)
{
    initial_.atoms = atoms_in;
    initial_.comment = comment;
}

bool NEBDriver::setFinalFromFile(const std::string& final_file, std::string* err) {
    if (initial_.atoms.empty()) {
        if (err) *err = "Error: initial structure is not set.";
        return false;
    }

    if (!final_.readXYZ(final_file)) {
        if (err) *err = "Error: failed to read final structure: " + final_file;
        return false;
    }

    if (!initial_.isCompatible(final_)) {
        if (err) *err = "Error: initial and final structures are not compatible (atom symbols/order mismatch).";
        return false;
    }

    // Keep API semantics consistent with setStructuresFromFiles(): if alignment is enabled,
    // we also align when the initial structure was provided from memory.
    if (use_alignment_) {
        std::cout << "Aligning structures using Fortran RMSD..." << std::endl;
        const double rmsd_before = initial_.calculateRMSD(final_);
        std::cout << "RMSD before alignment: " << std::fixed << std::setprecision(6) << rmsd_before << std::endl;

        const std::string temp_initial = util::makeTempFilePath("neb_initial_", ".xyz").string();
        const std::string temp_final = util::makeTempFilePath("neb_final_", ".xyz").string();

        std::error_code ec;
        auto cleanup = [&]() {
            ec.clear();
            std::filesystem::remove(temp_initial, ec);
            ec.clear();
            std::filesystem::remove(temp_final, ec);
        };

        if (!initial_.writeXYZ(temp_initial) || !final_.writeXYZ(temp_final)) {
            std::cerr << "Warning: cannot create temporary files for alignment; using original structures." << std::endl;
            cleanup();
            return true;
        }

        if (aligner_.alignAndReplace(temp_initial, temp_final)) {
            geom::Structure final_aligned;
            if (final_aligned.readXYZ(temp_final)) {
                const double rmsd_after = initial_.calculateRMSD(final_aligned);
                std::cout << "RMSD after alignment: " << std::fixed << std::setprecision(6) << rmsd_after << std::endl;

                if (rmsd_after < rmsd_before) {
                    final_ = final_aligned;
                    std::cout << "Alignment improved RMSD, using aligned final structure." << std::endl;
                } else {
                    std::cout << "Alignment did not improve RMSD, using original final structure." << std::endl;
                }
            } else {
                std::cerr << "Warning: cannot read aligned structure, using original." << std::endl;
            }
        } else {
            std::cerr << "Warning: alignment failed, using original structures." << std::endl;
        }

        cleanup();
    }

    return true;
}

void NEBDriver::collectEndpointArrays_(std::vector<std::string>& symbols,
                                       std::vector<double>& x0,
                                       std::vector<double>& x1) const
{
    const size_t n = initial_.size();
    symbols.resize(n);
    x0.assign(3 * n, 0.0);
    x1.assign(3 * n, 0.0);

    for (size_t i = 0; i < n; ++i) {
        symbols[i] = initial_.atoms[i].symbol;
        x0[3 * i + 0] = initial_.atoms[i].x;
        x0[3 * i + 1] = initial_.atoms[i].y;
        x0[3 * i + 2] = initial_.atoms[i].z;

        x1[3 * i + 0] = final_.atoms[i].x;
        x1[3 * i + 1] = final_.atoms[i].y;
        x1[3 * i + 2] = final_.atoms[i].z;
    }
}

void NEBDriver::loadImagesFromCartesian_(const std::vector<std::string>& symbols,
                                         const std::vector<std::vector<double>>& imgs_xyz,
                                         PathMethod actual_method,
                                         PathMethod requested_method)
{
    images_.clear();
    images_.resize(num_images_);

    for (int img = 0; img < num_images_; ++img) {
        images_[img].atoms.resize(symbols.size());
        images_[img].comment = makeDirectImageComment(actual_method, requested_method, img + 1);
        for (size_t i = 0; i < symbols.size(); ++i) {
            images_[img].atoms[i].symbol = symbols[i];
            images_[img].atoms[i].x = imgs_xyz[img][3 * i + 0];
            images_[img].atoms[i].y = imgs_xyz[img][3 * i + 1];
            images_[img].atoms[i].z = imgs_xyz[img][3 * i + 2];
        }
    }
}

NEBDriver::PathBuildResult NEBDriver::generateLICPath_(PathMethod requested_method,
                                                       const std::string& intro_message)
{
    if (!intro_message.empty()) {
        std::cout << intro_message << std::endl;
    }
    std::cout << "  Generating path using Cartesian linear interpolation (LIC)..." << std::endl;

    PathBuildResult result;
    result.path_generated = true;
    result.actual_method = PathMethod::LIC;
    result.attempted_methods = {PathMethod::LIC};
    result.detail = (requested_method == PathMethod::LIC)
        ? "Path generated with LIC."
        : ("Requested " + toString(requested_method) + " was not realized; LIC generated the path.");

    images_.clear();
    images_.resize(num_images_);

    for (int img = 0; img < num_images_; ++img) {
        const double factor = static_cast<double>(img + 1) / (num_images_ + 1);

        images_[img].atoms.resize(initial_.size());
        images_[img].comment = makeDirectImageComment(PathMethod::LIC, requested_method, img + 1);

        for (size_t i = 0; i < initial_.size(); ++i) {
            images_[img].atoms[i].symbol = initial_.atoms[i].symbol;
            images_[img].atoms[i].x = initial_.atoms[i].x + factor * (final_.atoms[i].x - initial_.atoms[i].x);
            images_[img].atoms[i].y = initial_.atoms[i].y + factor * (final_.atoms[i].y - initial_.atoms[i].y);
            images_[img].atoms[i].z = initial_.atoms[i].z + factor * (final_.atoms[i].z - initial_.atoms[i].z);
        }
    }

    return result;
}

NEBDriver::PathBuildResult NEBDriver::generateLIICPath_(std::string* err) {
    if (err) err->clear();
    std::cout << "  Generating path using linear interpolation in internal coordinates (LIIC)..." << std::endl;

    PathBuildResult result;
    result.attempted_methods.push_back(PathMethod::LIIC);

    std::vector<std::string> symbols;
    std::vector<double> x0;
    std::vector<double> x1;
    collectEndpointArrays_(symbols, x0, x1);

    std::vector<std::vector<double>> imgs_xyz;
    std::string liic_reason;

    try {
        bool ok = ICInterp::interpolate_liic(symbols, x0, x1, num_images_, imgs_xyz, liic_options_, &liic_reason);
        if (ok) {
            loadImagesFromCartesian_(symbols, imgs_xyz, PathMethod::LIIC, PathMethod::LIIC);
            result.path_generated = true;
            result.actual_method = PathMethod::LIIC;
            result.detail = "Path generated with LIIC.";
            return result;
        }
        std::cerr << "  LIIC failed: " << liic_reason << std::endl;
    } catch (const ChemData::UnknownElementError& e) {
        result.fatal = true;
        result.detail = e.what();
        if (err) *err = e.what();
        return result;
    } catch (const std::exception& e) {
        liic_reason = e.what();
        std::cerr << "  LIIC failed with exception: " << e.what() << std::endl;
    }

    if (enable_liic_to_dm_fallback_) {
        result.attempted_methods.push_back(PathMethod::DM);
        std::cout << "  Attempting distance-matrix fallback (DM)..." << std::endl;

        std::string dm_reason;
        try {
            bool ok_dm = ICInterp::interpolate_dm(symbols, x0, x1, num_images_, imgs_xyz, dm_options_, &dm_reason);
            if (ok_dm) {
                loadImagesFromCartesian_(symbols, imgs_xyz, PathMethod::DM, PathMethod::LIIC);
                result.path_generated = true;
                result.actual_method = PathMethod::DM;
                result.detail = decorateFailureReason("LIIC failed and DM generated the path", liic_reason) + ".";
                return result;
            }
            std::cerr << "  DM fallback failed: " << dm_reason << std::endl;
            liic_reason = decorateFailureReason(liic_reason.empty() ? "LIIC failure" : liic_reason,
                                                "DM fallback also failed: " + dm_reason);
        } catch (const ChemData::UnknownElementError& e) {
            result.fatal = true;
            result.detail = e.what();
            if (err) *err = e.what();
            return result;
        } catch (const std::exception& e) {
            std::cerr << "  DM fallback failed with exception: " << e.what() << std::endl;
            liic_reason = decorateFailureReason(liic_reason.empty() ? "LIIC failure" : liic_reason,
                                                std::string("DM fallback exception: ") + e.what());
        }
    } else {
        std::cout << "  LIIC->DM fallback disabled; skipping DM." << std::endl;
    }

    PathBuildResult lic_result = generateLICPath_(PathMethod::LIIC,
                                                  "  Falling back to Cartesian linear interpolation (LIC)...");
    lic_result.attempted_methods = result.attempted_methods;
    lic_result.attempted_methods.push_back(PathMethod::LIC);
    lic_result.detail = decorateFailureReason("LIIC could not be realized; LIC generated the path", liic_reason) + ".";
    return lic_result;
}

NEBDriver::PathBuildResult NEBDriver::generateDMPath_(std::string* err) {
    if (err) err->clear();
    std::cout << "  Generating path using distance-matrix interpolation (DM)..." << std::endl;

    PathBuildResult result;
    result.attempted_methods.push_back(PathMethod::DM);

    std::vector<std::string> symbols;
    std::vector<double> x0;
    std::vector<double> x1;
    collectEndpointArrays_(symbols, x0, x1);

    std::vector<std::vector<double>> imgs_xyz;
    std::string dm_reason;

    try {
        bool ok = ICInterp::interpolate_dm(symbols, x0, x1, num_images_, imgs_xyz, dm_options_, &dm_reason);
        if (ok) {
            loadImagesFromCartesian_(symbols, imgs_xyz, PathMethod::DM, PathMethod::DM);
            result.path_generated = true;
            result.actual_method = PathMethod::DM;
            result.detail = "Path generated with DM.";
            return result;
        }
        std::cerr << "  DM failed: " << dm_reason << std::endl;
    } catch (const ChemData::UnknownElementError& e) {
        result.fatal = true;
        result.detail = e.what();
        if (err) *err = e.what();
        return result;
    } catch (const std::exception& e) {
        dm_reason = e.what();
        std::cerr << "  DM failed with exception: " << e.what() << std::endl;
    }

    PathBuildResult lic_result = generateLICPath_(PathMethod::DM,
                                                  "  Falling back to Cartesian linear interpolation (LIC)...");
    lic_result.attempted_methods = result.attempted_methods;
    lic_result.attempted_methods.push_back(PathMethod::LIC);
    lic_result.detail = decorateFailureReason("DM could not be realized; LIC generated the path", dm_reason) + ".";
    return lic_result;
}

bool NEBDriver::optimizeNEB_(std::string* err) {
    std::cout << "Performing Nudged Elastic Band (NEB) optimization..." << std::endl;

    const int natoms = static_cast<int>(initial_.size());
    if (natoms <= 0) {
        if (err) *err = "Error: no atoms in input structure.";
        return false;
    }

    const double spring_constant = use_external_engine_ ? engine_cfg_.spring_constant : 1.0;

    std::unique_ptr<NEBEngine> engine;
    if (use_external_engine_) {
        engine = std::make_unique<ExternalEngine>(engine_cfg_);
    } else {
        engine = std::make_unique<DistancePenaltyEngine>(initial_, final_);
    }

    std::vector<std::vector<double>> real_forces; // [img][3N]
    std::vector<std::vector<double>> neb_forces;

    for (int cycle = 0; cycle < max_iterations_; ++cycle) {
        double max_force = 0.0;

        std::string engine_err;
        if (!engine->compute_forces(images_, cycle, real_forces, &engine_err)) {
            if (err) *err = "Error: engine failed: " + engine_err;
            return false;
        }

        if (real_forces.size() != static_cast<size_t>(num_images_)) {
            if (err) {
                *err = "Error: engine returned wrong number of images (got " + std::to_string(real_forces.size()) +
                       ", expected " + std::to_string(num_images_) + ")";
            }
            return false;
        }

        std::string core_err;
        if (!assembleNEBForces(initial_, final_, images_, real_forces, spring_constant,
                               neb_forces, &max_force, &core_err)) {
            if (err) *err = "Error: NEB core failed: " + core_err;
            return false;
        }

        // Update images
        for (int img = 0; img < num_images_; ++img) {
            for (int i = 0; i < natoms; ++i) {
                images_[img].atoms[i].x += step_init_ * neb_forces[img][3 * i + 0];
                images_[img].atoms[i].y += step_init_ * neb_forces[img][3 * i + 1];
                images_[img].atoms[i].z += step_init_ * neb_forces[img][3 * i + 2];
            }
        }

        if (cycle % 100 == 0 || cycle < 10) {
            std::cout << "  NEB cycle " << cycle << ", max force = "
                      << std::scientific << std::setprecision(4) << max_force << std::endl;
        }

        if (max_force < convergence_threshold_) {
            std::cout << "  NEB converged after " << cycle << " cycles!" << std::endl;
            break;
        }
    }

    return true;
}

void NEBDriver::updateNEBImageComments_(Method requested_mode,
                                        const PathBuildResult& init_result)
{
    const PathMethod requested_method = requestedPathMethod(requested_mode);
    for (int img = 0; img < num_images_; ++img) {
        images_[img].comment = makeNEBImageComment(requested_mode,
                                                   init_result.actual_method,
                                                   requested_method,
                                                   img + 1);
    }
}

void NEBDriver::resetRunReport_(Method method) {
    last_run_report_ = RunReport{};
    last_run_report_.requested_mode = method;
    last_run_report_.requested_path_method = requestedPathMethod(method);
    last_run_report_.optimized_with_neb = modeUsesNEB(method);
}

void NEBDriver::finalizeRunReport_(Method method,
                                   const PathBuildResult& result,
                                   bool optimized_with_neb,
                                   const std::string& detail_override)
{
    last_run_report_.requested_mode = method;
    last_run_report_.requested_path_method = requestedPathMethod(method);
    last_run_report_.actual_path_method = result.actual_method;
    last_run_report_.attempted_path_methods = result.attempted_methods;
    last_run_report_.optimized_with_neb = optimized_with_neb;
    last_run_report_.path_generated = result.path_generated && !result.fatal;
    last_run_report_.used_fallback =
        (result.actual_method != PathMethod::NONE && result.actual_method != last_run_report_.requested_path_method);
    last_run_report_.detail = detail_override.empty() ? result.detail : detail_override;

    if (result.fatal || !last_run_report_.path_generated) {
        last_run_report_.status = RunStatus::FATAL_ERROR;
    } else if (last_run_report_.used_fallback) {
        last_run_report_.status = RunStatus::GENERATED_WITH_FALLBACK;
    } else {
        last_run_report_.status = RunStatus::GENERATED_EXACT;
    }
}

bool NEBDriver::run(Method method, std::string* err) {
    images_.clear();
    resetRunReport_(method);
    if (err) err->clear();

    PathBuildResult path_result;

    switch (method) {
        case Method::LIC:
            path_result = generateLICPath_(PathMethod::LIC);
            finalizeRunReport_(method, path_result, false);
            return last_run_report_.path_generated;

        case Method::LIIC:
            path_result = generateLIICPath_(err);
            finalizeRunReport_(method, path_result, false);
            if (!last_run_report_.path_generated && err && err->empty()) {
                *err = last_run_report_.detail;
            }
            return last_run_report_.path_generated;

        case Method::DM:
            path_result = generateDMPath_(err);
            finalizeRunReport_(method, path_result, false);
            if (!last_run_report_.path_generated && err && err->empty()) {
                *err = last_run_report_.detail;
            }
            return last_run_report_.path_generated;

        case Method::NEB:
        case Method::NEB_LIIC: {
            path_result = (method == Method::NEB)
                ? generateLICPath_(PathMethod::LIC)
                : generateLIICPath_(err);

            if (path_result.fatal || !path_result.path_generated) {
                finalizeRunReport_(method, path_result, true);
                if (err && err->empty()) {
                    *err = last_run_report_.detail;
                }
                return false;
            }

            if (!optimizeNEB_(err)) {
                PathBuildResult failed = path_result;
                failed.path_generated = false;
                failed.fatal = true;
                std::string detail = toString(method) + " failed during NEB optimization";
                if (err && !err->empty()) {
                    detail += ": " + *err;
                }
                finalizeRunReport_(method, failed, true, detail);
                return false;
            }

            updateNEBImageComments_(method, path_result);

            std::ostringstream detail;
            detail << toString(method) << " completed with " << toString(path_result.actual_method)
                   << " initialization";
            if (path_result.actual_method != requestedPathMethod(method)) {
                detail << " (fallback from " << toString(requestedPathMethod(method)) << ")";
            }
            detail << ".";

            finalizeRunReport_(method, path_result, true, detail.str());
            if (err) err->clear();
            return true;
        }

        default:
            if (err) *err = "Error: unknown method.";
            path_result.fatal = true;
            path_result.detail = "Error: unknown method.";
            finalizeRunReport_(method, path_result, false);
            return false;
    }
}

bool NEBDriver::writeResults(const std::string& prefix, bool multiframe) const {
    if (multiframe) {
        return writeMultiframeXYZ_(prefix);
    }
    return writeSeparateXYZ_(prefix);
}

bool NEBDriver::writeSeparateXYZ_(const std::string& prefix) const {
    std::string filename = prefix + "00.xyz";
    if (!initial_.writeXYZ(filename, "Initial structure")) return false;
    std::cout << "Wrote " << filename << std::endl;

    for (int i = 0; i < num_images_; ++i) {
        filename = prefix + (i + 1 < 10 ? "0" : "") + std::to_string(i + 1) + ".xyz";
        if (!images_[i].writeXYZ(filename)) return false;
        std::cout << "Wrote " << filename << std::endl;
    }

    filename = prefix + (num_images_ + 1 < 10 ? "0" : "") + std::to_string(num_images_ + 1) + ".xyz";
    if (!final_.writeXYZ(filename, "Final structure")) return false;
    std::cout << "Wrote " << filename << std::endl;

    return true;
}

bool NEBDriver::writeMultiframeXYZ_(const std::string& prefix) const {
    const std::string filename = prefix + "trajectory.xyz";
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot create file " << filename << std::endl;
        return false;
    }

    auto write_frame = [&](const geom::Structure& s, const std::string& cmt) {
        file << s.size() << "\n";
        file << cmt << "\n";
        file << std::fixed << std::setprecision(8);
        for (const auto& atom : s.atoms) {
            file << std::setw(2) << std::left << atom.symbol << " "
                 << std::setw(15) << atom.x << " "
                 << std::setw(15) << atom.y << " "
                 << std::setw(15) << atom.z << "\n";
        }
    };

    write_frame(initial_, "Initial structure");
    for (int i = 0; i < num_images_; ++i) {
        write_frame(images_[i], images_[i].comment);
    }
    write_frame(final_, "Final structure");

    std::cout << "Wrote " << filename << std::endl;
    return true;
}

} // namespace neb
