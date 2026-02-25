#include "neb_driver.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <system_error>
#include <vector>

#include "util.h"
#include "zmat/covalent_radii.h"

namespace neb {

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

void NEBDriver::setDMFallback(bool enable) {
    enable_dm_fallback_ = enable;
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

void NEBDriver::performLIC_() {
    std::cout << "  Initializing intermediate images using Cartesian linear interpolation (LIC)..." << std::endl;

    images_.clear();
    images_.resize(num_images_);

    for (int img = 0; img < num_images_; ++img) {
        const double factor = static_cast<double>(img + 1) / (num_images_ + 1);

        images_[img].atoms.resize(initial_.size());
        images_[img].comment = "LIC intermediate image " + std::to_string(img + 1);

        for (size_t i = 0; i < initial_.size(); ++i) {
            images_[img].atoms[i].symbol = initial_.atoms[i].symbol;
            images_[img].atoms[i].x = initial_.atoms[i].x + factor * (final_.atoms[i].x - initial_.atoms[i].x);
            images_[img].atoms[i].y = initial_.atoms[i].y + factor * (final_.atoms[i].y - initial_.atoms[i].y);
            images_[img].atoms[i].z = initial_.atoms[i].z + factor * (final_.atoms[i].z - initial_.atoms[i].z);
        }
    }
}

bool NEBDriver::performLIIC_(std::string* err) {
    std::cout << "  Initializing intermediate images using internal-coordinate interpolation (LIIC)..." << std::endl;

    const size_t n = initial_.size();
    std::vector<std::string> symbols(n);
    std::vector<double> x0(3 * n, 0.0), x1(3 * n, 0.0);

    for (size_t i = 0; i < n; ++i) {
        symbols[i] = initial_.atoms[i].symbol;
        x0[3 * i + 0] = initial_.atoms[i].x;
        x0[3 * i + 1] = initial_.atoms[i].y;
        x0[3 * i + 2] = initial_.atoms[i].z;

        x1[3 * i + 0] = final_.atoms[i].x;
        x1[3 * i + 1] = final_.atoms[i].y;
        x1[3 * i + 2] = final_.atoms[i].z;
    }

    std::vector<std::vector<double>> imgs_xyz;
    std::string local_err;

    try {
        bool ok = ICInterp::interpolate_liic(symbols, x0, x1, num_images_, imgs_xyz, liic_options_, &local_err);

        if (!ok) {
            std::cerr << "  LIIC failed: " << local_err << std::endl;

            if (enable_dm_fallback_) {
                std::cout << "  Falling back to distance-matrix interpolation (DM)..." << std::endl;
                bool ok_dm = ICInterp::interpolate_dm(symbols, x0, x1, num_images_, imgs_xyz, dm_options_, &local_err);
                if (!ok_dm) {
                    std::cerr << "  DM fallback failed: " << local_err << std::endl;
                    std::cout << "  Falling back to Cartesian LIC..." << std::endl;
                    performLIC_();
                    return false;
                }

                images_.clear();
                images_.resize(num_images_);
                for (int img = 0; img < num_images_; ++img) {
                    images_[img].atoms.resize(n);
                    images_[img].comment = "DM (fallback) intermediate image " + std::to_string(img + 1);
                    for (size_t i = 0; i < n; ++i) {
                        images_[img].atoms[i].symbol = symbols[i];
                        images_[img].atoms[i].x = imgs_xyz[img][3 * i + 0];
                        images_[img].atoms[i].y = imgs_xyz[img][3 * i + 1];
                        images_[img].atoms[i].z = imgs_xyz[img][3 * i + 2];
                    }
                }
                return true;
            }

            std::cout << "  DM fallback disabled. Falling back to Cartesian LIC..." << std::endl;
            performLIC_();
            return false;
        }

        images_.clear();
        images_.resize(num_images_);
        for (int img = 0; img < num_images_; ++img) {
            images_[img].atoms.resize(n);
            images_[img].comment = "LIIC intermediate image " + std::to_string(img + 1);
            for (size_t i = 0; i < n; ++i) {
                images_[img].atoms[i].symbol = symbols[i];
                images_[img].atoms[i].x = imgs_xyz[img][3 * i + 0];
                images_[img].atoms[i].y = imgs_xyz[img][3 * i + 1];
                images_[img].atoms[i].z = imgs_xyz[img][3 * i + 2];
            }
        }

        return true;

    } catch (const ChemData::UnknownElementError& e) {
        if (err) *err = e.what();
        // Fatal: do NOT fall back to LIC/DM, because the input is invalid.
        return false;
    } catch (const std::exception& e) {
        // Non-fatal: preserve legacy behavior (fall back)
        std::cerr << "  LIIC failed with exception: " << e.what() << std::endl;
        if (enable_dm_fallback_) {
            std::cout << "  Falling back to distance-matrix interpolation (DM)..." << std::endl;
            try {
                bool ok_dm = ICInterp::interpolate_dm(symbols, x0, x1, num_images_, imgs_xyz, dm_options_, &local_err);
                if (!ok_dm) {
                    std::cerr << "  DM fallback failed: " << local_err << std::endl;
                    std::cout << "  Falling back to Cartesian LIC..." << std::endl;
                    performLIC_();
                    return false;
                }

                images_.clear();
                images_.resize(num_images_);
                for (int img = 0; img < num_images_; ++img) {
                    images_[img].atoms.resize(n);
                    images_[img].comment = "DM (fallback) intermediate image " + std::to_string(img + 1);
                    for (size_t i = 0; i < n; ++i) {
                        images_[img].atoms[i].symbol = symbols[i];
                        images_[img].atoms[i].x = imgs_xyz[img][3 * i + 0];
                        images_[img].atoms[i].y = imgs_xyz[img][3 * i + 1];
                        images_[img].atoms[i].z = imgs_xyz[img][3 * i + 2];
                    }
                }
                return true;
            } catch (const ChemData::UnknownElementError& ee) {
                if (err) *err = ee.what();
                return false;
            } catch (const std::exception& ee) {
                std::cerr << "  DM fallback failed with exception: " << ee.what() << std::endl;
                std::cout << "  Falling back to Cartesian LIC..." << std::endl;
                performLIC_();
                return false;
            }
        }

        std::cout << "  DM fallback disabled. Falling back to Cartesian LIC..." << std::endl;
        performLIC_();
        return false;
    }
}

bool NEBDriver::performDM_(std::string* err) {
    std::cout << "  Initializing intermediate images using distance-matrix interpolation (DM)..." << std::endl;

    const size_t n = initial_.size();
    std::vector<std::string> symbols(n);
    std::vector<double> x0(3 * n, 0.0), x1(3 * n, 0.0);

    for (size_t i = 0; i < n; ++i) {
        symbols[i] = initial_.atoms[i].symbol;
        x0[3 * i + 0] = initial_.atoms[i].x;
        x0[3 * i + 1] = initial_.atoms[i].y;
        x0[3 * i + 2] = initial_.atoms[i].z;

        x1[3 * i + 0] = final_.atoms[i].x;
        x1[3 * i + 1] = final_.atoms[i].y;
        x1[3 * i + 2] = final_.atoms[i].z;
    }

    std::vector<std::vector<double>> imgs_xyz;
    std::string local_err;

    try {
        bool ok = ICInterp::interpolate_dm(symbols, x0, x1, num_images_, imgs_xyz, dm_options_, &local_err);
        if (!ok) {
            std::cerr << "  DM failed: " << local_err << std::endl;
            std::cout << "  Falling back to Cartesian LIC..." << std::endl;
            performLIC_();
            return false;
        }

        images_.clear();
        images_.resize(num_images_);
        for (int img = 0; img < num_images_; ++img) {
            images_[img].atoms.resize(n);
            images_[img].comment = "DM intermediate image " + std::to_string(img + 1);
            for (size_t i = 0; i < n; ++i) {
                images_[img].atoms[i].symbol = symbols[i];
                images_[img].atoms[i].x = imgs_xyz[img][3 * i + 0];
                images_[img].atoms[i].y = imgs_xyz[img][3 * i + 1];
                images_[img].atoms[i].z = imgs_xyz[img][3 * i + 2];
            }
        }

        return true;

    } catch (const ChemData::UnknownElementError& e) {
        if (err) *err = e.what();
        // Fatal: invalid input.
        return false;
    } catch (const std::exception& e) {
        std::cerr << "  DM failed with exception: " << e.what() << std::endl;
        std::cout << "  Falling back to Cartesian LIC..." << std::endl;
        performLIC_();
        return false;
    }
}

bool NEBDriver::performNEB_(bool init_with_liic, std::string* err) {
    std::cout << "Performing Nudged Elastic Band (NEB) optimization..." << std::endl;

    // Initialize images
    if (init_with_liic) {
        std::string liic_err;
        bool ok = performLIIC_(&liic_err);
        if (!ok && !liic_err.empty()) {
            // Fatal (unknown element): abort NEB.
            if (err) *err = liic_err;
            return false;
        }
        // Non-fatal LIIC failure may have fallen back to LIC; continue.
    } else {
        performLIC_();
    }

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

bool NEBDriver::run(Method method, std::string* err) {
    // Clear any previous path
    images_.clear();

    switch (method) {
        case Method::LIC:
            performLIC_();
            return true;

        case Method::LIIC: {
            std::string local;
            (void)performLIIC_(&local);
            if (!local.empty()) {
                if (err) *err = local;
                // Fatal only if local is set (currently only unknown element)
                return false;
            }
            return true;
        }

        case Method::DM: {
            std::string local;
            (void)performDM_(&local);
            if (!local.empty()) {
                if (err) *err = local;
                return false;
            }
            return true;
        }

        case Method::NEB:
            return performNEB_(false, err);

        case Method::NEB_LIIC:
            return performNEB_(true, err);

        default:
            if (err) *err = "Error: unknown method.";
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
