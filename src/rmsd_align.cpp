#include "rmsd_align.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits.h>
#include <system_error>
#include <vector>

#include <sys/wait.h>
#include <unistd.h>

#include "util.h"

namespace rmsd {

std::string findCalcRMSDExecutable(const std::string& default_path) {
    auto isExecutable = [](const std::string& path) -> bool {
        return !path.empty() && access(path.c_str(), X_OK) == 0;
    };

    const std::string default_name = std::filesystem::path(default_path).filename().string();
    std::vector<std::string> backend_names;
    if (!default_name.empty()) {
        backend_names.push_back(default_name);
    }
    for (const std::string& candidate : {std::string("calc_rmsd_xyz"), std::string("allign_lapack"), std::string("allign_eigen")}) {
        if (std::find(backend_names.begin(), backend_names.end(), candidate) == backend_names.end()) {
            backend_names.push_back(candidate);
        }
    }

    // 0) Explicit default path, if it already points to an executable.
    if (isExecutable(default_path)) {
        return default_path;
    }

    // 1) PATH lookup in fallback order.
    for (const std::string& backend : backend_names) {
        const std::string command = "which " + backend + " 2>/dev/null";
        FILE* pipe = popen(command.c_str(), "r");
        if (pipe != nullptr) {
            char buffer[PATH_MAX];
            if (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
                pclose(pipe);
                std::string path(buffer);
                if (!path.empty() && path.back() == '\n') {
                    path.pop_back();
                }
                if (isExecutable(path)) {
                    return path;
                }
            } else {
                pclose(pipe);
            }
        }
    }

    // 2) Same directory as current executable (Linux), also in fallback order.
    {
        char exe_path[PATH_MAX];
        ssize_t len = readlink("/proc/self/exe", exe_path, sizeof(exe_path) - 1);
        if (len != -1) {
            exe_path[len] = '\0';
            std::string exe_dir(exe_path);
            size_t last_slash = exe_dir.find_last_of('/');
            if (last_slash != std::string::npos) {
                exe_dir = exe_dir.substr(0, last_slash + 1);
                for (const std::string& backend : backend_names) {
                    const std::string candidate = exe_dir + backend;
                    if (isExecutable(candidate)) {
                        return candidate;
                    }
                }
            }
        }
    }

    // 3) Final fallback keeps the historical default for caller-side logging/errors.
    return default_path;
}

std::filesystem::path FortranRMSDAligner::inferAlignedOutputPath(const std::filesystem::path& mobile_file) {
    const std::string ext = mobile_file.extension().string();

    if (ext == ".xyz") {
        return mobile_file.parent_path() / (mobile_file.stem().string() + "_new.xyz");
    }
    if (ext == ".gjf") {
        return mobile_file.parent_path() / (mobile_file.stem().string() + "_new.gjf");
    }

    // Fallback (should not happen because calc_rmsd_xyz rejects unknown extensions)
    return mobile_file.parent_path() / (mobile_file.filename().string() + "_new" + ext);
}

FortranRMSDAligner::FortranRMSDAligner(std::string exec_path)
    : rmsd_executable_(std::move(exec_path)) {}

bool FortranRMSDAligner::alignStructures(const std::string& reference_file, const std::string& mobile_file) {
    // Quote paths to support spaces
    const std::string command = util::shellQuote(rmsd_executable_) + " " +
                                util::shellQuote(reference_file) + " " +
                                util::shellQuote(mobile_file);

    std::cout << "  Running RMSD alignment backend: " << command << std::endl;

    util::CommandResult result = util::runCommand(command);
    if (!(result.exited && result.exit_code == 0)) {
        std::cerr << "  Error: RMSD alignment backend failed (" << util::formatCommandFailure(result) << ")" << std::endl;
        return false;
    }

    const std::filesystem::path aligned_path = inferAlignedOutputPath(std::filesystem::path(mobile_file));
    if (!std::filesystem::exists(aligned_path)) {
        std::cerr << "  Error: Expected aligned file " << aligned_path.string() << " was not created" << std::endl;
        return false;
    }

    std::cout << "  Alignment completed. Aligned structure saved as: " << aligned_path.string() << std::endl;
    return true;
}

bool FortranRMSDAligner::alignAndReplace(const std::string& reference_file, const std::string& mobile_file) {
    if (!alignStructures(reference_file, mobile_file)) {
        return false;
    }

    const std::filesystem::path aligned_path = inferAlignedOutputPath(std::filesystem::path(mobile_file));

    // Copy aligned file back to original (avoid external cp; handles spaces)
    std::error_code ec;
    std::filesystem::copy_file(aligned_path, mobile_file,
                               std::filesystem::copy_options::overwrite_existing, ec);
    if (ec) {
        std::cerr << "  Error: Failed to replace original file: " << ec.message() << std::endl;
        return false;
    }

    // Remove temporary aligned file (non-fatal if it fails)
    ec.clear();
    std::filesystem::remove(aligned_path, ec);

    return true;
}

void FortranRMSDAligner::setExecutablePath(const std::string& path) {
    rmsd_executable_ = path;
}

} // namespace rmsd
