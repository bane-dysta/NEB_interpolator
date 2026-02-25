#include "rmsd_align.h"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits.h>
#include <system_error>

#include <sys/wait.h>
#include <unistd.h>

#include "util.h"

namespace rmsd {

std::string findCalcRMSDExecutable(const std::string& default_path) {
    // 1) PATH
    {
        const std::string command = "which calc_rmsd_xyz 2>/dev/null";
        FILE* pipe = popen(command.c_str(), "r");
        if (pipe != nullptr) {
            char buffer[PATH_MAX];
            if (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
                pclose(pipe);
                std::string path(buffer);
                if (!path.empty() && path.back() == '\n') {
                    path.pop_back();
                }
                if (!path.empty()) {
                    return path;
                }
            } else {
                pclose(pipe);
            }
        }
    }

    // 2) Same directory as current executable (Linux)
    {
        char exe_path[PATH_MAX];
        ssize_t len = readlink("/proc/self/exe", exe_path, sizeof(exe_path) - 1);
        if (len != -1) {
            exe_path[len] = '\0';
            std::string exe_dir(exe_path);
            size_t last_slash = exe_dir.find_last_of('/');
            if (last_slash != std::string::npos) {
                exe_dir = exe_dir.substr(0, last_slash + 1);
                std::string candidate = exe_dir + "calc_rmsd_xyz";
                if (access(candidate.c_str(), X_OK) == 0) {
                    return candidate;
                }
            }
        }
    }

    // 3) Fallback
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

    std::cout << "  Running Fortran RMSD alignment: " << command << std::endl;

    util::CommandResult result = util::runCommand(command);
    if (!(result.exited && result.exit_code == 0)) {
        std::cerr << "  Error: Fortran RMSD alignment failed (" << util::formatCommandFailure(result) << ")" << std::endl;
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
