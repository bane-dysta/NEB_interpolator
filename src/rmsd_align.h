#ifndef RMSD_ALIGN_H
#define RMSD_ALIGN_H

#include <filesystem>
#include <string>

namespace rmsd {

// Locate an available alignment+RMSD backend executable.
//
// Fallback order:
//   1) calc_rmsd_xyz
//   2) allign_lapack
//   3) allign_eigen
//
// For each backend name, search order is:
//   a) PATH (via `which`)
//   b) Directory containing the current executable (Linux: /proc/self/exe)
//   c) Explicit default_path (checked first if non-empty)
std::string findCalcRMSDExecutable(const std::string& default_path = "./calc_rmsd_xyz");

class FortranRMSDAligner {
public:
    explicit FortranRMSDAligner(std::string exec_path = "./calc_rmsd_xyz");

    // Align mobile_file onto reference_file using calc_rmsd_xyz.
    // Does NOT overwrite the mobile file.
    // Returns true if the aligned output file was generated.
    bool alignStructures(const std::string& reference_file, const std::string& mobile_file);

    // Align and replace the mobile file with the aligned geometry.
    bool alignAndReplace(const std::string& reference_file, const std::string& mobile_file);

    void setExecutablePath(const std::string& path);

private:
    std::string rmsd_executable_;

    // calc_rmsd_xyz supports .xyz and .gjf and writes <name>_new.<ext>
    static std::filesystem::path inferAlignedOutputPath(const std::filesystem::path& mobile_file);
};

} // namespace rmsd

#endif // RMSD_ALIGN_H
