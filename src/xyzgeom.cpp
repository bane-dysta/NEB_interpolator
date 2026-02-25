#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <set>
#include <numeric>
#include <cstring>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <system_error>
#include "geometry.h"
#include "util.h"
#include "neb_driver.h"
#include "rmsd_align.h"

using Atom = geom::Atom;


class MoleculeAnalyzer {
private:
    std::vector<Atom> atoms;
    std::string currentFile;
    bool useBohr = false; // false = Angstrom, true = Bohr
    const double BOHR_TO_ANG = 0.529177211;
    
    // Parse Gview-style index selection string.
    //
    // Two variants are provided:
    //   - parseIndicesOrdered(): preserves the user input order (needed when the
    //     order has semantic meaning, e.g. angle/dihedral/vector direction).
    //   - parseIndicesSet(): set semantics (sorted + unique), useful for
    //     "select a group of atoms" operations.
    //
    // Supports both full-width and half-width commas.
    std::vector<int> parseIndicesOrdered(const std::string& str) {
        std::vector<int> indices;
        std::string cleaned = str;

        // Replace full-width comma with half-width
        size_t pos = 0;
        while ((pos = cleaned.find("，")) != std::string::npos) {
            cleaned.replace(pos, 3, ","); // Full-width comma is 3 bytes
        }

        auto is_space = [](unsigned char c) { return std::isspace(c) != 0; };

        auto parse_int = [](const std::string& s, int& out) -> bool {
            try {
                size_t idx = 0;
                int v = std::stoi(s, &idx);
                if (idx != s.size()) return false;
                out = v;
                return true;
            } catch (...) {
                return false;
            }
        };

        std::stringstream ss(cleaned);
        std::string token;

        while (std::getline(ss, token, ',')) {
            // Remove whitespace
            token.erase(std::remove_if(token.begin(), token.end(), is_space), token.end());
            if (token.empty()) continue;

            size_t dashPos = token.find('-');
            if (dashPos != std::string::npos) {
                // Range format: start-end
                int start = 0, end = 0;
                if (!parse_int(token.substr(0, dashPos), start) || !parse_int(token.substr(dashPos + 1), end)) {
                    continue;
                }
                if (start > end) std::swap(start, end);
                for (int i = start; i <= end; i++) {
                    if (i > 0 && static_cast<size_t>(i) <= atoms.size()) {
                        indices.push_back(i - 1); // Convert to 0-based
                    }
                }
            } else {
                // Single number
                int idx = 0;
                if (!parse_int(token, idx)) continue;
                if (idx > 0 && static_cast<size_t>(idx) <= atoms.size()) {
                    indices.push_back(idx - 1); // Convert to 0-based
                }
            }
        }

        return indices;
    }

    std::set<int> parseIndicesSet(const std::string& str) {
        const std::vector<int> ordered = parseIndicesOrdered(str);
        return std::set<int>(ordered.begin(), ordered.end());
    }

    static bool hasDuplicates(const std::vector<int>& v) {
        std::set<int> s(v.begin(), v.end());
        return s.size() != v.size();
    }

    bool swapAtomsBySelection_(const std::string& selection, std::string* err) {
        auto indices = parseIndicesSet(selection);
        if (indices.size() != 2) {
            if (err) *err = "Error: Please select exactly 2 atoms.";
            return false;
        }

        auto it = indices.begin();
        const int i = *it;
        ++it;
        const int j = *it;

        if (i < 0 || j < 0 || static_cast<size_t>(i) >= atoms.size() || static_cast<size_t>(j) >= atoms.size()) {
            if (err) *err = "Error: Atom index out of range.";
            return false;
        }

        std::swap(atoms[i], atoms[j]);
        return true;
    }

    bool mirrorThroughPlaneBySelection_(const std::string& selection, std::string* err) {
        auto indices = parseIndicesSet(selection);
        if (indices.size() != 3) {
            if (err) *err = "Error: Please select exactly 3 atoms.";
            return false;
        }

        std::vector<int> idx(indices.begin(), indices.end());
        if (idx[0] < 0 || idx[1] < 0 || idx[2] < 0 ||
            static_cast<size_t>(idx[0]) >= atoms.size() ||
            static_cast<size_t>(idx[1]) >= atoms.size() ||
            static_cast<size_t>(idx[2]) >= atoms.size()) {
            if (err) *err = "Error: Atom index out of range.";
            return false;
        }

        // Calculate plane equation ax + by + cz + d = 0
        double v1x = atoms[idx[1]].x - atoms[idx[0]].x;
        double v1y = atoms[idx[1]].y - atoms[idx[0]].y;
        double v1z = atoms[idx[1]].z - atoms[idx[0]].z;

        double v2x = atoms[idx[2]].x - atoms[idx[0]].x;
        double v2y = atoms[idx[2]].y - atoms[idx[0]].y;
        double v2z = atoms[idx[2]].z - atoms[idx[0]].z;

        double a, b, c;
        crossProduct(v1x, v1y, v1z, v2x, v2y, v2z, a, b, c);
        normalize(a, b, c);

        const double norm2 = a * a + b * b + c * c;
        if (norm2 < 1e-20) {
            if (err) *err = "Error: Invalid plane (three points are collinear or duplicated).";
            return false;
        }

        double d = -(a * atoms[idx[0]].x + b * atoms[idx[0]].y + c * atoms[idx[0]].z);

        // Mirror all atoms
        for (auto& atom : atoms) {
            double dist = a * atom.x + b * atom.y + c * atom.z + d;
            atom.x = atom.x - 2 * dist * a;
            atom.y = atom.y - 2 * dist * b;
            atom.z = atom.z - 2 * dist * c;
        }

        return true;
    }
    
    double getDistance(const Atom& a1, const Atom& a2) {
        double dx = a2.x - a1.x;
        double dy = a2.y - a1.y;
        double dz = a2.z - a1.z;
        double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
        return useBohr ? dist / BOHR_TO_ANG : dist;
    }
    
    void getVector(const Atom& a1, const Atom& a2, double& vx, double& vy, double& vz) {
        vx = a2.x - a1.x;
        vy = a2.y - a1.y;
        vz = a2.z - a1.z;
        if (useBohr) {
            vx /= BOHR_TO_ANG;
            vy /= BOHR_TO_ANG;
            vz /= BOHR_TO_ANG;
        }
    }
    
    double dotProduct(double x1, double y1, double z1, double x2, double y2, double z2) {
        return x1*x2 + y1*y2 + z1*z2;
    }
    
    void crossProduct(double x1, double y1, double z1, double x2, double y2, double z2,
                      double& rx, double& ry, double& rz) {
        rx = y1*z2 - z1*y2;
        ry = z1*x2 - x1*z2;
        rz = x1*y2 - y1*x2;
    }
    
    void normalize(double& x, double& y, double& z) {
        double mag = std::sqrt(x*x + y*y + z*z);
        if (mag > 1e-10) {
            x /= mag; y /= mag; z /= mag;
        }
    }
    
public:
    bool loadXYZ(const std::string& filename) {
        geom::Structure s;
        if (!s.readXYZ(filename)) {
            return false;
        }

        atoms = std::move(s.atoms);
        currentFile = filename;
        return true;
    }
    
    void saveXYZ(const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Cannot create file " << filename << std::endl;
            return;
        }
        
        file << atoms.size() << std::endl;
        file << "Generated by Molecular Geometry Analyzer" << std::endl;
        
        for (const auto& atom : atoms) {
            file << std::setw(4) << atom.symbol 
                 << std::fixed << std::setprecision(6)
                 << std::setw(12) << atom.x
                 << std::setw(12) << atom.y
                 << std::setw(12) << atom.z << std::endl;
        }
        
        file.close();
        std::cout << "File saved: " << filename << std::endl;
    }
    
    void printXYZ() {
        std::cout << atoms.size() << std::endl;
        std::cout << "Current structure" << std::endl;
        for (const auto& atom : atoms) {
            std::cout << std::setw(4) << atom.symbol 
                     << std::fixed << std::setprecision(6)
                     << std::setw(12) << atom.x
                     << std::setw(12) << atom.y
                     << std::setw(12) << atom.z << std::endl;
        }
    }
    
    // Function -1: Print current structure to screen
    void printCurrentStructure() {
        printXYZ();
    }
    
    // Function 13: Align with second XYZ file using calc_rmsd_xyz
    void alignWithSecondXYZ() {
        std::cout << "Enter reference XYZ filename to align with: ";
        std::string refFile;
        std::getline(std::cin, refFile);

        // Check if reference file exists
        std::ifstream check(refFile);
        if (!check.is_open()) {
            std::cerr << "Error: Cannot open reference file " << refFile << std::endl;
            return;
        }
        check.close();

        // Save current structure to a unique temp file
        const std::string tempFile = util::makeTempFilePath("xyzgeom_mobile_", ".xyz").string();
        saveXYZ(tempFile);

        // Use the shared RMSD alignment module (same lookup + invocation behavior as neb_interpolator)
        const std::string rmsd_exec = rmsd::findCalcRMSDExecutable();
        rmsd::FortranRMSDAligner aligner(rmsd_exec);

        std::cout << "Aligning with calc_rmsd_xyz: " << rmsd_exec << std::endl;

        if (!aligner.alignAndReplace(refFile, tempFile)) {
            std::cerr << "Error: Alignment failed" << std::endl;
            std::error_code ec;
            std::filesystem::remove(tempFile, ec);
            (void)ec;
            return;
        }

        // Load aligned structure back into memory
        const std::string oldFile = currentFile;
        if (!loadXYZ(tempFile)) {
            std::cerr << "Error: Cannot load aligned structure" << std::endl;
            std::error_code ec;
            std::filesystem::remove(tempFile, ec);
            (void)ec;
            return;
        }
        currentFile = oldFile;

        std::cout << "Structure aligned successfully and loaded into memory." << std::endl;

        // Clean up temp file
        std::error_code ec;
        std::filesystem::remove(tempFile, ec);
        (void)ec;
    }

    // Function 14: Path interpolation with second XYZ file (LIC / LIIC)
    void performNEBWithSecondXYZ() {
        std::cout << "Enter final XYZ filename for interpolation: ";
        std::string finalFile;
        std::getline(std::cin, finalFile);
        
        // Check if final file exists
        std::ifstream check(finalFile);
        if (!check.is_open()) {
            std::cerr << "Error: Cannot open final file " << finalFile << std::endl;
            return;
        }
        check.close();
        
        std::cout << "Enter number of intermediate images (default 5): ";
        std::string input;
        std::getline(std::cin, input);
        int numImages = 5;
        if (!input.empty()) {
            try {
                numImages = std::stoi(input);
            } catch (...) {
                std::cerr << "Warning: Invalid number of images, using default 5." << std::endl;
                numImages = 5;
            }
        }
        if (numImages <= 0) {
            std::cerr << "Warning: Number of images must be positive, using default 5." << std::endl;
            numImages = 5;
        }
        
        std::cout << "Choose method LIC or LIIC (default: LIC): ";
        std::getline(std::cin, input);

        // Normalize input: lowercase + strip whitespace
        std::string m = input;
        m.erase(std::remove_if(m.begin(), m.end(), [](unsigned char c){ return std::isspace(c) != 0; }), m.end());
        std::transform(m.begin(), m.end(), m.begin(), [](unsigned char c){ return static_cast<char>(std::tolower(c)); });

        bool useLIIC = (m == "liic" || m == "i");
        
        std::cout << "Enter output prefix (default: neb_): ";
        std::string prefix;
        std::getline(std::cin, prefix);
        if (prefix.empty()) prefix = "neb_";
        
        // Shared NEB / interpolation driver (same implementation as neb_interpolator)
        neb::NEBDriver driver(numImages);
        driver.setRMSDExecutable(rmsd::findCalcRMSDExecutable());

        // Convert current atoms to driver format (geom::Atom)
        std::vector<geom::Atom> init_atoms;
        init_atoms.reserve(atoms.size());
        for (const auto& atom : atoms) {
            init_atoms.emplace_back(atom.symbol, atom.x, atom.y, atom.z);
        }
        driver.setInitialFromMemory(init_atoms, "Initial structure (xyzgeom)");

        std::string err;
        if (!driver.setFinalFromFile(finalFile, &err)) {
            std::cerr << err << std::endl;
            return;
        }

        if (useLIIC) {
            if (!driver.run(neb::Method::LIIC, &err)) {
                std::cerr << (err.empty() ? "Error: LIIC failed" : err) << std::endl;
                return;
            }
        } else {
            if (!driver.run(neb::Method::LIC, &err)) {
                std::cerr << (err.empty() ? "Error: LIC failed" : err) << std::endl;
                return;
            }
        }

        if (!driver.writeResults(prefix, /*multiframe=*/false)) {
            std::cerr << "Error: Failed to write results" << std::endl;
            return;
        }

        std::cout << "Interpolation completed successfully!" << std::endl;
    }
    
    // Function 1: Calculate distance and vector between 2 atoms
    void calculateDistanceVector() {
        std::cout << "Enter indices of 2 atoms (i,j; vector points i -> j): ";
        std::string input;
        std::getline(std::cin, input);
        
        auto idx = parseIndicesOrdered(input);
        if (idx.size() != 2) {
            std::cout << "Error: Please select exactly 2 atoms." << std::endl;
            return;
        }

        if (hasDuplicates(idx)) {
            std::cout << "Error: Duplicate atom indices are not allowed." << std::endl;
            return;
        }
        
        double dist = getDistance(atoms[idx[0]], atoms[idx[1]]);
        
        double vx, vy, vz;
        getVector(atoms[idx[0]], atoms[idx[1]], vx, vy, vz);
        
        std::string unit = useBohr ? "Bohr" : "Angstrom";
        std::cout << "\nDistance: " << std::fixed << std::setprecision(6) 
                  << dist << " " << unit << std::endl;
        std::cout << "Vector: (" << vx << ", " << vy << ", " << vz << ") " << unit << std::endl;
    }
    
    // Function 2: Swap two atoms
    void swapAtoms() {
        std::cout << "Enter indices of 2 atoms to swap (Gview format): ";
        std::string input;
        std::getline(std::cin, input);

        std::string err;
        if (!swapAtomsBySelection_(input, &err)) {
            std::cout << err << std::endl;
            return;
        }
        std::cout << "Atoms swapped successfully." << std::endl;
    }
    
    // Function 3: Calculate bond angle and plane normal
    void calculateBondAngle() {
        std::cout << "Enter indices of 3 atoms (i,j,k; j is the vertex): ";
        std::string input;
        std::getline(std::cin, input);
        
        auto idx = parseIndicesOrdered(input);
        if (idx.size() != 3) {
            std::cout << "Error: Please select exactly 3 atoms." << std::endl;
            return;
        }

        if (hasDuplicates(idx)) {
            std::cout << "Error: Duplicate atom indices are not allowed." << std::endl;
            return;
        }
        
        // Calculate vectors
        double v1x = atoms[idx[0]].x - atoms[idx[1]].x;
        double v1y = atoms[idx[0]].y - atoms[idx[1]].y;
        double v1z = atoms[idx[0]].z - atoms[idx[1]].z;
        
        double v2x = atoms[idx[2]].x - atoms[idx[1]].x;
        double v2y = atoms[idx[2]].y - atoms[idx[1]].y;
        double v2z = atoms[idx[2]].z - atoms[idx[1]].z;
        
        // Calculate angle
        double dot = dotProduct(v1x, v1y, v1z, v2x, v2y, v2z);
        double mag1 = std::sqrt(v1x*v1x + v1y*v1y + v1z*v1z);
        double mag2 = std::sqrt(v2x*v2x + v2y*v2y + v2z*v2z);
        const double denom = mag1 * mag2;
        if (denom < 1e-12) {
            std::cout << "Error: Cannot compute angle (degenerate geometry; zero-length vector)." << std::endl;
            return;
        }
        double cosv = dot / denom;
        cosv = std::clamp(cosv, -1.0, 1.0);
        double angle = std::acos(cosv) * 180.0 / M_PI;
        
        // Calculate plane normal
        double nx, ny, nz;
        crossProduct(v1x, v1y, v1z, v2x, v2y, v2z, nx, ny, nz);
        normalize(nx, ny, nz);
        
        std::cout << "\nBond angle: " << std::fixed << std::setprecision(2) 
                  << angle << " degrees" << std::endl;
        std::cout << "Plane normal vector: (" << std::setprecision(6) 
                  << nx << ", " << ny << ", " << nz << ")" << std::endl;
    }
    
    // Function 4: Calculate angle between two planes
    void calculatePlaneAngle() {
        std::cout << "Enter indices of first 3 atoms (Gview format): ";
        std::string input1;
        std::getline(std::cin, input1);
        
        auto indices1 = parseIndicesSet(input1);
        if (indices1.size() != 3) {
            std::cout << "Error: Please select exactly 3 atoms for first plane." << std::endl;
            return;
        }
        
        std::cout << "Enter indices of second 3 atoms (Gview format): ";
        std::string input2;
        std::getline(std::cin, input2);
        
        auto indices2 = parseIndicesSet(input2);
        if (indices2.size() != 3) {
            std::cout << "Error: Please select exactly 3 atoms for second plane." << std::endl;
            return;
        }
        
        std::vector<int> idx1(indices1.begin(), indices1.end());
        std::vector<int> idx2(indices2.begin(), indices2.end());
        
        // Calculate normal for plane 1
        double v1x = atoms[idx1[1]].x - atoms[idx1[0]].x;
        double v1y = atoms[idx1[1]].y - atoms[idx1[0]].y;
        double v1z = atoms[idx1[1]].z - atoms[idx1[0]].z;
        
        double v2x = atoms[idx1[2]].x - atoms[idx1[0]].x;
        double v2y = atoms[idx1[2]].y - atoms[idx1[0]].y;
        double v2z = atoms[idx1[2]].z - atoms[idx1[0]].z;
        
        double n1x, n1y, n1z;
        crossProduct(v1x, v1y, v1z, v2x, v2y, v2z, n1x, n1y, n1z);
        normalize(n1x, n1y, n1z);
        const double n1norm2 = n1x*n1x + n1y*n1y + n1z*n1z;
        if (n1norm2 < 1e-20) {
            std::cout << "Error: First plane is invalid (points are collinear or duplicated)." << std::endl;
            return;
        }
        
        // Calculate normal for plane 2
        v1x = atoms[idx2[1]].x - atoms[idx2[0]].x;
        v1y = atoms[idx2[1]].y - atoms[idx2[0]].y;
        v1z = atoms[idx2[1]].z - atoms[idx2[0]].z;
        
        v2x = atoms[idx2[2]].x - atoms[idx2[0]].x;
        v2y = atoms[idx2[2]].y - atoms[idx2[0]].y;
        v2z = atoms[idx2[2]].z - atoms[idx2[0]].z;
        
        double n2x, n2y, n2z;
        crossProduct(v1x, v1y, v1z, v2x, v2y, v2z, n2x, n2y, n2z);
        normalize(n2x, n2y, n2z);
        const double n2norm2 = n2x*n2x + n2y*n2y + n2z*n2z;
        if (n2norm2 < 1e-20) {
            std::cout << "Error: Second plane is invalid (points are collinear or duplicated)." << std::endl;
            return;
        }
        
        // Calculate (acute) angle between planes: invariant to normal direction.
        double dot = dotProduct(n1x, n1y, n1z, n2x, n2y, n2z);
        dot = std::clamp(dot, -1.0, 1.0);
        double angle = std::acos(std::fabs(dot)) * 180.0 / M_PI;
        
        std::cout << "\nAngle between planes: " << std::fixed << std::setprecision(2) 
                  << angle << " degrees" << std::endl;
    }
    
    // Function 5: Mirror molecule through plane
    void mirrorThroughPlane() {
        std::cout << "Enter indices of 3 atoms defining the plane (Gview format): ";
        std::string input;
        std::getline(std::cin, input);

        std::string err;
        if (!mirrorThroughPlaneBySelection_(input, &err)) {
            std::cout << err << std::endl;
            return;
        }
        std::cout << "Molecule mirrored successfully." << std::endl;
    }
    
    // Function 6: Calculate dihedral angle
    void calculateDihedralAngle() {
        std::cout << "Enter indices of 4 atoms (i,j,k,l in order): ";
        std::string input;
        std::getline(std::cin, input);
        
        auto idx = parseIndicesOrdered(input);
        if (idx.size() != 4) {
            std::cout << "Error: Please select exactly 4 atoms." << std::endl;
            return;
        }

        if (hasDuplicates(idx)) {
            std::cout << "Error: Duplicate atom indices are not allowed." << std::endl;
            return;
        }
        
        // Vectors
        double v1x = atoms[idx[1]].x - atoms[idx[0]].x;
        double v1y = atoms[idx[1]].y - atoms[idx[0]].y;
        double v1z = atoms[idx[1]].z - atoms[idx[0]].z;
        
        double v2x = atoms[idx[2]].x - atoms[idx[1]].x;
        double v2y = atoms[idx[2]].y - atoms[idx[1]].y;
        double v2z = atoms[idx[2]].z - atoms[idx[1]].z;
        
        double v3x = atoms[idx[3]].x - atoms[idx[2]].x;
        double v3y = atoms[idx[3]].y - atoms[idx[2]].y;
        double v3z = atoms[idx[3]].z - atoms[idx[2]].z;
        
        // Normal vectors
        double n1x, n1y, n1z;
        crossProduct(v1x, v1y, v1z, v2x, v2y, v2z, n1x, n1y, n1z);
        
        double n2x, n2y, n2z;
        crossProduct(v2x, v2y, v2z, v3x, v3y, v3z, n2x, n2y, n2z);
        
        // Calculate dihedral angle
        double x = dotProduct(n1x, n1y, n1z, n2x, n2y, n2z);
        double y = std::sqrt(n1x*n1x + n1y*n1y + n1z*n1z) * std::sqrt(n2x*n2x + n2y*n2y + n2z*n2z);

        if (y < 1e-12) {
            std::cout << "Error: Cannot compute dihedral (degenerate geometry; plane normal is near zero)." << std::endl;
            return;
        }

        double cosv = x / y;
        cosv = std::clamp(cosv, -1.0, 1.0);
        double angle = std::acos(cosv) * 180.0 / M_PI;
        
        // Determine sign
        double cx, cy, cz;
        crossProduct(n1x, n1y, n1z, n2x, n2y, n2z, cx, cy, cz);
        if (dotProduct(cx, cy, cz, v2x, v2y, v2z) < 0) {
            angle = -angle;
        }
        
        std::cout << "\nDihedral angle: " << std::fixed << std::setprecision(2) 
                  << angle << " degrees" << std::endl;
    }
    
    // Function 7: Calculate geometric center
    void calculateGeometricCenter() {
        std::cout << "Enter indices of atoms (Gview format): ";
        std::string input;
        std::getline(std::cin, input);
        
        auto indices = parseIndicesSet(input);
        if (indices.empty()) {
            std::cout << "Error: Please select at least one atom." << std::endl;
            return;
        }
        
        double cx = 0, cy = 0, cz = 0;
        for (int idx : indices) {
            cx += atoms[idx].x;
            cy += atoms[idx].y;
            cz += atoms[idx].z;
        }
        
        cx /= indices.size();
        cy /= indices.size();
        cz /= indices.size();
        
        if (useBohr) {
            cx /= BOHR_TO_ANG;
            cy /= BOHR_TO_ANG;
            cz /= BOHR_TO_ANG;
        }
        
        std::string unit = useBohr ? "Bohr" : "Angstrom";
        std::cout << "\nGeometric center: (" << std::fixed << std::setprecision(6) 
                  << cx << ", " << cy << ", " << cz << ") " << unit << std::endl;
    }
    
    // Function 8: Export selected atoms
    void exportSelectedAtoms() {
        std::cout << "Enter indices of atoms to export (Gview format): ";
        std::string input;
        std::getline(std::cin, input);
        
        auto indices = parseIndicesSet(input);
        if (indices.empty()) {
            std::cout << "Error: Please select at least one atom." << std::endl;
            return;
        }
        
        std::cout << "Enter output filename: ";
        std::string filename;
        std::getline(std::cin, filename);
        
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Cannot create file " << filename << std::endl;
            return;
        }
        
        file << indices.size() << std::endl;
        file << "Selected atoms from " << currentFile << std::endl;
        
        for (int idx : indices) {
            file << std::setw(4) << atoms[idx].symbol 
                 << std::fixed << std::setprecision(6)
                 << std::setw(12) << atoms[idx].x
                 << std::setw(12) << atoms[idx].y
                 << std::setw(12) << atoms[idx].z << std::endl;
        }
        
        file.close();
        std::cout << "File saved: " << filename << std::endl;
    }
    
    // Function 9: Find atoms within radius
    void findAtomsWithinRadius() {
        std::cout << "Enter index of central atom: ";
        std::string input;
        std::getline(std::cin, input);
        
        auto indices = parseIndicesSet(input);
        if (indices.size() != 1) {
            std::cout << "Error: Please select exactly 1 atom." << std::endl;
            return;
        }
        
        int centerIdx = *indices.begin();
        
        std::cout << "Enter radius (in " << (useBohr ? "Bohr" : "Angstrom") << "): ";
        double radius;
        std::cin >> radius;
        std::cin.ignore();
        
        if (useBohr) {
            radius *= BOHR_TO_ANG; // Convert to Angstrom for calculation
        }
        
        std::set<int> nearbyAtoms;
        for (size_t i = 0; i < atoms.size(); i++) {
            if (static_cast<int>(i) != centerIdx) {
                double dist = std::sqrt(
                    std::pow(atoms[i].x - atoms[centerIdx].x, 2) +
                    std::pow(atoms[i].y - atoms[centerIdx].y, 2) +
                    std::pow(atoms[i].z - atoms[centerIdx].z, 2)
                );
                
                if (dist <= radius) {
                    nearbyAtoms.insert(static_cast<int>(i));
                }
            }
        }
        
        if (nearbyAtoms.empty()) {
            std::cout << "No atoms found within radius." << std::endl;
            return;
        }
        
        // Output Gview-style indices
        std::cout << "\nAtoms within radius (Gview indices): ";
        bool first = true;
        int start = -1, end = -1;
        for (int idx : nearbyAtoms) {
            if (start == -1) {
                start = end = idx + 1; // Convert to 1-based
            } else if (idx + 1 == end + 1) {
                end = idx + 1;
            } else {
                if (!first) std::cout << ",";
                if (start == end) {
                    std::cout << start;
                } else if (end == start + 1) {
                    std::cout << start << "," << end;
                } else {
                    std::cout << start << "-" << end;
                }
                first = false;
                start = end = idx + 1;
            }
        }
        if (!first) std::cout << ",";
        if (start == end) {
            std::cout << start;
        } else if (end == start + 1) {
            std::cout << start << "," << end;
        } else {
            std::cout << start << "-" << end;
        }
        std::cout << std::endl;
        
        // Save to file
        std::cout << "Enter output filename: ";
        std::string filename;
        std::getline(std::cin, filename);
        
        std::ofstream file(filename);
        file << nearbyAtoms.size() + 1 << std::endl; // +1 for center atom
        file << "Atoms within " << radius << " Angstrom of atom " << (centerIdx + 1) << std::endl;
        
        // Write center atom first
        file << std::setw(4) << atoms[centerIdx].symbol 
             << std::fixed << std::setprecision(6)
             << std::setw(12) << atoms[centerIdx].x
             << std::setw(12) << atoms[centerIdx].y
             << std::setw(12) << atoms[centerIdx].z << std::endl;
        
        // Write nearby atoms
        for (int idx : nearbyAtoms) {
            file << std::setw(4) << atoms[idx].symbol 
                 << std::fixed << std::setprecision(6)
                 << std::setw(12) << atoms[idx].x
                 << std::setw(12) << atoms[idx].y
                 << std::setw(12) << atoms[idx].z << std::endl;
        }
        
        file.close();
        std::cout << "File saved: " << filename << std::endl;
    }
    
    // Function 10: Export current structure
    void exportCurrentStructure() {
        std::cout << "Enter output filename: ";
        std::string filename;
        std::getline(std::cin, filename);
        saveXYZ(filename);
    }
    
    // Function 11: Load new XYZ file
    void loadNewFile() {
        std::cout << "Enter XYZ filename to load: ";
        std::string filename;
        std::getline(std::cin, filename);
        
        if (loadXYZ(filename)) {
            std::cout << "File loaded successfully. " << atoms.size() << " atoms." << std::endl;
        }
    }
    
    // Function 12: Toggle units
    void toggleUnits() {
        useBohr = !useBohr;
        std::cout << "Distance unit changed to: " << (useBohr ? "Bohr" : "Angstrom") << std::endl;
    }
    
    void showMenu() {
        std::cout << "\n========== Molecular Geometry Analyzer ==========" << std::endl;
        std::cout << "Current file: " << currentFile << std::endl;
        std::cout << "Atoms loaded: " << atoms.size() << std::endl;
        std::cout << "Distance unit: " << (useBohr ? "Bohr" : "Angstrom") << std::endl;
        std::cout << "\n-1. Print current structure to screen" << std::endl;
        std::cout << "1.  Calculate distance and vector between 2 atoms" << std::endl;
        std::cout << "2.  Swap two atoms" << std::endl;
        std::cout << "3.  Calculate bond angle and plane normal" << std::endl;
        std::cout << "4.  Calculate angle between two planes" << std::endl;
        std::cout << "5.  Mirror molecule through plane" << std::endl;
        std::cout << "6.  Calculate dihedral angle" << std::endl;
        std::cout << "7.  Calculate geometric center" << std::endl;
        std::cout << "8.  Export selected atoms to XYZ" << std::endl;
        std::cout << "9.  Find atoms within radius" << std::endl;
        std::cout << "10. Export current structure to XYZ" << std::endl;
        std::cout << "11. Load new XYZ file" << std::endl;
        std::cout << "12. Toggle distance unit (Bohr/Angstrom)" << std::endl;
        std::cout << "13. Align with second XYZ file (RMSD)" << std::endl;
        std::cout << "14. Interpolation with second XYZ file (LIC/LIIC)" << std::endl;
        std::cout << "0.  Exit" << std::endl;
        std::cout << "\nEnter choice: ";
    }
    
    void runCommandLine(int argc, char* argv[]) {
        // Command line mode for specific operations
        if (argc < 3) return;
        
        std::string operation = argv[2];
        
        if (operation == "--swap" && argc >= 4) {
            std::string err;
            if (!swapAtomsBySelection_(argv[3], &err)) {
                std::cerr << err << std::endl;
                return;
            }
            printXYZ();
        } else if (operation == "--mirror" && argc >= 4) {
            std::string err;
            if (!mirrorThroughPlaneBySelection_(argv[3], &err)) {
                std::cerr << err << std::endl;
                return;
            }
            printXYZ();
        } else if (operation == "--print") {
            printXYZ();
        }
    }
    
    void runInteractive() {
        while (true) {
            showMenu();
            
            int choice;
            std::cin >> choice;
            std::cin.ignore(); // Clear newline
            
            if (atoms.empty() && choice != 0 && choice != 11) {
                std::cout << "Error: No structure loaded. Please load an XYZ file first." << std::endl;
                continue;
            }
            
            switch (choice) {
                case -1: printCurrentStructure(); break;
                case 1: calculateDistanceVector(); break;
                case 2: swapAtoms(); break;
                case 3: calculateBondAngle(); break;
                case 4: calculatePlaneAngle(); break;
                case 5: mirrorThroughPlane(); break;
                case 6: calculateDihedralAngle(); break;
                case 7: calculateGeometricCenter(); break;
                case 8: exportSelectedAtoms(); break;
                case 9: findAtomsWithinRadius(); break;
                case 10: exportCurrentStructure(); break;
                case 11: loadNewFile(); break;
                case 12: toggleUnits(); break;
                case 13: alignWithSecondXYZ(); break;
                case 14: performNEBWithSecondXYZ(); break;
                case 0:
                    std::cout << "Exiting..." << std::endl;
                    return;
                default:
                    std::cout << "Invalid choice. Please try again." << std::endl;
            }
        }
    }
};

int main(int argc, char* argv[]) {
    MoleculeAnalyzer analyzer;
    
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <xyz_file> [options]" << std::endl;
        std::cout << "\nOptions for command-line processing:" << std::endl;
        std::cout << "  --swap <indices>    Swap two atoms and output new XYZ" << std::endl;
        std::cout << "  --mirror <indices>  Mirror through plane defined by 3 atoms" << std::endl;
        std::cout << "  --print             Print current structure to screen" << std::endl;
        std::cout << "\nExample: " << argv[0] << " molecule.xyz --swap 1,3" << std::endl;
        std::cout << "\nIf no options provided, enters interactive mode." << std::endl;
        return 1;
    }
    
    if (!analyzer.loadXYZ(argv[1])) {
        return 1;
    }
    
    if (argc > 2) {
        // Command line mode
        analyzer.runCommandLine(argc, argv);
    } else {
        // Interactive mode
        std::cout << "Loaded " << argv[1] << " successfully." << std::endl;
        analyzer.runInteractive();
    }
    
    return 0;
}