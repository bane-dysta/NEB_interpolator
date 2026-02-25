#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>

namespace geom {

struct Atom {
    std::string symbol;
    double x, y, z;

    Atom() : symbol(""), x(0.0), y(0.0), z(0.0) {}
    Atom(const std::string& s, double x_, double y_, double z_) : symbol(s), x(x_), y(y_), z(z_) {}
};

static inline std::string rstripCR(std::string s) {
    if (!s.empty() && s.back() == '\r') s.pop_back();
    return s;
}

class Structure {
public:
    std::vector<Atom> atoms;
    std::string comment;

    bool readXYZ(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Cannot open file " << filename << std::endl;
            return false;
        }

        int natoms;
        if (!(file >> natoms)) {
            std::cerr << "Error: Cannot read number of atoms from " << filename << std::endl;
            return false;
        }

        // Consume rest of the first line completely (CRLF safe)
        file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        // Read (possibly empty) comment line
        std::getline(file, comment);
        comment = rstripCR(comment);

        atoms.clear();
        atoms.reserve(natoms);

        for (int i = 0; i < natoms; ++i) {
            std::string sym;
            double x, y, z;
            if (!(file >> sym >> x >> y >> z)) {
                std::cerr << "Error: Cannot read atom " << i + 1 << " from " << filename << std::endl;
                return false;
            }
            atoms.emplace_back(sym, x, y, z);
        }

        return atoms.size() == static_cast<size_t>(natoms);
    }

    bool writeXYZ(const std::string& filename, const std::string& custom_comment = "") const {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Cannot create file " << filename << std::endl;
            return false;
        }

        file << atoms.size() << "\n";
        file << (custom_comment.empty() ? comment : custom_comment) << "\n";

        file << std::fixed << std::setprecision(8);
        for (const auto& atom : atoms) {
            file << std::setw(2) << std::left << atom.symbol << " "
                 << std::setw(15) << atom.x << " "
                 << std::setw(15) << atom.y << " "
                 << std::setw(15) << atom.z << "\n";
        }
        return true;
    }

    size_t size() const { return atoms.size(); }

    bool isCompatible(const Structure& other) const {
        if (atoms.size() != other.atoms.size()) return false;
        for (size_t i = 0; i < atoms.size(); ++i) {
            if (atoms[i].symbol != other.atoms[i].symbol) return false;
        }
        return true;
    }

    // RMSD without alignment
    double calculateRMSD(const Structure& other) const {
        if (atoms.size() != other.atoms.size()) return -1.0;

        double sum = 0.0;
        for (size_t i = 0; i < atoms.size(); ++i) {
            double dx = atoms[i].x - other.atoms[i].x;
            double dy = atoms[i].y - other.atoms[i].y;
            double dz = atoms[i].z - other.atoms[i].z;
            sum += dx * dx + dy * dy + dz * dz;
        }
        return std::sqrt(sum / atoms.size());
    }
};

} // namespace geom

#endif // GEOMETRY_H
