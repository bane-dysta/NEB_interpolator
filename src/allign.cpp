#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

struct AtomRecord {
    std::string element;
    Eigen::Vector3d coord = Eigen::Vector3d::Zero();
    std::string raw_line;
    bool has_atomic_number = false;
    int atomic_number = 0;
};

struct Molecule {
    std::vector<AtomRecord> atoms;
    std::vector<std::string> header_lines;
    std::vector<std::string> footer_lines;
    std::string xyz_comment;
};

struct AlignmentResult {
    double rmsd = 0.0;
    Eigen::Vector3d trans1 = Eigen::Vector3d::Zero();
    Eigen::Vector3d trans2 = Eigen::Vector3d::Zero();
    Eigen::Matrix3d rotation = Eigen::Matrix3d::Identity();
};

[[noreturn]] void fail(const std::string& msg) {
    throw std::runtime_error(msg);
}

std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return s;
}

std::string get_extension_lower(const std::string& filename) {
    return to_lower(fs::path(filename).extension().string());
}

std::vector<std::string> split_ws(const std::string& line) {
    std::istringstream iss(line);
    std::vector<std::string> tokens;
    std::string token;
    while (iss >> token) {
        tokens.push_back(token);
    }
    return tokens;
}

std::pair<int, int> parse_range(const std::string& s, const std::string& label) {
    auto pos = s.find('-');
    if (pos == std::string::npos) {
        fail("ERROR in subroutine rmsd: error " + label + ", no '-' found.\n" + label + " = " + s);
    }
    int a = 0, b = 0;
    try {
        a = std::stoi(s.substr(0, pos));
        b = std::stoi(s.substr(pos + 1));
    } catch (...) {
        fail("ERROR in subroutine rmsd: invalid " + label + ".\n" + label + " = " + s);
    }
    return {a, b};
}

AtomRecord parse_xyz_atom_line(const std::string& line) {
    auto tokens = split_ws(line);
    if (tokens.size() < 4) {
        fail("ERROR in subroutine read_coor_from_xyz: invalid xyz atom line.\nline = " + line);
    }
    AtomRecord atom;
    atom.element = tokens[0];
    atom.coord = Eigen::Vector3d(std::stod(tokens[1]), std::stod(tokens[2]), std::stod(tokens[3]));
    atom.raw_line = line;
    return atom;
}

AtomRecord parse_gjf_atom_line(const std::string& line) {
    auto tokens = split_ws(line);
    if (tokens.size() < 4) {
        fail("error in subroutine read_coor_from_gjf: wrong content in gjf.\nline = " + line);
    }
    AtomRecord atom;
    atom.element = tokens[0];
    atom.raw_line = line;
    if (tokens.size() <= 4) {
        atom.coord = Eigen::Vector3d(std::stod(tokens[1]), std::stod(tokens[2]), std::stod(tokens[3]));
    } else {
        atom.has_atomic_number = true;
        atom.atomic_number = std::stoi(tokens[1]);
        atom.coord = Eigen::Vector3d(std::stod(tokens[2]), std::stod(tokens[3]), std::stod(tokens[4]));
    }
    return atom;
}

Molecule read_xyz(const std::string& filename) {
    std::ifstream fin(filename);
    if (!fin) {
        fail("ERROR in subroutine rmsd: file " + filename + " does not exist!");
    }

    Molecule mol;
    std::string line;
    if (!std::getline(fin, line)) {
        fail("ERROR in subroutine read_natom_from_xyz: invalid natom.\nproblematic file: " + filename);
    }

    int natom = 0;
    try {
        natom = std::stoi(split_ws(line).at(0));
    } catch (...) {
        fail("ERROR in subroutine read_natom_from_xyz: invalid natom.\nproblematic file: " + filename + "\nnatom = " + line);
    }
    if (natom <= 0) {
        fail("ERROR in subroutine read_natom_from_xyz: invalid natom.\nproblematic file: " + filename + "\nnatom = " + std::to_string(natom));
    }

    if (!std::getline(fin, mol.xyz_comment)) {
        mol.xyz_comment.clear();
    }

    mol.atoms.reserve(static_cast<size_t>(natom));
    for (int i = 0; i < natom; ++i) {
        if (!std::getline(fin, line)) {
            fail("ERROR in subroutine read_coor_from_xyz: natom mismatch.\nproblematic file: " + filename);
        }
        mol.atoms.push_back(parse_xyz_atom_line(line));
    }
    return mol;
}

Molecule read_gjf(const std::string& filename) {
    std::ifstream fin(filename);
    if (!fin) {
        fail("ERROR in subroutine rmsd: file " + filename + " does not exist!");
    }

    Molecule mol;
    std::string line;
    int blank_count = 0;
    while (std::getline(fin, line)) {
        mol.header_lines.push_back(line);
        if (line.empty()) {
            ++blank_count;
        }
        if (blank_count == 2) {
            break;
        }
    }
    if (blank_count != 2) {
        fail("error in subroutine read_coor_from_gjf: wrong content in gjf.\nproblematic file: " + filename);
    }

    if (!std::getline(fin, line)) {
        fail("error in subroutine read_coor_from_gjf: wrong content in gjf.\nproblematic file: " + filename);
    }
    mol.header_lines.push_back(line); // charge and multiplicity

    while (std::getline(fin, line)) {
        if (line.empty()) {
            mol.footer_lines.push_back(line);
            break;
        }
        mol.atoms.push_back(parse_gjf_atom_line(line));
    }

    while (std::getline(fin, line)) {
        mol.footer_lines.push_back(line);
    }

    return mol;
}

Molecule read_molecule(const std::string& filename, const std::string& ext) {
    if (ext == ".xyz") {
        return read_xyz(filename);
    }
    if (ext == ".gjf") {
        return read_gjf(filename);
    }
    fail("ERROR: unsupported file format: " + filename);
}

Eigen::MatrixXd coords_from_range(const Molecule& mol, int begin1, int end1) {
    const int n = end1 - begin1 + 1;
    Eigen::MatrixXd coords(3, n);
    for (int i = 0; i < n; ++i) {
        coords.col(i) = mol.atoms.at(static_cast<size_t>(begin1 - 1 + i)).coord;
    }
    return coords;
}

Eigen::Vector3d move_coor_to_center(Eigen::MatrixXd& coords) {
    Eigen::Vector3d trans = -coords.rowwise().sum() / static_cast<double>(coords.cols());
    coords.colwise() += trans;
    return trans;
}

AlignmentResult rmsd_align(Eigen::MatrixXd& coor1, const Eigen::MatrixXd& coor2) {
    if (coor1.rows() != 3 || coor2.rows() != 3 || coor1.cols() != coor2.cols()) {
        fail("ERROR in subroutine rmsd: coordinate shape mismatch.");
    }

    AlignmentResult result;
    result.trans1 = move_coor_to_center(coor1);

    Eigen::MatrixXd centered2 = coor2;
    result.trans2 = move_coor_to_center(centered2);

    Eigen::Matrix3d H = coor1 * centered2.transpose();
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3d U = svd.matrixU();
    Eigen::Matrix3d V = svd.matrixV();

    Eigen::Matrix3d D = Eigen::Matrix3d::Identity();
    if ((U * V.transpose()).determinant() < 0.0) {
        D(2, 2) = -1.0;
    }

    result.rotation = U * D * V.transpose();
    coor1 = result.rotation.transpose() * coor1;

    double sumsq = 0.0;
    for (int i = 0; i < coor1.cols(); ++i) {
        sumsq += (coor1.col(i) - centered2.col(i)).squaredNorm();
    }
    result.rmsd = std::sqrt(sumsq / static_cast<double>(coor1.cols()));
    return result;
}

std::string output_filename(const std::string& filename, const std::string& ext) {
    fs::path p(filename);
    return (p.parent_path() / (p.stem().string() + "_new" + ext)).string();
}

void write_xyz(const std::string& filename, const Molecule& mol) {
    std::ofstream fout(filename);
    if (!fout) {
        fail("ERROR: failed to write file " + filename);
    }
    fout << mol.atoms.size() << "\n";
    fout << mol.xyz_comment << "\n";
    fout << std::fixed << std::setprecision(8);
    for (const auto& atom : mol.atoms) {
        fout << std::setw(2) << std::left << atom.element << ' ' << std::right
             << std::setw(16) << atom.coord.x()
             << std::setw(16) << atom.coord.y()
             << std::setw(16) << atom.coord.z() << '\n';
    }
}

void write_gjf(const std::string& filename, const Molecule& mol) {
    std::ofstream fout(filename);
    if (!fout) {
        fail("ERROR: failed to write file " + filename);
    }
    for (const auto& line : mol.header_lines) {
        fout << line << '\n';
    }
    fout << std::fixed << std::setprecision(8);
    for (const auto& atom : mol.atoms) {
        fout << std::setw(2) << std::left << atom.element << ' ' << std::right
             << std::setw(16) << atom.coord.x()
             << std::setw(16) << atom.coord.y()
             << std::setw(16) << atom.coord.z() << '\n';
    }
    for (const auto& line : mol.footer_lines) {
        fout << line << '\n';
    }
}

void write_molecule(const std::string& src_filename, const std::string& ext, const Molecule& mol) {
    const std::string out = output_filename(src_filename, ext);
    if (ext == ".xyz") {
        write_xyz(out, mol);
    } else if (ext == ".gjf") {
        write_gjf(out, mol);
    } else {
        fail("ERROR: unsupported output format: " + ext);
    }
}

int main(int argc, char* argv[]) {
    try {
        if (!(argc == 3 || argc == 5)) {
            std::cerr << "\n ERROR in subroutine rmsd: wrong command line arguments!\n";
            std::cerr << " Format: ./calc_rmsd_xyz_eigen file1 file2 [range1] [range2]\n";
            std::cerr << "\n Example 1: ./calc_rmsd_xyz_eigen a.gjf b.gjf\n";
            std::cerr << " Example 2: ./calc_rmsd_xyz_eigen a.xyz b.xyz\n";
            std::cerr << " Example 3: ./calc_rmsd_xyz_eigen a.gjf b.xyz 1-5 2-6\n";
            return 1;
        }

        std::string fname1 = argv[1];
        std::string fname2 = argv[2];
        std::string ext1 = get_extension_lower(fname1);
        std::string ext2 = get_extension_lower(fname2);

        if (ext1 != ".gjf" && ext1 != ".xyz") {
            fail("ERROR: unsupported file format for file1: " + fname1 + "\nSupported formats: .gjf, .xyz");
        }
        if (ext2 != ".gjf" && ext2 != ".xyz") {
            fail("ERROR: unsupported file format for file2: " + fname2 + "\nSupported formats: .gjf, .xyz");
        }
        if (!fs::exists(fname1)) {
            fail("ERROR in subroutine rmsd: file " + fname1 + " does not exist!");
        }
        if (!fs::exists(fname2)) {
            fail("ERROR in subroutine rmsd: file " + fname2 + " does not exist!");
        }

        int idx1 = 0, idx2 = 0, idx3 = 0, idx4 = 0;
        if (argc == 5) {
            auto r1 = parse_range(argv[3], "range1");
            auto r2 = parse_range(argv[4], "range2");
            idx1 = r1.first;
            idx2 = r1.second;
            idx3 = r2.first;
            idx4 = r2.second;
        }

        if (idx1 < 0 || idx2 < 0 || idx3 < 0 || idx4 < 0 || idx1 > idx2 || idx3 > idx4) {
            fail("ERROR in subroutine rmsd: input range not valid.");
        }
        if ((idx2 - idx1) != (idx4 - idx3)) {
            fail("ERROR in subroutine rmsd: number of atoms in two ranges are not equal.");
        }

        // Keep coordinates in file1 fixed; transform file2.
        Molecule fixed = read_molecule(fname1, ext1);
        Molecule moving = read_molecule(fname2, ext2);

        if (argc == 3) {
            idx1 = 1;
            idx2 = static_cast<int>(fixed.atoms.size());
            idx3 = 1;
            idx4 = static_cast<int>(moving.atoms.size());
        }

        if (static_cast<int>(fixed.atoms.size()) < idx2) {
            fail("ERROR in subroutine rmsd_wrapper: error range1.");
        }
        if (static_cast<int>(moving.atoms.size()) < idx4) {
            fail("ERROR in subroutine rmsd_wrapper: error range2.");
        }

        const int natom = idx2 - idx1 + 1;
        if (natom != idx4 - idx3 + 1) {
            fail("ERROR in subroutine rmsd_wrapper: number of atoms in two ranges are not equal.");
        }

        Eigen::MatrixXd coor1 = coords_from_range(moving, idx3, idx4);
        Eigen::MatrixXd coor2 = coords_from_range(fixed, idx1, idx2);
        AlignmentResult res = rmsd_align(coor1, coor2);

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "RMSD = " << std::setw(12) << res.rmsd << '\n';

        for (auto& atom : moving.atoms) {
            atom.coord += res.trans1;
            atom.coord = res.rotation.transpose() * atom.coord;
            atom.coord -= res.trans2;
        }

        write_molecule(fname2, ext2, moving);
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << ex.what() << '\n';
        return 1;
    }
}
