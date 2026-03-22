#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

extern "C" {
void dgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k,
            const double* alpha, const double* A, const int* lda,
            const double* B, const int* ldb, const double* beta,
            double* C, const int* ldc);

void dgesvd_(const char* jobu, const char* jobvt, const int* m, const int* n,
             double* A, const int* lda, double* S, double* U, const int* ldu,
             double* VT, const int* ldvt, double* work, const int* lwork, int* info);
}

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    Vec3& operator+=(const Vec3& other) {
        x += other.x; y += other.y; z += other.z; return *this;
    }
    Vec3& operator-=(const Vec3& other) {
        x -= other.x; y -= other.y; z -= other.z; return *this;
    }
};

Vec3 operator+(Vec3 a, const Vec3& b) { return a += b; }
Vec3 operator-(Vec3 a, const Vec3& b) { return a -= b; }
Vec3 operator/(Vec3 a, double s) { a.x /= s; a.y /= s; a.z /= s; return a; }

struct Mat3 {
    // column-major, compatible with Fortran LAPACK/BLAS
    double a[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    static Mat3 identity() {
        Mat3 m;
        m(0,0) = m(1,1) = m(2,2) = 1.0;
        return m;
    }

    double& operator()(int r, int c) { return a[c * 3 + r]; }
    double operator()(int r, int c) const { return a[c * 3 + r]; }
};

Vec3 matvec(const Mat3& M, const Vec3& v) {
    return {
        M(0,0) * v.x + M(0,1) * v.y + M(0,2) * v.z,
        M(1,0) * v.x + M(1,1) * v.y + M(1,2) * v.z,
        M(2,0) * v.x + M(2,1) * v.y + M(2,2) * v.z
    };
}

Mat3 transpose(const Mat3& M) {
    Mat3 T;
    for (int r = 0; r < 3; ++r) for (int c = 0; c < 3; ++c) T(r, c) = M(c, r);
    return T;
}

Mat3 matmul(const Mat3& A, const Mat3& B) {
    Mat3 C;
    const char n = 'N';
    const int m = 3, k = 3, l = 3;
    const double alpha = 1.0, beta = 0.0;
    dgemm_(&n, &n, &m, &l, &k, &alpha, A.a, &m, B.a, &k, &beta, C.a, &m);
    return C;
}

double det3(const Mat3& M) {
    return M(0,0) * (M(1,1) * M(2,2) - M(1,2) * M(2,1))
         - M(0,1) * (M(1,0) * M(2,2) - M(1,2) * M(2,0))
         + M(0,2) * (M(1,0) * M(2,1) - M(1,1) * M(2,0));
}

struct AtomRecord {
    std::string element;
    Vec3 coord;
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
    Vec3 trans1;
    Vec3 trans2;
    Mat3 rotation = Mat3::identity();
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
    while (iss >> token) tokens.push_back(token);
    return tokens;
}

std::pair<int, int> parse_range(const std::string& s, const std::string& label) {
    auto pos = s.find('-');
    if (pos == std::string::npos) fail("ERROR in subroutine rmsd: error " + label + ", no '-' found.\n" + label + " = " + s);
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
    if (tokens.size() < 4) fail("ERROR in subroutine read_coor_from_xyz: invalid xyz atom line.\nline = " + line);
    AtomRecord atom;
    atom.element = tokens[0];
    atom.coord = {std::stod(tokens[1]), std::stod(tokens[2]), std::stod(tokens[3])};
    atom.raw_line = line;
    return atom;
}

AtomRecord parse_gjf_atom_line(const std::string& line) {
    auto tokens = split_ws(line);
    if (tokens.size() < 4) fail("error in subroutine read_coor_from_gjf: wrong content in gjf.\nline = " + line);
    AtomRecord atom;
    atom.element = tokens[0];
    atom.raw_line = line;
    if (tokens.size() <= 4) {
        atom.coord = {std::stod(tokens[1]), std::stod(tokens[2]), std::stod(tokens[3])};
    } else {
        atom.has_atomic_number = true;
        atom.atomic_number = std::stoi(tokens[1]);
        atom.coord = {std::stod(tokens[2]), std::stod(tokens[3]), std::stod(tokens[4])};
    }
    return atom;
}

Molecule read_xyz(const std::string& filename) {
    std::ifstream fin(filename);
    if (!fin) fail("ERROR in subroutine rmsd: file " + filename + " does not exist!");

    Molecule mol;
    std::string line;
    if (!std::getline(fin, line)) fail("ERROR in subroutine read_natom_from_xyz: invalid natom.\nproblematic file: " + filename);

    int natom = 0;
    try {
        natom = std::stoi(split_ws(line).at(0));
    } catch (...) {
        fail("ERROR in subroutine read_natom_from_xyz: invalid natom.\nproblematic file: " + filename + "\nnatom = " + line);
    }
    if (natom <= 0) fail("ERROR in subroutine read_natom_from_xyz: invalid natom.\nproblematic file: " + filename + "\nnatom = " + std::to_string(natom));

    if (!std::getline(fin, mol.xyz_comment)) mol.xyz_comment.clear();
    mol.atoms.reserve(static_cast<size_t>(natom));
    for (int i = 0; i < natom; ++i) {
        if (!std::getline(fin, line)) fail("ERROR in subroutine read_coor_from_xyz: natom mismatch.\nproblematic file: " + filename);
        mol.atoms.push_back(parse_xyz_atom_line(line));
    }
    return mol;
}

Molecule read_gjf(const std::string& filename) {
    std::ifstream fin(filename);
    if (!fin) fail("ERROR in subroutine rmsd: file " + filename + " does not exist!");

    Molecule mol;
    std::string line;
    int blank_count = 0;
    while (std::getline(fin, line)) {
        mol.header_lines.push_back(line);
        if (line.empty()) ++blank_count;
        if (blank_count == 2) break;
    }
    if (blank_count != 2) fail("error in subroutine read_coor_from_gjf: wrong content in gjf.\nproblematic file: " + filename);

    if (!std::getline(fin, line)) fail("error in subroutine read_coor_from_gjf: wrong content in gjf.\nproblematic file: " + filename);
    mol.header_lines.push_back(line);

    while (std::getline(fin, line)) {
        if (line.empty()) {
            mol.footer_lines.push_back(line);
            break;
        }
        mol.atoms.push_back(parse_gjf_atom_line(line));
    }
    while (std::getline(fin, line)) mol.footer_lines.push_back(line);
    return mol;
}

Molecule read_molecule(const std::string& filename, const std::string& ext) {
    if (ext == ".xyz") return read_xyz(filename);
    if (ext == ".gjf") return read_gjf(filename);
    fail("ERROR: unsupported file format: " + filename);
}

std::vector<double> coords_from_range(const Molecule& mol, int begin1, int end1) {
    const int n = end1 - begin1 + 1;
    std::vector<double> coords(static_cast<size_t>(3 * n));
    for (int i = 0; i < n; ++i) {
        const auto& c = mol.atoms.at(static_cast<size_t>(begin1 - 1 + i)).coord;
        coords[static_cast<size_t>(i * 3 + 0)] = c.x;
        coords[static_cast<size_t>(i * 3 + 1)] = c.y;
        coords[static_cast<size_t>(i * 3 + 2)] = c.z;
    }
    return coords;
}

Vec3 move_coor_to_center(std::vector<double>& coords) {
    const int n = static_cast<int>(coords.size() / 3);
    Vec3 sum{};
    for (int i = 0; i < n; ++i) {
        sum.x += coords[static_cast<size_t>(i * 3 + 0)];
        sum.y += coords[static_cast<size_t>(i * 3 + 1)];
        sum.z += coords[static_cast<size_t>(i * 3 + 2)];
    }
    Vec3 trans{-sum.x / n, -sum.y / n, -sum.z / n};
    for (int i = 0; i < n; ++i) {
        coords[static_cast<size_t>(i * 3 + 0)] += trans.x;
        coords[static_cast<size_t>(i * 3 + 1)] += trans.y;
        coords[static_cast<size_t>(i * 3 + 2)] += trans.z;
    }
    return trans;
}

AlignmentResult rmsd_align(std::vector<double>& coor1, const std::vector<double>& coor2_in) {
    if (coor1.size() != coor2_in.size() || coor1.size() % 3 != 0) fail("ERROR in subroutine rmsd: coordinate shape mismatch.");
    const int natom = static_cast<int>(coor1.size() / 3);

    AlignmentResult result;
    result.trans1 = move_coor_to_center(coor1);

    std::vector<double> coor2 = coor2_in;
    result.trans2 = move_coor_to_center(coor2);

    Mat3 H;
    {
        const char transa = 'N', transb = 'T';
        const int m = 3, n = 3, k = natom;
        const double alpha = 1.0, beta = 0.0;
        const int lda = 3, ldb = 3, ldc = 3;
        dgemm_(&transa, &transb, &m, &n, &k, &alpha, coor1.data(), &lda, coor2.data(), &ldb, &beta, H.a, &ldc);
    }

    double s[3] = {0.0, 0.0, 0.0};
    Mat3 U, VT;
    int m = 3, n = 3, lda = 3, ldu = 3, ldvt = 3, info = 0;
    int lwork = -1;
    double wkopt = 0.0;
    std::vector<double> Hcopy(std::begin(H.a), std::end(H.a));
    const char jobu = 'A', jobvt = 'A';
    dgesvd_(&jobu, &jobvt, &m, &n, Hcopy.data(), &lda, s, U.a, &ldu, VT.a, &ldvt, &wkopt, &lwork, &info);
    if (info != 0) fail("ERROR in subroutine rmsd: LAPACK dgesvd workspace query failed. info = " + std::to_string(info));
    lwork = static_cast<int>(wkopt);
    std::vector<double> work(static_cast<size_t>(std::max(1, lwork)));
    Hcopy.assign(std::begin(H.a), std::end(H.a));
    dgesvd_(&jobu, &jobvt, &m, &n, Hcopy.data(), &lda, s, U.a, &ldu, VT.a, &ldvt, work.data(), &lwork, &info);
    if (info != 0) fail("ERROR in subroutine rmsd: LAPACK dgesvd error info /= 0. info = " + std::to_string(info));

    Mat3 UVT = matmul(U, VT);
    Mat3 D = Mat3::identity();
    if (det3(UVT) < 0.0) D(2, 2) = -1.0;
    result.rotation = matmul(matmul(U, D), VT);
    const Mat3 RT = transpose(result.rotation);

    for (int i = 0; i < natom; ++i) {
        Vec3 p{coor1[static_cast<size_t>(i * 3 + 0)], coor1[static_cast<size_t>(i * 3 + 1)], coor1[static_cast<size_t>(i * 3 + 2)]};
        Vec3 q = matvec(RT, p);
        coor1[static_cast<size_t>(i * 3 + 0)] = q.x;
        coor1[static_cast<size_t>(i * 3 + 1)] = q.y;
        coor1[static_cast<size_t>(i * 3 + 2)] = q.z;
    }

    double sumsq = 0.0;
    for (int i = 0; i < natom; ++i) {
        const double dx = coor1[static_cast<size_t>(i * 3 + 0)] - coor2[static_cast<size_t>(i * 3 + 0)];
        const double dy = coor1[static_cast<size_t>(i * 3 + 1)] - coor2[static_cast<size_t>(i * 3 + 1)];
        const double dz = coor1[static_cast<size_t>(i * 3 + 2)] - coor2[static_cast<size_t>(i * 3 + 3 - 1)];
        sumsq += dx * dx + dy * dy + dz * dz;
    }
    result.rmsd = std::sqrt(sumsq / natom);
    return result;
}

std::string output_filename(const std::string& filename, const std::string& ext) {
    fs::path p(filename);
    return (p.parent_path() / (p.stem().string() + "_new" + ext)).string();
}

void write_xyz(const std::string& filename, const Molecule& mol) {
    std::ofstream fout(filename);
    if (!fout) fail("ERROR: failed to write file " + filename);
    fout << mol.atoms.size() << "\n";
    fout << mol.xyz_comment << "\n";
    fout << std::fixed << std::setprecision(8);
    for (const auto& atom : mol.atoms) {
        fout << std::setw(2) << std::left << atom.element << ' ' << std::right
             << std::setw(16) << atom.coord.x
             << std::setw(16) << atom.coord.y
             << std::setw(16) << atom.coord.z << '\n';
    }
}

void write_gjf(const std::string& filename, const Molecule& mol) {
    std::ofstream fout(filename);
    if (!fout) fail("ERROR: failed to write file " + filename);
    for (const auto& line : mol.header_lines) fout << line << '\n';
    fout << std::fixed << std::setprecision(8);
    for (const auto& atom : mol.atoms) {
        fout << std::setw(2) << std::left << atom.element << ' ' << std::right
             << std::setw(16) << atom.coord.x
             << std::setw(16) << atom.coord.y
             << std::setw(16) << atom.coord.z << '\n';
    }
    for (const auto& line : mol.footer_lines) fout << line << '\n';
}

void write_molecule(const std::string& src_filename, const std::string& ext, const Molecule& mol) {
    const std::string out = output_filename(src_filename, ext);
    if (ext == ".xyz") write_xyz(out, mol);
    else if (ext == ".gjf") write_gjf(out, mol);
    else fail("ERROR: unsupported output format: " + ext);
}

int main(int argc, char* argv[]) {
    try {
        if (!(argc == 3 || argc == 5)) {
            std::cerr << "\n ERROR in subroutine rmsd: wrong command line arguments!\n";
            std::cerr << " Format: ./calc_rmsd_xyz_lapack file1 file2 [range1] [range2]\n";
            std::cerr << "\n Example 1: ./calc_rmsd_xyz_lapack a.gjf b.gjf\n";
            std::cerr << " Example 2: ./calc_rmsd_xyz_lapack a.xyz b.xyz\n";
            std::cerr << " Example 3: ./calc_rmsd_xyz_lapack a.gjf b.xyz 1-5 2-6\n";
            return 1;
        }

        std::string fname1 = argv[1];
        std::string fname2 = argv[2];
        std::string ext1 = get_extension_lower(fname1);
        std::string ext2 = get_extension_lower(fname2);

        if (ext1 != ".gjf" && ext1 != ".xyz") fail("ERROR: unsupported file format for file1: " + fname1 + "\nSupported formats: .gjf, .xyz");
        if (ext2 != ".gjf" && ext2 != ".xyz") fail("ERROR: unsupported file format for file2: " + fname2 + "\nSupported formats: .gjf, .xyz");
        if (!fs::exists(fname1)) fail("ERROR in subroutine rmsd: file " + fname1 + " does not exist!");
        if (!fs::exists(fname2)) fail("ERROR in subroutine rmsd: file " + fname2 + " does not exist!");

        int idx1 = 0, idx2 = 0, idx3 = 0, idx4 = 0;
        if (argc == 5) {
            auto r1 = parse_range(argv[3], "range1");
            auto r2 = parse_range(argv[4], "range2");
            idx1 = r1.first; idx2 = r1.second; idx3 = r2.first; idx4 = r2.second;
        }

        if (idx1 < 0 || idx2 < 0 || idx3 < 0 || idx4 < 0 || idx1 > idx2 || idx3 > idx4) fail("ERROR in subroutine rmsd: input range not valid.");
        if ((idx2 - idx1) != (idx4 - idx3)) fail("ERROR in subroutine rmsd: number of atoms in two ranges are not equal.");

        Molecule fixed = read_molecule(fname1, ext1);
        Molecule moving = read_molecule(fname2, ext2);

        if (argc == 3) {
            idx1 = 1; idx2 = static_cast<int>(fixed.atoms.size());
            idx3 = 1; idx4 = static_cast<int>(moving.atoms.size());
        }

        if (static_cast<int>(fixed.atoms.size()) < idx2) fail("ERROR in subroutine rmsd_wrapper: error range1.");
        if (static_cast<int>(moving.atoms.size()) < idx4) fail("ERROR in subroutine rmsd_wrapper: error range2.");

        auto coor1 = coords_from_range(moving, idx3, idx4);
        auto coor2 = coords_from_range(fixed, idx1, idx2);
        AlignmentResult res = rmsd_align(coor1, coor2);

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "RMSD = " << std::setw(12) << res.rmsd << '\n';

        Mat3 RT = transpose(res.rotation);
        for (auto& atom : moving.atoms) {
            atom.coord += res.trans1;
            atom.coord = matvec(RT, atom.coord);
            atom.coord -= res.trans2;
        }

        write_molecule(fname2, ext2, moving);
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << ex.what() << '\n';
        return 1;
    }
}
