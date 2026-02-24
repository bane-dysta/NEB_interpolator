#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <map>
#include <filesystem>
#include <system_error>
#include <sys/wait.h>
#include <unistd.h>
#include <limits.h>

#include "internal_ic.h"

// 查找 calc_rmsd_xyz 可执行文件的路径
std::string findCalcRMSDExecutable(const std::string& default_path = "./calc_rmsd_xyz") {
    // 1. 首先检查是否在 PATH 中
    std::string command = "which calc_rmsd_xyz 2>/dev/null";
    FILE* pipe = popen(command.c_str(), "r");
    if (pipe != nullptr) {
        char buffer[PATH_MAX];
        if (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
            pclose(pipe);
            std::string path(buffer);
            // 移除末尾的换行符
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
    
    // 2. 检查 neb_interpolator 所在目录
    char exe_path[PATH_MAX];
    ssize_t len = readlink("/proc/self/exe", exe_path, sizeof(exe_path) - 1);
    if (len != -1) {
        exe_path[len] = '\0';
        std::string exe_dir(exe_path);
        size_t last_slash = exe_dir.find_last_of('/');
        if (last_slash != std::string::npos) {
            exe_dir = exe_dir.substr(0, last_slash + 1);
            std::string candidate = exe_dir + "calc_rmsd_xyz";
            // 检查文件是否存在且可执行
            if (access(candidate.c_str(), X_OK) == 0) {
                return candidate;
            }
        }
    }
    
    // 3. 如果都找不到，返回默认路径
    return default_path;
}

struct Atom {
    std::string symbol;
    double x, y, z;
    
    Atom() : symbol(""), x(0.0), y(0.0), z(0.0) {}
    Atom(const std::string& s, double x_, double y_, double z_) 
        : symbol(s), x(x_), y(y_), z(z_) {}
};

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
        
        file.ignore();
        std::getline(file, comment);
        
        atoms.clear();
        atoms.reserve(natoms);
        
        for (int i = 0; i < natoms; ++i) {
            std::string symbol;
            double x, y, z;
            if (!(file >> symbol >> x >> y >> z)) {
                std::cerr << "Error: Cannot read atom " << i+1 << " from " << filename << std::endl;
                return false;
            }
            atoms.emplace_back(symbol, x, y, z);
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
    
    // 计算RMSD
    double calculateRMSD(const Structure& other) const {
        if (atoms.size() != other.atoms.size()) return -1.0;
        
        double sum = 0.0;
        for (size_t i = 0; i < atoms.size(); ++i) {
            double dx = atoms[i].x - other.atoms[i].x;
            double dy = atoms[i].y - other.atoms[i].y;
            double dz = atoms[i].z - other.atoms[i].z;
            sum += dx*dx + dy*dy + dz*dz;
        }
        return std::sqrt(sum / atoms.size());
    }
};

class FortranRMSDAligner {
private:
    std::string rmsd_executable;
    
public:
    FortranRMSDAligner(const std::string& exec_path = "./calc_rmsd_xyz") 
        : rmsd_executable(exec_path) {}
    
    // 使用Fortran程序进行对齐
    bool alignStructures(const std::string& reference_file, const std::string& mobile_file) {
        // 构建命令
        std::string command = rmsd_executable + " " + reference_file + " " + mobile_file;
        
        std::cout << "  Running Fortran RMSD alignment: " << command << std::endl;
        
        // 执行Fortran程序
        int result = system(command.c_str());
        
        if (result != 0) {
            std::cerr << "  Error: Fortran RMSD alignment failed with exit code " << result << std::endl;
            return false;
        }
        
        // 检查输出文件是否存在
        std::string aligned_file = mobile_file;
        size_t pos = aligned_file.find(".xyz");
        if (pos != std::string::npos) {
            aligned_file = aligned_file.substr(0, pos) + "_new.xyz";
        } else {
            aligned_file += "_new.xyz";
        }
        
        std::ifstream check(aligned_file);
        if (!check.is_open()) {
            std::cerr << "  Error: Expected aligned file " << aligned_file << " was not created" << std::endl;
            return false;
        }
        check.close();
        
        std::cout << "  Alignment completed. Aligned structure saved as: " << aligned_file << std::endl;
        return true;
    }
    
    // 对齐并替换原文件
    bool alignAndReplace(const std::string& reference_file, const std::string& mobile_file) {
        if (!alignStructures(reference_file, mobile_file)) {
            return false;
        }
        
        // 生成对齐后文件名
        std::string aligned_file = mobile_file;
        size_t pos = aligned_file.find(".xyz");
        if (pos != std::string::npos) {
            aligned_file = aligned_file.substr(0, pos) + "_new.xyz";
        } else {
            aligned_file += "_new.xyz";
        }
        
        // 将对齐后的文件复制回原文件（避免依赖外部 cp 命令，且能正确处理空格路径）
        std::error_code ec;
        std::filesystem::copy_file(aligned_file, mobile_file,
                                   std::filesystem::copy_options::overwrite_existing, ec);
        if (ec) {
            std::cerr << "  Error: Failed to replace original file: " << ec.message() << std::endl;
            return false;
        }

        // 删除临时文件（避免依赖外部 rm 命令）
        ec.clear();
        std::filesystem::remove(aligned_file, ec);
        // Ignore removal errors (not fatal)

        return true;
    }
    
    void setExecutablePath(const std::string& path) {
        rmsd_executable = path;
    }
};

// ============================================================
// External engine mode
// ============================================================
//
// This program can optionally call an external "engine" every NEB cycle
// to obtain physical gradients/forces for each intermediate image.
//
// Interface (fixed-format text files):
// - Input:  blocks of XYZ (type=xyz)
// - Output: blocks of gradients (type=gradient) or forces (type=force)
//
// Example output-gradient file (two images):
//   3
//   image=1 units=AU type=gradient
//   O   0.01  -0.02  0.00
//   H  -0.01   0.01  0.00
//   H   0.00   0.01  0.00
//   3
//   image=2 units=AU type=gradient
//   ...
//

static inline std::string trim(const std::string& s) {
    const char* ws = " \t\r\n";
    const size_t b = s.find_first_not_of(ws);
    if (b == std::string::npos) return "";
    const size_t e = s.find_last_not_of(ws);
    return s.substr(b, e - b + 1);
}

// Very small shell-quoting helper for /bin/sh -c
static inline std::string shellQuote(const std::string& s) {
    // wrap in single quotes; escape existing single quotes
    std::string out;
    out.reserve(s.size() + 2);
    out.push_back('\'');
    for (char c : s) {
        if (c == '\'') {
            out += "'\\''"; // close, escape, reopen
        } else {
            out.push_back(c);
        }
    }
    out.push_back('\'');
    return out;
}

static inline std::map<std::string, std::string> parseHeaderKV(const std::string& header_line) {
    std::map<std::string, std::string> kv;
    std::istringstream iss(header_line);
    std::string tok;
    while (iss >> tok) {
        const size_t eq = tok.find('=');
        if (eq == std::string::npos) continue;
        std::string k = tok.substr(0, eq);
        std::string v = tok.substr(eq + 1);
        kv[trim(k)] = trim(v);
    }
    return kv;
}

struct ExternalEngineConfig {
    bool enabled = false;
    std::string cmd;                         // command line
    std::string in_file = "neb_engine_in.dat";
    std::string out_file = "neb_engine_out.dat";
    std::string units = "Angstrom";          // only a label written to input/output headers
    bool keep_cycle_files = false;           // keep _cycle#### files
    bool output_is_force = false;            // default: output is gradient; if true, interpret as force
    double spring_constant = 1.0;            // NEB spring constant
    int run_every = 1;                       // run engine every N cycles
};

static inline std::string addCycleSuffix(const std::string& filename, int cycle) {
    // insert _cycle#### before the last '.' (or append)
    std::ostringstream oss;
    oss << "_cycle" << std::setw(4) << std::setfill('0') << cycle;
    const std::string suf = oss.str();

    const size_t dot = filename.find_last_of('.');
    if (dot == std::string::npos) {
        return filename + suf;
    }
    return filename.substr(0, dot) + suf + filename.substr(dot);
}

static bool writeEngineInputXYZ(const std::string& path,
                                const std::vector<Structure>& images,
                                const std::string& units,
                                int cycle,
                                std::string* err) {
    std::ofstream out(path);
    if (!out.is_open()) {
        if (err) *err = "Cannot create engine input file: " + path;
        return false;
    }
    out << std::fixed << std::setprecision(8);
    for (size_t img = 0; img < images.size(); ++img) {
        const Structure& s = images[img];
        out << s.atoms.size() << "\n";
        out << "image=" << (img + 1) << " units=" << units << " type=xyz cycle=" << cycle << "\n";
        for (const auto& atom : s.atoms) {
            out << std::setw(2) << std::left << atom.symbol << " "
                << std::setw(15) << atom.x << " "
                << std::setw(15) << atom.y << " "
                << std::setw(15) << atom.z << "\n";
        }
    }
    return true;
}

struct EngineVectorBlock {
    int natoms = 0;
    int image_index = -1; // 1-based
    std::string units;
    std::string type;
    std::vector<std::string> symbols;
    std::vector<double> vec; // length 3*natoms
};

static bool readEngineVectorBlocks(const std::string& path,
                                   std::vector<EngineVectorBlock>& blocks,
                                   std::string* err) {
    std::ifstream in(path);
    if (!in.is_open()) {
        if (err) *err = "Cannot open engine output file: " + path;
        return false;
    }

    blocks.clear();
    while (true) {
        int natoms = 0;
        if (!(in >> natoms)) {
            break; // EOF
        }
        std::string dummy;
        std::getline(in, dummy); // consume rest of line

        std::string header;
        if (!std::getline(in, header)) {
            if (err) *err = "Unexpected EOF while reading header line in: " + path;
            return false;
        }
        header = trim(header);

        EngineVectorBlock b;
        b.natoms = natoms;
        auto kv = parseHeaderKV(header);
        if (kv.count("image")) {
            try {
                b.image_index = std::stoi(kv["image"]);
            } catch (...) {
                b.image_index = -1;
            }
        }
        if (kv.count("units")) b.units = kv["units"];
        if (kv.count("type")) b.type = kv["type"];

        b.symbols.resize(natoms);
        b.vec.resize(3 * static_cast<size_t>(natoms), 0.0);
        for (int i = 0; i < natoms; ++i) {
            std::string sym;
            double a = 0.0, bb = 0.0, c = 0.0;
            if (!(in >> sym >> a >> bb >> c)) {
                if (err) *err = "Unexpected EOF while reading vectors in: " + path;
                return false;
            }
            b.symbols[i] = sym;
            b.vec[3 * i + 0] = a;
            b.vec[3 * i + 1] = bb;
            b.vec[3 * i + 2] = c;
        }
        blocks.push_back(std::move(b));
    }

    if (blocks.empty()) {
        if (err) *err = "Engine output file is empty or unreadable: " + path;
        return false;
    }
    return true;
}

static std::string replaceAll(std::string s, const std::string& from, const std::string& to) {
    if (from.empty()) return s;
    size_t pos = 0;
    while ((pos = s.find(from, pos)) != std::string::npos) {
        s.replace(pos, from.size(), to);
        pos += to.size();
    }
    return s;
}

static bool runExternalEngineOnce(const ExternalEngineConfig& cfg,
                                  int cycle,
                                  const std::string& in_path,
                                  const std::string& out_path,
                                  int natoms,
                                  int nimages,
                                  std::string* err) {
    if (cfg.cmd.empty()) {
        if (err) *err = "External engine enabled but --engine-cmd is empty";
        return false;
    }

    // Provide context via environment variables too.
    setenv("NEB_ENGINE_IN", in_path.c_str(), 1);
    setenv("NEB_ENGINE_OUT", out_path.c_str(), 1);
    {
        std::ostringstream oss;
        oss << cycle;
        setenv("NEB_CYCLE", oss.str().c_str(), 1);
    }
    {
        std::ostringstream oss;
        oss << natoms;
        setenv("NEB_NATOMS", oss.str().c_str(), 1);
    }
    {
        std::ostringstream oss;
        oss << nimages;
        setenv("NEB_NIMAGES", oss.str().c_str(), 1);
    }

    // Build command line
    std::string cmdline = cfg.cmd;
    const bool has_in = (cmdline.find("{in}") != std::string::npos);
    const bool has_out = (cmdline.find("{out}") != std::string::npos);
    const bool has_cycle = (cmdline.find("{cycle}") != std::string::npos);

    cmdline = replaceAll(cmdline, "{in}", shellQuote(in_path));
    cmdline = replaceAll(cmdline, "{out}", shellQuote(out_path));
    if (has_cycle) {
        cmdline = replaceAll(cmdline, "{cycle}", std::to_string(cycle));
    }

    // If user didn't use placeholders, append <in> <out>
    if (!has_in && !has_out) {
        cmdline += " " + shellQuote(in_path) + " " + shellQuote(out_path);
    }

    std::cout << "  [engine] cycle " << cycle << ": " << cmdline << std::endl;
    int ret = system(cmdline.c_str());
    if (ret != 0) {
        if (err) {
            std::ostringstream oss;
            oss << "External engine command failed (exit code " << ret << ")";
            *err = oss.str();
        }
        return false;
    }

    std::ifstream chk(out_path);
    if (!chk.is_open()) {
        if (err) *err = "External engine finished but did not create output file: " + out_path;
        return false;
    }
    return true;
}

static bool assembleEngineVectors(const std::vector<EngineVectorBlock>& blocks,
                                  const std::vector<Structure>& images,
                                  int num_images,
                                  std::vector<std::vector<double>>& vectors,
                                  std::string* err) {
    if (num_images <= 0) {
        if (err) *err = "num_images must be positive";
        return false;
    }
    if (images.size() != static_cast<size_t>(num_images)) {
        if (err) *err = "Internal error: images size mismatch";
        return false;
    }

    const int natoms = static_cast<int>(images[0].size());
    vectors.assign(num_images, std::vector<double>(3 * static_cast<size_t>(natoms), 0.0));
    std::vector<bool> filled(num_images, false);

    int next_seq = 1;
    for (const auto& b : blocks) {
        int img_idx = b.image_index;
        if (img_idx < 1 || img_idx > num_images) {
            // Assign sequentially to the next available image slot
            while (next_seq <= num_images && filled[next_seq - 1]) {
                ++next_seq;
            }
            if (next_seq > num_images) {
                if (err) *err = "Too many blocks in engine output";
                return false;
            }
            img_idx = next_seq;
        }

        const int slot = img_idx - 1;
        if (filled[slot]) {
            if (err) {
                std::ostringstream oss;
                oss << "Duplicate engine block for image=" << img_idx;
                *err = oss.str();
            }
            return false;
        }

        if (b.natoms != natoms) {
            if (err) {
                std::ostringstream oss;
                oss << "Engine block natoms mismatch for image=" << img_idx
                    << " (got " << b.natoms << ", expected " << natoms << ")";
                *err = oss.str();
            }
            return false;
        }

        // Validate atom symbols to avoid silent reorder issues
        for (int i = 0; i < natoms; ++i) {
            if (i < static_cast<int>(b.symbols.size())) {
                if (b.symbols[i] != images[slot].atoms[i].symbol) {
                    if (err) {
                        std::ostringstream oss;
                        oss << "Atom symbol mismatch in engine output for image=" << img_idx
                            << " atom=" << (i + 1)
                            << " (got '" << b.symbols[i] << "', expected '" << images[slot].atoms[i].symbol << "')";
                        *err = oss.str();
                    }
                    return false;
                }
            }
        }

        if (b.vec.size() != vectors[slot].size()) {
            if (err) {
                std::ostringstream oss;
                oss << "Engine vector length mismatch for image=" << img_idx;
                *err = oss.str();
            }
            return false;
        }

        vectors[slot] = b.vec;
        filled[slot] = true;
    }

    for (int i = 0; i < num_images; ++i) {
        if (!filled[i]) {
            if (err) {
                std::ostringstream oss;
                oss << "Missing engine output block for image=" << (i + 1);
                *err = oss.str();
            }
            return false;
        }
    }

    return true;
}

class NEBInterpolator {
private:
    Structure initial, final;
    std::vector<Structure> images;
    int num_images;
    double step_init;
    double convergence_threshold;
    int max_iterations;
    bool use_alignment;
    FortranRMSDAligner aligner;

    // External engine mode (optional, used by NEB)
    ExternalEngineConfig engine_cfg;
    bool use_external_engine = false;
    std::vector<std::vector<double>> last_engine_vec; // [img][3N]
    bool have_last_engine_vec = false;

    // IIC / DM interpolation options (dependency-free, GSM-inspired)
    ICInterp::IICOptions iic_options;
    ICInterp::DMOptions dm_options;
    bool enable_dm_fallback;

    static double distance(const Atom& a1, const Atom& a2) {
        double dx = a1.x - a2.x;
        double dy = a1.y - a2.y;
        double dz = a1.z - a2.z;
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }
    
    void computeDistanceMatrix(const Structure& structure, std::vector<std::vector<double>>& dist_matrix) {
        size_t n = structure.size();
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                if (i == j) {
                    dist_matrix[i][j] = 1e10;
                } else {
                    dist_matrix[i][j] = distance(structure.atoms[i], structure.atoms[j]);
                }
            }
        }
    }
    
public:
    NEBInterpolator(int num_img = 5, double step = 0.0001, double conv_thresh = 0.01, 
                   int max_iter = 10000, bool align = true, const std::string& rmsd_exec = "./calc_rmsd_xyz") 
        : num_images(num_img), step_init(step), convergence_threshold(conv_thresh), 
          max_iterations(max_iter), use_alignment(align), aligner(rmsd_exec), enable_dm_fallback(true) {}

    void setExternalEngine(const ExternalEngineConfig& cfg) {
        engine_cfg = cfg;
        use_external_engine = cfg.enabled;
        have_last_engine_vec = false;
        last_engine_vec.clear();
    }
    
    bool setStructures(const std::string& initial_file, const std::string& final_file) {
        if (!initial.readXYZ(initial_file)) return false;
        if (!final.readXYZ(final_file)) return false;
        
        if (!initial.isCompatible(final)) {
            std::cerr << "Error: Initial and final structures are not compatible" << std::endl;
            return false;
        }
        
        // 如果启用对齐，使用Fortran程序进行对齐
        if (use_alignment) {
            std::cout << "Aligning structures using Fortran RMSD..." << std::endl;
            double rmsd_before = initial.calculateRMSD(final);
            std::cout << "RMSD before alignment: " << std::fixed << std::setprecision(6) << rmsd_before << std::endl;
            
            // 创建临时文件用于对齐
            std::string temp_final = "temp_final.xyz";
            if (!final.writeXYZ(temp_final)) {
                std::cerr << "Error: Cannot create temporary file for alignment" << std::endl;
                return false;
            }
            
            // 执行对齐
            if (aligner.alignAndReplace(initial_file, temp_final)) {
                // 读取对齐后的结构
                Structure final_aligned;
                if (final_aligned.readXYZ(temp_final)) {
                    double rmsd_after = initial.calculateRMSD(final_aligned);
                    std::cout << "RMSD after alignment: " << std::fixed << std::setprecision(6) << rmsd_after << std::endl;
                    
                    if (rmsd_after < rmsd_before) {
                        final = final_aligned;
                        std::cout << "Alignment improved RMSD, using aligned final structure." << std::endl;
                    } else {
                        std::cout << "Alignment did not improve RMSD, using original final structure." << std::endl;
                    }
                } else {
                    std::cerr << "Warning: Cannot read aligned structure, using original" << std::endl;
                }
            } else {
                std::cerr << "Warning: Alignment failed, using original structures" << std::endl;
            }
            
            // 清理临时文件
            std::error_code ec;
            std::filesystem::remove(temp_final, ec);
            // Ignore removal errors
        }
        
        return true;
    }
    
    void performLIIC() {
        std::cout << "  Initializing intermediate images..." << std::endl;
        
        images.clear();
        images.resize(num_images);
        
        for (int img = 0; img < num_images; ++img) {
            double factor = static_cast<double>(img + 1) / (num_images + 1);
            
            images[img].atoms.resize(initial.size());
            images[img].comment = "LIIC intermediate image " + std::to_string(img + 1);
            
            for (size_t i = 0; i < initial.size(); ++i) {
                images[img].atoms[i].symbol = initial.atoms[i].symbol;
                images[img].atoms[i].x = initial.atoms[i].x + factor * (final.atoms[i].x - initial.atoms[i].x);
                images[img].atoms[i].y = initial.atoms[i].y + factor * (final.atoms[i].y - initial.atoms[i].y);
                images[img].atoms[i].z = initial.atoms[i].z + factor * (final.atoms[i].z - initial.atoms[i].z);
            }
        }
    }


    // GSM-inspired internal-coordinate interpolation (IIC) with optional DM fallback.
    // Returns true if IIC succeeded (or DM fallback succeeded). If both fail, falls back to LIIC and returns false.
    bool performIIC() {
        std::cout << "  Initializing intermediate images using internal-coordinate interpolation (IIC)..." << std::endl;

        const size_t n = initial.size();
        std::vector<std::string> symbols(n);
        std::vector<double> x0(3*n, 0.0), x1(3*n, 0.0);

        for (size_t i = 0; i < n; ++i) {
            symbols[i] = initial.atoms[i].symbol;
            x0[3*i+0] = initial.atoms[i].x;
            x0[3*i+1] = initial.atoms[i].y;
            x0[3*i+2] = initial.atoms[i].z;

            x1[3*i+0] = final.atoms[i].x;
            x1[3*i+1] = final.atoms[i].y;
            x1[3*i+2] = final.atoms[i].z;
        }

        std::vector<std::vector<double>> imgs_xyz;
        std::string err;
        bool ok = ICInterp::interpolate_iic(symbols, x0, x1, num_images, imgs_xyz, iic_options, &err);

        if (!ok) {
            std::cerr << "  IIC failed: " << err << std::endl;

            if (enable_dm_fallback) {
                std::cout << "  Falling back to distance-matrix interpolation (DM)..." << std::endl;
                bool ok_dm = ICInterp::interpolate_dm(symbols, x0, x1, num_images, imgs_xyz, dm_options, &err);
                if (!ok_dm) {
                    std::cerr << "  DM fallback failed: " << err << std::endl;
                    std::cout << "  Falling back to Cartesian LIIC..." << std::endl;
                    performLIIC();
                    return false;
                }
                // Fill images from DM result
                images.clear();
                images.resize(num_images);
                for (int img = 0; img < num_images; ++img) {
                    images[img].atoms.resize(n);
                    images[img].comment = "DM (fallback) intermediate image " + std::to_string(img + 1);
                    for (size_t i = 0; i < n; ++i) {
                        images[img].atoms[i].symbol = symbols[i];
                        images[img].atoms[i].x = imgs_xyz[img][3*i+0];
                        images[img].atoms[i].y = imgs_xyz[img][3*i+1];
                        images[img].atoms[i].z = imgs_xyz[img][3*i+2];
                    }
                }
                return true;
            } else {
                std::cout << "  DM fallback disabled. Falling back to Cartesian LIIC..." << std::endl;
                performLIIC();
                return false;
            }
        }

        // Fill images from IIC result
        images.clear();
        images.resize(num_images);
        for (int img = 0; img < num_images; ++img) {
            images[img].atoms.resize(n);
            images[img].comment = "IIC intermediate image " + std::to_string(img + 1);
            for (size_t i = 0; i < n; ++i) {
                images[img].atoms[i].symbol = symbols[i];
                images[img].atoms[i].x = imgs_xyz[img][3*i+0];
                images[img].atoms[i].y = imgs_xyz[img][3*i+1];
                images[img].atoms[i].z = imgs_xyz[img][3*i+2];
            }
        }
        return true;
    }

    // Distance-matrix interpolation (DM). Returns true if succeeded, else falls back to LIIC and returns false.
    bool performDM() {
        std::cout << "  Initializing intermediate images using distance-matrix interpolation (DM)..." << std::endl;

        const size_t n = initial.size();
        std::vector<std::string> symbols(n);
        std::vector<double> x0(3*n, 0.0), x1(3*n, 0.0);

        for (size_t i = 0; i < n; ++i) {
            symbols[i] = initial.atoms[i].symbol;
            x0[3*i+0] = initial.atoms[i].x;
            x0[3*i+1] = initial.atoms[i].y;
            x0[3*i+2] = initial.atoms[i].z;

            x1[3*i+0] = final.atoms[i].x;
            x1[3*i+1] = final.atoms[i].y;
            x1[3*i+2] = final.atoms[i].z;
        }

        std::vector<std::vector<double>> imgs_xyz;
        std::string err;
        bool ok = ICInterp::interpolate_dm(symbols, x0, x1, num_images, imgs_xyz, dm_options, &err);

        if (!ok) {
            std::cerr << "  DM failed: " << err << std::endl;
            std::cout << "  Falling back to Cartesian LIIC..." << std::endl;
            performLIIC();
            return false;
        }

        images.clear();
        images.resize(num_images);
        for (int img = 0; img < num_images; ++img) {
            images[img].atoms.resize(n);
            images[img].comment = "DM intermediate image " + std::to_string(img + 1);
            for (size_t i = 0; i < n; ++i) {
                images[img].atoms[i].symbol = symbols[i];
                images[img].atoms[i].x = imgs_xyz[img][3*i+0];
                images[img].atoms[i].y = imgs_xyz[img][3*i+1];
                images[img].atoms[i].z = imgs_xyz[img][3*i+2];
            }
        }
        return true;
    }

    bool performNEB(bool init_with_iic = false) {
        std::cout << "Performing Nudged Elastic Band (NEB) interpolation..." << std::endl;
        
        if (init_with_iic) {
            performIIC();
        } else {
            performLIIC();
        }
        
        // Spring constant: keep old default (1.0) for internal mode,
        // but allow overriding in external-engine mode.
        const double spring_constant = use_external_engine ? engine_cfg.spring_constant : 1.0;
        
        for (int cycle = 0; cycle < max_iterations; ++cycle) {
            double max_force = 0.0;
            std::vector<Structure> forces(num_images);

            for (int img = 0; img < num_images; ++img) {
                forces[img].atoms.resize(images[img].size());
                for (size_t i = 0; i < images[img].size(); ++i) {
                    forces[img].atoms[i] = Atom(images[img].atoms[i].symbol, 0.0, 0.0, 0.0);
                }
            }

            if (!use_external_engine) {
                // ============================================================
                // Original (internal) NEB-like force field (distance preservation)
                // ============================================================
                for (int img = 0; img < num_images; ++img) {
                    Structure& current = images[img];
                    Structure* prev = (img == 0) ? &initial : &images[img-1];
                    Structure* next = (img == num_images-1) ? &final : &images[img+1];

                    // 计算切线向量（legacy: per-atom)
                    std::vector<double> tangent_x(current.size(), 0.0);
                    std::vector<double> tangent_y(current.size(), 0.0);
                    std::vector<double> tangent_z(current.size(), 0.0);

                    for (size_t i = 0; i < current.size(); ++i) {
                        double dx_prev = current.atoms[i].x - prev->atoms[i].x;
                        double dy_prev = current.atoms[i].y - prev->atoms[i].y;
                        double dz_prev = current.atoms[i].z - prev->atoms[i].z;

                        double dx_next = next->atoms[i].x - current.atoms[i].x;
                        double dy_next = next->atoms[i].y - current.atoms[i].y;
                        double dz_next = next->atoms[i].z - current.atoms[i].z;

                        // Tangent direction should follow the path: R_{i+1} - R_{i-1}
                        // Here dx_next = R_{i+1}-R_i and dx_prev = R_i-R_{i-1}, so
                        // R_{i+1}-R_{i-1} = dx_next + dx_prev.
                        tangent_x[i] = dx_next + dx_prev;
                        tangent_y[i] = dy_next + dy_prev;
                        tangent_z[i] = dz_next + dz_prev;

                        double norm = std::sqrt(tangent_x[i]*tangent_x[i] + tangent_y[i]*tangent_y[i] + tangent_z[i]*tangent_z[i]);
                        if (norm > 1e-12) {
                            tangent_x[i] /= norm;
                            tangent_y[i] /= norm;
                            tangent_z[i] /= norm;
                        }
                    }

                    // 弹簧力计算
                    for (size_t i = 0; i < current.size(); ++i) {
                        double spring_force_x = spring_constant * ((next->atoms[i].x - current.atoms[i].x) - (current.atoms[i].x - prev->atoms[i].x));
                        double spring_force_y = spring_constant * ((next->atoms[i].y - current.atoms[i].y) - (current.atoms[i].y - prev->atoms[i].y));
                        double spring_force_z = spring_constant * ((next->atoms[i].z - current.atoms[i].z) - (current.atoms[i].z - prev->atoms[i].z));

                        double proj = spring_force_x * tangent_x[i] + spring_force_y * tangent_y[i] + spring_force_z * tangent_z[i];

                        forces[img].atoms[i].x += proj * tangent_x[i];
                        forces[img].atoms[i].y += proj * tangent_y[i];
                        forces[img].atoms[i].z += proj * tangent_z[i];
                    }

                    // 垂直力计算（简化的距离保持）
                    size_t n = current.size();
                    std::vector<std::vector<double>> current_dist(n, std::vector<double>(n));
                    std::vector<std::vector<double>> target_dist(n, std::vector<double>(n));

                    computeDistanceMatrix(current, current_dist);

                    double factor = static_cast<double>(img + 1) / (num_images + 1);
                    std::vector<std::vector<double>> initial_dist(n, std::vector<double>(n));
                    std::vector<std::vector<double>> final_dist(n, std::vector<double>(n));
                    computeDistanceMatrix(initial, initial_dist);
                    computeDistanceMatrix(final, final_dist);

                    for (size_t i = 0; i < n; ++i) {
                        for (size_t j = 0; j < n; ++j) {
                            target_dist[i][j] = initial_dist[i][j] + factor * (final_dist[i][j] - initial_dist[i][j]);
                        }
                    }

                    for (size_t i = 0; i < n; ++i) {
                        double perp_force_x = 0.0, perp_force_y = 0.0, perp_force_z = 0.0;

                        for (size_t j = 0; j < n; ++j) {
                            if (i != j && current_dist[i][j] > 1e-12) {
                                double dist_diff = target_dist[i][j] - current_dist[i][j];
                                double pos_diff_x = current.atoms[i].x - current.atoms[j].x;
                                double pos_diff_y = current.atoms[i].y - current.atoms[j].y;
                                double pos_diff_z = current.atoms[i].z - current.atoms[j].z;

                                double force_magnitude = 2.0 * dist_diff / current_dist[i][j];

                                perp_force_x += force_magnitude * pos_diff_x;
                                perp_force_y += force_magnitude * pos_diff_y;
                                perp_force_z += force_magnitude * pos_diff_z;
                            }
                        }

                        double proj = perp_force_x * tangent_x[i] + perp_force_y * tangent_y[i] + perp_force_z * tangent_z[i];

                        forces[img].atoms[i].x += perp_force_x - proj * tangent_x[i];
                        forces[img].atoms[i].y += perp_force_y - proj * tangent_y[i];
                        forces[img].atoms[i].z += perp_force_z - proj * tangent_z[i];
                    }

                    for (size_t i = 0; i < current.size(); ++i) {
                        double force_mag = std::sqrt(forces[img].atoms[i].x * forces[img].atoms[i].x +
                                                     forces[img].atoms[i].y * forces[img].atoms[i].y +
                                                     forces[img].atoms[i].z * forces[img].atoms[i].z);
                        max_force = std::max(max_force, force_mag);
                    }
                }
            } else {
                // ============================================================
                // External engine NEB: F = F_real_perp + F_spring_parallel
                // ============================================================

                // 1) Run engine (possibly every N cycles) to get gradients/forces
                const int natoms = static_cast<int>(initial.size());
                const bool need_run = (!have_last_engine_vec) || (engine_cfg.run_every <= 1) || (cycle % engine_cfg.run_every == 0);
                if (need_run) {
                    std::string in_path = engine_cfg.in_file;
                    std::string out_path = engine_cfg.out_file;
                    if (engine_cfg.keep_cycle_files) {
                        in_path = addCycleSuffix(in_path, cycle);
                        out_path = addCycleSuffix(out_path, cycle);
                    }

                    std::string err;
                    if (!writeEngineInputXYZ(in_path, images, engine_cfg.units, cycle, &err)) {
                        std::cerr << "  [engine] Error: " << err << std::endl;
                        return false;
                    }

                    if (!runExternalEngineOnce(engine_cfg, cycle, in_path, out_path, natoms, num_images, &err)) {
                        std::cerr << "  [engine] Error: " << err << std::endl;
                        return false;
                    }

                    std::vector<EngineVectorBlock> blocks;
                    if (!readEngineVectorBlocks(out_path, blocks, &err)) {
                        std::cerr << "  [engine] Error: " << err << std::endl;
                        return false;
                    }

                    if (!assembleEngineVectors(blocks, images, num_images, last_engine_vec, &err)) {
                        std::cerr << "  [engine] Error: " << err << std::endl;
                        return false;
                    }
                    have_last_engine_vec = true;
                }

                // 2) NEB force assembly using last_engine_vec
                for (int img = 0; img < num_images; ++img) {
                    Structure& current = images[img];
                    const Structure* prev = (img == 0) ? &initial : &images[img - 1];
                    const Structure* next = (img == num_images - 1) ? &final : &images[img + 1];

                    const size_t n = current.size();
                    std::vector<double> t(3 * n, 0.0);

                    // Tangent vector t = R_{i+1} - R_{i-1}
                    for (size_t i = 0; i < n; ++i) {
                        t[3 * i + 0] = next->atoms[i].x - prev->atoms[i].x;
                        t[3 * i + 1] = next->atoms[i].y - prev->atoms[i].y;
                        t[3 * i + 2] = next->atoms[i].z - prev->atoms[i].z;
                    }

                    double tnorm2 = 0.0;
                    for (double v : t) tnorm2 += v * v;
                    const double tnorm = std::sqrt(tnorm2);
                    std::vector<double> t_unit(3 * n, 0.0);
                    if (tnorm > 1e-16) {
                        for (size_t k = 0; k < t_unit.size(); ++k) t_unit[k] = t[k] / tnorm;
                    }

                    // Spring force along tangent (standard NEB spacing term)
                    double d_next2 = 0.0, d_prev2 = 0.0;
                    for (size_t i = 0; i < n; ++i) {
                        const double dxn = next->atoms[i].x - current.atoms[i].x;
                        const double dyn = next->atoms[i].y - current.atoms[i].y;
                        const double dzn = next->atoms[i].z - current.atoms[i].z;
                        d_next2 += dxn * dxn + dyn * dyn + dzn * dzn;

                        const double dxp = current.atoms[i].x - prev->atoms[i].x;
                        const double dyp = current.atoms[i].y - prev->atoms[i].y;
                        const double dzp = current.atoms[i].z - prev->atoms[i].z;
                        d_prev2 += dxp * dxp + dyp * dyp + dzp * dzp;
                    }
                    const double d_next = std::sqrt(d_next2);
                    const double d_prev = std::sqrt(d_prev2);
                    const double spring_scalar = (tnorm > 1e-16) ? (spring_constant * (d_next - d_prev)) : 0.0;

                    // Real force from engine output (gradient -> force = -grad; force -> as-is)
                    const std::vector<double>& vec = last_engine_vec[img];
                    double dot_ft = 0.0;
                    for (size_t i = 0; i < n; ++i) {
                        const double vx = vec[3 * i + 0];
                        const double vy = vec[3 * i + 1];
                        const double vz = vec[3 * i + 2];

                        const double fx = engine_cfg.output_is_force ? vx : -vx;
                        const double fy = engine_cfg.output_is_force ? vy : -vy;
                        const double fz = engine_cfg.output_is_force ? vz : -vz;

                        dot_ft += fx * t_unit[3 * i + 0] + fy * t_unit[3 * i + 1] + fz * t_unit[3 * i + 2];
                    }

                    for (size_t i = 0; i < n; ++i) {
                        const double vx = vec[3 * i + 0];
                        const double vy = vec[3 * i + 1];
                        const double vz = vec[3 * i + 2];

                        // Real force
                        const double fx = engine_cfg.output_is_force ? vx : -vx;
                        const double fy = engine_cfg.output_is_force ? vy : -vy;
                        const double fz = engine_cfg.output_is_force ? vz : -vz;

                        // Perpendicular component: F_perp = F - (F·t) t
                        const double tx = t_unit[3 * i + 0];
                        const double ty = t_unit[3 * i + 1];
                        const double tz = t_unit[3 * i + 2];
                        const double fx_perp = fx - dot_ft * tx;
                        const double fy_perp = fy - dot_ft * ty;
                        const double fz_perp = fz - dot_ft * tz;

                        // Add spring force along tangent
                        const double fx_tot = fx_perp + spring_scalar * tx;
                        const double fy_tot = fy_perp + spring_scalar * ty;
                        const double fz_tot = fz_perp + spring_scalar * tz;

                        forces[img].atoms[i].x = fx_tot;
                        forces[img].atoms[i].y = fy_tot;
                        forces[img].atoms[i].z = fz_tot;

                        const double force_mag = std::sqrt(fx_tot * fx_tot + fy_tot * fy_tot + fz_tot * fz_tot);
                        max_force = std::max(max_force, force_mag);
                    }
                }
            }
            
            // 更新位置
            for (int img = 0; img < num_images; ++img) {
                for (size_t i = 0; i < images[img].size(); ++i) {
                    images[img].atoms[i].x += step_init * forces[img].atoms[i].x;
                    images[img].atoms[i].y += step_init * forces[img].atoms[i].y;
                    images[img].atoms[i].z += step_init * forces[img].atoms[i].z;
                }
            }
            
            if (cycle % 100 == 0 || cycle < 10) {
                std::cout << "  NEB cycle " << cycle << ", max force = " << std::scientific << std::setprecision(4) << max_force << std::endl;
            }
            
            if (max_force < convergence_threshold) {
                std::cout << "  NEB converged after " << cycle << " cycles!" << std::endl;
                break;
            }
        }

        return true;
    }
    
    bool writeResults(const std::string& prefix = "", bool multiframe = false) {
        if (multiframe) {
            return writeMultiframeXYZ(prefix);
        } else {
            return writeSeparateXYZ(prefix);
        }
    }
    
    void setAlignment(bool enable) {
        use_alignment = enable;
    }
    
    void setRMSDExecutable(const std::string& path) {
        aligner.setExecutablePath(path);
    }
    
private:
    bool writeSeparateXYZ(const std::string& prefix = "") {
        std::string filename = prefix + "00.xyz";
        if (!initial.writeXYZ(filename, "Initial structure")) return false;
        std::cout << "Wrote " << filename << std::endl;
        
        for (int i = 0; i < num_images; ++i) {
            filename = prefix + (i + 1 < 10 ? "0" : "") + std::to_string(i + 1) + ".xyz";
            if (!images[i].writeXYZ(filename)) return false;
            std::cout << "Wrote " << filename << std::endl;
        }
        
        filename = prefix + (num_images + 1 < 10 ? "0" : "") + std::to_string(num_images + 1) + ".xyz";
        if (!final.writeXYZ(filename, "Final structure")) return false;
        std::cout << "Wrote " << filename << std::endl;
        
        return true;
    }
    
    bool writeMultiframeXYZ(const std::string& prefix = "") {
        std::string filename = prefix + "trajectory.xyz";
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Cannot create file " << filename << std::endl;
            return false;
        }
        
        file << std::fixed << std::setprecision(8);
        
        file << initial.atoms.size() << "\n";
        file << "Frame 0: Initial structure\n";
        for (const auto& atom : initial.atoms) {
            file << std::setw(2) << std::left << atom.symbol << " "
                 << std::setw(15) << atom.x << " "
                 << std::setw(15) << atom.y << " "
                 << std::setw(15) << atom.z << "\n";
        }
        
        for (int i = 0; i < num_images; ++i) {
            file << images[i].atoms.size() << "\n";
            file << "Frame " << (i + 1) << ": Intermediate image " << (i + 1) << "\n";
            for (const auto& atom : images[i].atoms) {
                file << std::setw(2) << std::left << atom.symbol << " "
                     << std::setw(15) << atom.x << " "
                     << std::setw(15) << atom.y << " "
                     << std::setw(15) << atom.z << "\n";
            }
        }
        
        file << final.atoms.size() << "\n";
        file << "Frame " << (num_images + 1) << ": Final structure\n";
        for (const auto& atom : final.atoms) {
            file << std::setw(2) << std::left << atom.symbol << " "
                 << std::setw(15) << atom.x << " "
                 << std::setw(15) << atom.y << " "
                 << std::setw(15) << atom.z << "\n";
        }
        
        std::cout << "Wrote multiframe trajectory to " << filename << std::endl;
        std::cout << "Total frames: " << (num_images + 2) << std::endl;
        
        return true;
    }

public:
    void setIICOptions(const ICInterp::IICOptions& opt) { iic_options = opt; }
    void setDMOptions(const ICInterp::DMOptions& opt) { dm_options = opt; }
    void setDMFallback(bool enable) { enable_dm_fallback = enable; }

    void setParameters(int num_img, double step, double conv_thresh, int max_iter) {
        num_images = num_img;
        step_init = step;
        convergence_threshold = conv_thresh;
        max_iterations = max_iter;
    }
};

void printUsage(const char* program_name) {
    std::cout << "Usage: " << program_name << " [options] <initial.xyz> <final.xyz>\n"
              << "\nOptions:\n"
              << "  -n, --nimages NUM        Number of intermediate images (default: 5)\n"
              << "  -m, --method METHOD      Method: liic | iic | dm | neb | neb-iic (default: neb)\n"
              << "  -p, --prefix PREFIX      Output filename prefix (default: empty)\n"
              << "  -o, --output MODE        Output mode: separate or multiframe (default: separate)\n"
              << "  -s, --step STEP          Initial step size for NEB optimization (default: 0.0001)\n"
              << "  -c, --conv THRESHOLD     Convergence threshold (default: 0.01)\n"
              << "  -i, --maxiter ITER       Maximum iterations for NEB (default: 10000)\n"
              << "  -a, --align              Enable structure alignment using Fortran RMSD (default: enabled)\n"
              << "  --no-align               Disable structure alignment\n"
              << "  -r, --rmsd-exec PATH     Path to Fortran RMSD executable (default: ./calc_rmsd_xyz)\n"
              << "  -h, --help               Show this help message\n"
              << "\nIIC options (used by -m iic or -m neb-iic):\n"
              << "  --bond-factor F          Bond cutoff factor (default: 1.25)\n"
              << "  --fd-step H              Finite-difference step for B matrix in IIC (default: 1e-4)\n"
              << "  --ev-thresh T            Eigenvalue threshold for DLC selection (default: 1e-3)\n"
              << "  --iic-maxiter N          Max back-transform iterations (default: 50)\n"
              << "  --iic-tol T              RMS primitive residual tolerance (default: 1e-4)\n"
              << "  --iic-damp D             Damping added to eigenvalues (default: 1e-8)\n"
              << "  --iic-max-step S         Max cartesian step per iteration (default: 0.20)\n"
              << "\nDM options (used as fallback when IIC fails, or -m dm):\n"
              << "  --dm-maxiter N           Max DM iterations (default: 800)\n"
              << "  --dm-step S              DM gradient step size (default: 5e-3)\n"
              << "  --dm-tol T               RMS distance-error tolerance (default: 1e-3)\n"
              << "  --dm-max-step S          Max cartesian step per iteration (default: 0.20)\n"
              << "  --no-dm-fallback         Disable DM fallback (fallback directly to LIIC)\n"
              << "\nExternal engine options (for NEB / NEB-IIC only):\n"
              << "  --engine-cmd CMD         Enable external-engine mode; run CMD each NEB cycle\n"
              << "                           If CMD contains {in}/{out}/{cycle}, they will be replaced.\n"
              << "                           Otherwise, the program appends: <infile> <outfile>\n"
              << "  --engine-in FILE         Engine input filename (default: neb_engine_in.dat)\n"
              << "  --engine-out FILE        Engine output filename (default: neb_engine_out.dat)\n"
              << "  --engine-units U         Units label written in engine I/O headers (default: Angstrom)\n"
              << "  --engine-vector TYPE     Interpret engine output as: gradient | force (default: gradient)\n"
              << "  --engine-spring K        Spring constant used in external-engine NEB (default: 1.0)\n"
              << "  --engine-every N         Run engine every N cycles (default: 1)\n"
              << "  --engine-keep-files      Keep per-cycle engine I/O files (adds _cycle#### suffix)\n"
              << "\nOutput modes:\n"
              << "  separate     Generate separate XYZ files (00.xyz, 01.xyz, ...)\n"
              << "  multiframe   Generate single trajectory.xyz with all frames\n"
              << "\nAlignment:\n"
              << "  This program can use a Fortran RMSD calculator for precise structure alignment.\n"
              << "  Make sure the calc_rmsd_xyz executable is compiled and accessible.\n"
              << "\nExamples:\n"
              << "  " << program_name << " -n 10 -m neb -p reaction_ initial.xyz final.xyz\n"
              << "  " << program_name << " -n 10 -m neb-iic --bond-factor 1.2 --fd-step 1e-4 --ev-thresh 1e-3 initial.xyz final.xyz\n"
              << "  " << program_name << " --no-align -o multiframe -n 5 -m iic initial.xyz final.xyz\n"
              << "  " << program_name << " -m dm -n 8 initial.xyz final.xyz\n"
              << "  " << program_name << " -m neb -n 8 --engine-cmd \"python engine.py {in} {out}\" initial.xyz final.xyz\n";
}

int main(int argc, char* argv[]) {
    std::string initial_file, final_file, prefix = "";
    std::string method = "neb";
    std::string output_mode = "separate";
    std::string rmsd_executable = findCalcRMSDExecutable();
    int num_images = 5;
    double step_size = 0.0001;
    double conv_threshold = 0.01;
    int max_iterations = 10000;
    bool use_alignment = true;

    // IIC/DM options (dependency-free)
    ICInterp::IICOptions iic_opt;
    ICInterp::DMOptions dm_opt;
    bool dm_fallback = true;

    // External engine options (optional)
    ExternalEngineConfig engine_cfg;
    
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printUsage(argv[0]);
            return 0;
        } else if ((strcmp(argv[i], "-n") == 0 || strcmp(argv[i], "--nimages") == 0) && i + 1 < argc) {
            num_images = std::atoi(argv[++i]);
            if (num_images <= 0) {
                std::cerr << "Error: Number of images must be positive" << std::endl;
                return 1;
            }
        } else if ((strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "--method") == 0) && i + 1 < argc) {
            method = argv[++i];
            if (method != "liic" && method != "iic" && method != "dm" && method != "neb" && method != "neb-iic") {
                std::cerr << "Error: Method must be \'liic\', \'iic\', \'dm\', \'neb\', or \'neb-iic\'" << std::endl;
                return 1;
            }
        } else if ((strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "--prefix") == 0) && i + 1 < argc) {
            prefix = argv[++i];
        } else if ((strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) && i + 1 < argc) {
            output_mode = argv[++i];
            if (output_mode != "separate" && output_mode != "multiframe") {
                std::cerr << "Error: Output mode must be 'separate' or 'multiframe'" << std::endl;
                return 1;
            }
        } else if ((strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--step") == 0) && i + 1 < argc) {
            step_size = std::atof(argv[++i]);
            if (step_size <= 0) {
                std::cerr << "Error: Step size must be positive" << std::endl;
                return 1;
            }
        } else if ((strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--conv") == 0) && i + 1 < argc) {
            conv_threshold = std::atof(argv[++i]);
            if (conv_threshold <= 0) {
                std::cerr << "Error: Convergence threshold must be positive" << std::endl;
                return 1;
            }
        } else if ((strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--maxiter") == 0) && i + 1 < argc) {
            max_iterations = std::atoi(argv[++i]);
            if (max_iterations <= 0) {
                std::cerr << "Error: Maximum iterations must be positive" << std::endl;
                return 1;
            }
        } else if (strcmp(argv[i], "-a") == 0 || strcmp(argv[i], "--align") == 0) {
            use_alignment = true;
        } else if (strcmp(argv[i], "--no-align") == 0) {
            use_alignment = false;
        } else if ((strcmp(argv[i], "-r") == 0 || strcmp(argv[i], "--rmsd-exec") == 0) && i + 1 < argc) {
            rmsd_executable = argv[++i];
        } else if (strcmp(argv[i], "--bond-factor") == 0 && i + 1 < argc) {
            iic_opt.bond_factor = std::atof(argv[++i]);
            if (iic_opt.bond_factor <= 0.0) {
                std::cerr << "Error: --bond-factor must be positive" << std::endl;
                return 1;
            }
        } else if (strcmp(argv[i], "--fd-step") == 0 && i + 1 < argc) {
            iic_opt.fd_step = std::atof(argv[++i]);
            if (iic_opt.fd_step <= 0.0) {
                std::cerr << "Error: --fd-step must be positive" << std::endl;
                return 1;
            }
        } else if (strcmp(argv[i], "--ev-thresh") == 0 && i + 1 < argc) {
            iic_opt.ev_thresh = std::atof(argv[++i]);
            if (iic_opt.ev_thresh < 0.0) {
                std::cerr << "Error: --ev-thresh must be non-negative" << std::endl;
                return 1;
            }
        } else if (strcmp(argv[i], "--iic-maxiter") == 0 && i + 1 < argc) {
            iic_opt.max_iter = std::atoi(argv[++i]);
            if (iic_opt.max_iter <= 0) {
                std::cerr << "Error: --iic-maxiter must be positive" << std::endl;
                return 1;
            }
        } else if (strcmp(argv[i], "--iic-tol") == 0 && i + 1 < argc) {
            iic_opt.tol = std::atof(argv[++i]);
            if (iic_opt.tol <= 0.0) {
                std::cerr << "Error: --iic-tol must be positive" << std::endl;
                return 1;
            }
        } else if (strcmp(argv[i], "--iic-damp") == 0 && i + 1 < argc) {
            iic_opt.damp = std::atof(argv[++i]);
            if (iic_opt.damp < 0.0) {
                std::cerr << "Error: --iic-damp must be non-negative" << std::endl;
                return 1;
            }
        } else if (strcmp(argv[i], "--iic-max-step") == 0 && i + 1 < argc) {
            iic_opt.max_cart_step = std::atof(argv[++i]);
            if (iic_opt.max_cart_step <= 0.0) {
                std::cerr << "Error: --iic-max-step must be positive" << std::endl;
                return 1;
            }
        } else if (strcmp(argv[i], "--dm-maxiter") == 0 && i + 1 < argc) {
            dm_opt.max_iter = std::atoi(argv[++i]);
            if (dm_opt.max_iter <= 0) {
                std::cerr << "Error: --dm-maxiter must be positive" << std::endl;
                return 1;
            }
        } else if (strcmp(argv[i], "--dm-step") == 0 && i + 1 < argc) {
            dm_opt.step = std::atof(argv[++i]);
            if (dm_opt.step <= 0.0) {
                std::cerr << "Error: --dm-step must be positive" << std::endl;
                return 1;
            }
        } else if (strcmp(argv[i], "--dm-tol") == 0 && i + 1 < argc) {
            dm_opt.tol = std::atof(argv[++i]);
            if (dm_opt.tol <= 0.0) {
                std::cerr << "Error: --dm-tol must be positive" << std::endl;
                return 1;
            }
        } else if (strcmp(argv[i], "--dm-max-step") == 0 && i + 1 < argc) {
            dm_opt.max_cart_step = std::atof(argv[++i]);
            if (dm_opt.max_cart_step <= 0.0) {
                std::cerr << "Error: --dm-max-step must be positive" << std::endl;
                return 1;
            }
        } else if (strcmp(argv[i], "--no-dm-fallback") == 0) {
            dm_fallback = false;
        } else if ((strcmp(argv[i], "--engine-cmd") == 0 || strcmp(argv[i], "--engine") == 0 || strcmp(argv[i], "--external-engine") == 0) && i + 1 < argc) {
            engine_cfg.enabled = true;
            engine_cfg.cmd = argv[++i];
        } else if (strcmp(argv[i], "--engine-in") == 0 && i + 1 < argc) {
            engine_cfg.in_file = argv[++i];
        } else if (strcmp(argv[i], "--engine-out") == 0 && i + 1 < argc) {
            engine_cfg.out_file = argv[++i];
        } else if (strcmp(argv[i], "--engine-units") == 0 && i + 1 < argc) {
            engine_cfg.units = argv[++i];
        } else if (strcmp(argv[i], "--engine-vector") == 0 && i + 1 < argc) {
            std::string v = argv[++i];
            if (v == "gradient" || v == "grad") {
                engine_cfg.output_is_force = false;
            } else if (v == "force" || v == "forces") {
                engine_cfg.output_is_force = true;
            } else {
                std::cerr << "Error: --engine-vector must be 'gradient' or 'force'" << std::endl;
                return 1;
            }
        } else if (strcmp(argv[i], "--engine-spring") == 0 && i + 1 < argc) {
            engine_cfg.spring_constant = std::atof(argv[++i]);
            if (engine_cfg.spring_constant <= 0.0) {
                std::cerr << "Error: --engine-spring must be positive" << std::endl;
                return 1;
            }
        } else if (strcmp(argv[i], "--engine-every") == 0 && i + 1 < argc) {
            engine_cfg.run_every = std::atoi(argv[++i]);
            if (engine_cfg.run_every <= 0) {
                std::cerr << "Error: --engine-every must be positive" << std::endl;
                return 1;
            }
        } else if (strcmp(argv[i], "--engine-keep-files") == 0) {
            engine_cfg.keep_cycle_files = true;
        } else if (argv[i][0] != '-') {
            if (initial_file.empty()) {
                initial_file = argv[i];
            } else if (final_file.empty()) {
                final_file = argv[i];
            } else {
                std::cerr << "Error: Too many input files specified" << std::endl;
                return 1;
            }
        } else {
            std::cerr << "Error: Unknown option " << argv[i] << std::endl;
            printUsage(argv[0]);
            return 1;
        }
    }
    
    if (initial_file.empty() || final_file.empty()) {
        std::cerr << "Error: Must specify both initial and final structure files" << std::endl;
        printUsage(argv[0]);
        return 1;
    }
    
    std::cout << "NEB/LIIC/IIC/DM Interpolation Program with Fortran RMSD Alignment\n"
              << "=========================================================\n"
              << "Method: " << method << "\n"
              << "Number of images: " << num_images << "\n"
              << "Structure alignment: " << (use_alignment ? "enabled (Fortran RMSD)" : "disabled") << "\n"
              << "RMSD executable: " << rmsd_executable << "\n"
              << "Initial structure: " << initial_file << "\n"
              << "Final structure: " << final_file << "\n"
              << "Output mode: " << output_mode << "\n"
              << "Output prefix: " << (prefix.empty() ? "(none)" : prefix) << "\n"
              << "External engine: " << (engine_cfg.enabled ? "enabled" : "disabled") << "\n"
              << std::endl;

    if (engine_cfg.enabled) {
        std::cout << "  Engine cmd: " << engine_cfg.cmd << "\n"
                  << "  Engine I/O: in='" << engine_cfg.in_file << "', out='" << engine_cfg.out_file << "'\n"
                  << "  Engine units label: " << engine_cfg.units << "\n"
                  << "  Output vector: " << (engine_cfg.output_is_force ? "force" : "gradient") << "\n"
                  << "  Spring constant (k): " << engine_cfg.spring_constant << "\n"
                  << "  Run every: " << engine_cfg.run_every << " cycle(s)\n"
                  << "  Keep per-cycle files: " << (engine_cfg.keep_cycle_files ? "yes" : "no") << "\n"
                  << std::endl;

        if (method != "neb" && method != "neb-iic") {
            std::cout << "  Note: external engine is only used by NEB methods; current method '" << method
                      << "' will ignore it.\n" << std::endl;
        }
    }

    if (method == "iic" || method == "neb-iic") {
        std::cout << "IIC options: bond_factor=" << iic_opt.bond_factor
                  << ", fd_step=" << iic_opt.fd_step
                  << ", ev_thresh=" << iic_opt.ev_thresh
                  << ", iic_maxiter=" << iic_opt.max_iter
                  << ", iic_tol=" << iic_opt.tol
                  << ", iic_damp=" << iic_opt.damp
                  << ", iic_max_step=" << iic_opt.max_cart_step
                  << "\n";
        std::cout << "DM fallback: " << (dm_fallback ? "enabled" : "disabled") << "\n";
    }
    if (method == "dm") {
        std::cout << "DM options: dm_maxiter=" << dm_opt.max_iter
                  << ", dm_step=" << dm_opt.step
                  << ", dm_tol=" << dm_opt.tol
                  << ", dm_max_step=" << dm_opt.max_cart_step
                  << "\n";
    }

    NEBInterpolator interpolator(num_images, step_size, conv_threshold, max_iterations, use_alignment, rmsd_executable);

    interpolator.setIICOptions(iic_opt);
    interpolator.setDMOptions(dm_opt);
    interpolator.setDMFallback(dm_fallback);
    interpolator.setExternalEngine(engine_cfg);
    
    if (!interpolator.setStructures(initial_file, final_file)) {
        return 1;
    }
    
    if (method == "liic") {
        interpolator.performLIIC();
    } else if (method == "iic") {
        interpolator.performIIC();
    } else if (method == "dm") {
        interpolator.performDM();
    } else if (method == "neb") {
        if (!interpolator.performNEB(false)) {
            std::cerr << "Error: NEB optimization failed" << std::endl;
            return 1;
        }
    } else if (method == "neb-iic") {
        if (!interpolator.performNEB(true)) {
            std::cerr << "Error: NEB optimization failed" << std::endl;
            return 1;
        }
    } else {
        std::cerr << "Error: Unknown method: " << method << std::endl;
        return 1;
    }
    
    bool multiframe = (output_mode == "multiframe");
    if (!interpolator.writeResults(prefix, multiframe)) {
        std::cerr << "Error: Failed to write output files" << std::endl;
        return 1;
    }
    
    std::cout << "\nInterpolation completed successfully!" << std::endl;
    return 0;
}