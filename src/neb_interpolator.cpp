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
        
        // 将对齐后的文件复制回原文件
        std::string copy_command = "cp " + aligned_file + " " + mobile_file;
        int result = system(copy_command.c_str());
        
        if (result != 0) {
            std::cerr << "  Error: Failed to replace original file" << std::endl;
            return false;
        }
        
        // 删除临时文件
        std::string remove_command = "rm " + aligned_file;
        int ret = system(remove_command.c_str());
        (void)ret;
        
        return true;
    }
    
    void setExecutablePath(const std::string& path) {
        rmsd_executable = path;
    }
};

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
            int ret = system(("rm -f " + temp_final).c_str());
            (void)ret;
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

    void performNEB(bool init_with_iic = false) {
        std::cout << "Performing Nudged Elastic Band (NEB) interpolation..." << std::endl;
        
        if (init_with_iic) {
            performIIC();
        } else {
            performLIIC();
        }
        
        double spring_constant = 1.0;
        
        for (int cycle = 0; cycle < max_iterations; ++cycle) {
            double max_force = 0.0;
            std::vector<Structure> forces(num_images);
            
            for (int img = 0; img < num_images; ++img) {
                forces[img].atoms.resize(images[img].size());
                for (size_t i = 0; i < images[img].size(); ++i) {
                    forces[img].atoms[i] = Atom(images[img].atoms[i].symbol, 0.0, 0.0, 0.0);
                }
            }
            
            // NEB力计算
            for (int img = 0; img < num_images; ++img) {
                Structure& current = images[img];
                Structure* prev = (img == 0) ? &initial : &images[img-1];
                Structure* next = (img == num_images-1) ? &final : &images[img+1];
                
                // 计算切线向量
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
                    
                    tangent_x[i] = dx_next - dx_prev;
                    tangent_y[i] = dy_next - dy_prev;
                    tangent_z[i] = dz_next - dz_prev;
                    
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
              << "  " << program_name << " -m dm -n 8 initial.xyz final.xyz\n";
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
              << std::endl;

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
        interpolator.performNEB(false);
    } else if (method == "neb-iic") {
        interpolator.performNEB(true);
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