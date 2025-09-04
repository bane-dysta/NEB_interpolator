#ifndef NEB_INTERPOLATOR_H
#define NEB_INTERPOLATOR_H

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <iostream>

struct NEBAtom {
    std::string symbol;
    double x, y, z;
    
    NEBAtom() : symbol(""), x(0.0), y(0.0), z(0.0) {}
    NEBAtom(const std::string& s, double x_, double y_, double z_) 
        : symbol(s), x(x_), y(y_), z(z_) {}
};

class NEBStructure {
public:
    std::vector<NEBAtom> atoms;
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
    
    bool isCompatible(const NEBStructure& other) const {
        if (atoms.size() != other.atoms.size()) return false;
        for (size_t i = 0; i < atoms.size(); ++i) {
            if (atoms[i].symbol != other.atoms[i].symbol) return false;
        }
        return true;
    }
    
    double calculateRMSD(const NEBStructure& other) const {
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

class SimpleNEBInterpolator {
private:
    NEBStructure initial, final;
    std::vector<NEBStructure> images;
    int num_images;
    double step_init;
    double convergence_threshold;
    int max_iterations;
    
    static double distance(const NEBAtom& a1, const NEBAtom& a2) {
        double dx = a1.x - a2.x;
        double dy = a1.y - a2.y;
        double dz = a1.z - a2.z;
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }
    
    void computeDistanceMatrix(const NEBStructure& structure, std::vector<std::vector<double>>& dist_matrix) {
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
    SimpleNEBInterpolator(int num_img = 5, double step = 0.0001, double conv_thresh = 0.01, int max_iter = 10000) 
        : num_images(num_img), step_init(step), convergence_threshold(conv_thresh), max_iterations(max_iter) {}
    
    void setInitialFromMemory(const std::vector<NEBAtom>& atoms_in) {
        initial.atoms = atoms_in;
        initial.comment = "Initial structure from memory";
    }
    
    bool setFinalFromFile(const std::string& final_file) {
        if (!final.readXYZ(final_file)) return false;
        
        if (!initial.isCompatible(final)) {
            std::cerr << "Error: Initial and final structures are not compatible" << std::endl;
            return false;
        }
        
        return true;
    }
    
    void performLIIC() {
        std::cout << "  Initializing LIIC intermediate images..." << std::endl;
        
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
    
    void performNEB() {
        std::cout << "Performing Nudged Elastic Band (NEB) interpolation..." << std::endl;
        
        performLIIC();  // Start from LIIC
        
        double spring_constant = 1.0;
        
        for (int cycle = 0; cycle < max_iterations; ++cycle) {
            double max_force = 0.0;
            std::vector<NEBStructure> forces(num_images);
            
            for (int img = 0; img < num_images; ++img) {
                forces[img].atoms.resize(images[img].size());
                for (size_t i = 0; i < images[img].size(); ++i) {
                    forces[img].atoms[i] = NEBAtom(images[img].atoms[i].symbol, 0.0, 0.0, 0.0);
                }
            }
            
            // NEB force calculation
            for (int img = 0; img < num_images; ++img) {
                NEBStructure& current = images[img];
                NEBStructure* prev = (img == 0) ? &initial : &images[img-1];
                NEBStructure* next = (img == num_images-1) ? &final : &images[img+1];
                
                // Calculate tangent vector
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
                
                // Spring force calculation
                for (size_t i = 0; i < current.size(); ++i) {
                    double spring_force_x = spring_constant * ((next->atoms[i].x - current.atoms[i].x) - (current.atoms[i].x - prev->atoms[i].x));
                    double spring_force_y = spring_constant * ((next->atoms[i].y - current.atoms[i].y) - (current.atoms[i].y - prev->atoms[i].y));
                    double spring_force_z = spring_constant * ((next->atoms[i].z - current.atoms[i].z) - (current.atoms[i].z - prev->atoms[i].z));
                    
                    double proj = spring_force_x * tangent_x[i] + spring_force_y * tangent_y[i] + spring_force_z * tangent_z[i];
                    
                    forces[img].atoms[i].x += proj * tangent_x[i];
                    forces[img].atoms[i].y += proj * tangent_y[i]; 
                    forces[img].atoms[i].z += proj * tangent_z[i];
                }
                
                // Perpendicular force calculation (simplified distance preservation)
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
            
            // Update positions
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
    
    bool writeResults(const std::string& prefix = "") {
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
    
    void setParameters(int num_img, double step, double conv_thresh, int max_iter) {
        num_images = num_img;
        step_init = step;
        convergence_threshold = conv_thresh;
        max_iterations = max_iter;
    }
};

#endif // NEB_INTERPOLATOR_H