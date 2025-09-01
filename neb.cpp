#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <cstring>

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
        
        file.ignore(); // Skip newline after number
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
        
        if (atoms.size() != static_cast<size_t>(natoms)) {
            std::cerr << "Error: Expected " << natoms << " atoms, got " << atoms.size() << std::endl;
            return false;
        }
        
        return true;
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
        if (atoms.size() != other.atoms.size()) {
            return false;
        }
        
        for (size_t i = 0; i < atoms.size(); ++i) {
            if (atoms[i].symbol != other.atoms[i].symbol) {
                return false;
            }
        }
        return true;
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
                    dist_matrix[i][j] = 1e10; // Large number to avoid division by zero
                } else {
                    dist_matrix[i][j] = distance(structure.atoms[i], structure.atoms[j]);
                }
            }
        }
    }
    
    double computeObjectiveFunction(const Structure& current, const std::vector<std::vector<double>>& target_dist) {
        size_t n = current.size();
        std::vector<std::vector<double>> current_dist(n, std::vector<double>(n));
        computeDistanceMatrix(current, current_dist);
        
        double objective = 0.0;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i + 1; j < n; ++j) {
                double diff = target_dist[i][j] - current_dist[i][j];
                objective += diff * diff / (current_dist[i][j] * current_dist[i][j] * current_dist[i][j] * current_dist[i][j]);
            }
        }
        return objective;
    }
    
    void computeGradient(Structure& structure, const std::vector<std::vector<double>>& target_dist, double step) {
        size_t n = structure.size();
        std::vector<std::vector<double>> current_dist(n, std::vector<double>(n));
        computeDistanceMatrix(structure, current_dist);
        
        for (size_t i = 0; i < n; ++i) {
            double grad_x = 0.0, grad_y = 0.0, grad_z = 0.0;
            
            for (size_t j = 0; j < n; ++j) {
                if (i != j) {
                    double dist_diff = target_dist[i][j] - current_dist[i][j];
                    double pos_diff_x = structure.atoms[i].x - structure.atoms[j].x;
                    double pos_diff_y = structure.atoms[i].y - structure.atoms[j].y;
                    double pos_diff_z = structure.atoms[i].z - structure.atoms[j].z;
                    
                    double factor = 2.0 * dist_diff * (2.0 * target_dist[i][j] - current_dist[i][j]) / 
                                   (current_dist[i][j] * current_dist[i][j] * current_dist[i][j] * current_dist[i][j] * current_dist[i][j]);
                    
                    grad_x += factor * pos_diff_x;
                    grad_y += factor * pos_diff_y;
                    grad_z += factor * pos_diff_z;
                }
            }
            
            structure.atoms[i].x += step * grad_x;
            structure.atoms[i].y += step * grad_y;
            structure.atoms[i].z += step * grad_z;
        }
    }
    
public:
    NEBInterpolator(int num_img = 5, double step = 0.0001, double conv_thresh = 0.01, int max_iter = 10000) 
        : num_images(num_img), step_init(step), convergence_threshold(conv_thresh), max_iterations(max_iter) {}
    
    bool setStructures(const std::string& initial_file, const std::string& final_file) {
        if (!initial.readXYZ(initial_file)) {
            return false;
        }
        if (!final.readXYZ(final_file)) {
            return false;
        }
        
        if (!initial.isCompatible(final)) {
            std::cerr << "Error: Initial and final structures are not compatible" << std::endl;
            return false;
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
    
    void performNEB() {
        std::cout << "Performing Nudged Elastic Band (NEB) interpolation..." << std::endl;
        
        // Start with LIIC
        performLIIC();
        
        double spring_constant = 1.0;
        
        // NEB optimization with spring forces and perpendicular forces
        for (int cycle = 0; cycle < max_iterations; ++cycle) {
            double max_force = 0.0;
            std::vector<Structure> forces(num_images);
            
            // Initialize force structures
            for (int img = 0; img < num_images; ++img) {
                forces[img].atoms.resize(images[img].size());
                for (size_t i = 0; i < images[img].size(); ++i) {
                    forces[img].atoms[i] = Atom(images[img].atoms[i].symbol, 0.0, 0.0, 0.0);
                }
            }
            
            // Compute forces for each image
            for (int img = 0; img < num_images; ++img) {
                
                // Spring forces (parallel to the band)
                Structure& current = images[img];
                Structure* prev = (img == 0) ? &initial : &images[img-1];
                Structure* next = (img == num_images-1) ? &final : &images[img+1];
                
                // Compute tangent vector
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
                    
                    // Simplified tangent (could be improved with energy considerations)
                    tangent_x[i] = dx_next - dx_prev;
                    tangent_y[i] = dy_next - dy_prev;
                    tangent_z[i] = dz_next - dz_prev;
                    
                    // Normalize tangent
                    double norm = std::sqrt(tangent_x[i]*tangent_x[i] + tangent_y[i]*tangent_y[i] + tangent_z[i]*tangent_z[i]);
                    if (norm > 1e-12) {
                        tangent_x[i] /= norm;
                        tangent_y[i] /= norm;
                        tangent_z[i] /= norm;
                    }
                }
                
                // Spring forces along the band
                for (size_t i = 0; i < current.size(); ++i) {
                    double spring_force_x = spring_constant * ((next->atoms[i].x - current.atoms[i].x) - (current.atoms[i].x - prev->atoms[i].x));
                    double spring_force_y = spring_constant * ((next->atoms[i].y - current.atoms[i].y) - (current.atoms[i].y - prev->atoms[i].y));
                    double spring_force_z = spring_constant * ((next->atoms[i].z - current.atoms[i].z) - (current.atoms[i].z - prev->atoms[i].z));
                    
                    // Project spring force onto tangent
                    double proj = spring_force_x * tangent_x[i] + spring_force_y * tangent_y[i] + spring_force_z * tangent_z[i];
                    
                    forces[img].atoms[i].x += proj * tangent_x[i];
                    forces[img].atoms[i].y += proj * tangent_y[i]; 
                    forces[img].atoms[i].z += proj * tangent_z[i];
                }
                
                // Perpendicular forces (simplified distance-based potential)
                size_t n = current.size();
                std::vector<std::vector<double>> current_dist(n, std::vector<double>(n));
                std::vector<std::vector<double>> target_dist(n, std::vector<double>(n));
                
                computeDistanceMatrix(current, current_dist);
                
                // Target distances by linear interpolation
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
                
                // Compute perpendicular forces
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
                    
                    // Project out the tangent component (make it perpendicular)
                    double proj = perp_force_x * tangent_x[i] + perp_force_y * tangent_y[i] + perp_force_z * tangent_z[i];
                    
                    forces[img].atoms[i].x += perp_force_x - proj * tangent_x[i];
                    forces[img].atoms[i].y += perp_force_y - proj * tangent_y[i];
                    forces[img].atoms[i].z += perp_force_z - proj * tangent_z[i];
                }
                
                // Track maximum force
                for (size_t i = 0; i < current.size(); ++i) {
                    double force_mag = std::sqrt(forces[img].atoms[i].x * forces[img].atoms[i].x + 
                                               forces[img].atoms[i].y * forces[img].atoms[i].y + 
                                               forces[img].atoms[i].z * forces[img].atoms[i].z);
                    max_force = std::max(max_force, force_mag);
                }
            }
            
            // Update positions based on forces
            for (int img = 0; img < num_images; ++img) {
                for (size_t i = 0; i < images[img].size(); ++i) {
                    images[img].atoms[i].x += step_init * forces[img].atoms[i].x;
                    images[img].atoms[i].y += step_init * forces[img].atoms[i].y;
                    images[img].atoms[i].z += step_init * forces[img].atoms[i].z;
                }
            }
            
            if (cycle % 100 == 0 || cycle < 10) {
                std::cout << "  NEB cycle " << cycle << ", max force = " << max_force << std::endl;
            }
            
            // Check convergence
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
    
private:
    bool writeSeparateXYZ(const std::string& prefix = "") {
        // Write initial structure
        std::string filename = prefix + "00.xyz";
        if (!initial.writeXYZ(filename, "Initial structure")) {
            return false;
        }
        std::cout << "Wrote " << filename << std::endl;
        
        // Write intermediate images
        for (int i = 0; i < num_images; ++i) {
            filename = prefix + (i + 1 < 10 ? "0" : "") + std::to_string(i + 1) + ".xyz";
            if (!images[i].writeXYZ(filename)) {
                return false;
            }
            std::cout << "Wrote " << filename << std::endl;
        }
        
        // Write final structure
        filename = prefix + (num_images + 1 < 10 ? "0" : "") + std::to_string(num_images + 1) + ".xyz";
        if (!final.writeXYZ(filename, "Final structure")) {
            return false;
        }
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
        
        // Write initial structure
        file << initial.atoms.size() << "\n";
        file << "Frame 0: Initial structure\n";
        for (const auto& atom : initial.atoms) {
            file << std::setw(2) << std::left << atom.symbol << " "
                 << std::setw(15) << atom.x << " "
                 << std::setw(15) << atom.y << " "
                 << std::setw(15) << atom.z << "\n";
        }
        
        // Write intermediate images
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
        
        // Write final structure
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
              << "  -n, --nimages NUM      Number of intermediate images (default: 5)\n"
              << "  -m, --method METHOD    Interpolation method: liic or neb (default: neb)\n"
              << "  -p, --prefix PREFIX    Output filename prefix (default: empty)\n"
              << "  -o, --output MODE      Output mode: separate or multiframe (default: separate)\n"
              << "  -s, --step STEP        Initial step size for NEB optimization (default: 0.0001)\n"
              << "  -c, --conv THRESHOLD   Convergence threshold (default: 0.01)\n"
              << "  -i, --maxiter ITER     Maximum iterations for NEB (default: 10000)\n"
              << "  -h, --help             Show this help message\n"
              << "\nOutput modes:\n"
              << "  separate     Generate separate XYZ files (00.xyz, 01.xyz, ...)\n"
              << "  multiframe   Generate single trajectory.xyz with all frames\n"
              << "\nExamples:\n"
              << "  " << program_name << " -n 10 -m neb -p reaction_ initial.xyz final.xyz\n"
              << "  " << program_name << " -o multiframe -n 5 -m liic initial.xyz final.xyz\n";
}

int main(int argc, char* argv[]) {
    std::string initial_file, final_file, prefix = "";
    std::string method = "neb";
    std::string output_mode = "separate";
    int num_images = 5;
    double step_size = 0.0001;
    double conv_threshold = 0.01;
    int max_iterations = 10000;
    
    // Parse command line arguments
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
            if (method != "liic" && method != "neb") {
                std::cerr << "Error: Method must be 'liic' or 'neb'" << std::endl;
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
    
    std::cout << "NEB/LIIC Interpolation Program\n"
              << "==============================\n"
              << "Method: " << (method == "neb" ? "NEB" : "LIIC") << "\n"
              << "Number of images: " << num_images << "\n"
              << "Initial structure: " << initial_file << "\n"
              << "Final structure: " << final_file << "\n"
              << "Output mode: " << output_mode << "\n"
              << "Output prefix: " << (prefix.empty() ? "(none)" : prefix) << "\n"
              << std::endl;
    
    NEBInterpolator interpolator(num_images, step_size, conv_threshold, max_iterations);
    
    if (!interpolator.setStructures(initial_file, final_file)) {
        return 1;
    }
    
    if (method == "liic") {
        interpolator.performLIIC();
    } else {
        interpolator.performNEB();
    }
    
    bool multiframe = (output_mode == "multiframe");
    if (!interpolator.writeResults(prefix, multiframe)) {
        std::cerr << "Error: Failed to write output files" << std::endl;
        return 1;
    }
    
    std::cout << "\nInterpolation completed successfully!" << std::endl;
    return 0;
}