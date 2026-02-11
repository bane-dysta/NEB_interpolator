#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <limits>
#include <sstream>
#include <iomanip>

#include "atom.h"
#include "covalent_radii.h"

class ZMatrixGenerator {
private:
    // 常数定义
    static constexpr double BOHR_TO_ANG = 0.52917721067;
    static constexpr double ANG_TO_BOHR = 1.0 / BOHR_TO_ANG;
    static constexpr double PI = 3.14159265358979;
    static constexpr double BOND_CRITERION = 1.15;
    static constexpr double ANGLE_TOLERANCE = 0.5;
    
    std::vector<Atom> atoms_;
    std::vector<std::vector<int>> zmatrix_;  // [atom_index][0=bond_ref, 1=angle_ref, 2=dihedral_ref]
    std::vector<std::vector<int>> neighbors_;  // 邻居列表
    
    std::string input_file_;
    std::string output_file_;

public:
    ZMatrixGenerator(const std::string& input_file) : input_file_(input_file) {
        // 生成输出文件名
        size_t dot_pos = input_file_.find_last_of('.');
        if (dot_pos != std::string::npos) {
            output_file_ = input_file_.substr(0, dot_pos) + ".zmat";
        } else {
            output_file_ = input_file_ + ".zmat";
        }
    }
    
    void generate() {
        readXYZFile();
        generateNeighborList();
        generateZMatrix();
        outputZMatrix();
    }

private:
    void readXYZFile() {
        std::ifstream file(input_file_);
        if (!file.is_open()) {
            throw std::runtime_error("Error: Cannot open file " + input_file_);
        }
        
        int natom;
        file >> natom;
        
        // 跳过注释行
        std::string line;
        std::getline(file, line);  // 跳过数字后的换行
        std::getline(file, line);  // 跳过注释行
        
        atoms_.reserve(natom);
        
        for (int i = 0; i < natom; ++i) {
            std::string element;
            double x, y, z;
            file >> element >> x >> y >> z;
            
            // 转换为bohr单位
            x *= ANG_TO_BOHR;
            y *= ANG_TO_BOHR;
            z *= ANG_TO_BOHR;
            
            int atomic_num = ChemData::CovalentRadii::getAtomicNumber(element);
            atoms_.emplace_back(element, x, y, z, atomic_num);
        }
        
        file.close();
    }
    
    double atomDistance(int i, int j) const {
        const Atom& atom_i = atoms_[i];
        const Atom& atom_j = atoms_[j];
        
        double dx = atom_i.x - atom_j.x;
        double dy = atom_i.y - atom_j.y;
        double dz = atom_i.z - atom_j.z;
        
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }
    
    double atomAngle(int i, int j, int k) const {
        // 计算角度 i-j-k (j为顶点)
        const Atom& atom_i = atoms_[i];
        const Atom& atom_j = atoms_[j];
        const Atom& atom_k = atoms_[k];
        
        // 向量 j->i 和 j->k
        double v1x = atom_i.x - atom_j.x;
        double v1y = atom_i.y - atom_j.y;
        double v1z = atom_i.z - atom_j.z;
        
        double v2x = atom_k.x - atom_j.x;
        double v2y = atom_k.y - atom_j.y;
        double v2z = atom_k.z - atom_j.z;
        
        double dot_prod = v1x*v2x + v1y*v2y + v1z*v2z;
        double mag1 = std::sqrt(v1x*v1x + v1y*v1y + v1z*v1z);
        double mag2 = std::sqrt(v2x*v2x + v2y*v2y + v2z*v2z);
        
        if (mag1 < 1e-10 || mag2 < 1e-10) {
            return 0.0;
        }
        
        double cos_angle = std::max(-1.0, std::min(1.0, dot_prod / (mag1 * mag2)));
        return std::acos(cos_angle) * 180.0 / PI;
    }
    
    double dihedralAngle(int i, int j, int k, int l) const {
        // 计算二面角 i-j-k-l
        const Atom& atom_i = atoms_[i];
        const Atom& atom_j = atoms_[j];
        const Atom& atom_k = atoms_[k];
        const Atom& atom_l = atoms_[l];
        
        // 向量
        double v1x = atom_i.x - atom_j.x;
        double v1y = atom_i.y - atom_j.y;
        double v1z = atom_i.z - atom_j.z;
        
        double v2x = atom_k.x - atom_j.x;
        double v2y = atom_k.y - atom_j.y;
        double v2z = atom_k.z - atom_j.z;
        
        double v3x = atom_l.x - atom_k.x;
        double v3y = atom_l.y - atom_k.y;
        double v3z = atom_l.z - atom_k.z;
        
        // 叉积得到法向量
        double n1x = v1y*v2z - v1z*v2y;
        double n1y = v1z*v2x - v1x*v2z;
        double n1z = v1x*v2y - v1y*v2x;
        
        double n2x = v2y*v3z - v2z*v3y;
        double n2y = v2z*v3x - v2x*v3z;
        double n2z = v2x*v3y - v2y*v3x;
        
        double dot_prod = n1x*n2x + n1y*n2y + n1z*n2z;
        double mag1 = std::sqrt(n1x*n1x + n1y*n1y + n1z*n1z);
        double mag2 = std::sqrt(n2x*n2x + n2y*n2y + n2z*n2z);
        
        if (mag1 < 1e-10 || mag2 < 1e-10) {
            return 0.0;
        }
        
        // 叉积用于确定符号
        double cross_prod = (n1y*n2z - n1z*n2y)*v2x + 
                           (n1z*n2x - n1x*n2z)*v2y + 
                           (n1x*n2y - n1y*n2x)*v2z;
        
        return std::atan2(cross_prod, dot_prod) * 180.0 / PI;
    }
    
    void generateNeighborList() {
        size_t natom = atoms_.size();
        neighbors_.resize(natom);
        
        for (size_t i = 0; i < natom; ++i) {
            for (size_t j = 0; j < natom; ++j) {
                if (i == j) continue;
                
                double dist = atomDistance(i, j);
                double sum_cov = ChemData::CovalentRadii::getRadius(atoms_[i].atomic_number) +
                                ChemData::CovalentRadii::getRadius(atoms_[j].atomic_number);
                
                if (dist < BOND_CRITERION * sum_cov) {
                    neighbors_[i].push_back(j);
                }
            }
        }
    }
    
    void generateZMatrix() {
        size_t natom = atoms_.size();
        zmatrix_.resize(natom, std::vector<int>(3, -1));  // -1表示未设置
        
        for (size_t iatm = 1; iatm < natom; ++iatm) {  // 从第2个原子开始
            // 寻找键连接参考
            int i1 = findBondReference(iatm);
            zmatrix_[iatm][0] = i1;
            
            if (iatm >= 2) {  // 从第3个原子开始需要角度参考
                int i2 = findAngleReference(iatm, i1);
                zmatrix_[iatm][1] = i2;
                
                if (iatm >= 3) {  // 从第4个原子开始需要二面角参考
                    int i3 = findDihedralReference(iatm, i1, i2);
                    zmatrix_[iatm][2] = i3;
                }
            }
        }
    }
    
    int findBondReference(int iatm) {
        // 优先寻找键连接的原子
        for (int neighbor : neighbors_[iatm]) {
            if (neighbor < iatm) {  // 只考虑之前的原子
                return neighbor;
            }
        }
        
        // 如果没有键连接，寻找最近的原子
        double min_dist = std::numeric_limits<double>::max();
        int closest = -1;
        
        for (int j = 0; j < iatm; ++j) {
            double dist = atomDistance(iatm, j);
            if (dist < min_dist) {
                min_dist = dist;
                closest = j;
            }
        }
        
        return closest;
    }
    
    int findAngleReference(int iatm, int i1) {
        // 优先寻找与i1键连接的原子
        for (int neighbor : neighbors_[i1]) {
            if (neighbor < iatm && neighbor != i1) {
                if (atoms_.size() > 3) {
                    double angle = atomAngle(iatm, i1, neighbor);
                    if (angle < ANGLE_TOLERANCE || angle > (180.0 - ANGLE_TOLERANCE)) {
                        continue;  // 避免接近线性的角度
                    }
                }
                return neighbor;
            }
        }
        
        // 如果没有找到，寻找距离i1最近的原子
        double min_dist = std::numeric_limits<double>::max();
        int closest = -1;
        
        for (int j = 0; j < iatm; ++j) {
            if (j == i1) continue;
            
            double dist = atomDistance(i1, j);
            if (dist < min_dist) {
                if (atoms_.size() > 3) {
                    double angle = atomAngle(iatm, i1, j);
                    if (angle < ANGLE_TOLERANCE || angle > (180.0 - ANGLE_TOLERANCE)) {
                        continue;
                    }
                }
                min_dist = dist;
                closest = j;
            }
        }
        
        if (closest == -1) {
            throw std::runtime_error("Error: Cannot find suitable angle reference for atom " + 
                                   std::to_string(iatm + 1));
        }
        
        return closest;
    }
    
    int findDihedralReference(int iatm, int i1, int i2) {
        // 优先寻找与i2键连接的原子
        for (int neighbor : neighbors_[i2]) {
            if (neighbor < iatm && neighbor != i1 && neighbor != i2) {
                double angle = atomAngle(i1, i2, neighbor);
                if (angle < ANGLE_TOLERANCE || angle > (180.0 - ANGLE_TOLERANCE)) {
                    continue;
                }
                return neighbor;
            }
        }
        
        // 如果没有找到，寻找距离i2最近的原子
        double min_dist = std::numeric_limits<double>::max();
        int closest = -1;
        
        for (int j = 0; j < iatm; ++j) {
            if (j == i1 || j == i2) continue;
            
            double dist = atomDistance(i2, j);
            if (dist < min_dist) {
                double angle = atomAngle(i1, i2, j);
                if (angle < ANGLE_TOLERANCE || angle > (180.0 - ANGLE_TOLERANCE)) {
                    continue;
                }
                min_dist = dist;
                closest = j;
            }
        }
        
        if (closest == -1) {
            throw std::runtime_error("Error: Cannot find suitable dihedral reference for atom " + 
                                   std::to_string(iatm + 1));
        }
        
        return closest;
    }
    
    void outputZMatrix() {
        std::ofstream file(output_file_);
        if (!file.is_open()) {
            throw std::runtime_error("Error: Cannot create output file " + output_file_);
        }
        
        file << "# Z-Matrix format\n";
        file << "# Atom  Bond_ref  Bond(A)    Angle_ref  Angle(deg)  Dihedral_ref  Dihedral(deg)\n";
        
        // 第一个原子
        file << std::left << std::setw(4) << atoms_[0].name << "\n";
        
        // 第二个原子
        if (atoms_.size() >= 2) {
            int i1 = zmatrix_[1][0];
            double bond_dist = atomDistance(1, i1) * BOHR_TO_ANG;
            
            file << std::left << std::setw(4) << atoms_[1].name
                 << std::right << std::setw(6) << (i1 + 1)
                 << std::setw(12) << std::fixed << std::setprecision(6) << bond_dist
                 << "\n";
        }
        
        // 第三个及以后的原子
        for (size_t i = 2; i < atoms_.size(); ++i) {
            int i1 = zmatrix_[i][0];
            int i2 = zmatrix_[i][1];
            
            double bond_dist = atomDistance(i, i1) * BOHR_TO_ANG;
            double bond_angle = atomAngle(i, i1, i2);
            
            file << std::left << std::setw(4) << atoms_[i].name
                 << std::right << std::setw(6) << (i1 + 1)
                 << std::setw(12) << std::fixed << std::setprecision(6) << bond_dist
                 << std::setw(8) << (i2 + 1)
                 << std::setw(10) << std::fixed << std::setprecision(4) << bond_angle;
            
            if (i >= 3) {
                int i3 = zmatrix_[i][2];
                double dihedral = dihedralAngle(i, i1, i2, i3);
                
                file << std::setw(8) << (i3 + 1)
                     << std::setw(12) << std::fixed << std::setprecision(4) << dihedral;
            }
            
            file << "\n";
        }
        
        file.close();
        std::cout << "Z-matrix generated and saved to: " << output_file_ << std::endl;
    }
};

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Error: No input file specified\n";
        std::cerr << "Usage: " << argv[0] << " input.xyz\n";
        return 1;
    }
    
    try {
        ZMatrixGenerator generator(argv[1]);
        generator.generate();
        std::cout << "Z-matrix generation completed successfully!" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}