#ifndef COVALENT_RADII_H
#define COVALENT_RADII_H

#include <unordered_map>
#include <string>

namespace ChemData {

// 共价半径数据集类型
enum class RadiiType {
    DALTON,    // Dalton Trans., 2008, 2832-2838 (默认)
    SURESH     // J. Phys. Chem. A 2001, 105, 5940-5944
};

class CovalentRadii {
public:
    // 获取共价半径 (单位: Ångström)
    static double getRadius(int atomic_number, RadiiType type = RadiiType::DALTON);
    static double getRadius(const std::string& element, RadiiType type = RadiiType::DALTON);
    
    // 获取原子序数
    static int getAtomicNumber(const std::string& element);
    
    // 获取元素符号
    static std::string getElementSymbol(int atomic_number);
    
    // 检查元素是否支持
    static bool isSupported(int atomic_number, RadiiType type = RadiiType::DALTON);
    static bool isSupported(const std::string& element, RadiiType type = RadiiType::DALTON);
    
    // 获取数据集信息
    static std::string getDatasetInfo(RadiiType type);
    static size_t getSupportedElementCount(RadiiType type);

private:
    // 两套共价半径数据 (Ångström)
    static const std::unordered_map<int, double> radii_data_;        // Dalton数据集
    static const std::unordered_map<int, double> suresh_radii_data_; // Suresh数据集
    static const std::unordered_map<std::string, int> element_to_atomic_number_;
    static const std::unordered_map<int, std::string> atomic_number_to_element_;
    
    // 根据类型获取相应的数据集
    static const std::unordered_map<int, double>& getDataset(RadiiType type);
};

// 内联实现
inline double CovalentRadii::getRadius(int atomic_number, RadiiType type) {
    const auto& dataset = getDataset(type);
    auto it = dataset.find(atomic_number);
    return (it != dataset.end()) ? it->second : 0.0;
}

inline double CovalentRadii::getRadius(const std::string& element, RadiiType type) {
    return getRadius(getAtomicNumber(element), type);
}

inline int CovalentRadii::getAtomicNumber(const std::string& element) {
    auto it = element_to_atomic_number_.find(element);
    return (it != element_to_atomic_number_.end()) ? it->second : 6; // 默认返回C
}

inline std::string CovalentRadii::getElementSymbol(int atomic_number) {
    auto it = atomic_number_to_element_.find(atomic_number);
    return (it != atomic_number_to_element_.end()) ? it->second : "C";
}

inline bool CovalentRadii::isSupported(int atomic_number, RadiiType type) {
    const auto& dataset = getDataset(type);
    return dataset.find(atomic_number) != dataset.end();
}

inline bool CovalentRadii::isSupported(const std::string& element, RadiiType type) {
    return isSupported(getAtomicNumber(element), type);
}

inline const std::unordered_map<int, double>& CovalentRadii::getDataset(RadiiType type) {
    return (type == RadiiType::SURESH) ? suresh_radii_data_ : radii_data_;
}

inline std::string CovalentRadii::getDatasetInfo(RadiiType type) {
    switch (type) {
        case RadiiType::DALTON:
            return "Dalton Trans., 2008, 2832-2838";
        case RadiiType::SURESH:
            return "J. Phys. Chem. A 2001, 105, 5940-5944";
        default:
            return "Unknown dataset";
    }
}

inline size_t CovalentRadii::getSupportedElementCount(RadiiType type) {
    const auto& dataset = getDataset(type);
    return dataset.size();
}

} // namespace ChemData

#endif // COVALENT_RADII_H