#ifndef COVALENT_RADII_H
#define COVALENT_RADII_H

#include <unordered_map>
#include <string>
#include <algorithm>
#include <cctype>
#include <stdexcept>

namespace ChemData {

// Raised when an element symbol (e.g. from an XYZ file) cannot be recognized.
class UnknownElementError : public std::runtime_error {
public:
    explicit UnknownElementError(const std::string& msg)
        : std::runtime_error(msg) {}
};

// Raised when an atomic number cannot be mapped back to an element symbol.
class UnknownAtomicNumberError : public std::runtime_error {
public:
    explicit UnknownAtomicNumberError(const std::string& msg)
        : std::runtime_error(msg) {}
};

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
    // Trim whitespace (defensive)
    std::string sym = element;
    auto not_space = [](unsigned char c) { return std::isspace(c) == 0; };
    sym.erase(sym.begin(), std::find_if(sym.begin(), sym.end(), not_space));
    sym.erase(std::find_if(sym.rbegin(), sym.rend(), not_space).base(), sym.end());

    if (sym.empty()) {
        throw UnknownElementError("Error: empty element symbol");
    }

    // Normalize capitalization: first letter uppercase, remaining letters lowercase.
    sym[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(sym[0])));
    for (size_t i = 1; i < sym.size(); ++i) {
        sym[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(sym[i])));
    }

    auto it = element_to_atomic_number_.find(sym);
    if (it == element_to_atomic_number_.end()) {
        throw UnknownElementError("Error: unknown element symbol '" + element + "' (normalized to '" + sym + "')");
    }
    return it->second;
}

inline std::string CovalentRadii::getElementSymbol(int atomic_number) {
    auto it = atomic_number_to_element_.find(atomic_number);
    if (it == atomic_number_to_element_.end()) {
        throw UnknownAtomicNumberError("Error: unknown atomic number " + std::to_string(atomic_number));
    }
    return it->second;
}

inline bool CovalentRadii::isSupported(int atomic_number, RadiiType type) {
    const auto& dataset = getDataset(type);
    return dataset.find(atomic_number) != dataset.end();
}

inline bool CovalentRadii::isSupported(const std::string& element, RadiiType type) {
    try {
        return isSupported(getAtomicNumber(element), type);
    } catch (...) {
        return false;
    }
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