#ifndef ATOM_H
#define ATOM_H

#include <string>

struct Atom {
    std::string name;
    double x, y, z;         // 坐标 (bohr)
    int atomic_number;
    
    Atom() : x(0.0), y(0.0), z(0.0), atomic_number(1) {}
    
    Atom(const std::string& element, double x_coord, double y_coord, double z_coord, int atomic_num)
        : name(element), x(x_coord), y(y_coord), z(z_coord), atomic_number(atomic_num) {}
};

#endif // ATOM_H