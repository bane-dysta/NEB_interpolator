#ifndef NEB_CORE_H
#define NEB_CORE_H

#include <string>
#include <vector>

#include "geometry.h"

namespace neb {

// Assemble NEB forces for each intermediate image.
//
// This is a stateless (pure) function: it does NOT call any external engine,
// does NOT do file I/O, and does NOT update coordinates. It only assembles
// the NEB force from the provided real forces.
//
// Convention:
//   forces[img] is a flattened 3N vector: [x1,y1,z1,x2,y2,z2,...]
//
// NEB formula (standard, simplified):
//   F_i = F_real_i^\perp + F_spring_i^\parallel
//       = (F_real_i - (F_real_i·tau_i) tau_i) + k (|R_{i+1}-R_i| - |R_i-R_{i-1}|) tau_i
//
// Tangent tau_i is computed as a single global 3N vector and normalized:
//   tau_i = normalize(R_{i+1} - R_{i-1})
//
// Returns true on success; on failure returns false and sets *err (if provided).
bool assembleNEBForces(const geom::Structure& initial,
                       const geom::Structure& final,
                       const std::vector<geom::Structure>& images,
                       const std::vector<std::vector<double>>& real_forces,
                       double spring_constant,
                       std::vector<std::vector<double>>& neb_forces,
                       double* max_force,
                       std::string* err);

} // namespace neb

#endif // NEB_CORE_H
