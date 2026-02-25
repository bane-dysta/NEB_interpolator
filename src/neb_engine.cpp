#include "neb_engine.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>

#include "util.h"

namespace neb {

namespace {

// ============================================================
// External engine I/O contract
// ============================================================
//
// The external engine interface is fixed-format text files.
//
// Input file: blocks of XYZ geometries (type=xyz)
//   <natoms>
//   image=<i> units=<U> type=xyz cycle=<cycle>
//   <sym> <x> <y> <z>
//   ...
//
// Output file: blocks of gradients or forces (type=gradient or type=force)
//   <natoms>
//   image=<i> units=<U> type=gradient
//   <sym> <vx> <vy> <vz>
//   ...
//
// The output block must be provided for each image (1..nimages). If the engine
// omits image=, blocks are assigned sequentially.
//
// Notes:
// - Units are treated as labels only (no conversion).
// - If config.output_is_force == false, vectors are interpreted as gradients
//   and converted to forces by force = -grad.

static inline std::string trim(const std::string& s) {
    const char* ws = " \t\r\n";
    const size_t b = s.find_first_not_of(ws);
    if (b == std::string::npos) return "";
    const size_t e = s.find_last_not_of(ws);
    return s.substr(b, e - b + 1);
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
                                const std::vector<geom::Structure>& images,
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
        const geom::Structure& s = images[img];
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

    if (cycle == 0) {
        if (has_in && !has_out) {
            std::cerr << "Warning: --engine-cmd contains {in} but not {out}; appending output path at the end." << std::endl;
        } else if (!has_in && has_out) {
            std::cerr << "Warning: --engine-cmd contains {out} but not {in}; appending input path at the end." << std::endl;
        }
    }

    cmdline = replaceAll(cmdline, "{in}", util::shellQuote(in_path));
    cmdline = replaceAll(cmdline, "{out}", util::shellQuote(out_path));
    if (has_cycle) {
        cmdline = replaceAll(cmdline, "{cycle}", std::to_string(cycle));
    }

    // If user didn't provide placeholders, append missing args.
    // This is intentionally forgiving: if the user includes only {in} or only
    // {out}, we still add the missing counterpart so engines that expect both
    // paths keep working.
    if (!has_in) {
        cmdline += " " + util::shellQuote(in_path);
    }
    if (!has_out) {
        cmdline += " " + util::shellQuote(out_path);
    }

    std::cout << "  [engine] cycle " << cycle << ": " << cmdline << std::endl;

    util::CommandResult res = util::runCommand(cmdline);
    if (!(res.exited && res.exit_code == 0)) {
        if (err) {
            std::ostringstream oss;
            oss << "External engine command failed (" << util::formatCommandFailure(res) << ")";
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
                                  const std::vector<geom::Structure>& images,
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

} // namespace

// ============================================================
// ExternalEngine
// ============================================================

ExternalEngine::ExternalEngine(const ExternalEngineConfig& cfg) : cfg_(cfg) {}

bool ExternalEngine::compute_forces(const std::vector<geom::Structure>& images,
                                    int cycle,
                                    std::vector<std::vector<double>>& forces,
                                    std::string* err) {
    if (images.empty()) {
        if (err) *err = "ExternalEngine: no images provided";
        return false;
    }

    const int nimages = static_cast<int>(images.size());
    const int natoms = static_cast<int>(images[0].size());
    if (natoms <= 0) {
        if (err) *err = "ExternalEngine: natoms is zero";
        return false;
    }
    for (int img = 1; img < nimages; ++img) {
        if (static_cast<int>(images[img].size()) != natoms) {
            if (err) *err = "ExternalEngine: inconsistent natoms across images";
            return false;
        }
    }

    if (!have_cache_ && cfg_.run_every > 1) {
        std::cerr << "Warning: --engine-run-every " << cfg_.run_every
                  << " will reuse the previous cycle's forces for intermediate cycles. "
                  << "This is an approximation and may harm convergence." << std::endl;
    }

    const bool need_run = (!have_cache_) || (cfg_.run_every <= 1) || (cycle % cfg_.run_every == 0);
    if (need_run) {
        std::string in_path = cfg_.in_file;
        std::string out_path = cfg_.out_file;
        if (cfg_.keep_cycle_files) {
            in_path = addCycleSuffix(in_path, cycle);
            out_path = addCycleSuffix(out_path, cycle);
        }

        std::string local_err;
        if (!writeEngineInputXYZ(in_path, images, cfg_.units, cycle, &local_err)) {
            if (err) *err = local_err;
            return false;
        }
        if (!runExternalEngineOnce(cfg_, cycle, in_path, out_path, natoms, nimages, &local_err)) {
            if (err) *err = local_err;
            return false;
        }

        std::vector<EngineVectorBlock> blocks;
        if (!readEngineVectorBlocks(out_path, blocks, &local_err)) {
            if (err) *err = local_err;
            return false;
        }

        std::vector<std::vector<double>> vec;
        if (!assembleEngineVectors(blocks, images, nimages, vec, &local_err)) {
            if (err) *err = local_err;
            return false;
        }

        // Convert to real forces (gradient -> force = -grad; force -> as-is)
        forces.assign(nimages, std::vector<double>(3 * static_cast<size_t>(natoms), 0.0));
        for (int img = 0; img < nimages; ++img) {
            if (vec[img].size() != forces[img].size()) {
                if (err) *err = "ExternalEngine: vector length mismatch";
                return false;
            }
            if (cfg_.output_is_force) {
                forces[img] = vec[img];
            } else {
                for (size_t k = 0; k < vec[img].size(); ++k) {
                    forces[img][k] = -vec[img][k];
                }
            }
        }

        last_forces_ = forces;
        have_cache_ = true;
    } else {
        forces = last_forces_;
    }

    return true;
}

// ============================================================
// DistancePenaltyEngine
// ============================================================

static inline double dist_atom(const geom::Atom& a, const geom::Atom& b) {
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

DistancePenaltyEngine::DistancePenaltyEngine(const geom::Structure& initial, const geom::Structure& final)
    : natoms_(static_cast<int>(initial.size())) {
    const size_t n = static_cast<size_t>(natoms_);
    d0_.assign(n * n, 0.0);
    d1_.assign(n * n, 0.0);

    for (int i = 0; i < natoms_; ++i) {
        for (int j = i + 1; j < natoms_; ++j) {
            const double di = dist_atom(initial.atoms[i], initial.atoms[j]);
            const double df = dist_atom(final.atoms[i], final.atoms[j]);
            d0_[n * i + j] = di;
            d0_[n * j + i] = di;
            d1_[n * i + j] = df;
            d1_[n * j + i] = df;
        }
    }
}

bool DistancePenaltyEngine::compute_forces(const std::vector<geom::Structure>& images,
                                          int /*cycle*/,
                                          std::vector<std::vector<double>>& forces,
                                          std::string* err) {
    const int num_images = static_cast<int>(images.size());
    if (natoms_ <= 0) {
        if (err) *err = "DistancePenaltyEngine: natoms is zero";
        return false;
    }

    if (num_images == 0) {
        forces.clear();
        return true;
    }

    const size_t n = static_cast<size_t>(natoms_);
    forces.assign(num_images, std::vector<double>(3 * n, 0.0));

    for (int img = 0; img < num_images; ++img) {
        const geom::Structure& s = images[img];
        if (static_cast<int>(s.size()) != natoms_) {
            if (err) *err = "DistancePenaltyEngine: natoms mismatch in image";
            return false;
        }

        const double factor = static_cast<double>(img + 1) / static_cast<double>(num_images + 1);

        // Pairwise penalty forces (NO spring, NO projection here)
        for (int i = 0; i < natoms_; ++i) {
            for (int j = i + 1; j < natoms_; ++j) {
                const double dx = s.atoms[i].x - s.atoms[j].x;
                const double dy = s.atoms[i].y - s.atoms[j].y;
                const double dz = s.atoms[i].z - s.atoms[j].z;
                const double dij = std::sqrt(dx * dx + dy * dy + dz * dz);
                if (dij < 1e-12) continue;

                const double dtarget = d0_[n * i + j] + factor * (d1_[n * i + j] - d0_[n * i + j]);
                const double diff = dtarget - dij; // (target - current)
                const double scale = 2.0 * diff / dij;

                // Force on i: +scale * (ri - rj); on j: opposite
                forces[img][3 * i + 0] += scale * dx;
                forces[img][3 * i + 1] += scale * dy;
                forces[img][3 * i + 2] += scale * dz;

                forces[img][3 * j + 0] -= scale * dx;
                forces[img][3 * j + 1] -= scale * dy;
                forces[img][3 * j + 2] -= scale * dz;
            }
        }
    }

    return true;
}

} // namespace neb
