#include <cstring>
#include <iostream>
#include <string>

#include "neb_driver.h"
#include "rmsd_align.h"

static void printUsage(const char* program_name) {
    std::cout
        << "Usage: " << program_name << " [options] <initial.xyz> <final.xyz>\n"
        << "\nOptions:\n"
        << "  -n, --nimages NUM        Number of intermediate images (default: 5)\n"
        << "  -m, --method METHOD      Method: lic | liic | dm | neb | neb-liic (default: neb)\n"
        << "  -p, --prefix PREFIX      Output filename prefix (default: empty)\n"
        << "  -o, --output MODE        Output mode: separate or multiframe (default: separate)\n"
        << "  -s, --step STEP          Step size for NEB optimization update (default: 0.0001)\n"
        << "  -c, --conv THRESHOLD     Convergence threshold (default: 0.01)\n"
        << "  -i, --maxiter ITER       Maximum iterations for NEB (default: 10000)\n"
        << "  -a, --align              Enable structure alignment using calc_rmsd_xyz (default: enabled)\n"
        << "  --no-align               Disable structure alignment\n"
        << "  -r, --rmsd-exec PATH     Path to calc_rmsd_xyz (default: auto-detect)\n"
        << "  -h, --help               Show this help message\n"
        << "\nLIIC options (used by -m liic or -m neb-liic):\n"
        << "  --bond-factor F          Bond cutoff factor (default: 1.25)\n"
        << "  --fd-step H              Finite-difference step for B matrix (default: 1e-4)\n"
        << "  --ev-thresh T            Eigenvalue threshold for DLC selection (default: 1e-3)\n"
        << "  --liic-maxiter N         Max back-transform iterations (default: 50)\n"
        << "  --liic-tol T             RMS primitive residual tolerance (default: 1e-4)\n"
        << "  --liic-damp D            Damping added to eigenvalues (default: 1e-8)\n"
        << "  --liic-max-step S        Max cartesian step per iteration (default: 0.20)\n"
        << "  --liic-verbose V         LIIC verbosity: 0=silent, 1=per-image, 2=per-iter (default: 0)\n"
        << "\nDM options (used as fallback when LIIC fails, or -m dm):\n"
        << "  --dm-maxiter N           Max DM iterations (default: 800)\n"
        << "  --dm-step S              DM gradient step size (default: 5e-3)\n"
        << "  --dm-tol T               RMS distance-error tolerance (default: 1e-3)\n"
        << "  --dm-max-step S          Max cartesian step per iteration (default: 0.20)\n"
        << "  --no-dm-fallback         Disable DM fallback (fallback directly to LIC)\n"
        << "\nExternal engine options (for NEB / NEB-LIIC only):\n"
        << "  (Note) If external engine is NOT enabled, NEB/NEB-LIIC uses DistancePenaltyEngine (virtual distance-penalty forces).\n"
        << "  --engine-cmd CMD         Enable external-engine mode; run CMD each NEB cycle\n"
        << "                           If CMD contains {in}/{out}/{cycle}, they will be replaced.\n"
        << "                           Otherwise, the program appends: <infile> <outfile>\n"
        << "  --engine-in FILE         Engine input filename (default: neb_engine_in.dat)\n"
        << "  --engine-out FILE        Engine output filename (default: neb_engine_out.dat)\n"
        << "  --engine-units U         Units label written in engine I/O headers (default: Angstrom)\n"
        << "  --engine-output-is-force Interpret engine output as force (default: gradient)\n"
        << "  --engine-spring K        Spring constant used in external-engine NEB (default: 1.0)\n"
        << "  --engine-run-every N     Run engine every N cycles (default: 1). N>1 reuses the previous cycle's forces/gradients on skipped cycles.\n"
        << "  --engine-keep-cycle-files Keep per-cycle engine I/O files (adds _cycle#### suffix)\n"
        << "\nOutput modes:\n"
        << "  separate     Generate separate XYZ files (00.xyz, 01.xyz, ...)\n"
        << "  multiframe   Generate a single trajectory.xyz with all frames\n"
        << "\nExamples:\n"
        << "  " << program_name << " -n 10 -m neb -p reaction_ initial.xyz final.xyz\n"
        << "  " << program_name << " -n 10 -m neb-liic --bond-factor 1.2 --fd-step 1e-4 --ev-thresh 1e-3 initial.xyz final.xyz\n"
        << "  " << program_name << " --no-align -o multiframe -n 5 -m liic initial.xyz final.xyz\n"
        << "  " << program_name << " -m dm -n 8 initial.xyz final.xyz\n"
        << "  " << program_name << " -m neb -n 8 --engine-cmd \"python engine.py {in} {out}\" initial.xyz final.xyz\n";
}

static bool isValidMethod(const std::string& m) {
    // Canonical
    if (m == "lic" || m == "liic" || m == "dm" || m == "neb" || m == "neb-liic") return true;
    // Hidden legacy aliases (kept for personal backward-compatibility)
    if (m == "iic" || m == "neb-iic") return true;
    return false;
}

int main(int argc, char* argv[]) {
    std::string initial_file;
    std::string final_file;
    std::string prefix;
    std::string method = "neb";
    std::string output_mode = "separate";

    std::string rmsd_executable = rmsd::findCalcRMSDExecutable();

    int num_images = 5;
    double step_size = 0.0001;
    double conv_threshold = 0.01;
    int max_iterations = 10000;
    bool use_alignment = true;

    // LIIC/DM options (dependency-free)
    ICInterp::LIICOptions liic_opt;
    ICInterp::DMOptions dm_opt;
    bool dm_fallback = true;

    // External engine options (optional)
    neb::ExternalEngineConfig engine_cfg;

    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "-h") == 0 || std::strcmp(argv[i], "--help") == 0) {
            printUsage(argv[0]);
            return 0;
        }

        if ((std::strcmp(argv[i], "-n") == 0 || std::strcmp(argv[i], "--nimages") == 0) && i + 1 < argc) {
            num_images = std::atoi(argv[++i]);
            if (num_images <= 0) {
                std::cerr << "Error: number of images must be positive" << std::endl;
                return 1;
            }
            continue;
        }

        if ((std::strcmp(argv[i], "-m") == 0 || std::strcmp(argv[i], "--method") == 0) && i + 1 < argc) {
            method = argv[++i];
            if (!isValidMethod(method)) {
                std::cerr << "Error: method must be one of: lic, liic, dm, neb, neb-liic" << std::endl;
                return 1;
            }
            continue;
        }

        if ((std::strcmp(argv[i], "-p") == 0 || std::strcmp(argv[i], "--prefix") == 0) && i + 1 < argc) {
            prefix = argv[++i];
            continue;
        }

        if ((std::strcmp(argv[i], "-o") == 0 || std::strcmp(argv[i], "--output") == 0) && i + 1 < argc) {
            output_mode = argv[++i];
            if (output_mode != "separate" && output_mode != "multiframe") {
                std::cerr << "Error: output mode must be 'separate' or 'multiframe'" << std::endl;
                return 1;
            }
            continue;
        }

        if ((std::strcmp(argv[i], "-s") == 0 || std::strcmp(argv[i], "--step") == 0) && i + 1 < argc) {
            step_size = std::atof(argv[++i]);
            if (step_size <= 0.0) {
                std::cerr << "Error: step size must be positive" << std::endl;
                return 1;
            }
            continue;
        }

        if ((std::strcmp(argv[i], "-c") == 0 || std::strcmp(argv[i], "--conv") == 0) && i + 1 < argc) {
            conv_threshold = std::atof(argv[++i]);
            if (conv_threshold <= 0.0) {
                std::cerr << "Error: convergence threshold must be positive" << std::endl;
                return 1;
            }
            continue;
        }

        if ((std::strcmp(argv[i], "-i") == 0 || std::strcmp(argv[i], "--maxiter") == 0) && i + 1 < argc) {
            max_iterations = std::atoi(argv[++i]);
            if (max_iterations <= 0) {
                std::cerr << "Error: maximum iterations must be positive" << std::endl;
                return 1;
            }
            continue;
        }

        if (std::strcmp(argv[i], "-a") == 0 || std::strcmp(argv[i], "--align") == 0) {
            use_alignment = true;
            continue;
        }

        if (std::strcmp(argv[i], "--no-align") == 0) {
            use_alignment = false;
            continue;
        }

        if ((std::strcmp(argv[i], "-r") == 0 || std::strcmp(argv[i], "--rmsd-exec") == 0) && i + 1 < argc) {
            rmsd_executable = argv[++i];
            continue;
        }

        // --- LIIC options (new canonical flags) ---
        if (std::strcmp(argv[i], "--bond-factor") == 0 && i + 1 < argc) {
            liic_opt.bond_factor = std::atof(argv[++i]);
            if (liic_opt.bond_factor <= 0.0) {
                std::cerr << "Error: --bond-factor must be positive" << std::endl;
                return 1;
            }
            continue;
        }

        if (std::strcmp(argv[i], "--fd-step") == 0 && i + 1 < argc) {
            liic_opt.fd_step = std::atof(argv[++i]);
            if (liic_opt.fd_step <= 0.0) {
                std::cerr << "Error: --fd-step must be positive" << std::endl;
                return 1;
            }
            continue;
        }

        if (std::strcmp(argv[i], "--ev-thresh") == 0 && i + 1 < argc) {
            liic_opt.ev_thresh = std::atof(argv[++i]);
            if (liic_opt.ev_thresh < 0.0) {
                std::cerr << "Error: --ev-thresh must be non-negative" << std::endl;
                return 1;
            }
            continue;
        }

        if (std::strcmp(argv[i], "--liic-maxiter") == 0 && i + 1 < argc) {
            liic_opt.max_iter = std::atoi(argv[++i]);
            if (liic_opt.max_iter <= 0) {
                std::cerr << "Error: --liic-maxiter must be positive" << std::endl;
                return 1;
            }
            continue;
        }

        if (std::strcmp(argv[i], "--liic-tol") == 0 && i + 1 < argc) {
            liic_opt.tol = std::atof(argv[++i]);
            if (liic_opt.tol <= 0.0) {
                std::cerr << "Error: --liic-tol must be positive" << std::endl;
                return 1;
            }
            continue;
        }

        if (std::strcmp(argv[i], "--liic-damp") == 0 && i + 1 < argc) {
            liic_opt.damp = std::atof(argv[++i]);
            if (liic_opt.damp < 0.0) {
                std::cerr << "Error: --liic-damp must be non-negative" << std::endl;
                return 1;
            }
            continue;
        }

        if (std::strcmp(argv[i], "--liic-max-step") == 0 && i + 1 < argc) {
            liic_opt.max_cart_step = std::atof(argv[++i]);
            if (liic_opt.max_cart_step <= 0.0) {
                std::cerr << "Error: --liic-max-step must be positive" << std::endl;
                return 1;
            }
            continue;
        }

        if (std::strcmp(argv[i], "--liic-verbose") == 0 && i + 1 < argc) {
            liic_opt.verbose = std::atoi(argv[++i]);
            if (liic_opt.verbose < 0) liic_opt.verbose = 0;
            continue;
        }

        // --- Legacy LIIC option aliases (hidden; do not document) ---
        if (std::strcmp(argv[i], "--iic-maxiter") == 0 && i + 1 < argc) {
            liic_opt.max_iter = std::atoi(argv[++i]);
            continue;
        }
        if (std::strcmp(argv[i], "--iic-tol") == 0 && i + 1 < argc) {
            liic_opt.tol = std::atof(argv[++i]);
            continue;
        }
        if (std::strcmp(argv[i], "--iic-damp") == 0 && i + 1 < argc) {
            liic_opt.damp = std::atof(argv[++i]);
            continue;
        }
        if (std::strcmp(argv[i], "--iic-max-step") == 0 && i + 1 < argc) {
            liic_opt.max_cart_step = std::atof(argv[++i]);
            continue;
        }

        // --- DM options ---
        if (std::strcmp(argv[i], "--dm-maxiter") == 0 && i + 1 < argc) {
            dm_opt.max_iter = std::atoi(argv[++i]);
            if (dm_opt.max_iter <= 0) {
                std::cerr << "Error: --dm-maxiter must be positive" << std::endl;
                return 1;
            }
            continue;
        }

        if (std::strcmp(argv[i], "--dm-step") == 0 && i + 1 < argc) {
            dm_opt.step = std::atof(argv[++i]);
            if (dm_opt.step <= 0.0) {
                std::cerr << "Error: --dm-step must be positive" << std::endl;
                return 1;
            }
            continue;
        }

        if (std::strcmp(argv[i], "--dm-tol") == 0 && i + 1 < argc) {
            dm_opt.tol = std::atof(argv[++i]);
            if (dm_opt.tol <= 0.0) {
                std::cerr << "Error: --dm-tol must be positive" << std::endl;
                return 1;
            }
            continue;
        }

        if (std::strcmp(argv[i], "--dm-max-step") == 0 && i + 1 < argc) {
            dm_opt.max_cart_step = std::atof(argv[++i]);
            if (dm_opt.max_cart_step <= 0.0) {
                std::cerr << "Error: --dm-max-step must be positive" << std::endl;
                return 1;
            }
            continue;
        }

        if (std::strcmp(argv[i], "--no-dm-fallback") == 0) {
            dm_fallback = false;
            continue;
        }

        // --- External engine options ---
        if ((std::strcmp(argv[i], "--engine-cmd") == 0 || std::strcmp(argv[i], "--engine") == 0 ||
             std::strcmp(argv[i], "--external-engine") == 0) &&
            i + 1 < argc) {
            engine_cfg.enabled = true;
            engine_cfg.cmd = argv[++i];
            continue;
        }

        if (std::strcmp(argv[i], "--engine-in") == 0 && i + 1 < argc) {
            engine_cfg.in_file = argv[++i];
            continue;
        }

        if (std::strcmp(argv[i], "--engine-out") == 0 && i + 1 < argc) {
            engine_cfg.out_file = argv[++i];
            continue;
        }

        if (std::strcmp(argv[i], "--engine-units") == 0 && i + 1 < argc) {
            engine_cfg.units = argv[++i];
            continue;
        }

        if (std::strcmp(argv[i], "--engine-output-is-force") == 0) {
            engine_cfg.output_is_force = true;
            continue;
        }

        if (std::strcmp(argv[i], "--engine-spring") == 0 && i + 1 < argc) {
            engine_cfg.spring_constant = std::atof(argv[++i]);
            if (engine_cfg.spring_constant <= 0.0) {
                std::cerr << "Error: --engine-spring must be positive" << std::endl;
                return 1;
            }
            continue;
        }

        if (std::strcmp(argv[i], "--engine-run-every") == 0 && i + 1 < argc) {
            engine_cfg.run_every = std::atoi(argv[++i]);
            if (engine_cfg.run_every <= 0) {
                std::cerr << "Error: --engine-run-every must be positive" << std::endl;
                return 1;
            }
            continue;
        }

        if (std::strcmp(argv[i], "--engine-keep-cycle-files") == 0) {
            engine_cfg.keep_cycle_files = true;
            continue;
        }

        // Hidden convenience aliases (kept for personal scripts)
        if (std::strcmp(argv[i], "--engine-vector") == 0 && i + 1 < argc) {
            const std::string v = argv[++i];
            if (v == "gradient" || v == "grad") {
                engine_cfg.output_is_force = false;
            } else if (v == "force" || v == "forces") {
                engine_cfg.output_is_force = true;
            } else {
                std::cerr << "Error: --engine-vector must be 'gradient' or 'force'" << std::endl;
                return 1;
            }
            continue;
        }

        if (std::strcmp(argv[i], "--engine-every") == 0 && i + 1 < argc) {
            engine_cfg.run_every = std::atoi(argv[++i]);
            continue;
        }

        if (std::strcmp(argv[i], "--engine-keep-files") == 0) {
            engine_cfg.keep_cycle_files = true;
            continue;
        }

        // --- Positional args ---
        if (argv[i][0] != '-') {
            if (initial_file.empty()) {
                initial_file = argv[i];
                continue;
            }
            if (final_file.empty()) {
                final_file = argv[i];
                continue;
            }
            std::cerr << "Error: too many input files specified" << std::endl;
            return 1;
        }

        std::cerr << "Error: unknown option " << argv[i] << std::endl;
        printUsage(argv[0]);
        return 1;
    }

    if (initial_file.empty() || final_file.empty()) {
        std::cerr << "Error: must specify both initial and final structure files" << std::endl;
        printUsage(argv[0]);
        return 1;
    }

    // Canonicalize legacy method names (do not advertise)
    if (method == "iic") {
        std::cerr << "Warning: method 'iic' is deprecated; using 'liic'." << std::endl;
        method = "liic";
    }
    if (method == "neb-iic") {
        std::cerr << "Warning: method 'neb-iic' is deprecated; using 'neb-liic'." << std::endl;
        method = "neb-liic";
    }

    neb::Method method_enum = neb::Method::NEB;
    if (method == "lic") {
        method_enum = neb::Method::LIC;
    } else if (method == "liic") {
        method_enum = neb::Method::LIIC;
    } else if (method == "dm") {
        method_enum = neb::Method::DM;
    } else if (method == "neb") {
        method_enum = neb::Method::NEB;
    } else if (method == "neb-liic") {
        method_enum = neb::Method::NEB_LIIC;
    } else {
        std::cerr << "Error: unknown method '" << method << "'" << std::endl;
        return 1;
    }

    std::cout
        << "NEB / Path interpolation tool\n"
        << "=============================\n"
        << "Method: " << method << "\n"
        << "Intermediate images: " << num_images << "\n"
        << "Alignment: " << (use_alignment ? "enabled (calc_rmsd_xyz)" : "disabled") << "\n"
        << "RMSD executable: " << rmsd_executable << "\n"
        << "Initial structure: " << initial_file << "\n"
        << "Final structure: " << final_file << "\n"
        << "Output mode: " << output_mode << "\n"
        << "Output prefix: " << (prefix.empty() ? "(none)" : prefix) << "\n"
        << "External engine: " << (engine_cfg.enabled ? "enabled" : "disabled") << "\n"
        << std::endl;

    if (engine_cfg.enabled) {
        std::cout
            << "  Engine cmd: " << engine_cfg.cmd << "\n"
            << "  Engine I/O: in='" << engine_cfg.in_file << "', out='" << engine_cfg.out_file << "'\n"
            << "  Engine units label: " << engine_cfg.units << "\n"
            << "  Output vector: " << (engine_cfg.output_is_force ? "force" : "gradient") << "\n"
            << "  Spring constant (k): " << engine_cfg.spring_constant << "\n"
            << "  Run every: " << engine_cfg.run_every << " cycle(s)\n"
            << "  Keep per-cycle files: " << (engine_cfg.keep_cycle_files ? "yes" : "no") << "\n"
            << std::endl;

        if (method != "neb" && method != "neb-liic") {
            std::cout
                << "  Note: external engine is only used by NEB methods; current method '" << method
                << "' will ignore it.\n"
                << std::endl;
        }
    }

    if (method == "liic" || method == "neb-liic") {
        std::cout
            << "LIIC options: bond_factor=" << liic_opt.bond_factor
            << ", fd_step=" << liic_opt.fd_step
            << ", ev_thresh=" << liic_opt.ev_thresh
            << ", liic_maxiter=" << liic_opt.max_iter
            << ", liic_tol=" << liic_opt.tol
            << ", liic_damp=" << liic_opt.damp
            << ", liic_max_step=" << liic_opt.max_cart_step
            << ", liic_verbose=" << liic_opt.verbose
            << "\n";
        std::cout << "DM fallback: " << (dm_fallback ? "enabled" : "disabled") << "\n";
    }

    if (method == "dm") {
        std::cout
            << "DM options: dm_maxiter=" << dm_opt.max_iter
            << ", dm_step=" << dm_opt.step
            << ", dm_tol=" << dm_opt.tol
            << ", dm_max_step=" << dm_opt.max_cart_step
            << "\n";
    }

    neb::NEBDriver driver(num_images, step_size, conv_threshold, max_iterations, use_alignment, rmsd_executable);
    driver.setLIICOptions(liic_opt);
    driver.setDMOptions(dm_opt);
    driver.setDMFallback(dm_fallback);
    driver.setExternalEngine(engine_cfg);

    std::string err;
    if (!driver.setStructuresFromFiles(initial_file, final_file, &err)) {
        std::cerr << err << std::endl;
        return 1;
    }

    if (!driver.run(method_enum, &err)) {
        if (!err.empty()) {
            std::cerr << err << std::endl;
        } else {
            std::cerr << "Error: failed to generate path." << std::endl;
        }
        return 1;
    }

    const bool multiframe = (output_mode == "multiframe");
    if (!driver.writeResults(prefix, multiframe)) {
        std::cerr << "Error: failed to write output files" << std::endl;
        return 1;
    }

    std::cout << "\nInterpolation completed successfully!" << std::endl;
    return 0;
}
