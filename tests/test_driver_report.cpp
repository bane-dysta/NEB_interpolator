#include <iostream>
#include <string>
#include <vector>

#include "neb_driver.h"

static void require(bool cond, const std::string& msg) {
    if (!cond) {
        std::cerr << "[FAIL] " << msg << std::endl;
        std::exit(1);
    }
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::cerr << "Usage: test_driver_report <bondshift_init.xyz> <bondshift_final.xyz> <water.xyz> <water_rotated.xyz>" << std::endl;
        return 2;
    }

    const std::string bondshift_init = argv[1];
    const std::string bondshift_final = argv[2];
    const std::string water_init = argv[3];
    const std::string water_final = argv[4];
    std::string err;

    // Exact LIC path: no fallback.
    {
        neb::NEBDriver driver(3, 0.0001, 0.01, 100, false, "./calc_rmsd_xyz");
        require(driver.setStructuresFromFiles(bondshift_init, bondshift_final, &err), err);
        require(driver.run(neb::Method::LIC, &err), "LIC run failed: " + err);
        const auto& report = driver.lastRunReport();
        require(report.status == neb::RunStatus::GENERATED_EXACT, "LIC should report GENERATED_EXACT");
        require(report.actual_path_method == neb::PathMethod::LIC, "LIC actual path method mismatch");
        require(report.attempted_path_methods.size() == 1 && report.attempted_path_methods[0] == neb::PathMethod::LIC,
                "LIC attempt chain mismatch");
        require(driver.lastRunUsedRequestedPathMethod(), "LIC should realize the requested path method exactly");
    }

    // Force LIIC to fail so the report must expose DM fallback.
    {
        neb::NEBDriver driver(3, 0.0001, 0.01, 100, false, "./calc_rmsd_xyz");
        ICInterp::LIICOptions liic_opt;
        liic_opt.ev_thresh = 1e9; // ensure no DLC modes survive
        driver.setLIICOptions(liic_opt);
        require(driver.setStructuresFromFiles(bondshift_init, bondshift_final, &err), err);
        require(driver.run(neb::Method::LIIC, &err), "LIIC run should still generate a fallback path: " + err);
        const auto& report = driver.lastRunReport();
        require(report.status == neb::RunStatus::GENERATED_WITH_FALLBACK,
                "LIIC fallback should report GENERATED_WITH_FALLBACK");
        require(report.requested_path_method == neb::PathMethod::LIIC, "LIIC requested path method mismatch");
        require(report.actual_path_method == neb::PathMethod::DM, "LIIC fallback should land on DM");
        require(report.attempted_path_methods.size() == 2, "LIIC fallback chain length mismatch");
        require(report.attempted_path_methods[0] == neb::PathMethod::LIIC, "LIIC fallback chain should start with LIIC");
        require(report.attempted_path_methods[1] == neb::PathMethod::DM, "LIIC fallback chain should end with DM");
        require(!driver.lastRunUsedRequestedPathMethod(), "LIIC fallback should not count as exact realization");
    }

    // Force DM to fail so the report must expose LIC fallback.
    {
        neb::NEBDriver driver(3, 0.0001, 0.01, 100, false, "./calc_rmsd_xyz");
        ICInterp::DMOptions dm_opt;
        dm_opt.max_iter = 1;
        dm_opt.tol = 1e-12;
        driver.setDMOptions(dm_opt);
        require(driver.setStructuresFromFiles(water_init, water_final, &err), err);
        require(driver.run(neb::Method::DM, &err), "DM run should still generate a fallback path: " + err);
        const auto& report = driver.lastRunReport();
        require(report.status == neb::RunStatus::GENERATED_WITH_FALLBACK,
                "DM fallback should report GENERATED_WITH_FALLBACK");
        require(report.requested_path_method == neb::PathMethod::DM, "DM requested path method mismatch");
        require(report.actual_path_method == neb::PathMethod::LIC, "DM fallback should land on LIC");
        require(report.attempted_path_methods.size() == 2, "DM fallback chain length mismatch");
        require(report.attempted_path_methods[0] == neb::PathMethod::DM, "DM fallback chain should start with DM");
        require(report.attempted_path_methods[1] == neb::PathMethod::LIC, "DM fallback chain should end with LIC");
        require(!driver.lastRunUsedRequestedPathMethod(), "DM fallback should not count as exact realization");
    }

    std::cout << "All driver report checks passed." << std::endl;
    return 0;
}
