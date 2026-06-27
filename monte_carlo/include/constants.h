#pragma once

#include <numbers>

namespace cn
{
    constexpr int table_dim = 25;
    constexpr int nparmx = 23, npermx = 23;
    constexpr double pi = std::numbers::pi_v<double>;
    constexpr double c = 2.99792458e10; //cm/s
    constexpr double bk = 1.380658e-16; // erg/k
    constexpr double bmag = 0.5, totnn = 1e9, dke = 0.16, ec = 4.8e-10;
    constexpr double pr_mass = 1.67e-24, db = 0.1;
    constexpr int n_report = 5'000'000;

    constexpr double root_finding_tol = 1e-6; //fortran is 1e-5
}