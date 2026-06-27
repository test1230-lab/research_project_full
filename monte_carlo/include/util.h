#pragma once

#include <array>


//idk where to put these structs/functions

inline double sqr(double x)
{
    return x*x;
}

using vec3 = std::array<double, 3>;

struct NeutralInfo
{
    //--- from config ---
    double mass_amu;
    double rel_abundance;
    double rm;
    double ar;
    double br;
    double eps;//eV

    //derived
    double mu = 0.0;
    double beta = 0.0;
    double chi_cutoff = 0.0;
    double mass = 0.0;
    double log_eps_over_bk = 0.0;
};

struct RCS_Result
{
    double rcb, rmag, b0;
};
