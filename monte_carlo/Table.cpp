#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#include "Table.h"
#include "SimConfig.h"
#include "util.h"
#include "scattering.h"


void Table::init_table()
{
    /*  Initialize tables to determine bcutoff.  One table is created for
    each neutral constituent.
    element TABLE(I,J,K):  I = index of neutral constituent
                            J = tangent alpha
                            K = b, the impact parameter          */

    for (int i = 0; i < n_neutrals; i++)
    {
        for (int j = 0; j < cn::table_dim; j++)
        {
            const double tan_a = std::pow(10.0, cn::dke*j - 1.84);
            RCS_Result res = scatter::rcs(true, tan_a, 0.0);
            bs[i, j] = std::max(res.b0, 0.5);
            for (int k = 0; k < cn::table_dim; k++)
            {
                const double b = cn::db*k + bs[i, j];
                res = scatter::rcs(false, tan_a, b);

                const double sa = scatter::scattering_angle<25>(res.rmag, res.rcb, b, tan_a);
                table[i, j, k] = std::abs(sa);
            }
        }
    }
}

void Table::init_rates(const SimConfig& cfg)
{
    /*  Determine kmax, the largest collision-rate kernel over all neutrals,
    sampled at the maximum relative speed gmax and a low reference speed gl.
    alfa is the mean free path scale derived from kmax.  */

    const auto& d = cfg.get_derived();
    const auto& neutrals = cfg.get_neutrals();
    const double gmax = d.gmax;

    constexpr double gl = 10.0;
    for (int i = 0; i < n_neutrals; i++)
    {
        const NeutralInfo& nu = neutrals[i];

        const double ta_max = 0.5*nu.mu*gmax*gmax/nu.eps;
        const double bcomxn = bcutoff(ta_max, i, nu.chi_cutoff);
        const double bco_max = bcomxn*nu.rm;
        const double kpmax = gmax*cn::pi*bco_max*bco_max;

        const double tal = 0.5*nu.mu*gl*gl/nu.eps;
        const double bcxn = bcutoff(tal, i, nu.chi_cutoff);
        const double bcax = bcxn*nu.rm;
        const double kpl = gl*cn::pi*bcax*bcax;

        const double gmax_factor = 0.5*nu.mu*gmax*gmax/cn::bk;
        kmax = std::max({kpmax, kpl, kmax, 2.0*gmax*sqr(nu.ar - nu.br*std::log10(gmax_factor))});
    }

    kmax *= 1.01;
    alfa = ((cn::ec*cn::bmag)/(d.migrm*cn::c))/(cn::totnn*kmax);
}

double Table::bintp(double chic, int ineut, int imin2, int indx) const
{
    const double x = std::abs(chic);
    const double x0 = table[ineut, imin2, indx];       // TABLE(.,IMIN2,INDX)
    const double x1 = table[ineut, imin2, indx + 1];   // TABLE(.,.,INDX+1)
    const double xm1 = table[ineut, imin2, indx - 1];   // TABLE(.,.,INDX-1)

    const double fofx = (x - x0 )*(x - x1)/((xm1 - x0 )*(xm1 - x1))*cn::db*(indx - 1)
                      + (x - xm1)*(x - x1)/((x0  - xm1)*(x0  - x1))*cn::db*indx
                      + (x - xm1)*(x - x0)/((x1  - xm1)*(x1  - x0))*cn::db*(indx + 1);

    return fofx + bs[ineut, imin2];
}

double Table::bcutoff(double tan_a, int ineut, double chic) const
{
    // https://mathworld.wolfram.com/NaturalLogarithmof10.html
    constexpr double ln10 = 2.302585092994045684;

    const double abschic = std::abs(chic);
    const double log10tan_a = std::log10(tan_a);

    // 0-indexed tan-alpha row for this tangent alpha (Fortran's IMIN1 minus 1)
    const int imin1 = static_cast<int>(1.0 + 1.84/cn::dke + log10tan_a/cn::dke) - 1;

    int imin = cn::table_dim - 1;                  // 24  (Fortran default 25)
    if (imin1 <= cn::table_dim - 1) 
    {
        imin = imin1; // Fortran IMIN1 <= 25
    }

    if (imin1 < 0)
    {
        imin = 0; // Fortran IMIN1 <= 0
    }      

    int imin2 = imin;

    // first column (descending b) whose angle exceeds |chic|; returns -1 (and
    // warns) if the cutoff is off the small-b end of the row
    auto find_indx = [&](int row) -> int {
        int k;
        for (k = cn::table_dim - 2; k >= 0; --k)
        {
            if (row >= 0 && (row <= cn::table_dim - 1) && (table[ineut, row, k] > abschic))
            {
                break;
            }
        }

        if (k < 0)
        {
            std::cerr << " bcf-- chi cutoff too large, not in table\n";
        }

        return k;
    };

    // bintp touches indx-1..indx+1, so indx must be in [1,23] and row in [0,24].
    // The Fortran read past the table on these error paths and STOPped afterward;
    // we throw before any out-of-bounds access (same net outcome).
    auto guarded_bintp = [&](int row, int col) -> double {
        if (row < 0 || row > cn::table_dim - 1)
        {
            throw std::runtime_error("bcf: imin2 out of range = " + std::to_string(row));
        }

        if (col < 1 || col > cn::table_dim - 2)
        {
            throw std::runtime_error("bcf: index out of range = " + std::to_string(col));
        }

        return bintp(chic, ineut, row, col);
    };

    // advance the tan-alpha row until the largest-b angle drops below |chic|
    while (imin2 <= cn::table_dim - 1 && (table[ineut, imin2, cn::table_dim - 1] >= abschic))
    {
        imin2++;
    }

    int indx = find_indx(imin2);
    double bcutof;

    if (imin1 < 0)
    {
        // polarization-dominated: extrapolate analytically (1/4 power)
        bcutof = guarded_bintp(imin2, indx);

        //bcutof *= std::pow(std::pow(10.0, -1.84 + imin2*cn::dke)/tan_a, 0.25);
        bcutof *= std::exp(0.25*ln10*(imin2*cn::dke - log10tan_a - 1.84));
    }
    else if (imin1 < cn::table_dim - 1)                // Fortran IMIN1 < 25
    {
        bcutof = guarded_bintp(imin2, indx);
        if (imin2 != imin1)
        {
            // row was advanced (column out of range originally) -> analytic interp
            //bcutof *= std::pow(std::pow(10.0, -1.84 + imin2*cn::dke)/tan_a, 0.25);

            bcutof *= std::exp(0.25*ln10*(imin2*cn::dke - log10tan_a - 1.84));
        }
        else
        {
            // two-point interpolation across adjacent tan-alpha (energy) rows
            imin2++;
            indx = find_indx(imin2);
            const double bcuto1 = guarded_bintp(imin2, indx);
            const double p = (1.84 + log10tan_a)/cn::dke - imin1;
            bcutof = (1.0 - p)*bcutof + p*bcuto1;
        }
    }
    else // Fortran IMIN1 >= 25
    {
        // 1/r^12 regime: scaling law (1/12 power)
        bcutof = guarded_bintp(imin2, indx);
        if (imin1 > cn::table_dim - 1) // Fortran IMIN1 > 25
        {
            //bcutof *= std::pow(std::pow(10.0, -1.84 + 24.0*cn::dke)/tan_a, 1.0/12.0);
            bcutof *= std::exp((1.0/12.0)*ln10*(24.0*cn::dke - log10tan_a - 1.84));
        }
    }

    if (imin2 > cn::table_dim - 1 || imin2 < 0)
    {
        throw std::out_of_range(
            "Table::bcutoff: tan-alpha row index imin2 = " + std::to_string(imin2)
            + " is outside the valid range [0, " + std::to_string(cn::table_dim - 1) + "]");
    }

    if (indx > cn::table_dim - 2 || indx < 1)
    {
        throw std::out_of_range(
            "Table::bcutoff: impact-parameter column index indx = " + std::to_string(indx)
            + " is outside the valid range [1, " + std::to_string(cn::table_dim - 2) + "]");
    }

    return bcutof;
}

void Table::init_bcutoff_lookup(const SimConfig& cfg)
{
    const auto& neutrals = cfg.get_neutrals();

    for (int ineut = 0; ineut < n_neutrals; ineut++)
    {
        const NeutralInfo& nu = neutrals[ineut];

        for (int ita = 0; ita < bcon_ta_dim; ita++)
        {
            const double log_ta = chi_log_ta_min + bcon_dlog_ta*ita;
            const double ta = std::pow(10.0, log_ta);
            bcon_table[ineut, ita] = bcutoff(ta, ineut, nu.chi_cutoff);
        }
    }
}

double Table::bcutoff_lookup(double tan_a, double log_ta, int ineut, double chic) const
{
    if (!std::isfinite(log_ta) || log_ta < chi_log_ta_min || log_ta > chi_log_ta_max)
    {
        return bcutoff(tan_a, ineut, chic);
    }

    const double x = (log_ta - chi_log_ta_min)*bcon_inv_dlog_ta;
    const int i0 = static_cast<int>(x);
    if (i0 >= bcon_ta_dim - 1)
    {
        return bcon_table[ineut, bcon_ta_dim - 1];
    }

    const int i1 = i0 + 1;
    const double fx = x - i0;
    const double b0 = bcon_table[ineut, i0];
    const double b1 = bcon_table[ineut, i1];
    const double rel_jump = std::abs(b1 - b0)/std::max(std::max(std::abs(b0), std::abs(b1)), 1e-300);
    if (rel_jump > 5e-3)
    {
        return bcutoff(tan_a, ineut, chic);
    }

    return (1.0 - fx)*b0 + fx*b1;
}

void Table::init_chi_table(const SimConfig& cfg)
{
    // Tabulate the polarization scattering angle chi(ta, bn) once, on a
    // (log10(ta), u) grid where u = bn/bcon in [0,1]. Built with the exact
    // rcs + Gauss quadrature; the Monte Carlo loop then only interpolates.
    // Must run after init_table() (uses bcutoff, which needs table/bs).
    const auto& neutrals = cfg.get_neutrals();

    for (int ineut = 0; ineut < n_neutrals; ineut++)
    {
        const NeutralInfo& nu = neutrals[ineut];

        for (int ita = 0; ita < chi_ta_dim; ita++)
        {
            const double log_ta = chi_log_ta_min + chi_dlog_ta*ita;
            const double ta = std::pow(10.0, log_ta);
            const double bcon = bcutoff(ta, ineut, nu.chi_cutoff);

            for (int iu = 0; iu < chi_u_dim; iu++)
            {
                if (iu == 0)
                {
                    // limit as bn -> 0: head-on, chi -> pi. avoids degenerate
                    // root-finding at exactly b = 0.
                    chi_table[ineut, ita, iu] = cn::pi;
                    continue;
                }

                const double u = static_cast<double>(iu)/(chi_u_dim - 1);
                const double bn = u*bcon;

                const RCS_Result res = scatter::rcs(false, ta, bn);
                // Store |chi|: scattering_angle returns the signed deflection,
                // which is large-and-negative for attractive/orbiting impact
                // parameters. Only the magnitude is physical here because collis()
                // randomizes the azimuth (the Fortran TABLE likewise stores
                // ABS(CHI)). Storing the signed value would (a) cancel under
                // interpolation across the sign change and (b) be clamped away in
                // chi_lookup, zeroing every attractive collision.
                chi_table[ineut, ita, iu] = std::abs(scatter::scattering_angle<20>(res.rmag, res.rcb, bn, ta));
            }
        }
    }
}

double Table::chi_lookup(double tan_a, int ineut, double bn, double bcon) const
{
    if (tan_a <= 0.0 || bcon <= 0.0)
    {
        return 0.0;
    }

    return chi_lookup_u(std::log10(tan_a), ineut, bn/bcon);
}

double Table::chi_lookup_u(double log_ta, int ineut, double u) const
{
    if (!std::isfinite(log_ta))
    {
        return 0.0;
    }

    // fractional coordinate along log10(ta), clamped into the grid
    double x = (log_ta - chi_log_ta_min)*chi_inv_dlog_ta;
    x = std::clamp(x, 0.0, static_cast<double>(chi_ta_dim - 1));
    const int i0 = static_cast<int>(x);
    const int i1 = std::min(i0 + 1, chi_ta_dim - 1);
    const double fx = x - i0;

    // fractional coordinate along u = bn/bcon in [0,1]
    double y = std::clamp(u, 0.0, 1.0)*chi_inv_du;
    y = std::clamp(y, 0.0, static_cast<double>(chi_u_dim - 1));
    const int j0 = static_cast<int>(y);
    const int j1 = std::min(j0 + 1, chi_u_dim - 1);
    const double fy = y - j0;

    const double c00 = chi_table[ineut, i0, j0];
    const double c10 = chi_table[ineut, i1, j0];
    const double c01 = chi_table[ineut, i0, j1];
    const double c11 = chi_table[ineut, i1, j1];

    const double c0 = (1.0 - fx)*c00 + fx*c10;
    const double c1 = (1.0 - fx)*c01 + fx*c11;

    // chi_table holds |chi| (>= 0), which can exceed pi in the orbiting regime,
    // so only floor at 0 (guard interpolation undershoot) -- do NOT clamp to pi,
    // which would discard large-angle/orbiting deflections.
    return std::max((1.0 - fy)*c0 + fy*c1, 0.0);
}
