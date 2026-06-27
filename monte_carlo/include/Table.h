#pragma once

#include "constants.h"
#include "SimConfig.h"
#include "mdarray.h"

class Table
{
public:
    Table(const SimConfig& cfg)
    : n_neutrals(cfg.get_n_neutrals()),
      table(n_neutrals, cn::table_dim, cn::table_dim),
      bs(n_neutrals, cn::table_dim),
      chi_table(n_neutrals, chi_ta_dim, chi_u_dim),
      bcon_table(n_neutrals, bcon_ta_dim)
    {
        init_table();
        init_bcutoff_lookup(cfg);
        init_chi_table(cfg);   // uses bcutoff(), so must run after init_table()
        init_rates(cfg);
    }

    // cutoff impact parameter for neutral ineut at tangent-alpha tan_a,
    // for the scattering-angle cutoff chic
    double bcutoff(double tan_a, int ineut, double chic) const;
    double bcutoff_lookup(double tan_a, double log_ta, int ineut, double chic) const;

    // bilinear lookup of the polarization scattering angle on a dedicated
    // (log10(ta), u = bn/bcon) grid built once at construction. bcon is the
    // dimensionless cutoff from bcutoff(); u in [0,1] always covers the
    // polarization range, so no exact fallback is needed.
    double chi_lookup(double tan_a, int ineut, double bn, double bcon) const;
    double chi_lookup_u(double log_ta, int ineut, double u) const;

    // collision-rate kernel and time-step scale. these are the only derived
    // quantities that need the scattering-angle table (via bcutoff); everything
    // else config-derived lives in SimConfig::Derived.
    double get_kmax() const { return kmax; }
    double get_alfa() const { return alfa; }

private:
    // dedicated polarization scattering-angle table: axes are log10(ta) and
    // u = bn/bcon in [0,1]. Generous log10(ta) range + clamping in chi_lookup
    // keep out-of-range values safe (tighten the range once ta is measured).
    static constexpr int chi_ta_dim = 256;
    static constexpr int chi_u_dim = 512;
    static constexpr double chi_log_ta_min = -6.0;
    static constexpr double chi_log_ta_max =  4.0;
    static constexpr double chi_dlog_ta = (chi_log_ta_max - chi_log_ta_min) / (chi_ta_dim - 1);
    static constexpr double chi_inv_dlog_ta = 1.0 / chi_dlog_ta;
    static constexpr double chi_inv_du = chi_u_dim - 1;   // 1/du, du = 1/(dim-1)
    static constexpr int bcon_ta_dim = 2048;
    static constexpr double bcon_dlog_ta = (chi_log_ta_max - chi_log_ta_min) / (bcon_ta_dim - 1);
    static constexpr double bcon_inv_dlog_ta = 1.0 / bcon_dlog_ta;

    int n_neutrals;
    array3d<double> table;
    array2d<double> bs;
    array3d<double> chi_table;
    array2d<double> bcon_table;

    double kmax = 0.0;
    double alfa = 0.0;

    void init_table();
    void init_bcutoff_lookup(const SimConfig& cfg);
    void init_chi_table(const SimConfig& cfg);
    void init_rates(const SimConfig& cfg);

    double bintp(double chic, int ineut, int imin2, int indx) const;
};
