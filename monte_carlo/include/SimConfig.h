#pragma once

#include <string>
#include <vector>
#include <print>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <stdexcept>
#include <format>

#include "toml.hpp"

#include "util.h"
#include "constants.h"

class SimConfig
{
public:
    // quantities derived from the raw config that depend only on config inputs
    // (not on the scattering-angle table). table-dependent quantities (kmax,
    // alfa) live in Table.
    struct Derived
    {
        double vd = 0.0;       // ExB drift speed
        double migrm = 0.0;    // ion mass [g]
        double tabun = 0.0;    // total relative abundance of all neutrals
        double dlvpar = 0.0;   // v-parallel bin width
        double dlvper = 0.0;   // v-perp bin width
        double scvpar = 0.0;   // 1/dlvpar
        double scvper = 0.0;   // 1/dlvper
        double gmax = 0.0;     // max relative speed (last neutral, as in fortran)
    };

    SimConfig() = default;

    void init(const std::string& cfg_filename)
    {
        read_config(cfg_filename);
        compute_derived();
    }

    void print_config_summary()
    {
        std::print("====== Monte-Carlo collision simulation ======\n");
        std::print("     output prefix:   \t{}\n", out_filename);
        std::print("     neutral temp:    \t{}[K]\n", neutral_temp);
        std::print("     electric field:  \t{}[mV/m]\n", electric_field);
        std::print("     ion mass:        \t{}[amu]\n", ion_mass);
        std::print("     neutrals:        \t{}\n", neutrals.size());
        std::print("     total collisions:\t{:.1f}M\n", total_collisions/1e6);
        std::print("     threads:         \t{}\n",  n_threads);
        std::print("==============================================\n");
    }

    int get_nthreads() const { return n_threads; }
    uint64_t get_total_collisions() const { return total_collisions; }
    double get_neutral_temp() const { return neutral_temp; }
    double get_electric_field() const { return electric_field; }
    double get_ion_mass() const { return ion_mass; }
    std::size_t get_n_neutrals() const  { return neutrals.size(); }
    const std::string& get_out_filename() const { return out_filename; }

    double get_dlvpar() const {return der.dlvpar; }
    double get_dlvper() const {return der.dlvper; }

    // neutrals with their derived fields (eps in erg, mass, mu, beta, chi_cutoff)
    // filled in by compute_derived().
    const std::vector<NeutralInfo>& get_neutrals() const { return neutrals; }
    const Derived& get_derived() const { return der; }

private:
    std::string out_filename;
    double neutral_temp; // neutral temp in kelvin
    double electric_field; // [mV/m]
    double ion_mass; // ion mass in amu
    std::vector<NeutralInfo> neutrals;
    uint64_t total_collisions;
    int n_threads;

    Derived der{};

    void read_config(const std::string& cfg_filename)
    {
        const auto data = toml::parse(cfg_filename);

        out_filename = toml::find_or<std::string>(data, "out_filename", "output.dat");
        neutral_temp = toml::find<double>(data, "neutral_temp");
        electric_field = toml::find<double>(data, "electric_field");
        total_collisions = toml::find<uint64_t>(data, "num_collisions");
        n_threads = toml::find<int>(data, "num_threads");
        ion_mass = toml::find<double>(data, "ion_mass");

        //checking of the config is valid
        if (neutral_temp <= 0.0)
        {
            const std::string s = std::format("Neutral Temp: {:.3e}", neutral_temp);
            throw std::runtime_error("Neutral temp <= zero " + s);
        }

        if (electric_field < 0.0)
        {
            const std::string s = std::format("Electric Field: {:.3e}", electric_field);
            throw std::runtime_error("Electric field magnitude < zero. " + s);
        }

        if (total_collisions == 0)
        {
            const std::string s = std::format("Total Collisions: {}", total_collisions);
            throw std::runtime_error("Total collisions must be greater than zero. " + s);
        }

        if (n_threads < 1)
        {
            const std::string s = std::format("Thread Count: {}", n_threads);
            throw std::runtime_error("Number of threads < one. " + s);
        }

        if (ion_mass <= 0.0)
        {
            const std::string s = std::format("Ion Mass: {:.6f}[amu]", ion_mass);
            throw std::runtime_error("Ion mass <= zero. " + s);
        }

        const auto& neutral_array = toml::find<toml::array>(data, "neutrals");

        for (const auto& neutral_info : neutral_array)
        {
            NeutralInfo ni;
            ni.mass_amu = toml::find<double>(neutral_info, "mass");
            ni.rel_abundance = toml::find<double>(neutral_info, "abundance");
            ni.rm = toml::find<double>(neutral_info, "rm");
            ni.ar = toml::find<double>(neutral_info, "ar");
            ni.br = toml::find<double>(neutral_info, "br");
            ni.eps = toml::find<double>(neutral_info, "epsilon");
            neutrals.push_back(ni);
        }
    }

    // fills in the derived NeutralInfo fields and the config-only Derived scalars.
    void compute_derived()
    {
        const double ev = electric_field*(1e-5)*(3.336e-3);
        der.vd = ev*cn::c/cn::bmag;
        der.migrm = ion_mass*cn::pr_mass;

        double tabun = 0.0;
        double sma = 0.0;
        for (NeutralInfo& n : neutrals)
        {
            n.eps *= 1.6e-12;
            n.log_eps_over_bk = std::log10(n.eps/cn::bk);
            n.mass = n.mass_amu*cn::pr_mass;
            n.mu = cn::pr_mass * n.mass_amu*ion_mass/(n.mass_amu + ion_mass);
            tabun += n.rel_abundance;
            n.beta = std::sqrt(cn::bk*neutral_temp/n.mass);
            sma += n.rel_abundance*n.mass;
        }

        der.tabun = tabun;

        // range of Vperp and Vpar
        const double avmngr = sma/tabun;

        const double a = 2.0*cn::bk*neutral_temp/der.migrm;
        const double b = (2.0/3.0)*der.vd*der.vd*avmngr/der.migrm;
        const double vpar = std::sqrt(a + b);

        der.dlvpar = 3.0*vpar/(cn::nparmx + 1);
        der.dlvper = der.dlvpar;

        // CHI cutoff value per neutral. gmax ends up holding the last neutral's
        // value, matching the fortran reference.
        double gmax = 0.0;
        for (NeutralInfo& n : neutrals)
        {
            const double dp2 = der.vd*der.vd/(2.0*n.beta*n.beta);
            const double vith2 = 2.0*n.beta*n.beta*(1.0 + 2.0*dp2/3.0)*n.mass_amu/ion_mass;
            const double gav = std::sqrt(vith2 + 2.0*n.beta*n.beta + der.vd*der.vd);
            gmax = 5.0*gav;
            const double dlvo2 = 0.5*std::min(der.dlvpar, der.dlvper);
            const double m2 = n.mass_amu/(ion_mass + n.mass_amu);
            const double v = dlvo2/gmax;
            n.chi_cutoff = v/m2;
        }

        der.gmax = gmax;
        der.scvpar = 1.0/der.dlvpar;
        der.scvper = 1.0/der.dlvper;
    }
};
