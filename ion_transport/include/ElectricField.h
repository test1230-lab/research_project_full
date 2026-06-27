#pragma once

#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <ranges>
#include <string_view>
#include <filesystem>
#include <algorithm>
#include <utility>
#include <cmath>
#include <array>
#include <optional>
#include <print>
#include <boost/math/interpolators/makima.hpp>
#include <boost/math/interpolators/pchip.hpp>

#include "mdarray.h"

using boost::math::interpolators::makima;
using boost::math::interpolators::pchip;

class ElectricField
{
public:
    ElectricField(std::string directory)
    {
        this->dir = directory;
        this->n_files = num_files();
        coeffs.resize(n_angles, n_files, n_cols);
        ion_thermal_speeds.resize(n_angles, n_files);
        read_coeffs_new_fmt();
        compute_interp_coeffs();
    }

    double compute_dist_discrete(int electric_field, int aspect_angle, double x) const;
    double compute_dist(double electric_field, int aspect_angle, double v) const;
    double compute_integral(double electric_field, int aspect_angle) const; //check

private:
    static constexpr int n_cols = 26, n_angles = 10;
    static constexpr std::string_view delim{"\t\t"};
    static constexpr int order = 5;

    static constexpr std::array<double, 2*n_cols> leg_a = [] {
        std::array<double, 2*n_cols> a{};
        for (int i = 1; i < 2*(n_cols - 1); i++) a[i] = (2.0*i + 1.0)/(i + 1.0);
        return a;
    }();

    static constexpr std::array<double, 2*n_cols> leg_b = [] {
        std::array<double, 2*n_cols> b{};
        for (int i = 1; i < 2*(n_cols - 1); i++) b[i] = double(i)/(i + 1.0);
        return b;
    }();

    using Spline1 = makima<std::vector<double>>;
    using Spline2 = pchip<std::vector<double>>;

    int n_files;
    std::filesystem::path dir;
    std::vector<double> e_field_vals;

    //"coeffs" is indexed as such: coeffs[aspect angle idx][e field idx][i]
    //"ion_thermal_speeds" is indexed as such: ion_thermal_speeds[aspect angle idx][e field idx]
    array3d<double> coeffs;
    array2d<double> ion_thermal_speeds;

    std::array<std::array<std::optional<Spline1>, n_cols>, n_angles> e_interp;
    std::array<std::optional<Spline2>, n_angles> ion_speed_interp;

    int num_files() const;
    void read_coeffs();
    void read_coeffs_new_fmt();
    int extract_number(const std::filesystem::path& p) const;
    void compute_interp_coeffs();
    double eval_legendre_series(const std::array<double, n_cols>& coeffs, double x) const;        
};