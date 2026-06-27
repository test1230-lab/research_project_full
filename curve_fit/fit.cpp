#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <span>
#include <print>

#include "./include/npy.hpp"

#include "./include/Chebyshev.h"
#include "./include/Legendre.h"
#include "./include/TailFit.h"

//TODO: these should not be hardcoded
constexpr double vmin_c = -3000.0;
constexpr double vmax_c = -vmin_c;
constexpr double v_trans = 500.0;
constexpr int cheb_deg = 45;
constexpr int leg_deg = 50;

double smoothstep(double x) 
{
    x = std::clamp(x, 0.0, 1.0);
    return x*x*x*(10.0 + x*(6.0*x - 15.0));
}

//arg is path of input file
int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cerr << "Incorrect number of args.\n The argument is the path of the input file.\n";
        return 1;
    }

    const std::string input_path{argv[1]};

    npy::npy_data npy_in = npy::read_npy<double>(input_path);

    const std::vector<unsigned long> shape = npy_in.shape;

    if (shape.size() != 2 || shape[1] != 2)
    {
        std::cerr << "Expected a 2D array with 2 columns.\n";
        return 1;
    }

    // The (i, j) indexing below assumes C (row-major) order. Bail out on a
    // Fortran-ordered array rather than silently reading it transposed.
    if (npy_in.fortran_order)
    {
        std::cerr << "Expected a C-ordered array, but the input is Fortran-ordered.\n";
        return 1;
    }

    const int rows = shape[0];
    const int cols = shape[1];

    std::vector<double> v(rows);
    std::vector<double> log_g(rows);

    for (int i = 0; i < rows; i++)
    {
        log_g[i] = std::log(std::max(1e-14, npy_in.data[i*cols + 1]));
        v[i] = npy_in.data[i*cols];
    }


    //const double vmin = *std::ranges::min_element(v);
    //const double vmax = *std::ranges::max_element(v);

    const double vmin = v.front();
    const double vmax = v.back();


    // find the [lo, hi) index range where vmin_c <= v < vmax_c
    //iterators
    const auto lo = std::ranges::lower_bound(v, vmin_c);   // first element >= vmin_c
    const auto hi = std::ranges::upper_bound(v, vmax_c);   // first element  > vmax_c

    const int offset = lo - v.begin(); // offset of central region
    const int count = hi - lo; // number of central elements

    std::span<const double> v_full{v};
    std::span<const double> vc = v_full.subspan(offset, count);
    std::span<const double> vl = v_full.subspan(0, offset);
    std::span<const double> vr = v_full.subspan(offset + count);

    // matching slices of log_g
    std::span<const double> log_g_full{log_g};
    std::span<const double> log_g_c = log_g_full.subspan(offset, count);
    std::span<const double> log_g_l = log_g_full.subspan(0, offset);
    std::span<const double> log_g_r = log_g_full.subspan(offset + count);

    //fit on log(g) in the central region
    Chebyshev cheb = Chebyshev::fit(vmin_c, vmax_c, cheb_deg, vc, log_g_c);

    //TODO: anchor so there is no discontinuty
    TailFit left_tail = TailFit::fit(vmin, vmin_c, vl, log_g_l);
    TailFit right_tail = TailFit::fit(vmax_c, vmax, vr, log_g_r);

    //maybe use a different container and filling method, idk
    std::vector<double> fitted(v.size());

    //stitch the tail fits with the central fit

    //v_trans is a velocity band width (m/s), not a sample count, so convert it
    //to index ranges via the velocity grid.
    const int n = static_cast<int>(v.size());
    const int lt_end = std::ranges::lower_bound(v, vmin_c + v_trans) - v.begin();
    const int rt_start = std::max<int>(lt_end, std::ranges::lower_bound(v, vmax_c - v_trans) - v.begin());

    //left tail fit
    for (int i = 0; i < offset; i++)
    {
        fitted[i] = left_tail(v[i]);
    }

    //transition region from left tail fit to center fit
    for (int i = offset; i < lt_end; i++)
    {
        const double s = smoothstep((v[i] - vmin_c)/v_trans);
        const double a = (1.0 - s)*left_tail(v[i]);
        const double b = s*cheb(v[i]);

        fitted[i] = a + b;
    }

    //central fit
    for (int i = lt_end; i < rt_start; i++)
    {
        fitted[i] = cheb(v[i]);
    }

    //transition from central fit to right tail fit
    for (int i = rt_start; i < offset + count; i++)
    {
        const double s = smoothstep((v[i] - (vmax_c - v_trans))/v_trans);
        const double a = (1.0 - s)*cheb(v[i]);
        const double b = s*right_tail(v[i]);

        fitted[i] = a + b;
    }

    //right tail fit
    for (int i = offset + count; i < n; i++)
    {
        fitted[i] = right_tail(v[i]);
    }

    //apply std::exp to each element to undo the log
    std::ranges::transform(fitted, fitted.begin(), [](double x){ return std::exp(x); });

    //final fit
    Legendre leg = Legendre::fit(vmin, vmax, leg_deg, v, fitted);
    std::vector<double> coeffs = leg.get_coeffs_stdvec();
  
    //print out the even coeffs
    for (int i = 0; i < coeffs.size(); i++)
    {
        if (i % 2 == 0) 
        { 
            std::print("Coefficient {}: {: .6e}\n", i, coeffs[i]); 
        }
    }

    return 0;
}