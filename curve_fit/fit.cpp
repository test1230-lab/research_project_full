#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <span>
#include <utility>
#include <print>
#include <stdexcept>

#include "./include/npy.hpp"

#include "./include/Chebyshev.h"
#include "./include/Legendre.h"
#include "./include/TailFit.h"

constexpr int cheb_deg = 45;
constexpr int leg_deg = 50;
constexpr double bk = 1.380649e-23; //https://physics.nist.gov/cgi-bin/cuu/Value?k
constexpr double kg_amu = 1.66053906892e-27; //https://physics.nist.gov/cgi-bin/cuu/Value?Rukg
constexpr double ion_mass = 16.0*kg_amu; //for atomic oxygen

struct FitParams
{
    double vmin_c, vmax_c, v_trans, tail_cutoff;
};

double sqr(double x)
{
    return x*x;
}

double smoothstep(double x) 
{
    x = std::clamp(x, 0.0, 1.0);
    return x*x*x*(10.0 + x*(6.0*x - 15.0));
}

std::pair<std::vector<double>, std::vector<double>> read_npy(const std::string& in)
{
    npy::npy_data npy_in = npy::read_npy<double>(in);

    const std::vector<unsigned long> shape = npy_in.shape;

    if (shape.size() != 2 || shape[1] != 2)
    {
        throw std::runtime_error("Expected a 2D array with 2 columns");
    }

    if (npy_in.fortran_order)
    {
        throw std::runtime_error("Expected a C-ordered array, but the input is Fortran-ordered");
    }

    const std::size_t rows = shape[0];
    const std::size_t cols = shape[1];

    if (rows == 0)
    {
        throw std::runtime_error("Expected at least one input row");
    }

    std::vector<double> v(rows);
    std::vector<double> g(rows);

    for (int i = 0; i < rows; i++)
    {
        g[i] = npy_in.data[i*cols + 1];
        v[i] = npy_in.data[i*cols];
    }

    //I want to avoid a copy
    return {std::move(v), std::move(g)};
}

//stddev is the g-weighted spread of v, i.e. vtherm/sqrt(2)
FitParams get_fit_params(double stddev)
{
    FitParams params;

    params.vmax_c = 1.75*stddev; //1.75 is just a number I chose
    params.vmin_c = -params.vmax_c;
    params.v_trans = stddev/8.0; //also just a number I chose
    params.tail_cutoff = stddev*4.0; //TODO: find best value for this

    return params;
}

double compute_v_therm(const std::vector<double>& v, const std::vector<double>& g)
{
    const int N = static_cast<int>(v.size());
    const double dv = v[1] - v[0];

    const double density = dv*std::accumulate(g.begin(), g.end(), 0.0);

    double vavg = 0.0;
    for (int i = 0; i < N; i++)
    {
        vavg += v[i]*g[i];
    }
    vavg *= dv/density;

    double temp = 0.0;
    for (int i = 0; i < N; i++)
    {
        temp += g[i]*sqr(v[i] - vavg);
    }
    temp *= dv*ion_mass/(bk*density);
    const double vtherm = std::sqrt(2.0*bk*temp/ion_mass);

    return vtherm;
}

//arg is path of input file
int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cerr << "Incorrect number of args.\n The argument is the path of the input file.\n";
        return 1;
    }

    std::pair<std::vector<double>, std::vector<double>> in;

    try
    {
        in = read_npy(argv[1]);
    }
    catch(const std::exception& e)
    {
        std::cerr << "error in reading file. Threw: " << e.what() << '\n';
        return 1;
    }

    const double vtherm = compute_v_therm(in.first, in.second);

    const FitParams fparams = get_fit_params(vtherm/std::sqrt(2.0));
    const double vmin_c = fparams.vmin_c;
    const double vmax_c = fparams.vmax_c;
    const double v_trans = fparams.v_trans;
    const double tail_cutoff = fparams.tail_cutoff;

    std::ranges::transform(in.second, in.second.begin(), [](double x){ return std::log(std::max(1e-16, x)); });

    //trying to avoid a copy
    const std::vector<double> v = std::move(in.first);
    const std::vector<double> log_g = std::move(in.second);

    const double vmin = v.front();
    const double vmax = v.back();

    // find the central index range where vmin_c <= v <= vmax_c
    //these are iterators
    const auto lo = std::ranges::lower_bound(v, vmin_c);   // first element >= vmin_c
    const auto hi = std::ranges::upper_bound(v, vmax_c);   // first element  > vmax_c

    //find tail cuttoffs
    const auto lt = std::ranges::lower_bound(v, -tail_cutoff);
    const auto rt = std::ranges::upper_bound(v, tail_cutoff);

    const std::size_t offset = static_cast<std::size_t>(lo - v.begin()); // offset of central region
    const std::size_t count = static_cast<std::size_t>(hi - lo); // number of central elements

    const std::size_t l_tail_start = static_cast<std::size_t>(lt - v.begin());
    const std::size_t l_tail_count = offset - l_tail_start; // -tail_cutoff up to the left seam
    const std::size_t r_tail_count = static_cast<std::size_t>(rt - v.begin()) - (offset + count); // right seam up to +tail_cutoff

    std::span<const double> v_full{v};
    std::span<const double> vc = v_full.subspan(offset, count);
    std::span<const double> vl = v_full.subspan(l_tail_start, l_tail_count);
    std::span<const double> vr = v_full.subspan(offset + count, r_tail_count);

    // matching slices of log_g
    std::span<const double> log_g_full{log_g};
    std::span<const double> log_g_c = log_g_full.subspan(offset, count);
    std::span<const double> log_g_l = log_g_full.subspan(l_tail_start, l_tail_count);
    std::span<const double> log_g_r = log_g_full.subspan(offset + count, r_tail_count);

    //fit on log(g) in the central region
    Chebyshev cheb = Chebyshev::fit(vmin_c, vmax_c, cheb_deg, vc, log_g_c);

    //tails anchored to the central fit's value and slope at the seams
    TailFit left_tail = TailFit::fit(vmin, vmin_c, vl, log_g_l, cheb);
    TailFit right_tail = TailFit::fit(vmax_c, vmax, vr, log_g_r, cheb);

    //stitch the tail fits with the central fit

    //v_trans is a velocity band width (m/s), not a sample count, so convert it
    //to index ranges via the velocity grid.
    const int n = static_cast<int>(v.size());
    const int lt_end = static_cast<int>(std::ranges::lower_bound(v, vmin_c + v_trans) - v.begin());
    const int rt_start = std::max(lt_end, static_cast<int>(std::ranges::lower_bound(v, vmax_c - v_trans) - v.begin()));

    std::vector<double> fitted(v.size());

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

    //undo log, then perform the final fit on the dist
    std::ranges::transform(fitted, fitted.begin(), [](double x){ return std::exp(x); });

    //final fit. the DAT consumers (ion_transport) evaluate the series at
    //x = v/(4*vtherm) and drop |x| > 1, so the coefficients must be Legendre
    //in exactly that variable: domain +-4*vtherm, fit only on points inside it
    const double leg_hw = 4.0*vtherm;

    //the projection is identically zero beyond its support and dist1d writes
    //the full support, so when the data stops short of +-4*vtherm extend it
    //with explicit zeros: the fit domain must be covered, and zero is the true
    //value out there. pad counts: smallest p with v.front() - p*dv <= -leg_hw
    //(mirrored on the right)
    const double dv_grid = v[1] - v[0];
    const int pad_l = std::max(0, static_cast<int>(std::ceil((v.front() + leg_hw)/dv_grid)));
    const int pad_r = std::max(0, static_cast<int>(std::ceil((leg_hw - v.back())/dv_grid)));

    if (pad_l + pad_r > 0)
    {
        std::print(std::cerr, "note: data [{:.1f}, {:.1f}] stops short of the fit domain +-{:.1f} (4*vtherm); zero-padded {}+{} bins\n",
                vmin, vmax, leg_hw, pad_l, pad_r);
    }

    std::vector<double> v_ext(n + pad_l + pad_r);
    std::vector<double> f_ext(n + pad_l + pad_r, 0.0); //zeros; data overwrites the middle

    for (int i = 0; i < pad_l; i++)
    {
        v_ext[i] = v.front() - (pad_l - i)*dv_grid;
    }

    std::ranges::copy(v, v_ext.begin() + pad_l);
    std::ranges::copy(fitted, f_ext.begin() + pad_l);

    for (int i = 0; i < pad_r; i++)
    {
        v_ext[pad_l + n + i] = v.back() + (i + 1)*dv_grid;
    }

    const std::size_t leg_off = static_cast<std::size_t>(std::ranges::lower_bound(v_ext, -leg_hw) - v_ext.begin());
    const std::size_t leg_cnt = static_cast<std::size_t>(std::ranges::upper_bound(v_ext, leg_hw) - v_ext.begin()) - leg_off;

    Legendre leg = Legendre::fit(-leg_hw, leg_hw, leg_deg,
            std::span<const double>{v_ext}.subspan(leg_off, leg_cnt),
            std::span<const double>{f_ext}.subspan(leg_off, leg_cnt));
    std::vector<double> coeffs = leg.get_coeffs_stdvec();

    std::print("{:.8f}\n", vtherm);
    for (int i = 0; i < coeffs.size(); i++)
    {
        if (i % 2 == 0)
        {
            std::print("{:.8e}\n", coeffs[i]);
        }
    }

    return 0;
}