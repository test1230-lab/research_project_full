#include <cstddef>
#include <iostream>
#include <vector>
#include <fstream>
#include <numbers>
#include <string>

#include <cmath>
#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <iomanip>

#include <boost/math/quadrature/gauss.hpp>
#include "./include/libInterpolate/Interpolate.hpp"
#include "./include/mdarray.h"
#include "./include/npy.hpp"

constexpr double pi = std::numbers::pi_v<double>;
constexpr double az_floor = -99.99;

constexpr double dv = 10.0;

constexpr int nphi_panels = 24;

bool is_equal_dbl(double a, double b, double eps = 1e-9)
{
    return std::abs(a - b) < eps;
}

double sqr(double x)
{
    return x*x;
}

struct Dist
{
    int npar, nper;
    double dlvpar, dlvper, cnorm;
    array2d<double> dist;
};

class InterpWrapper
{
public:
    InterpWrapper() = default;

    //(x0, dx) and (y0, dy): origin and spacing of each axis, so operator()
    //can turn a physical point into a cell index without a search.
    void set_data(const std::vector<double>& x, const std::vector<double>& y,
                  const std::vector<double>& z, const array2d<double>& dist,
                  double x0, double dx, double y0, double dy)
    {
        this->x0 = x0;
        this->dx = dx;
        this->y0 = y0;
        this->dy = dy;

        ni = static_cast<int>(dist.dim0());
        nj = static_cast<int>(dist.dim1());

        build_cubic_mask(dist);

        lin_interp.setData(x, y, z);
        cubic_interp.setData(x, y, z);
    }

    double operator()(double x, double y) const
    {
        const int i = cell_index(x, x0, dx, ni);
        const int j = cell_index(y, y0, dy, nj);

        return static_cast<bool>(cubic_mask[i, j]) ? cubic_interp(x, y) : lin_interp(x, y);
    }

private:
    _2D::BilinearInterpolator<double> lin_interp{};
    _2D::BicubicInterpolator<double> cubic_interp{};

    double x0{}, dx{1.0}, y0{}, dy{1.0};
    int ni{}, nj{};

    //per-cell, so dims are (ni-1) x (nj-1). char, not bool, to avoid
    //the std::vector<bool> specialization
    array2d<char> cubic_mask{};

    //physical coord -> lower-left node of its cell, clamped to a valid cell.
    static int cell_index(double v, double v0, double dv, int n)
    {
        const int idx = static_cast<int>(std::floor((v - v0)/dv));
        return std::clamp(idx, 0, n - 2);
    }

    //A cell is cubic-eligible only if its full 4x4 stencil ([i-1, i+2] x
    //[j-1, j+2]) is in bounds and floor-free; border cells fall back to linear.
    void build_cubic_mask(const array2d<double>& dist)
    {
        cubic_mask.resize(ni - 1, nj - 1);

        auto is_floor = [&](int i, int j)
        {
            return is_equal_dbl(dist[i, j], az_floor);
        };

        for (int i = 0; i < ni - 1; i++)
        {
            for (int j = 0; j < nj - 1; j++)
            {
                bool ok = (i - 1 >= 0) && (i + 2 < ni) &&
                          (j - 1 >= 0) && (j + 2 < nj);

                for (int ii = i - 1; ok && ii <= i + 2; ii++)
                {
                    for (int jj = j - 1; ok && jj <= j + 2; jj++)
                    {
                        if (is_floor(ii, jj))
                        {
                            ok = false;
                        }
                    }
                }

                cubic_mask[i, j] = ok ? 1u : 0u;
            }
        }
    }
};

//left col: v, right col: g
void write_1d_dist_npy(const std::string& filename, const std::vector<double>& v, const std::vector<double>& g)
{
    npy::npy_data_ptr<double> d;
    constexpr int n_cols = 2;
    const std::size_t N = v.size();

    std::vector<double> flat(n_cols*N);
    for (std::size_t i = 0; i < N; i++)
    {
        flat[i*n_cols] = v[i];
        flat[i*n_cols + 1] = g[i];
    }

    d.data_ptr = flat.data();
    d.shape = {N, n_cols};
    d.fortran_order = false;

    npy::write_npy(filename, d);
}

void write_1d_dist_csv(const std::string& path, const std::vector<double>& v, const std::vector<double>& g)
{
    std::ofstream file(path);

    file << std::scientific << std::setprecision(6);

    for (std::size_t i = 0; i < v.size(); i++)
    {
        file << v[i] << ',' << g[i];

        if (i < v.size() - 1)
        {
            file << '\n';
        }
    }
}

Dist read_log_dist_csv(const std::string& path)
{
    std::ifstream file(path);
    if (!file)
    {
        throw std::runtime_error("Could not open file");
    }

    Dist d;

    std::string line;
    std::getline(file, line);
    std::ranges::replace(line, ',', ' ');

    std::istringstream header(line);
    header >> d.npar >> d.nper >> d.dlvpar >> d.dlvper >> d.cnorm;

    //convert to SI
    d.dlvpar /= 100.0;
    d.dlvper /= 100.0;
    d.cnorm /= 100.0*100.0*100.0;

    d.dist.resize(d.npar, d.nper);

    //values are ln(f) in CGS; f -> SI shifts them by ln(100^3). convert here so
    //the returned Dist is fully SI (floor bins stay at the floor)
    const double log_conv = std::log(std::pow(100.0, 3));

    for (int i = 0; i < d.npar; i++)
    {
        std::getline(file, line);
        std::ranges::replace(line, ',', ' ');
        std::istringstream row(line);

        for (int j = 0; j < d.nper; j++)
        {
            row >> d.dist[i, j];

            if (d.dist[i, j] > az_floor)
            {
                d.dist[i, j] += log_conv;
            }
        }
    }

    return d;
}

using boost::math::quadrature::gauss;

template <class F>
double composite_gauss(F&& f, double a, double b, int panels)
{
    const double h = (b - a)/panels;
    double s = 0.0;

    for (int i = 0; i < panels; i++)
    {
        s += gauss<double, 15>::integrate(f, a + i*h, a + (i + 1)*h);
    }

    return s;
}

//creates 1d dist
template <class I>
double compute_1d_dist(const I& interp, double rmax, double alpha, double v1,
        double par_min, double par_max, double per_min, double per_max)
{
    constexpr double eps = 1e-12;

    const double ca = std::cos(alpha);
    const double sa = std::sin(alpha);

    //this is ugly
    auto radial = [&](double phi) -> double
    {
        const double cphi = std::cos(phi);
        const double sphi = std::sin(phi);

        //vpar(r) = v1*ca - A*r ; vper(r)^2 = (v1*sa + B*r)^2 + (C*r)^2
        const double A = sa*cphi;
        const double B = ca*cphi;
        const double C = sphi;

        //|vpar| <= par_max  ->  r in [ra_lo, ra_hi]
        double ra_lo = 0.0, ra_hi = rmax;
        if (std::abs(A) > eps)
        {
            const double r1 = (v1*ca - par_max)/A;
            const double r2 = (v1*ca + par_max)/A;
            ra_lo = std::min(r1, r2);
            ra_hi = std::max(r1, r2);
        }
        else if (std::abs(v1*ca) > par_max)
        {
            return 0.0;
        }

        //vper <= per_max  ->  a2*r^2 + b*r + c <= 0
        const double a2 = B*B + C*C;
        const double b = 2.0*v1*sa*B;
        const double c = sqr(v1*sa) - sqr(per_max);

        double rb_lo = 0.0, rb_hi = rmax;
        if (a2 > eps)
        {
            const double disc = b*b - 4.0*a2*c;
            if (disc < 0.0)
            {
                return 0.0;
            }
            const double sd = std::sqrt(disc);
            rb_lo = -(b + sd)/(2.0*a2);
            rb_hi = (sd - b)/(2.0*a2);
        }
        else if (c > 0.0)
        {
            return 0.0;
        }

        const double r_lo = std::max({0.0, ra_lo, rb_lo});
        const double r_hi = std::min({rmax, ra_hi, rb_hi});
        if (r_hi <= r_lo)
        {
            return 0.0;
        }

        auto fr = [&](double r) -> double
        {
            const double vpar = v1*ca - A*r;
            const double vper = std::sqrt(sqr(v1*sa + B*r) + sqr(C*r));

            const double vpar_c = std::clamp(std::abs(vpar), par_min, par_max);
            const double vper_c = std::clamp(vper,per_min, per_max);

            const double log_val = interp(vpar_c, vper_c);
            if (log_val <= az_floor)
            {
                return 0.0;
            }

            return r*std::exp(log_val);
        };

        return gauss<double, 20>::integrate(fr, r_lo, r_hi);
    };

    return composite_gauss(radial, 0.0, 2.0*pi, nphi_panels);
}


//args are input.csv output(dont include .csv) aspect_angle(deg)
int main(int argc, char* argv[])
{
    if (argc != 4)
    {
        std::cerr << "Invalid number of args. args are: input.csv output(dont include .csv) aspect_angle(deg)\n";
        return 1;
    }

    const std::string input_filename{argv[1]};
    const std::string output_filename{argv[2]};
    const double alpha = (pi/180.0)*std::stod(argv[3]);//convert to radians

    if(output_filename.ends_with(".csv"))
    {
        std::cerr << "output filename has an extension\n";
        return 1;
    }

    const Dist d = read_log_dist_csv(input_filename);

    std::vector<double> parv(d.npar);
    std::vector<double> perv(d.nper);

    for (int i = 0; i < d.npar; i++)
    {
        parv[i] = (i + 0.5)*d.dlvpar;
    }

    for (int i = 0; i < d.nper; i++)
    {
        perv[i] = (i + 0.5)*d.dlvper;
    }

    const double rmax = std::sqrt(parv.back()*parv.back() + perv.back()*perv.back());

    std::vector<double> xgrid;
    std::vector<double> ygrid;
    std::vector<double> zgrid;

    xgrid.reserve(d.dist.size());
    ygrid.reserve(d.dist.size());
    zgrid.reserve(d.dist.size());

    for (int i = 0; i < d.npar; i++)
    {
        for (int j = 0; j < d.nper; j++)
        {
            xgrid.push_back(parv[i]);
            ygrid.push_back(perv[j]);
            zgrid.push_back(d.dist[i, j]);
        }
    }

    InterpWrapper interp;
    interp.set_data(xgrid, ygrid, zgrid, d.dist, parv.front(), d.dlvpar, perv.front(), d.dlvper);

    const double par_min = parv.front();
    const double par_max = parv.back();
    const double per_min = perv.front();
    const double per_max = perv.back();

    //output grid: the projection's support is |v1| <= rmax, so size the grid
    //from the data instead of a hardcoded extent. this keeps the downstream
    //Legendre fit window (+-4*vtherm) inside the data at every E-field
    const double vmax = std::ceil(rmax/dv)*dv;
    const double vmin = -vmax;
    const int nv = static_cast<int>(std::round((vmax - vmin)/dv)) + 1;

    std::vector<double> v1(nv);
    std::vector<double> fv(nv);

    #pragma omp parallel for
    for (int i = 0; i < nv; i++)
    {
        const double v = vmin + dv*i;
        fv[i] = compute_1d_dist(interp, rmax, alpha, v, par_min, par_max, per_min, per_max);
        v1[i] = v;
    }

    std::cout << "\nwriting output to disk\n";
    write_1d_dist_csv(output_filename + ".csv", v1, fv);
    write_1d_dist_npy(output_filename + ".npy", v1, fv);
    std::cout << "Wrote files to disk.\n";

    return 0;
}