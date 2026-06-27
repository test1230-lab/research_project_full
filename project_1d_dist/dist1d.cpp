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

constexpr double vmin = -7500.0;
constexpr double vmax = 7500.0;
constexpr double dv = 10.0;
constexpr int nv = (vmax - vmin)/dv + 1;

constexpr int nphi_panels = 24;

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

// left col: v, right col: g
void write_1d_dist_npy(const std::string& filename, const std::vector<double>& v, const std::vector<double>& g)
{
    npy::npy_data_ptr<double> d;
	constexpr int n_cols = 2;
	const std::size_t N = v.size();

    std::vector<double> flat(n_cols*N);
    for (std::size_t i = 0; i < N; i++)
    {
        flat[i*n_cols + 0] = v[i];
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
    std::replace(line.begin(), line.end(), ',', ' ');

	std::istringstream header(line);
	header >> d.npar >> d.nper >> d.dlvpar >> d.dlvper >> d.cnorm;

	// Convert to SI
	d.dlvpar /= 100.0;
	d.dlvper /= 100.0;
	d.cnorm /= 100.0*100.0*100.0;

	d.dist.resize(d.npar, d.nper);

	for (int i = 0; i < d.npar; i++)
	{
		std::getline(file, line);
		std::replace(line.begin(), line.end(), ',', ' ');
		std::istringstream row(line);

		for (int j = 0; j < d.nper; j++)
		{
			row >> d.dist[i, j];
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

		// vpar(r) = v1*ca - A*r ; vper(r)^2 = (v1*sa + B*r)^2 + (C*r)^2
		const double A = sa*cphi;
		const double B = ca*cphi;
		const double C = sphi;

		// |vpar| <= par_max  ->  r in [ra_lo, ra_hi]
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

		// vper <= per_max  ->  a2*r^2 + b*r + c <= 0
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
			rb_lo = (-b - sd)/(2.0*a2);
			rb_hi = (-b + sd)/(2.0*a2);
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

	const double log_conv = std::log(std::pow(100.0, 3));

	const std::string input_filename{argv[1]};
	const std::string output_filename{argv[2]};
	const double alpha = (pi/180.0)*std::stod(argv[3]);//convert to radians

	if(output_filename.ends_with(".csv"))
	{
		std::cerr << "output filename has an extension\n";
		return 1;
	}

	Dist d = read_log_dist_csv(input_filename);

	for (int i = 0; i < d.npar; i++)
	{
		for (int j = 0; j < d.nper; j++)
		{
			const double x = d.dist[i, j];
			d.dist[i, j] = (x > az_floor) ? (x + log_conv) : az_floor;
		}
	}

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

	//TODO: I want to use bicubic for the center, bilinear for the edges
	_2D::BilinearInterpolator<double> interp;
	interp.setData(xgrid, ygrid, zgrid);

	const double par_min = parv.front();
	const double par_max = parv.back();
	const double per_min = perv.front();
	const double per_max = perv.back();

	std::vector<double> v1(nv);
	std::vector<double> fv1(nv);

	#pragma omp parallel for
	for (int i = 0; i < nv; i++)
	{
		/*if (i % 100 == 0)
		{
			const double percent = std::ceil(100.0*i/static_cast<double>(nv));

			std::cout << "\rProgress: " << std::fixed << std::setprecision(1)
					  << percent << "%" << std::flush;
		}*/

		const double v = vmin + dv*i;
		fv1[i] = compute_1d_dist(interp, rmax, alpha, v, par_min, par_max, per_min, per_max);
		v1[i] = v;
	}

	std::cout << "\nwriting output to disk\n";
	write_1d_dist_csv(output_filename + ".csv", v1, fv1);
	write_1d_dist_npy(output_filename + ".npy", v1, fv1);
	std::cout << "Wrote files to disk.\n";

	return 0;
}