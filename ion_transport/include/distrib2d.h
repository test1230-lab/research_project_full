#pragma once
#include <cmath>
#include <vector>
#include <numbers>

#include "mdarray.h"
#include "ElectricField.h"

class Distrib2D
{
public:
	Distrib2D(double v_min, double v_max, double v_step, double m, std::string coeff_path) : 
		vf_min(v_min), vf_max(v_max), dv(v_step), mass(m), ef(coeff_path)
	{
		nv = std::round((vf_max - vf_min) / dv) + 1;

		//populate vf arrays
		vf_vec.resize(nv);
		for (int i = 0; i < nv; i++)
		{
			vf_vec[i] = vf_min + dv*i;
		}
	}
	
	//moments at time steps
	struct Moments
	{
		Moments(int N) : density(N), avg_par_velocity(N),
						 ion_temp_par(N) , ion_temp_perp(N){}

		std::vector<double> density;
		std::vector<double> avg_par_velocity;
		std::vector<double> ion_temp_par;
		std::vector<double> ion_temp_perp;
	};

	array2d<double> get_f_vf_dist(double t, double dz) const;
	Moments get_moments(double t1, double t2, double dt, double dz) const;
	std::vector<double> calc_times(double t1, double t2, double dt) const;

private:
	double vf_min, vf_max, dv, mass;
	std::vector<double> vf_vec;
	int nv;
	ElectricField ef;

	static constexpr double pi = std::numbers::pi_v<double>;
	static constexpr double kb = 1.380649e-23; //https://physics.nist.gov/cgi-bin/cuu/Value?k
	static constexpr double acceleration = -5.0; //m s^-2
	static constexpr double e_field_ramp_rate = 1.0; //1 mV/m per second
	static constexpr double e_field_peak = 100.0; //  [mV/m]

	double electric_field(double t) const; // [mV/m]
	double boundary_density(double t) const; 

	double distribution_function(double v, double n, double t, int angle) const;

	double calc_dt(double vi, double vf, double t) const;
	double compute_vi_par(double vf_par, double t, double dz) const;
};