#include "./include/distrib2d.h"

double sqr(double x)
{
	return x*x;
}

// mV/m
double Distrib2D::electric_field(double t) const
{
	if (t <= 0)
	{
		return 0.0;
	}

    const double ramp = t*e_field_ramp_rate;
    if (ramp < e_field_peak)// rising edge
    {
    	return ramp;
    }

    return std::max(0.0, e_field_peak - (ramp - e_field_peak));   // falling edge, min of 0
}

//maybe precompute and interpolate
double Distrib2D::distribution_function(double v, double n, double t, int angle) const
{
	const double efield = electric_field(t);
	const double integral = ef.compute_integral(efield, angle);
	return n*ef.compute_dist(efield, angle, v)/integral;
}

double Distrib2D::boundary_density(double t) const
{
	return 1e12;
	
	/*const double a = 1e12;

	if (t < 0.0)
	{
		return a;
	}

	return a + (t/1000.0)*a;*/

	//return 1e12;

	/*if (t < 0.0)
	{
		return 1e12;
	}
	return std::clamp(1e12*std::pow(0.5, t/200.0), 0.5e12, 1e12);*/
}

double Distrib2D::calc_dt(double vi, double vf, double t) const
{
	return (vf - vi)/acceleration;
}

//should throw if sqrt arg is less than zero
double Distrib2D::compute_vi_par(double vf_par, double t, double dz) const
{
	return std::sqrt(vf_par*vf_par - 2.0*acceleration*dz);
}

//optimize!
array2d<double> Distrib2D::get_f_vf_dist(double t, double dz) const
{
	array2d<double> result(nv, nv);

	#pragma omp parallel for schedule(dynamic, 8)
	for (int i = 0; i < nv; i++)
	{
		const double vf_par = vf_vec[i];
		const double vi_par = compute_vi_par(vf_par, t, dz);
		const double dt = calc_dt(vi_par, vf_par, t); //same particle
		const double bnd_density = boundary_density(t - dt);
		const double f_vf_par = distribution_function(vi_par, bnd_density, t - dt, 0);

		for (int j = 0; j < nv; j++)
		{
			const double vf_perp = vf_vec[j];
			const double f_vf_perp = distribution_function(vf_perp, bnd_density, t - dt, 90);
			result[i, j] = f_vf_par*f_vf_perp/bnd_density;
		}
	}

	return result;
}

//utility function
std::vector<double> Distrib2D::calc_times(double t1, double t2, double dt) const
{
	const int N = std::round((t2 - t1) / dt);
	std::vector<double> times(N);

	for (int i = 0; i < N; i++)
	{
		times[i] = t1 + i*dt;
	}

	return times;
}

Distrib2D::Moments Distrib2D::get_moments(double t1, double t2, double dt, double dz) const
{
	const int N = std::round((t2 - t1) / dt);
	Moments m(N);

	#pragma omp parallel for schedule(dynamic, 8)
	for (int n = 0; n < N; n++)
	{
		const double t = t1 + dt*n;

		double m0_par = 0.0;
		double m1_par = 0.0;
		double m2_par = 0.0;

		double m0_perp = 0.0;
		double m2_perp = 0.0;

		double integral0 = 0.0;
		double integral1 = 0.0;

    	const double tau_perp = t - calc_dt(compute_vi_par(0.0, t, dz), 0.0, t);
    	const double bnd_density_perp = boundary_density(tau_perp);

		for (int i = 0; i < nv; i++)
		{
			const double vf_par = vf_vec[i];
			const double vi_par = compute_vi_par(vf_par, t, dz);
			const double tau = t - calc_dt(vi_par, vf_par, t);
			const double bnd_density = boundary_density(tau);
			const double f_par = distribution_function(vi_par, bnd_density, tau, 0);
			const double f_perp = distribution_function(vf_vec[i], bnd_density_perp, tau_perp, 90); //vi = vf for perp
			const double wi = 1.0 - 0.5*((i == 0) + (i == nv - 1));

			const double v = vf_par;

			m0_par += wi*f_par;
			m1_par += wi*f_par*v;
			m2_par += wi*f_par*v*v;

			m0_perp += wi*f_perp;
			m2_perp += wi*f_perp*v*v;

			const double s = bnd_density / dv;

			integral0 += wi*s*f_par/bnd_density;
			integral1 += wi*s*f_par*vf_par/bnd_density;
		}

		const double density = dv*dv*integral0;
		const double avg_vel_par = integral1/integral0; //dv*dv*integral1 / density simplified

		m.density[n] = density;
        m.avg_par_velocity[n] = avg_vel_par;
		
		const double a = mass/kb;
		m.ion_temp_par[n] = a*(m2_par/m0_par - sqr(m1_par/m0_par));
		m.ion_temp_perp[n] = a*m2_perp/m0_perp;
	}

	return m;
}
