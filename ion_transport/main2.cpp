#include <iostream>
#include <vector>
#include <string>
#include <print>

#include "./include/mdarray.h"
#include "./include/distrib2d.h"
#include "./include/npy.hpp" // https://github.com/llohse/libnpy


void  write_vec1d_npy(const std::string& filename, const std::vector<double>& vec)
{
    npy::npy_data_ptr<double> d;
    d.data_ptr = vec.data();
    d.shape = {vec.size()};
    d.fortran_order = false;

    npy::write_npy(filename, d);
}

void write_array2d_npy(const std::string& filename, const array2d<double>& arr)
{
    npy::npy_data_ptr<double> d;
    d.data_ptr = arr.data();
    d.shape = {arr.dim0(), arr.dim1()};
    d.fortran_order = false;

    npy::write_npy(filename, d);
}

//arg is coeff dir
int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cerr << "invalid arg count\n";
        return 1; 
    }

    const std::string dir(argv[1]);

    if (dir == "-h" || dir == "--help")
    {
        std::cerr << "The argument is for the folder that contains the coeffs\n";
        return 0;
    }

    const std::string moment_dir = "./moment_output/";
    const std::string grid_dir = "./output_2d/";

	const double dz = 100e3;
	const double mass = 2.66e-26;
    
    const double t1 = 0.0;
	const double t2 = 1000.0;
	const double dt = 2.0;

    const double vmin = -5500.0;
    const double vmax = -vmin;
    const double dv = 2.0;

    Distrib2D dist{vmin, vmax, dv, mass, dir};
    
    std::vector<double> times = dist.calc_times(t1, t2, dt);
    Distrib2D::Moments m = dist.get_moments(t1, t2, dt, dz);


    Distrib2D dist2{-3000.0, 3000.0, 10, mass, dir};

    std::vector<array2d<double>> dists(10);
    std::vector<std::string> titles(10);
    
    for (int i = 0; i < 9; i++)
    {
        const double t = 100.0 + 50.0*i;
        titles[i] =  std::to_string((int)t);
        dists[i] = dist2.get_f_vf_dist(t, dz);
    }

    std::print("writing outputs to disk\n");

    write_vec1d_npy(moment_dir + "time.npy", times);
    write_vec1d_npy(moment_dir + "density.npy", m.density);
    write_vec1d_npy(moment_dir + "avg_par_velocity.npy", m.avg_par_velocity);
    write_vec1d_npy(moment_dir + "ion_temp_par.npy", m.ion_temp_par);
    write_vec1d_npy(moment_dir + "ion_temp_perp.npy", m.ion_temp_perp);

    for (int i = 0; i < 9; i++)
    {
        const std::string filename = grid_dir + titles[i] + ".npy";
        write_array2d_npy(filename, dists[i]);
    }

	return 0;
}