#include <algorithm>
#include <iomanip>
#include <iostream>
#include <print>
#include <exception>
#include <vector>
#include <string>
#include <array>
#include <cmath>
#include <cstdint>
#include <atomic>
#include <fstream>
#include <chrono>
#include <random> 

#include "Simulation.h"
#include "Table.h"
#include "SimConfig.h"
#include "constants.h"
#include "npy.hpp"

//idk where to put these, util.h

template <int N, int M>
void write_array2d_npy(const std::string& filename, const std::array<std::array<double, M>, N>& arr)
{
    npy::npy_data_ptr<double> d;

    std::array<double, N*M> flat;
    for (int i = 0; i < N; i++)
    {
        std::copy(arr[i].begin(), arr[i].end(), flat.begin() + i*M);
    }

    d.data_ptr = flat.data();
    d.shape = {N, M};
    d.fortran_order = false;

    npy::write_npy(filename, d);
}

void write_dist_to_file(const std::string& path, const Simulation::Histogram& data, const SimConfig& cfg, double cnorm)
{
    std::ofstream file(path);

    //write header
    file << (cn::nparmx + 1) << ',' << (cn::npermx + 1) << ',' << cfg.get_dlvpar();
    file << ',' << cfg.get_dlvper() << ',' << cnorm << '\n';

    file << std::scientific << std::setprecision(6);

    //write dist
    for (int i = 0; i < cn::nparmx + 1; i++)
    {
        for (int j = 0; j < cn::npermx + 1; j++)
        {
            file << data[i][j];
            if (j != cn::npermx)
            {
                file << ',';
            }
        }

        if (i != cn::nparmx)
        {
            file << '\n';
        }
    }
}


// arg is config file filename
int main(int argc, char* argv[])
{
    constexpr uint64_t seed = 123'456'789'123ULL;
    //std::random_device rd;
    //const uint64_t seed = rd();

    if (argc != 2)
    {
        //I figured I would just stick with std::print since the next output uses it
        std::print(std::cerr, "Invalid arg count. the argument is the filename of the config\n");
        return 1;
    }

    auto t0 = std::chrono::steady_clock::now();

    SimConfig cfg;
    try
    {
        cfg.init(argv[1]);
    }
    catch (const std::exception& e)
    {
        std::print(std::cerr, "Failed to read config file/invalid config: {}\n{}\n", argv[1], e.what());
        return 1;
    }

    cfg.print_config_summary();

    const Table table{cfg};
    const int nthreads = cfg.get_nthreads();

    std::vector<Simulation> sims;
    sims.reserve(nthreads);

    std::mt19937_64 seeder{seed};
    for (int i = 0; i < nthreads; i++)
    {
        sims.push_back({cfg, table, seeder()});
    }

    auto t1 = std::chrono::steady_clock::now();
    uint64_t ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    uint64_t total_ms = ms;
    std::print("Initalized in {:.3f}[s]\n", ms/1000.0);

    t0 = std::chrono::steady_clock::now();

    std::atomic<uint64_t> collisions{0};

    #pragma omp parallel for num_threads(nthreads)
    for (int i = 0; i < nthreads; i++)
    {
        sims[i].run_sim(collisions, i);
    }

    // accumulate every thread's binned distribution into the master grid.
    //sum_hist takes z_master as a ref and writes to it
    Simulation::Histogram zm{}; //z_master
    for (const Simulation& s : sims)
    {
        s.sum_hist(zm);
    }

    double sum = 0.0;
    for (int ipar = 0; ipar <= cn::nparmx; ipar++)
    {
        for (int iper = 0; iper <= cn::npermx; iper++)
        {
            sum += zm[ipar][iper]*(iper + 0.5)*cfg.get_dlvper();
        }
    }

    const double cntnorm = 2.0*sum*cn::pi*cfg.get_dlvpar()*cfg.get_dlvper();

    constexpr double az_floor = -99.99;
    Simulation::Histogram az{};
    for (int ipar = 0; ipar <= cn::nparmx; ipar++)
    {
        for (int iper = 0; iper <= cn::npermx; iper++)
        {
            const double val = zm[ipar][iper];
            az[ipar][iper] = (val > 0.0) ? std::log(val/((iper + 0.5)*cntnorm)) : az_floor;
        }
    }

    t1 = std::chrono::steady_clock::now();
    ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

    std::print("Simulation and post-processing completed in {:.3f}[s]\n", ms/1000.0);
    std::print("Total time(excluding disk io for output): {:.3f}[s]\n", (total_ms + ms)/1000.0);

    const std::string& out = cfg.get_out_filename();
    write_array2d_npy<cn::nparmx + 1, cn::npermx + 1>(out + "av.npy", az);
    write_array2d_npy<cn::nparmx + 1, cn::npermx + 1>(out + "hist.npy", zm);
    write_dist_to_file(out + "_log_dist.csv", az, cfg, cntnorm);
    write_dist_to_file(out + "_lin_dist.csv", zm, cfg, cntnorm);

    std::print("\nResults have been written to disk.\n");

    return 0;
}