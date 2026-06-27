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
#include <numeric>

#include "Simulation.h"
#include "Table.h"
#include "SimConfig.h"
#include "constants.h"
#include "npy.hpp"

//idk where to put these

template <int N, int M>
void write_array2d_npy(const std::string& filename, const std::array<std::array<double, M>, N>& arr)
{
    npy::npy_data_ptr<double> d;

    std::array<double, N*M> flat;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            flat[j + i*M] = arr[i][j];
        }
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

    //write normalized 2d dist
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

    //initalize sim for each thread
    for (int i = 0; i < nthreads; i++)
    {
        //magic number is 64 bit golden ratio used by splitmix64
        const uint64_t tid = static_cast<uint64_t>(i);
        const uint64_t thread_seed = seed + 0x9E3779B97F4A7C15ULL*tid;
        sims.push_back({cfg, table, thread_seed});
    }

    auto t1 = std::chrono::steady_clock::now();
    int ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

    std::print("Initalized in {:.3f}s\n", ms/1000.0);

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

    //nested loop instead?
    //find normalization constant
    double sum = 0.0;
    for (int i = 0; i <= cn::npermx; i++)
    {
        const double w = (i + 0.5)*cfg.get_dlvper();
        sum += w*std::accumulate(zm[i].begin(), zm[i].end(), 0);
    }

    const double cntnorm = 2.0*sum*cn::pi*cfg.get_dlvpar()*cfg.get_dlvper();

    constexpr double az_floor = -99.99;
    Simulation::Histogram az{};
    for (int ipar = 0; ipar <= cn::nparmx; ipar++)
    {
        for (int iper = 0; iper <= cn::npermx; iper++)
        {
            if (zm[ipar][iper] > 0.0)
            {
                az[ipar][iper] = std::log(zm[ipar][iper]/((iper + 0.5)*cntnorm));
            }
            else
            {
                az[ipar][iper] = az_floor;
            }
        }
    }

    const std::string& out = cfg.get_out_filename();
    write_array2d_npy<cn::nparmx + 1, cn::npermx + 1>(out + "av.npy", az);
    write_array2d_npy<cn::nparmx + 1, cn::npermx + 1>(out + "hist.npy", zm);
    write_dist_to_file(out + "_log_dist.csv", az, cfg, cntnorm);

    std::print("\nResults have been written to disk.\n");

    return 0;
}
