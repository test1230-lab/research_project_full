#pragma once

#include <array>
#include <cstdint>
#include <atomic>

#include "constants.h"
#include "SimConfig.h"
#include "RandomGen.h"
#include "Table.h"
#include "util.h"

//one object per thread
class Simulation
{
public:
    using Histogram = std::array<std::array<double, cn::npermx + 1>, cn::nparmx + 1>;

    Simulation(const SimConfig& config, const Table& table, uint64_t seed)
        : cfg{config}, tbl{table}, rng{seed} {}

    void sum_hist(Histogram& other) const;
    void run_sim(std::atomic<uint64_t>& collisions, int tid);

private:
    const SimConfig& cfg;
    const Table& tbl;
    Histogram hist{};
    RandomGen rng;

    vec3 collis(const vec3& vi, const vec3& vn, const vec3& g, double gsq, double mn, double mi, double chi);
};
