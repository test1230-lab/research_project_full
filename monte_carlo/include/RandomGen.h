#pragma once

#include <array>
#include <cstdint>
#include <limits>
#include <boost/random/normal_distribution.hpp>

//xoshiro256++ Sebastiano Vigna
class RandomGen
{
public:
    using result_type = uint64_t;

    RandomGen(uint64_t seed)
    {
        x = seed;

        s[0] = splitmix_next();
        s[1] = splitmix_next();
        s[2] = splitmix_next();
        s[3] = splitmix_next();
    }

    uint64_t operator()()
    {
        return next();
    }

    //returns [0,1)
    double uniform()
    {
        return static_cast<double>(next() >> 11) * 0x1.0p-53;
    }

    //std dev = 1 mean = 0
    double normal_dist()
    {
        return normal(*this);
    }

    static constexpr uint64_t min()
    { 
        return std::numeric_limits<result_type>::min();
    }

    static constexpr uint64_t max()
    {
        return std::numeric_limits<result_type>::max();
    }

private:
    std::array<uint64_t, 4> s;
    uint64_t x;
    boost::random::normal_distribution<double> normal{0.0, 1.0};

    uint64_t rotl(uint64_t x, int k) const
    {
        return (x << k) | (x >> (64 - k));
    }
    
    uint64_t splitmix_next()
    {
        uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;   
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }

    uint64_t next()
    {
        const uint64_t result = rotl(s[0] + s[3], 23) + s[0];

        const uint64_t t = s[1] << 17;

        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];

        s[2] ^= t;

        s[3] = rotl(s[3], 45);

        return result;
    }
};