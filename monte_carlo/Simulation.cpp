#include <array>
#include <cmath>

#include "Simulation.h"
#include "Table.h"
#include "SimConfig.h"
#include "util.h"
#include "constants.h"

void Simulation::sum_hist(Histogram& other) const
{
    for (int i = 0; i < cn::nparmx + 1; i++)
    {
        for (int j = 0; j < cn::npermx + 1; j++)
        {
            other[i][j] += hist[i][j];
        }
    }
}

std::array<double, 3> Simulation::collis(const std::array<double, 3>& vi, const std::array<double, 3>& vn, 
                                        const std::array<double, 3>& g, double gsq, double mn, double mi, double chi)
{
    const double coschi = std::cos(chi);
    const double sinchi = std::sin(chi);

    const double normg = std::sqrt(gsq);
    const double mmi = mi/(mi + mn);
    const double mmn = mn/(mi + mn);

    //origin of center of mass coordinates
    const double ocmx = mmi*vi[0] + mmn*vn[0];
    const double ocmy = mmi*vi[1] + mmn*vn[1];
    const double ocmz = mmi*vi[2] + mmn*vn[2];

    //radius of the base of the scattering cone (norm of vector A)
    const double norma = std::abs(normg*sinchi);
    const double xlen = 2.0*norma;
    const double ylen = 2.0*norma;

    //gbar, i.e. vector
    const double gper = std::sqrt(g[1]*g[1] + g[2]*g[2]);
    if (gper == 0.0 || normg == 0.0)
    {
        return vi;
    }

    double x = 0.0;
    double y = 0.0;
    double r = 0.0;
    do
    {
        x = xlen*(rng.uniform() - 0.5);
        y = ylen*(rng.uniform() - 0.5);
        r = std::sqrt(x*x + y*y);
    } while (r > norma || r == 0.0);

    //renormalize x and y onto the circle
    x *= norma/r;
    y *= norma/r;

    const double aax = g[0]*coschi + x*gper/normg;
    const double aay = g[1]*coschi - x*g[1]*g[0]/(normg*gper) + y*g[2]/gper;
    const double aaz = g[2]*coschi - x*g[0]*g[2]/(normg*gper) - y*g[1]/gper;

    return {ocmx + mmn*aax, ocmy + mmn*aay, ocmz + mmn*aaz};
}


void Simulation::run_sim(std::atomic<uint64_t>& coll, int tid)
{
    const auto& d = cfg.get_derived();
    const auto& neutrals = cfg.get_neutrals();
    const std::size_t n_neutrals = neutrals.size();
    const uint64_t coll_tot = cfg.get_total_collisions();
    const int nthreads = cfg.get_nthreads();

    const uint64_t base = coll_tot / nthreads;
    const uint64_t rem = coll_tot % nthreads;
    const uint64_t n_collisions = base + (static_cast<uint64_t>(tid) < rem ? 1 : 0);

    const double alfa = tbl.get_alfa();
    const double kmax = tbl.get_kmax();

    std::array<double, 3> vaft{};
    std::array<double, 3> vbef{};
    std::array<double, 3> vneu{};
    std::array<double, 3> g{};

    double sf = 0.0;
    double cf = 1.0;

    double vper = -(d.vd + 1.0)*(rng.uniform() + 0.5);
    vaft[2] = rng.uniform()*100000.0;

    uint64_t reported = 0;

    //main loop
    for (uint64_t i = 0; i < n_collisions; i++)
    {
        const uint64_t done = i + 1;
        if (done % cn::n_report == 0 || done == n_collisions)
        {
            const uint64_t coll_now = (coll += done - reported);
            reported = done;

            #pragma omp critical(print_progress)
            {
                const double percent = 100.0 * coll_now/static_cast<double>(coll_tot);
                std::print(stdout, "\rprogress: {:.1f}M / {:.1f}M ({:.1f}%)     ", coll_now/1e6, coll_tot/1e6, percent);
                std::fflush(stdout);
            }
        }

        double delth = 0.0;
        int ineut = -1;
        double b2max = 0.0;
        double b2tot = 0.0;
        double b2 = 0.0;
        double b2ex = 0.0;
        double b2pol = 0.0;
        double gt2 = 0.0;
        double ta = 0.0;
        double log_ta = 0.0;
        double bcon = 0.0;

        while(true)
        {
            const double rnum = 1.0 - rng.uniform(); //now its (0,1]
            const double rdis = -std::log(rnum)*alfa;

            delth += rdis;

            //pick neutral
            double fabun = 0.0;
            const double pikneu = rng.uniform()*d.tabun;
            for (std::size_t j = 0; j < n_neutrals; j++)
            {
                ineut = j;
                fabun += neutrals[j].rel_abundance;
                if (fabun >= pikneu)
                {
                    break;
                }
            }

            const NeutralInfo& nu = neutrals[ineut];

            const double sd = std::sin(delth);
            const double cd = std::cos(delth);
            vbef[0] = vper*(sf*cd + cf*sd);
            vbef[1] = vper*(cf*cd - sf*sd);
            vbef[2] = vaft[2];

            vneu[0] = nu.beta*rng.normal_dist();
            vneu[1] = nu.beta*rng.normal_dist() - d.vd;
            vneu[2] = nu.beta*rng.normal_dist();

            g[0] = vbef[0] - vneu[0];
            g[1] = vbef[1] - vneu[1];
            g[2] = vbef[2] - vneu[2];

            gt2 = g[0]*g[0] + g[1]*g[1] + g[2]*g[2];
            const double gt = std::sqrt(gt2);
            const double energy = 0.5*nu.mu*gt2;
            ta = energy/nu.eps;
            log_ta = std::log10(ta);

            bcon = tbl.bcutoff_lookup(ta, log_ta, ineut, nu.chi_cutoff);
            const double bcut = bcon*nu.rm;
            b2pol = bcut*bcut;
            const double log_gfactor = log_ta + nu.log_eps_over_bk;
            b2ex = 2.0/cn::pi*sqr(nu.ar - nu.br*log_gfactor);
            b2tot = std::max(b2ex, b2pol);
            b2max = kmax/gt/cn::pi;
            b2 = rng.uniform()*b2max;

            if (b2 <= b2tot)
            {
                break;
            }
        }

        const int nv3 = static_cast<int>(std::abs(vbef[2])*d.scvpar);
        const int nvper = static_cast<int>(std::abs(vper)*d.scvper);
        if (nv3 <= cn::nparmx && nvper <= cn::npermx)
        {
            hist[nv3][nvper] += delth;
        }

        if (b2 < b2ex)
        {
            const double rand = rng.uniform();
            if (rand < 0.5)
            {
                std::swap(vbef, vneu);
                for (double& gi : g)
                {
                    gi = -gi;
                }
            }
        }

        if (b2 < b2pol)
        {
            const double b = std::sqrt(b2);
            const double u = (bcon > 0.0) ? b/(bcon*neutrals[ineut].rm) : 0.0;
            const double chipol = tbl.chi_lookup_u(log_ta, ineut, u);

            if (chipol != 0.0)
            {
                vaft = collis(vbef, vneu, g, gt2, neutrals[ineut].mass, d.migrm, chipol);
            }
            else
            {
                vaft = vbef;
            }
        }
        else
        {
            vaft = vbef;
        }

        vper = std::sqrt(vaft[0]*vaft[0] + vaft[1]*vaft[1]);
        if (vper != 0.0)
        {
            const double inv = 1.0/vper;
            sf = vaft[0]*inv;
            cf = vaft[1]*inv;
        }
        else
        {
            sf = 0.0;
            cf = 1.0;
        }
    }
}