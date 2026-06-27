#pragma once

#include <cmath>
#include <stdexcept>

#include <boost/math/tools/roots.hpp>
#include <boost/math/quadrature/gauss.hpp>

#include "constants.h"
#include "util.h"

namespace scatter
{
    inline bool abs_tol(double a, double b)
    {
        return std::abs(b - a) <= cn::root_finding_tol;
    }

    //multiplying several times is faster than std::pow unless -ffast-math is used
    inline double scattering_integrand(double rcb, double tan_a, double theta)
    {
        const double s = std::sin(theta);
        const double s2 = s*s;
        const double s4 = s2*s2;
        const double s6 = s4*s2;
        const double s8 = s6*s2;
        const double s10 = s8*s2;

        const double theta1 = s2;
        const double theta5 = s2 + s4 + s6 + s8 + s10;

        const double rcb4 = rcb*rcb*rcb*rcb;
        const double rcb12 = rcb4*rcb4*rcb4;

        const double x1 = 1.0 + theta5/(2.0*rcb12*tan_a);
        const double x2 = -1.5*theta1/(rcb4*tan_a);
        double z = x1 + x2;

        // z is the radicand; it vanishes at the classical turning point. Because
        // rcb comes from a root solve with finite tolerance (cn::root_finding_tol),
        // a quadrature node can land just past the turning point and drive z to a
        // tiny negative value. That is numerical noise relative to the magnitudes
        // being subtracted, so clamp it to the singularity instead of crashing;
        // only a radicand that is negative by a non-trivial fraction of those
        // magnitudes signals a genuine domain error.
        if (z < 0.0)
        {
            const double scale = std::abs(x1) + std::abs(x2);
            if (z < -1e-5 * scale)
            {
                throw std::domain_error("in function scattering_integrand: negative sqrt");
            }
            z = std::abs(z);
        }

        return 1.0/std::sqrt(z);
    }

    inline double find_rcb_root(double tan_a, double b)
    {
        const double ct = 0.5/tan_a;

        const double a1 = 1.0;
        const double a2 = -b*b;
        const double a3 = 3.0*ct;
        const double a7 = -ct;

        //this is the rcb polynomial
        auto rcb_poly = [a1, a2, a3, a7](double x) -> double
        {
            const double x2 = x*x;
            const double x4 = x2*x2;
            const double x5 = x4*x;
            const double x6 = x5*x;
            return a1*x6 + a2*x5 + a3*x4 + a7;
        };
        
        double xi = 0.0;
        double xf = 0.0;

        const double ap = 6.0*a1;
        const double bp = 5.0*a2;
        const double cp = 4.0*a3;

        const double cx = bp*bp;
        const double cy = -4.0*ap*cp;
        const double cxmy = cx + cy;

        //ugly
        //find interval in which root lies
        if (cxmy < 0)
        {
            xi = 0.0;
            xf = std::max(1.0, std::abs(a2) + std::abs(a7));
        }
        else
        {
            const double x = cy/cx;
            const double sqrtcxmy = std::sqrt(cxmy);
            double x1 = -(bp + sqrtcxmy)/(2.0*ap);
            if (std::abs(x) < 0.01)
            {
                x1 = -cp/bp;
            }
            const double x2 = (sqrtcxmy - bp)/(2.0*ap);
            const double x24 = x2*x2*x2*x2;
            const double x25 = x24*x2;
            const double x26 = x25*x2;

            const double trm_neg = a7 + a2*x25;
            const double trm_pos = a3*x24 + a1*x26;
            const double frcb = trm_pos + trm_neg;
            if (frcb >= 0.0)
            {
                const double ratio = frcb / (trm_pos - trm_neg);
                if (ratio < 1e-5)
                {
                    //found the root
                    xi = x2;
                    xf = x2;
                }
                else
                {
                    xi = 0.0;
                    xf = x1;
                }
            }
            else
            {
                xi = x2;
                xf = std::max(1.0, std::abs(a2) + std::abs(a7));
            }

        }

        if (xi == xf)
        {
            const double r = std::sqrt(std::max(xi, 0.0));
            return std::max(r, 1e-4);
        }

        std::uintmax_t max_iter = 200;
        auto res = boost::math::tools::toms748_solve(rcb_poly, xi, xf, abs_tol, max_iter);

        const double rcb =  std::sqrt((res.first + res.second)/2.0);
        return std::max(rcb, 1e-4);
    }

    inline double find_rc0_root(double tan_a)
    {
        const double a1 = 2.0;
        const double a3 = -3.0/tan_a;
        const double a7 = 5.0/tan_a;

        //this is the rc0 polynomial
        auto rc0_poly = [a1, a3, a7](double x) -> double
        {
            const double x2 = x*x;
            const double x4 = x2*x2;
            const double x6 = x4*x2;
            return a1*x6 + a3*x4 + a7;
        };

        double xi = 0.0;
        double xf = 0.0;

        const double x2 = std::sqrt((-2.0*a3/(3.0*a1)));
        const double x24 = x2*x2*x2*x2;
        const double x26 = x24*x2*x2;

        const double trm_neg = a3*x24;
        const double trm_pos = a1*x26 + a7;
        const double frc0 = trm_pos + trm_neg;

        //find interval in which root lies
        if (frc0 > 0.0)
        {
            const double ratio = frc0 / (trm_pos - trm_neg);
            if (ratio < 1e-5)
            {
                xi = x2;
                xf = x2;
            }
            else
            {
                xi = 1.0e-4;
                xf = 1.0e-4;
            }
        }
        else
        {
            xi = x2;
            xf = std::sqrt(-a3/a1);
        }
        
        if (xi == xf)
        {
            const double r = std::sqrt(std::max(xi, 0.0));
            return std::max(r, 1e-4);
        }

        std::uintmax_t max_iter = 200;
        auto res = boost::math::tools::toms748_solve(rc0_poly, xi, xf, abs_tol, max_iter);

        const double rc0 =  std::sqrt((res.first + res.second)/2.0);

        return std::max(rc0, 1e-4);
    }

    inline RCS_Result rcs(bool iflag, double tan_a, double b)
    {
        const double ct = 0.5/tan_a;
        const double rcb = find_rcb_root(tan_a, b);
        const double rc0 = find_rc0_root(tan_a);

        RCS_Result result{rcb, 0.5, 0.5};

        const double yxb = rcb / rc0;
        if (yxb < 0.98) 
        {
            result.rmag = std::asin(yxb);
        }

        if (iflag) 
        {
            const double rc04 = rc0*rc0*rc0*rc0;
            const double rc012 = rc04*rc04*rc04;

            const double arg = 1.0 + 3.0*ct/rc04 - ct/rc012;

            if (arg >= 0.0) 
            {
                result.b0 = rc0*std::sqrt(arg);
            }
        }

        return result;
    }


    template <unsigned N>
    inline double scattering_angle(double rm, double rcb, double b, double tan_a)
    {
        auto f = [=](double theta) -> double
        {
            return scattering_integrand(rcb, tan_a, theta);
        };

        using boost::math::quadrature::gauss;
        const double integral = gauss<double, N>::integrate(f, 0.0, rm)
                              + gauss<double, N>::integrate(f, rm, cn::pi/2.0);

        return cn::pi - 2.0*b/rcb*integral;
    }
}