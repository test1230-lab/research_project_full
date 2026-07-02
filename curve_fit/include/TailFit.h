#pragma once

#include <span>

#include <Eigen/Dense>

#include "Chebyshev.h"

using vec3 = Eigen::Vector3d;

//anchored quartic tail in centered log-space:
//  ell(x) = logf_s + s*m_s*x + q2*x^2 + q3*x^3 + q4*x^4,  x = (v - v_stitch)/s
//anchor (value logf_s, slope m_s) is fixed at the seam; only q = (q2,q3,q4) is fit
class TailFit
{
public:
    TailFit(double v_stitch, double scale, double logf_stitch, double dlogf_stitch, const vec3& q)
        : v_stitch(v_stitch), scale(scale), logf_stitch(logf_stitch), dlogf_stitch(dlogf_stitch), q(q) {}

    const vec3& get_coeffs() const { return q; }   // [q2, q3, q4]

    double operator()(double v) const;   //log-space value

    //anchor passed in explicitly; weights empty => unweighted
    static TailFit fit(double v_stitch, double scale, double logf_stitch, double dlogf_stitch, 
           std::span<const double> v, std::span<const double> y, std::span<const double> weights = {});

    //anchor read from the central Chebyshev fit at the seam
    static TailFit fit(double v_lo, double v_hi, std::span<const double> x,
           std::span<const double> y, const Chebyshev& cheb, std::span<const double> weights = {});

private:
    double v_stitch, scale, logf_stitch, dlogf_stitch;
    vec3 q;   // [q2, q3, q4]
};
