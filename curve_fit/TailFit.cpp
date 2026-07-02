#include "./include/TailFit.h"

#include <cmath>
#include <stdexcept>

double TailFit::operator()(double v) const
{
    const double x  = (v - v_stitch)/scale;
    const double x2 = x*x;
    const double x3 = x*x2;
    const double x4 = x2*x2;

    return logf_stitch + scale*dlogf_stitch*x + q(0)*x2 + q(1)*x3 + q(2)*x4;
}

TailFit TailFit::fit(double v_stitch, double scale, double logf_stitch, double dlogf_stitch,
        std::span<const double> v, std::span<const double> y, std::span<const double> weights)
{
    if (v.size() != y.size())
    {
        throw std::invalid_argument("TailFit::fit: v and y must have equal length.");
    }

    if (v.size() < 3)
    {
        throw std::invalid_argument("TailFit::fit: need >= 3 points to fit q2, q3, q4.");
    }

    if (scale <= 0.0)
    {
        throw std::invalid_argument("TailFit::fit: scale must be positive.");
    }

    const bool weighted = !weights.empty();
    if (weighted && weights.size() != v.size())
    {
        throw std::invalid_argument("TailFit::fit: weights must match v length (or be empty).");
    }

    const int n = v.size();
    eigen_mat A(n, 3);
    eigen_vec b(n);

    const double sms = scale*dlogf_stitch;   // s*m_s

    for (int i = 0; i < n; i++)
    {
        if (weighted && weights[i] < 0.0)
        {
            throw std::invalid_argument("TailFit::fit: weights must be non-negative.");
        }

        const double x  = (v[i] - v_stitch) / scale;
        const double x2 = x*x;
        const double r  = weighted ? std::sqrt(weights[i]) : 1.0;   //sqrt(w) row scaling

        A(i, 0) = r*x2;
        A(i, 1) = r*x2*x;
        A(i, 2) = r*x2*x2;

        //rhs with the anchor (value + slope) subtracted off
        b(i) = r*(y[i] - logf_stitch - sms*x);
    }

    //QR least-squares; avoid normal equations (squares the condition number)
    const vec3 q = A.colPivHouseholderQr().solve(b);

    return TailFit{v_stitch, scale, logf_stitch, dlogf_stitch, q};
}

TailFit TailFit::fit(double v_lo, double v_hi, std::span<const double> x,
        std::span<const double> y, const Chebyshev& cheb, std::span<const double> weights)
{
    //stitch at the inner edge (smaller |v|); scale = window width, so far edge -> |x| ~ 1
    const double v_stitch = (std::abs(v_lo) <= std::abs(v_hi)) ? v_lo : v_hi;
    const double scale = v_hi - v_lo;

    return fit(v_stitch, scale, cheb(v_stitch), cheb.first_deriv(v_stitch), x, y, weights);
}
