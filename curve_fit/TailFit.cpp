#include "./include/TailFit.h"

std::array<double, 5> TailFit::get_coeffs_arr() const
{
    std::array<double, 5> res{};
    
    for (int i = 0; i < 5; i++)
    {
        res[i] = coeffs(i);
    }

    return res;
}

double TailFit::operator()(double v) const
{
    const double u = map_to_domain(v, vmin, vmax);
    return tail_model(u, coeffs);
}

eigen_vec TailFit::operator()(const eigen_vec& v) const
{
    eigen_vec res{v.size()};

    //the compiler should optimize this, check later
    for (int i = 0; i < v.size(); i++)
    {
        const double u = map_to_domain(v(i), vmin, vmax);
        res(i) = tail_model(u, coeffs);
    }

    return res;
}

TailFit TailFit::fit(double vmin, double vmax, std::span<const double> x, std::span<const double> y)
{
    Eigen::Map<const eigen_vec> ex{x.data(), static_cast<Eigen::Index>(x.size())};
    Eigen::Map<const eigen_vec> ey{y.data(), static_cast<Eigen::Index>(y.size())};

    const eigen_mat A = create_v_matrix(ex, vmin, vmax);
    const eigen_vec coeffs = A.colPivHouseholderQr().solve(ey);

    return {vmin, vmax, coeffs};
}

double TailFit::tail_model(double x, const eigen_vec& c)
{
    const double x2 = x*x;
    const double x3 = x*x2;
    const double x4 = x2*x2;

    return c(0)*x4 + c(1)*x3 + c(2)*x2 + c(3)*x + c(4);
}

double TailFit::map_to_domain(double v, double vmin, double vmax)
{
    return 2.0*(v - vmin)/(vmax - vmin) - 1.0;
}

eigen_mat TailFit::create_v_matrix(Eigen::Ref<const eigen_vec> x, double vmin, double vmax)
{
    eigen_mat A(x.size(), 5); // 5 is degree + 1

    for (int i = 0; i < x.size(); i++)
    {
        const double u = map_to_domain(x(i), vmin, vmax);

        const double u2 = u*u;
        const double u3 = u*u2;
        const double u4 = u2*u2;

        A(i, 0) = u4;
        A(i, 1) = u3;
        A(i, 2) = u2;
        A(i, 3) = u;
        A(i, 4) = 1.0;
    }

    return A;
}