#include "./include/Legendre.h"

#include <stdexcept>

std::vector<double> Legendre::get_coeffs_stdvec() const
{
    std::vector<double> c(coeffs.data(), coeffs.data() + coeffs.size());
    return c;
}

//map from [vmin, vmax] to [-1, 1]
double Legendre::map_to_domain(double v, double vmin, double vmax)
{
    return 2.0*(v - vmin)/(vmax - vmin) - 1.0;
}

double Legendre::eval_scalar(double v) const
{
    const double u = map_to_domain(v, vmin, vmax);

    double y_kp1 = 0.0;
    double y_kp2 = 0.0;

    for (int k = degree; k >= 1; k--)
    {
        const double dk = static_cast<double>(k);
        const double y_k = ((2.0*dk + 1.0)/(dk + 1.0))*u*y_kp1 -
            ((dk + 1.0)/(dk + 2.0))*y_kp2 + coeffs(k);

        y_kp2 = y_kp1;
        y_kp1 = y_k;
    }
    return coeffs(0) + u*y_kp1 - y_kp2/2.0;
}

//Clenshaw’s algorithm
double Legendre::operator()(double v) const
{
    return eval_scalar(v);
}

//TODO: maybe use spans
eigen_vec Legendre::operator()(const eigen_vec& v) const
{
    eigen_vec res(v.size());

    for (int i = 0; i < v.size(); i++)
    {
        res(i) = eval_scalar(v(i));
    }

    return res;
}

eigen_mat Legendre::create_pv_matrix(Eigen::Ref<const eigen_vec> x, double vmin, double vmax, int deg)
{
    eigen_mat A(x.size(), deg + 1);

    for (int i = 0; i < x.size(); i++)
    {
        const double u = map_to_domain(x(i), vmin, vmax);

        A(i, 0) = 1.0;
        A(i, 1) = u;

        for (int l = 1; l < deg; l++)
        {
            const double w = static_cast<double>(l) + 1.0;
            A(i, l + 1) = ((2.0*l + 1.0)/w)*u*A(i, l) - (l/w)*A(i, l - 1);
        }
    }

    return A;
}

//factory
Legendre Legendre::fit(double vmin, double vmax, int deg, std::span<const double> x, std::span<const double> y)
{
    if (deg < 2)
    {
        throw std::invalid_argument("Degree of polynomial is less than two.");
    }

    Eigen::Map<const eigen_vec> ex{x.data(), static_cast<Eigen::Index>(x.size())};
    Eigen::Map<const eigen_vec> ey{y.data(), static_cast<Eigen::Index>(y.size())};

    const eigen_mat A = create_pv_matrix(ex, vmin, vmax, deg);
    const eigen_vec coeffs = A.colPivHouseholderQr().solve(ey);

    //construct obj and return in one line
    return {vmin, vmax, deg, coeffs};;
}