#include "./include/Chebyshev.h"

#include <stdexcept>

std::vector<double> Chebyshev::get_coeffs_stdvec() const
{
    std::vector<double> c(coeffs.data(), coeffs.data() + coeffs.size());
    return c;
}

//map from [vmin, vmax] to [-1, 1]
double Chebyshev::map_to_domain(double v, double vmin, double vmax)
{
    return 2.0*(v - vmin)/(vmax - vmin) - 1.0;
}

double Chebyshev::eval_scalar(double v) const
{
    const double u = map_to_domain(v, vmin, vmax);

    double b_i1 = 0.0;
    double b_i2 = 0.0;

    for (int i = coeffs.size() - 1; i >= 1; i--) 
    {
        const double b_i = 2.0*u*b_i1 - b_i2 + coeffs(i);
        b_i2 = b_i1;
        b_i1 = b_i;
    }

    return u*b_i1 - b_i2 + coeffs(0);
}

//d/dv of sum_k c_k T_k(u): carry T_k and T_k' through their recurrences together,
//then chain-rule by du/dv = 2/(vmax - vmin)
double Chebyshev::first_deriv(double v) const
{
    const double u = map_to_domain(v, vmin, vmax);
    const double du_dv = 2.0/(vmax - vmin);

    double t_prev = 1.0;   // T_0
    double t_curr = u;     // T_1
    double dt_prev = 0.0;  // T_0'
    double dt_curr = 1.0;  // T_1'

    double deriv = 0.0;//k=0 term is 0
    if (coeffs.size() > 1)
    {
        deriv += coeffs(1)*dt_curr;        //k=1
    }

    for (int k = 2; k < coeffs.size(); k++)
    {
        const double t_next = 2.0*u*t_curr - t_prev;
        const double dt_next = 2.0*t_curr + 2.0*u*dt_curr - dt_prev;

        deriv += coeffs(k)*dt_next;

        t_prev = t_curr;
        t_curr = t_next;
        dt_prev = dt_curr;
        dt_curr = dt_next;
    }

    return deriv*du_dv;
}

//Clenshaw’s algorithm
double Chebyshev::operator()(double v) const
{
    return eval_scalar(v);
}

//TODO: maybe use spans
eigen_vec Chebyshev::operator()(const eigen_vec& v) const
{
    eigen_vec res(v.size());

    for (int i = 0; i < v.size(); i++)
    {
        res(i) = eval_scalar(v(i));
    }

    return res;
}

eigen_mat Chebyshev::create_pv_matrix(Eigen::Ref<const eigen_vec> x, double vmin, double vmax, int deg)
{
    eigen_mat A(x.size(), deg + 1);

    for (int i = 0; i < x.size(); i++)
    {
        const double u = map_to_domain(x(i), vmin, vmax);

        A(i, 0) = 1.0; //T_0 = 1
        A(i, 1) = u; //T_1 = x

        for (int n = 1; n < deg; n++)
        {
            A(i, n + 1) = 2.0*u*A(i, n) - A(i, n - 1);
        }
    }

    return A;
}

//factory
Chebyshev Chebyshev::fit(double vmin, double vmax, int deg, std::span<const double> x, std::span<const double> y)
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
    return {vmin, vmax, deg, coeffs};
}