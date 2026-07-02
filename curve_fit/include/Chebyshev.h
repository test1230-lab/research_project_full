#pragma once

#include <vector>
#include <span>
#include <Eigen/Dense>

using eigen_vec = Eigen::VectorXd;
using eigen_mat = Eigen::MatrixXd;

//change some of the v vars to be named x instead
//degree is expected to be at least 2, will throw if it is < 2
class Chebyshev
{
public:
    Chebyshev(double lo, double hi, int d, const eigen_vec& c) //TODO: how should I pass c??
        : vmin(lo), vmax(hi), degree(d), coeffs(c){} //maybe std::move for coeffs

    const eigen_vec& get_coeffs() const { return coeffs; }
    int get_degree() const { return degree; }
    std::vector<double> get_coeffs_stdvec() const;
    double operator()(double v) const;
    eigen_vec operator()(const eigen_vec& v) const;
    double first_deriv(double v) const;

    //static Chebyshev fit(double xmin, double xmax, int deg, const eigen_vec& x, const eigen_vec& y);
    static Chebyshev fit(double xmin, double xmax, int deg, std::span<const double> x, std::span<const double> y);

private:
    double vmin, vmax;
    int degree;
    eigen_vec coeffs;

    static double map_to_domain(double v, double vmin, double vmax);
    static eigen_mat create_pv_matrix(Eigen::Ref<const eigen_vec> x, double vmin, double vmax, int deg);
    double eval_scalar(double v) const;

};
