#pragma once

#include <vector>
#include <span>
#include <Eigen/Dense>

using eigen_vec = Eigen::VectorXd;
using eigen_mat = Eigen::MatrixXd;

//degree is expected to be at least 2, will throw if it is < 2
class Legendre
{
public:
    Legendre(double lo, double hi, int d, const eigen_vec& c) 
        : vmin(lo), vmax(hi), degree(d), coeffs(c){} //maybe std::move for coeffs

    const eigen_vec& get_coeffs() const { return coeffs; }
    int get_degree() const { return degree; }
    std::vector<double> get_coeffs_stdvec() const;
    double operator()(double v) const;
    eigen_vec operator()(const eigen_vec& v) const;

    static Legendre fit(double xmin, double xmax, int deg, std::span<const double> x, std::span<const double> y);

private:
    double vmin, vmax;
    int degree;
    eigen_vec coeffs;

    static double map_to_domain(double v, double vmin, double vmax);
    static eigen_mat create_pv_matrix(Eigen::Ref<const eigen_vec> x, double vmin, double vmax, int deg);
    double eval_scalar(double v) const;
};