#pragma once

#include <span>
#include <array>

#include <Eigen/Dense>

using eigen_vec = Eigen::VectorXd;
using eigen_mat = Eigen::MatrixXd;

class TailFit
{
public:
    TailFit(double lo, double hi, const eigen_vec& c) : vmin(lo), vmax(hi), coeffs(c) {}
    
    const eigen_vec& get_coeffs() const { return coeffs; }
    std::array<double, 5> get_coeffs_arr() const;

    double operator()(double v) const;
    eigen_vec operator()(const eigen_vec& v) const;
    static TailFit fit(double vmin, double vmax, std::span<const double> x, std::span<const double> y);

private:
    double vmin, vmax;
    eigen_vec coeffs;

    static double tail_model(double x, const eigen_vec& c);
    static double map_to_domain(double v, double vmin, double vmax);
    static eigen_mat create_v_matrix(Eigen::Ref<const eigen_vec> x, double vmin, double vmax);
};