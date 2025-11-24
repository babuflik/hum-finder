#pragma once
#include <Eigen/Dense>
#include <vector>

struct Sig {
    Eigen::MatrixXd y;      // measurements
    Eigen::VectorXd t;      // time
    Eigen::MatrixXd u;      // input signal
    Eigen::MatrixXd x;      // state

    // Monte Carlo simulations
    std::vector<Eigen::MatrixXd> yMC;
    std::vector<Eigen::MatrixXd> xMC;

    // Cramer-Rao Lower Bound covariance matrix
    Eigen::MatrixXd Px;

    Sig() {}
    Sig(int ny, int N) : y(Eigen::MatrixXd::Zero(ny, N)), 
                          t(Eigen::VectorXd::Zero(N)),
                          u(Eigen::MatrixXd::Zero(0, N)),
                          x(Eigen::MatrixXd::Zero(0, N)),
                          Px(Eigen::MatrixXd::Zero(0,0)) {}
};
