#pragma once
#include <Eigen/Dense>
#include <vector>

struct Sig {
    Eigen::MatrixXd y;      // mätningar
    Eigen::VectorXd t;      // tid
    Eigen::MatrixXd u;      // insignal
    Eigen::MatrixXd x;      // tillstånd

    // MC-simuleringar
    std::vector<Eigen::MatrixXd> yMC;
    std::vector<Eigen::MatrixXd> xMC;

    // CRLB
    Eigen::MatrixXd Px;

    Sig() {}
    Sig(int ny, int N) : y(Eigen::MatrixXd::Zero(ny, N)), 
                          t(Eigen::VectorXd::Zero(N)),
                          u(Eigen::MatrixXd::Zero(0, N)),
                          x(Eigen::MatrixXd::Zero(0, N)),
                          Px(Eigen::MatrixXd::Zero(0,0)) {}
};
