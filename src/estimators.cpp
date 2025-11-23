#include "estimators.h"
#include <iostream>
#include <cmath>
#include <limits>

// ------------------------
// Least Squares
// ------------------------
std::tuple<Sig, SensorMod> ls(const SensorMod& s, const Sig& y) {
    if (s.pe.size() == 0) throw std::runtime_error("pe must be defined");

    int N = y.t.size();
    int nx = s.nn[0];
    int ny = s.nn[2];

    Eigen::MatrixXd x_est(nx, N);
    Eigen::MatrixXd y_est(ny, N);
    std::vector<Eigen::MatrixXd> Px_list, Py_list;

    for (int k=0; k<N; ++k) {
        Eigen::MatrixXd H(nx, nx); // numgrad approximation
        double eps = 1e-6;
        Eigen::VectorXd x0 = s.x0;
        Eigen::VectorXd f0 = s.h(y.t[k], x0, y.u.col(k), s.th);
        H.resize(ny, nx);
        for (int i=0; i<nx; ++i) {
            Eigen::VectorXd x1 = x0;
            x1[i] += eps;
            H.col(i) = (s.h(y.t[k], x1, y.u.col(k), s.th) - f0)/eps;
        }
        Eigen::MatrixXd PP = (H.transpose()*H).completeOrthogonalDecomposition().pseudoInverse();
        Eigen::VectorXd ff = H.transpose() * (y.y.col(k) - f0 + H*x0);
        x_est.col(k) = PP*ff;
        y_est.col(k) = f0;
        Px_list.push_back(PP*(H.transpose()*s.pe*H)*PP);
        Py_list.push_back(H*Px_list.back()*H.transpose());
    }

    Sig xhat;
    xhat.x = x_est;
    xhat.y = y_est;
    xhat.t = y.t;
    xhat.u = y.u;
    xhat.xMC = Px_list; // temporärt lagring av covariance

    SensorMod shat = s;
    if (N==1) {
        shat.x0 = x_est.col(0);
    }
    return {xhat, shat};
}

// ------------------------
// Weighted Least Squares
// ------------------------
std::tuple<Sig, SensorMod> wls(const SensorMod& s, const Sig& y) {
    if (s.pe.size() == 0) throw std::runtime_error("pe must be defined");

    int N = y.t.size();
    int nx = s.nn[0];
    int ny = s.nn[2];

    Eigen::MatrixXd x_est(nx, N);
    Eigen::MatrixXd y_est(ny, N);

    Eigen::MatrixXd R = s.pe;
    Eigen::MatrixXd Rinv = R.completeOrthogonalDecomposition().pseudoInverse();

    for (int k=0; k<N; ++k) {
        Eigen::MatrixXd H(nx, nx); 
        double eps = 1e-6;
        Eigen::VectorXd x0 = s.x0;
        Eigen::VectorXd f0 = s.h(y.t[k], x0, y.u.col(k), s.th);
        H.resize(ny, nx);
        for (int i=0; i<nx; ++i) {
            Eigen::VectorXd x1 = x0;
            x1[i] += eps;
            H.col(i) = (s.h(y.t[k], x1, y.u.col(k), s.th) - f0)/eps;
        }
        Eigen::MatrixXd PP = (H.transpose()*Rinv*H).completeOrthogonalDecomposition().pseudoInverse();
        Eigen::VectorXd ff = H.transpose()*Rinv*(y.y.col(k)-f0 + H*x0);
        x_est.col(k) = PP*ff;
        y_est.col(k) = f0;
    }

    Sig xhat;
    xhat.x = x_est;
    xhat.y = y_est;
    xhat.t = y.t;
    xhat.u = y.u;

    SensorMod shat = s;
    if (N==1) shat.x0 = x_est.col(0);
    return {xhat, shat};
}

// ------------------------
// LH1 (1D likelihood)
// ------------------------
std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> lh1(
    const SensorMod& s, const Sig& y, 
    const Eigen::VectorXd* x1, int ind)
{
    int nx = s.nn[0];
    if (ind >= nx) throw std::runtime_error("ind exceeds state dimension");
    Eigen::VectorXd grid = x1 ? *x1 : Eigen::VectorXd::LinSpaced(100, 0.0, 1.0);
    int N = grid.size();
    Eigen::VectorXd lh(N), px0 = Eigen::VectorXd::Ones(N), px(N);
    px = px0;

    for (int k=0; k<y.t.size(); ++k) {
        for (int n=0; n<N; ++n) {
            Eigen::VectorXd x = s.x0;
            x[ind] = grid[n];
            Eigen::VectorXd yhat = s.h(y.t[k], x, y.u.col(k), s.th);
            Eigen::VectorXd epsi = y.y.col(k) - yhat;
            double val = epsi.transpose()*s.pe.completeOrthogonalDecomposition().pseudoInverse()*epsi;
            lh[n] = std::exp(-0.5*val);
        }
        px = px.cwiseProduct(lh);
    }

    return {lh, grid, px};
}

// ------------------------
// LH2 (2D likelihood)
// ------------------------
std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd>
lh2(const SensorMod& s, const Sig& y,
    const Eigen::VectorXd* x1_ptr,
    const Eigen::VectorXd* x2_ptr,
    const std::vector<int>& ind)
{
    if (ind.size() != 2) throw std::runtime_error("ind must have size 2");
    int nx = s.nn[0];
    if (ind[0]>=nx || ind[1]>=nx) throw std::runtime_error("ind exceeds state dimension");

    Eigen::VectorXd x1 = x1_ptr ? *x1_ptr : Eigen::VectorXd::LinSpaced(30, s.x0[ind[0]]-0.5, s.x0[ind[0]]+0.5);
    Eigen::VectorXd x2 = x2_ptr ? *x2_ptr : Eigen::VectorXd::LinSpaced(30, s.x0[ind[1]]-0.5, s.x0[ind[1]]+0.5);

    Eigen::MatrixXd X1, X2;
    Eigen::MatrixXd lh(x1.size(), x2.size());
    Eigen::MatrixXd px(x1.size(), x2.size()), px0(px);

    for (int i=0; i<x1.size(); ++i) {
        for (int j=0; j<x2.size(); ++j) {
            Eigen::VectorXd xx = s.x0;
            xx[ind[0]] = x1[i];
            xx[ind[1]] = x2[j];
            double val = 0.0;
            for (int k=0; k<y.t.size(); ++k) {
                Eigen::VectorXd yhat = s.h(y.t[k], xx, y.u.col(k), s.th);
                Eigen::VectorXd epsi = y.y.col(k) - yhat;
                val += epsi.transpose()*s.pe.completeOrthogonalDecomposition().pseudoInverse()*epsi;
            }
            lh(i,j) = std::exp(-0.5*val);
            px(i,j) = lh(i,j); 
            px0(i,j) = 1.0;
        }
    }

    // Return meshgrid for plotting
    X1 = x1.replicate(1,x2.size());
    X2 = x2.transpose().replicate(x1.size(),1);

    return {lh, x1, x2, px, px0, X1, X2};
}
