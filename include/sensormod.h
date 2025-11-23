#pragma once
#include "sig.h"
#include <Eigen/Dense>
#include <functional>
#include <random>
#include <vector>
#include <string>
#include <cmath> // För std::pow, std::sqrt, M_PI

class SensorMod {
public:
    using DynFunc = std::function<Eigen::VectorXd(double,
                                                  const Eigen::VectorXd&,
                                                  const Eigen::VectorXd&,
                                                  const Eigen::VectorXd&)>;

    // -----------------------------------------------------
    // Befintliga medlemmar
    // -----------------------------------------------------
    Eigen::Vector4i nn;      // [nx, nu, ny, nth]

    Eigen::VectorXd x0;
    Eigen::VectorXd th;
    Eigen::MatrixXd pv; // process noise
    Eigen::MatrixXd pe; // measurement noise (ny x ny)
    DynFunc h;          // sensorfunktion
    double fs = 0.0;

    std::string name;
    std::string xlabel, ylabel, ulabel, thlabel;

    // Konstruktor (oförändrad)
    SensorMod(const DynFunc& hFunc, const Eigen::Vector4i& nn_)
        : h(hFunc), nn(nn_)
    {
        x0 = Eigen::VectorXd::Zero(nn[0]);
        th = Eigen::VectorXd::Zero(nn[3]);
        pv = Eigen::MatrixXd::Zero(nn[0], nn[0]);
        pe = Eigen::MatrixXd::Zero(nn[2], nn[2]);
    }
    
    // -----------------------------------------------------
    // 1. simulate (Original funktion)
    // -----------------------------------------------------
    Sig simulate(const Eigen::VectorXd& t,
                 const Eigen::MatrixXd* xOverride = nullptr,
                 const Eigen::MatrixXd* uOverride = nullptr,
                 int MC = 0)
    {
        int N = t.size();
        int nx = nn[0];
        int nu = nn[1];
        int ny = nn[2];

        Eigen::MatrixXd x = xOverride ? *xOverride : Eigen::MatrixXd::Zero(nx, N).colwise() + x0;
        Eigen::MatrixXd u = uOverride ? *uOverride : Eigen::MatrixXd::Zero(nu, N);
        Eigen::MatrixXd y(ny, N);

        std::random_device rd;
        std::mt19937 gen(rd());
        // OBS: Denna rad använder endast diagonalen (antar okorrelerat brus)
        Eigen::VectorXd pe_std = pe.diagonal().cwiseSqrt(); 

        for (int k = 0; k < N; ++k) {
            Eigen::VectorXd yk = h(t[k], x.col(k), u.col(k), th);
            for (int i = 0; i < ny; ++i) {
                std::normal_distribution<> d(0.0, pe_std[i]);
                yk[i] += d(gen);
            }
            y.col(k) = yk;
        }

        Sig z;
        z.y = y;
        z.t = t;
        z.x = x;
        z.u = u;

        if (MC > 0) {
            std::normal_distribution<> d_th(0.0, 1.0);
            for (int i = 0; i < MC; ++i) {
                Eigen::VectorXd thMC = th;
                for (int j = 0; j < nn[3]; ++j) thMC[j] += d_th(gen) * std::sqrt(pv(j,j));
                Eigen::MatrixXd yMC_i(ny, N);
                for (int k = 0; k < N; ++k) {
                    yMC_i.col(k) = h(t[k], x.col(k), u.col(k), thMC);
                }
                z.yMC.push_back(yMC_i);
                z.xMC.push_back(x);
            }
        }

        return z;
    }

    // -----------------------------------------------------
    // 2. likelihood_function (Ny funktion, motsvarar lh1/lh2)
    // -----------------------------------------------------
    /**
     * @brief Beräknar den totala Likelihood-funktionen L(x | y) över ett rutnät av tillståndspunkter.
     * @param y Mätdata som Sig-objekt.
     * @param xGrid Matris med tillståndspunkter att utvärdera (N_grid x nx).
     * @param state_indices Ej implementerad, men behålls för framtida kompatibilitet med lh1/lh2.
     * @return Eigen::VectorXd Vektor med sannolikhetsvärden (Likelihood) för varje gridpunkt.
     */
    Eigen::VectorXd likelihood_function(
        const Sig& y, 
        const Eigen::MatrixXd& xGrid,
        const std::vector<int>& state_indices = {} // {0} är default i C++
    )
    {
        int N_grid = xGrid.rows();
        int N_time = y.t.size();
        int ny = nn[2];
        int nx = nn[0];
        
        if (xGrid.cols() != nx) {
            throw std::runtime_error("likelihood_function: xGrid måste ha nx kolumner.");
        }
        if (pe.rows() != ny || pe.cols() != ny) {
            throw std::runtime_error("likelihood_function: Mätbruskovarians Pe är fel dimensionerad.");
        }

        Eigen::VectorXd L = Eigen::VectorXd::Ones(N_grid); // Initiera Likelihood till 1.0
        Eigen::VectorXd mu_v = Eigen::VectorXd::Zero(ny);  // Mätbrusets medelvärde är noll

        // Iterera över alla mätningstidpunkter (k)
        for (int k = 0; k < N_time; ++k) {
            Eigen::MatrixXd epsi_grid(N_grid, ny); // Mätfelsmatris (N_grid x ny)
            
            double t_k = y.t[k];
            // OBS: y.u kan vara tom (nu=0). Kontrollera om u ska användas eller ignoreras.
            Eigen::VectorXd u_k = (y.u.cols() > 0) ? y.u.col(k) : Eigen::VectorXd::Zero(nn[1]); 
            
            for (int i = 0; i < N_grid; ++i) {
                Eigen::VectorXd x_i = xGrid.row(i).transpose();
                Eigen::VectorXd y_hat_i = h(t_k, x_i, u_k, th);
                
                // Mätning y_k - Prediktion y_hat_i
                Eigen::VectorXd epsi_i = y.y.col(k) - y_hat_i;
                epsi_grid.row(i) = epsi_i.transpose();
            }

            // Beräkna PDF för felet p(epsi_k | Pe) = p(y_k | x)
            Eigen::VectorXd p_y_given_x = multivariate_normal_pdf(pe, mu_v, epsi_grid);
            
            // Multiplicera Likelihood: L_total = L_prev * p(y_k | x)
            L.array() *= p_y_given_x.array();
        }

        return L;
    }

private:
    // -----------------------------------------------------
    // 3. multivariate_normal_pdf (Ny intern hjälpfuntion)
    // -----------------------------------------------------
    Eigen::VectorXd multivariate_normal_pdf(
        const Eigen::MatrixXd& P, 
        const Eigen::VectorXd& mu, 
        const Eigen::MatrixXd& x) const
    {
        int nx = mu.size();
        int N = x.rows();
        
        if (P.rows() != nx || P.cols() != nx || x.cols() != nx) {
            // Detta ska fångas av anropande funktion, men är en sista utväg
            Eigen::VectorXd p(N);
            p.setZero();
            return p;
        }
        
        Eigen::VectorXd p(N);
        double detP = P.determinant();

        if (detP <= 0 || !P.isApprox(P.transpose())) {
            // Sätt sannolikheten till 0 om matrisen är singulär eller icke-symmetrisk
            p.setZero();
            return p;
        }

        // Konstant term: 1 / sqrt((2*pi)^nx * det(P))
        double constant = 1.0 / (std::pow(2.0 * M_PI, (double)nx / 2.0) * std::sqrt(detP));
        
        // Använder P.llt().solve() för numeriskt stabil lösning av P_inv * epsi
        Eigen::LLT<Eigen::MatrixXd> llt(P);
        if(llt.info() != Eigen::Success) {
            // Fallback om Cholesky-faktorisering misslyckas (icke-positivt definit)
            p.setZero();
            return p;
        }

        for (int i = 0; i < N; ++i) {
            Eigen::VectorXd epsi = x.row(i).transpose() - mu; 
            
            // Lösning av temp = P_inv * epsi (numeriskt stabilt)
            Eigen::VectorXd temp = llt.solve(epsi);
            
            // Kvadratisk form: V = epsi^T * P_inv * epsi
            double V = epsi.dot(temp); 

            // PDF: p = const * exp(-0.5 * V)
            p[i] = constant * std::exp(-0.5 * V);
        }

        return p;
    }
};