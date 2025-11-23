#pragma once
#include "sensormod.h" // Inkludera basklassen SensorMod
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <string>
#include <random>

// Hårdkodad check för kovariansmatris.
// Kräver ytterligare implementering för fullständig MATLAB iscov-logik.
inline bool is_valid_covariance(const Eigen::MatrixXd& P) {
    if (P.rows() != P.cols()) return false; // Måste vara kvadratisk
    if (!P.isApprox(P.transpose())) return false; // Måste vara symmetrisk
    
    // För positivt semidefiniteness (PSD), används Cholesky-faktorisering (LLT)
    // En striktare check än vad som visas här skulle krävas för fullständig MATLAB-paritet.
    Eigen::LLT<Eigen::MatrixXd> llt(P);
    return llt.info() == Eigen::Success;
}


class NDist : public SensorMod {
public:
    using Vector = Eigen::VectorXd;
    using Matrix = Eigen::MatrixXd;
    
    // ------------------------------------------------------------------
    // MATLAB properties
    // ------------------------------------------------------------------
    // Mu och P lagras i NDist, inte SensorMod
    Vector mu;
    Matrix P;

    int MC = 1000;          // Monte Carlo samples (default)
    long state;             // Random state (simulerar MATLAB's round(1e6*rand))
    std::vector<std::string> xlabel;

    // ------------------------------------------------------------------
    // Konstruktor (MATLAB ndist)
    // ------------------------------------------------------------------

    // 1. Tom konstruktor (ndist)
    NDist() : SensorMod(nullptr, Eigen::Vector4i::Zero()) {
        mu = Vector();
        P = Matrix();
        state = generate_random_state();
    }
    
    // 2. Fullständig konstruktor (ndist(mu, P))
    NDist(const Vector& mu_, const Matrix& P_) : SensorMod(nullptr, get_nn(mu_)) {
        initialize(mu_, P_);
        state = generate_random_state();
    }

    // 3. Kopieringskonstruktor (ndist(X) eller ndist(gmdist))
    // Vi använder en generisk kopieringsfunktion för att hantera ndist/gmdist-fall
    NDist(const NDist& other) : SensorMod(other) { // Anropa SensorMod kopieringskonstruktor
        // Kopiera ndist-specifika medlemmar
        mu = other.mu;
        P = other.P;
        MC = other.MC;
        xlabel = other.xlabel;
        state = generate_random_state(); // Ny random state vid kopiering
    }
    
    // OBS: Implementering för ndist(gmdist) kräver att gmdist är definierad
    // och har mean()/cov() metoder, vilket är utanför detta scope.

    // ------------------------------------------------------------------
    // Setters (Motsvarar MATLAB set.mu/set.P)
    // ------------------------------------------------------------------

    void set_mu(const Vector& new_mu) {
        if (!mu.isApprox(Vector()) && new_mu.size() != mu.size()) {
            throw std::runtime_error("NDist::set_mu: Cannot change the dimension of mu.");
        }
        mu = new_mu;
        // SensorMod::nn måste uppdateras om dimensionen ändras (om SensorMod används)
        // Detta är komplext i C++, så vi förlitar oss på konstruktorn.
    }

    void set_P(const Matrix& new_P) {
        if (!P.isApprox(Matrix()) && new_P.rows() != P.rows()) {
            throw std::runtime_error("NDist::set_P: Cannot change the dimension of P.");
        }
        if (!is_valid_covariance(new_P)) {
            throw std::runtime_error("NDist::set_P: P is not a valid covariance.");
        }
        P = new_P;
    }

    // ------------------------------------------------------------------
    // Indexering (Motsvarar MATLAB Xsubsref/X(i))
    // ------------------------------------------------------------------

    // Använder en funktion istället för operator() för att undvika komplexiteten med Eigen
    NDist subsref(const std::vector<int>& indices) const {
        // Kontrollera index
        for (int i : indices) {
            if (i < 1 || i > mu.size()) { // MATLAB-indexering startar vid 1
                throw std::runtime_error("NDist::subsref: Index out of range.");
            }
        }
        
        // Hämta delar av mu och P (C++ index = MATLAB index - 1)
        Vector new_mu(indices.size());
        Matrix new_P(indices.size(), indices.size());
        std::vector<std::string> new_xlabel(indices.size());

        for (size_t k = 0; k < indices.size(); ++k) {
            int i = indices[k] - 1; // Konvertera till nollbaserat index
            new_mu[k] = mu[i];
            new_xlabel[k] = xlabel[i];

            for (size_t j = 0; j < indices.size(); ++j) {
                int l = indices[j] - 1;
                new_P(k, j) = P(i, l);
            }
        }

        NDist Xi(new_mu, new_P);
        Xi.xlabel = new_xlabel;
        return Xi;
    }

private:
    long generate_random_state() const {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_real_distribution<> distrib(0.0, 1e6);
        return static_cast<long>(std::round(distrib(gen)));
    }

    // Hjälpfunktion för att sätta SensorMod::nn korrekt
    static Eigen::Vector4i get_nn(const Vector& mu_) {
        // nx = dimensionen av mu. nu, ny, nth är 0.
        return Eigen::Vector4i(mu_.size(), 0, 0, 0);
    }
    
    // Hjälpfunktion för att initialisera mu och P med felkontroll
    void initialize(const Vector& mu_, const Matrix& P_) {
        int n = mu_.size();
        
        if (P_.rows() != n || P_.cols() != n) {
            throw std::runtime_error("NDIST: P must be square with same dimension as mu.");
        }
        
        if (!P_.array().isFinite().all() || !mu_.array().isFinite().all()) {
            throw std::runtime_error("NDIST: mu and P must contain real, finite values.");
        }
        
        // Utför is_valid_covariance checken (simulerar MATLAB's iscov)
        if (!is_valid_covariance(P_)) {
            // Här skulle mer specifik felhantering ske, som i MATLAB switch-satsen
            throw std::runtime_error("NDIST: P is not a valid covariance matrix (e.g., not PSD or symmetric).");
        }
        
        mu = mu_;
        P = P_;
        
        // Sätt default xlabel
        xlabel.resize(n);
        for (int i = 0; i < n; ++i) {
            xlabel[i] = "x" + std::to_string(i + 1);
        }
    }
};