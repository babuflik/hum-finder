#pragma once
#include "sensormod.h" // Inkludera basklassen SensorMod
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <string>
#include <random>

// Simple check for covariance matrix validity.
// Requires additional implementation for full MATLAB iscov parity.
inline bool is_valid_covariance(const Eigen::MatrixXd& P) {
    if (P.rows() != P.cols()) return false; // Must be square
    if (!P.isApprox(P.transpose())) return false; // Must be symmetric
    
    // For positive semidefiniteness (PSD), use Cholesky factorization (LLT).
    // A stricter check would be required for full MATLAB iscov parity.
    Eigen::LLT<Eigen::MatrixXd> llt(P);
    return llt.info() == Eigen::Success;
}


class NDist : public SensorMod {
public:
    using Vector = Eigen::VectorXd;
    using Matrix = Eigen::MatrixXd;
    
    // ------------------------------------------------------------------
    // MATLAB-like properties
    // ------------------------------------------------------------------
    // mu and P are stored in NDist (not in SensorMod)
    Vector mu;
    Matrix P;

    int MC = 1000;          // Monte Carlo samples (default)
    long state;             // Random state (simulerar MATLAB's round(1e6*rand))
    std::vector<std::string> xlabel;

    // ------------------------------------------------------------------
    // Constructors (MATLAB ndist compatibility)
    // ------------------------------------------------------------------

    // 1. Tom konstruktor (ndist)
    NDist() : SensorMod(nullptr, Eigen::Vector4i::Zero()) {
        mu = Vector();
        P = Matrix();
        state = generate_random_state();
    }
    
    // 2. Full constructor (ndist(mu, P))
    NDist(const Vector& mu_, const Matrix& P_) : SensorMod(nullptr, get_nn(mu_)) {
        initialize(mu_, P_);
        state = generate_random_state();
    }

    // 3. Copy constructor (ndist(X) or ndist(gmdist))
    // We use a generic copy approach to handle ndist/gmdist cases
    NDist(const NDist& other) : SensorMod(other) { // Anropa SensorMod kopieringskonstruktor
        // Kopiera ndist-specifika medlemmar
        mu = other.mu;
        P = other.P;
        MC = other.MC;
        xlabel = other.xlabel;
        state = generate_random_state(); // Ny random state vid kopiering
    }
    
    // NOTE: Implementing ndist(gmdist) requires gmdist to be defined
    // and provide mean()/cov() methods; that is out of scope here.

    // ------------------------------------------------------------------
    // Setters (corresponding to MATLAB set.mu/set.P)
    // ------------------------------------------------------------------

    void set_mu(const Vector& new_mu) {
        if (!mu.isApprox(Vector()) && new_mu.size() != mu.size()) {
            throw std::runtime_error("NDist::set_mu: Cannot change the dimension of mu.");
        }
        mu = new_mu;
        // SensorMod::nn must be updated if dimension changes (when using SensorMod)
        // This is complex in C++; we rely on constructors to set it correctly.
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
    // Indexing (similar to MATLAB subsref X(i))
    // ------------------------------------------------------------------

    // Use a function instead of operator() to avoid Eigen-related complexity
    NDist subsref(const std::vector<int>& indices) const {
        // Kontrollera index
        for (int i : indices) {
            if (i < 1 || i > mu.size()) { // MATLAB-indexering startar vid 1
                throw std::runtime_error("NDist::subsref: Index out of range.");
            }
        }
        
        // Extract parts of mu and P (C++ index = MATLAB index - 1)
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

    // Helper to set SensorMod::nn correctly
    static Eigen::Vector4i get_nn(const Vector& mu_) {
        // nx = dimension of mu. nu, ny, nth are 0.
        return Eigen::Vector4i(mu_.size(), 0, 0, 0);
    }
    
    // Helper to initialize mu and P with validation
    void initialize(const Vector& mu_, const Matrix& P_) {
        int n = mu_.size();
        
        if (P_.rows() != n || P_.cols() != n) {
            throw std::runtime_error("NDIST: P must be square with same dimension as mu.");
        }
        
        if (!P_.array().isFinite().all() || !mu_.array().isFinite().all()) {
            throw std::runtime_error("NDIST: mu and P must contain real, finite values.");
        }
        
        // Perform is_valid_covariance check (simulates MATLAB's iscov)
        if (!is_valid_covariance(P_)) {
            // More specific error handling would be performed here, similar to MATLAB
            throw std::runtime_error("NDIST: P is not a valid covariance matrix (e.g., not PSD or symmetric).");
        }
        
        mu = mu_;
        P = P_;
        
        // Set default xlabel
        xlabel.resize(n);
        for (int i = 0; i < n; ++i) {
            xlabel[i] = "x" + std::to_string(i + 1);
        }
    }
};