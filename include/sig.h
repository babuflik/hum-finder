#pragma once
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <stdexcept>

/**
 * @brief SIG is the data class for signals (C++ implementation of MATLAB sig.m)
 * 
 * An object in the SIG class stores numerical data (y,x,u) from a system,
 * where y is the mandatory measurement, u is an optional input and
 * x is an optional state. It also stores uncertainties represented with
 * covariances or samples.
 * 
 * Fields (MATLAB compatible):
 * - y:      Primary data field (N x ny)
 * - t:      Time vector (N x 1)
 * - u:      Optional input (N x nu)
 * - x:      Optional state (N x nx)
 * - fs:     Sampling frequency (NaN for non-uniform data)
 * - Py:     Covariance for output (N x ny x ny)
 * - Px:     Covariance for state (N x nx x nx)
 * - yMC:    Monte Carlo realizations of output (MC x N x ny)
 * - xMC:    Monte Carlo realizations of state (MC x N x nx)
 * - MC:     Number of MC data
 * - ylabel: Labels for the output
 * - ulabel: Labels for the input
 * - xlabel: Labels for the state
 * - tlabel: Label for time
 * - nn:     Dimensions [nx, nu, ny]
 * - name:   Signal name
 * - desc:   Description
 * - marker: Discrete events marked in the signal
 * - markerlabel: Labels for events
 */
struct Sig {
    // Primary data fields
    Eigen::MatrixXd y;      // measurements (N x ny)
    Eigen::VectorXd t;      // time vector (N x 1)
    Eigen::MatrixXd u;      // input signal (N x nu)
    Eigen::MatrixXd x;      // state (N x nx)
    
    double fs = std::numeric_limits<double>::quiet_NaN(); // sampling frequency
    
    // Covariance matrices
    std::vector<Eigen::MatrixXd> Py; // Covariance for output (N matrices of ny x ny)
    std::vector<Eigen::MatrixXd> Px; // Covariance for state (N matrices of nx x nx)
    
    // Monte Carlo simulations
    std::vector<Eigen::MatrixXd> yMC; // MC realizations of output
    std::vector<Eigen::MatrixXd> xMC; // MC realizations of state
    int MC = 30; // Number of Monte Carlo samples (default 30 like MATLAB)
    
    // Labels
    std::vector<std::string> ylabel;
    std::vector<std::string> ulabel;
    std::vector<std::string> xlabel;
    std::string tlabel;
    
    // Dimensions [nx, nu, ny]
    Eigen::Vector3i nn = Eigen::Vector3i::Zero();
    
    // Metadata
    std::string name;
    std::string desc;
    std::vector<double> marker;
    std::string markerlabel;
    std::string moviefile;
    
    // User data (generic storage)
    void* userdata = nullptr;
    
    // Constructors
    
    /**
     * @brief Default constructor (empty signal)
     */
    Sig() {}
    
    /**
     * @brief Construct signal with y and sampling frequency (uniform sampling)
     * @param y_ Measurement data (N x ny)
     * @param fs_ Sampling frequency
     */
    Sig(const Eigen::MatrixXd& y_, double fs_) {
        int N = y_.rows();
        int ny = y_.cols();
        
        y = y_;
        fs = fs_;
        t = Eigen::VectorXd::LinSpaced(N, 0, (N-1) * fs_);
        u = Eigen::MatrixXd::Zero(N, 0);
        x = Eigen::MatrixXd::Zero(N, 0);
        
        nn << 0, 0, ny;
        init_labels();
    }
    
    /**
     * @brief Construct signal with y and time vector (non-uniform sampling)
     * @param y_ Measurement data (N x ny)
     * @param t_ Time vector (N x 1)
     */
    Sig(const Eigen::MatrixXd& y_, const Eigen::VectorXd& t_) {
        if (y_.rows() != t_.size()) {
            throw std::runtime_error("SIG: y and t must have same number of rows");
        }
        
        y = y_;
        t = t_;
        u = Eigen::MatrixXd::Zero(t_.size(), 0);
        x = Eigen::MatrixXd::Zero(t_.size(), 0);
        
        nn << 0, 0, y_.cols();
        check_uniform_sampling();
        init_labels();
    }
    
    /**
     * @brief Construct signal with y, t, and u
     * @param y_ Measurement data (N x ny)
     * @param t_ Time vector (N x 1) or scalar fs
     * @param u_ Input data (N x nu)
     */
    Sig(const Eigen::MatrixXd& y_, const Eigen::VectorXd& t_, const Eigen::MatrixXd& u_) {
        if (y_.rows() != t_.size() || y_.rows() != u_.rows()) {
            throw std::runtime_error("SIG: y, t, and u must have same number of rows");
        }
        
        y = y_;
        t = t_;
        u = u_;
        x = Eigen::MatrixXd::Zero(t_.size(), 0);
        
        nn << 0, u_.cols(), y_.cols();
        check_uniform_sampling();
        init_labels();
    }
    
    /**
     * @brief Construct signal with y, t, u, and x
     * @param y_ Measurement data (N x ny)
     * @param t_ Time vector (N x 1) or scalar fs
     * @param u_ Input data (N x nu)
     * @param x_ State data (N x nx)
     */
    Sig(const Eigen::MatrixXd& y_, const Eigen::VectorXd& t_, 
        const Eigen::MatrixXd& u_, const Eigen::MatrixXd& x_) {
        if (y_.rows() != t_.size() || y_.rows() != u_.rows() || y_.rows() != x_.rows()) {
            throw std::runtime_error("SIG: y, t, u, and x must have same number of rows");
        }
        
        y = y_;
        t = t_;
        u = u_;
        x = x_;
        
        nn << x_.cols(), u_.cols(), y_.cols();
        check_uniform_sampling();
        init_labels();
    }
    
    /**
     * @brief Get length of signal (number of time samples)
     */
    int length() const {
        return y.rows();
    }
    
    /**
     * @brief Get size along dimension (1=time, 2=y, 3=u, 4=x)
     */
    int size(int dim) const {
        switch(dim) {
            case 1: return y.rows();  // N (time samples)
            case 2: return y.cols();  // ny
            case 3: return u.cols();  // nu
            case 4: return x.cols();  // nx
            default: throw std::runtime_error("SIG::size: dim must be 1,2,3, or 4");
        }
    }
    
    /**
     * @brief Compute sample mean of y
     */
    Eigen::VectorXd mean_y() const {
        return y.colwise().mean().transpose();
    }
    
    /**
     * @brief Compute sample variance of y
     */
    Eigen::VectorXd var_y() const {
        Eigen::VectorXd m = mean_y();
        Eigen::MatrixXd centered = y.rowwise() - m.transpose();
        return (centered.array().square().colwise().sum() / (y.rows() - 1)).transpose();
    }
    
    /**
     * @brief Compute sample covariance matrix of y
     */
    Eigen::MatrixXd cov_y() const {
        Eigen::VectorXd m = mean_y();
        Eigen::MatrixXd centered = y.rowwise() - m.transpose();
        return (centered.transpose() * centered) / (y.rows() - 1);
    }
    
    /**
     * @brief Extract a subset of the signal by time range
     * @param t_start Start time
     * @param t_end End time
     */
    Sig extract(double t_start, double t_end) const {
        std::vector<int> indices;
        for (int i = 0; i < t.size(); ++i) {
            if (t(i) >= t_start && t(i) <= t_end) {
                indices.push_back(i);
            }
        }
        
        if (indices.empty()) {
            throw std::runtime_error("SIG::extract: no samples in specified range");
        }
        
        Sig subset;
        subset.y = Eigen::MatrixXd(indices.size(), y.cols());
        subset.t = Eigen::VectorXd(indices.size());
        subset.u = Eigen::MatrixXd(indices.size(), u.cols());
        subset.x = Eigen::MatrixXd(indices.size(), x.cols());
        
        for (size_t i = 0; i < indices.size(); ++i) {
            subset.y.row(i) = y.row(indices[i]);
            subset.t(i) = t(indices[i]);
            if (u.cols() > 0) subset.u.row(i) = u.row(indices[i]);
            if (x.cols() > 0) subset.x.row(i) = x.row(indices[i]);
        }
        
        subset.nn = nn;
        subset.ylabel = ylabel;
        subset.ulabel = ulabel;
        subset.xlabel = xlabel;
        subset.check_uniform_sampling();
        
        return subset;
    }
    
    /**
     * @brief Downsample signal by factor
     * @param factor Downsampling factor (keep every factor-th sample)
     */
    Sig downsample(int factor) const {
        if (factor < 1) {
            throw std::runtime_error("SIG::downsample: factor must be >= 1");
        }
        
        int N_new = (y.rows() + factor - 1) / factor;
        
        Sig down;
        down.y = Eigen::MatrixXd(N_new, y.cols());
        down.t = Eigen::VectorXd(N_new);
        down.u = Eigen::MatrixXd(N_new, u.cols());
        down.x = Eigen::MatrixXd(N_new, x.cols());
        
        for (int i = 0; i < N_new; ++i) {
            int idx = i * factor;
            if (idx < y.rows()) {
                down.y.row(i) = y.row(idx);
                down.t(i) = t(idx);
                if (u.cols() > 0) down.u.row(i) = u.row(idx);
                if (x.cols() > 0) down.x.row(i) = x.row(idx);
            }
        }
        
        down.nn = nn;
        down.ylabel = ylabel;
        down.ulabel = ulabel;
        down.xlabel = xlabel;
        down.fs = (!std::isnan(fs)) ? fs / factor : fs;
        
        return down;
    }
    
private:
    /**
     * @brief Initialize labels based on dimensions
     */
    void init_labels() {
        ylabel.clear();
        ulabel.clear();
        xlabel.clear();
        
        for (int i = 0; i < nn[2]; ++i) {
            ylabel.push_back("y" + std::to_string(i+1));
        }
        for (int i = 0; i < nn[1]; ++i) {
            ulabel.push_back("u" + std::to_string(i+1));
        }
        for (int i = 0; i < nn[0]; ++i) {
            xlabel.push_back("x" + std::to_string(i+1));
        }
    }
    
    /**
     * @brief Check if sampling is uniform and set fs accordingly
     */
    void check_uniform_sampling() {
        if (t.size() < 2) {
            fs = std::numeric_limits<double>::quiet_NaN();
            return;
        }
        
        Eigen::VectorXd dt = t.tail(t.size()-1) - t.head(t.size()-1);
        double max_diff = (dt.array() - dt.mean()).abs().maxCoeff();
        
        if (max_diff < 10 * std::numeric_limits<double>::epsilon()) {
            fs = dt.mean(); // Uniform sampling
        } else {
            fs = std::numeric_limits<double>::quiet_NaN(); // Non-uniform
        }
    }
};
