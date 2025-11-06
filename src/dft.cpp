#include "dft.h"
#include <cmath>

DFT::DFT(int size) : size_(size) {
    in_ = fftw_alloc_real(size_);
    out_ = fftw_alloc_complex(size_ / 2 + 1);
    plan_ = fftw_plan_dft_r2c_1d(size_, in_, out_, FFTW_ESTIMATE);
}

DFT::~DFT() {
    fftw_destroy_plan(plan_);
    fftw_free(in_);
    fftw_free(out_);
}

std::vector<std::complex<double>> DFT::compute(const std::vector<double>& input) {
    if (input.size() != size_t(size_)) {
        throw std::invalid_argument("Input size must match DFT size");
    }

    // Copy input data
    std::copy(input.begin(), input.end(), in_);

    // Execute FFT
    fftw_execute(plan_);

    // Copy results to output vector
    std::vector<std::complex<double>> result(size_ / 2 + 1);
    for (int i = 0; i < size_ / 2 + 1; i++) {
        result[i] = std::complex<double>(out_[i][0], out_[i][1]);
    }

    return result;
}

std::vector<double> DFT::getMagnitude(const std::vector<std::complex<double>>& dftResult) {
    std::vector<double> magnitude(dftResult.size());
    for (size_t i = 0; i < dftResult.size(); i++) {
        magnitude[i] = std::abs(dftResult[i]);
    }
    return magnitude;
}