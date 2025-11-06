#pragma once
#include <vector>
#include <complex>
#include <fftw3.h>

class DFT {
public:
    explicit DFT(int size);
    ~DFT();

    // Compute DFT for given input
    std::vector<std::complex<double>> compute(const std::vector<double>& input);
    
    // Get magnitude spectrum from DFT result
    static std::vector<double> getMagnitude(const std::vector<std::complex<double>>& dftResult);

private:
    int size_;
    fftw_plan plan_;
    double* in_;
    fftw_complex* out_;
};