#include "dft.h"
#include <cmath>
#include <vector>
#include <complex>

DFT::DFT(int size) : size(size) {
    data.resize(size);
    output.resize(size);
}

void DFT::setInput(const std::vector<double>& input) {
    if (input.size() != size) {
        throw std::invalid_argument("Input size must match DFT size.");
    }
    data = input;
}

void DFT::compute() {
    for (int k = 0; k < size; ++k) {
        std::complex<double> sum(0.0, 0.0);
        for (int n = 0; n < size; ++n) {
            double angle = 2.0 * M_PI * k * n / size;
            std::complex<double> w(cos(angle), -sin(angle));
            sum += data[n] * w;
        }
        output[k] = sum;
    }
}

std::vector<std::complex<double>> DFT::getOutput() const {
    return output;
}

std::vector<double> DFT::getMagnitude() const {
    std::vector<double> magnitudes(size);
    for (int i = 0; i < size; ++i) {
        magnitudes[i] = std::abs(output[i]);
    }
    return magnitudes;
}