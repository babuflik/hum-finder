#ifndef DFT_H
#define DFT_H

#include <vector>
#include <complex>

class DFT {
public:
    DFT();
    std::vector<std::complex<double>> compute(const std::vector<double>& inputSignal);
    std::vector<double> getMagnitude(const std::vector<std::complex<double>>& dftResult);
    std::vector<double> getFrequencyBins(int sampleRate, int numSamples);

private:
    void normalize(std::vector<std::complex<double>>& dftResult);
};

#endif // DFT_H