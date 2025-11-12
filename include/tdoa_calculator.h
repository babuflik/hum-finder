#ifndef TDOA_CALCULATOR_H
#define TDOA_CALCULATOR_H

#include <vector>
#include <complex>
#include <array>
#include <Eigen/Dense>
#include <fftw3.h>
#include "microphone.h"

// Definiera TDOA-konstanter
static constexpr int MIC_COUNT_TDOA = 4;
static constexpr int MEASUREMENT_DIM_TDOA = 6;
static constexpr double SAMPLE_RATE = 44100.0; // Standard samplingsfrekvens (Hz)
static constexpr size_t FFT_SIZE = Microphone::BUFFER_SIZE; // Använder samma storlek som Microphone bufferten

// Typdefinitioner för att förenkla koden
using ComplexVector = std::vector<std::complex<double>>;
using RealVector = std::vector<double>;
using TdoaVector = Eigen::Matrix<double, MEASUREMENT_DIM_TDOA, 1>;

class TdoaCalculator {
public:
    TdoaCalculator();

    // Huvudmetod för att beräkna TDOA-vektorn från råa buffertar
    TdoaVector calculate(
        const std::array<RealVector, MIC_COUNT_TDOA>& audioBuffers
    );

private:
    // --- Signalbehandlingssteg (MÅSTE IMPLEMENTERAS MED FFT-BIBLIOTEK) ---
    
    void recursiveFft(ComplexVector& a, bool inverse);

    // Använder ett Butterworth-filter för att isolera 100-300 Hz
    RealVector bandpassFilter(const RealVector& signal);

    // Utför FFT. I ett verkligt system använder man t.ex. FFTW här.
    ComplexVector fft(const RealVector& signal);
    
    // Utför IFFT. I ett verkligt system använder man t.ex. FFTW här.
    RealVector ifft(const ComplexVector& spectrum);
    
    // Beräknar TDOA i samples för ett par
    double calculateGccPhat(const RealVector& signal1, const RealVector& signal2);
};

#endif // TDOA_CALCULATOR_H