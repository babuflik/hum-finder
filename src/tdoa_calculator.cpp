#include "tdoa_calculator.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <iostream>

// Standard C++ header för PI (saknas ibland i vissa kompilatorer)
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

TdoaCalculator::TdoaCalculator() {
    // Kontrollera att FFT_SIZE är en potens av 2 (Krävs för den rekursiva FFT:n)
    if ((FFT_SIZE > 0) && ((FFT_SIZE & (FFT_SIZE - 1)) != 0)) {
        std::cerr << "VARNING: FFT_SIZE bör vara en potens av 2 för den aktuella FFT-implementationen!" << std::endl;
    }
}

// -----------------------------------------------------------------
// COOLEY-TUKEY REKURSIV FFT/IFFT IMPLEMENTATION
// -----------------------------------------------------------------

// Privata rekursiva funktionen för FFT/IFFT
void TdoaCalculator::recursiveFft(ComplexVector& a, bool inverse) {
    int n = a.size();
    if (n <= 1) return;

    // Dela upp i jämna och udda delar
    ComplexVector a0(n / 2), a1(n / 2);
    for (int i = 0; i < n / 2; i++) {
        a0[i] = a[2 * i];
        a1[i] = a[2 * i + 1];
    }

    // Rekursiva anrop
    recursiveFft(a0, inverse);
    recursiveFft(a1, inverse);

    // Kombinera resultaten
    double angle_sign = inverse ? 2 * M_PI : -2 * M_PI;
    std::complex<double> w(1), wn = std::exp(std::complex<double>(0, angle_sign / n));

    for (int k = 0; k < n / 2; k++) {
        a[k] = a0[k] + w * a1[k];
        a[k + n / 2] = a0[k] - w * a1[k];
        w *= wn;
    }
}

// Utför FFT
ComplexVector TdoaCalculator::fft(const RealVector& signal) {
    if (signal.size() != FFT_SIZE) {
        std::cerr << "FEL: Signalstorlek (" << signal.size() << ") matchar inte FFT_SIZE (" << FFT_SIZE << ")." << std::endl;
        return ComplexVector(FFT_SIZE, 0.0);
    }
    
    ComplexVector result(signal.begin(), signal.end());
    recursiveFft(result, false);
    return result;
}

// Utför IFFT
RealVector TdoaCalculator::ifft(const ComplexVector& spectrum) {
    if (spectrum.size() != FFT_SIZE) {
        std::cerr << "FEL: Spektrumstorlek (" << spectrum.size() << ") matchar inte FFT_SIZE (" << FFT_SIZE << ")." << std::endl;
        return RealVector(FFT_SIZE, 0.0);
    }

    ComplexVector complex_result = spectrum;
    
    // 1. Invers FFT (använd omvänd tecken)
    recursiveFft(complex_result, true);
    
    // 2. Normalisering
    RealVector real_result(FFT_SIZE);
    for (size_t i = 0; i < FFT_SIZE; ++i) {
        // Dela med N och ta realdelen (FFT/IFFT skalan)
        real_result[i] = complex_result[i].real() / FFT_SIZE;
    }

    return real_result;
}


// -----------------------------------------------------------------
// BANDPASS FILTER IMPLEMENTATION (Frequency Domain)
// -----------------------------------------------------------------

/*
 * Implementerar ett Bandpass-filter i frekvensdomänen.
 * Isolerar den önskade frekvensen 100 Hz till 300 Hz.
 * * F_res = SAMPLE_RATE / FFT_SIZE = 44100 / 1024 = 43.066 Hz/bin
 * * Cutoff låg: 100 Hz -> Bin 3 (129.19 Hz)
 * Cutoff hög: 300 Hz -> Bin 7 (301.46 Hz)
 * * Range: [Bin 3, Bin 7] och dess speglade motsvarighet [N - 7, N - 3].
 */
RealVector TdoaCalculator::bandpassFilter(const RealVector& signal) {
    if (signal.size() != FFT_SIZE) {
        // Bör inte hända om anropet kommer från calculateGccPhat, men en säkerhetskontroll
        return signal;
    }

    // 1. FFT
    ComplexVector spectrum = fft(signal);

    // Frekvensgränser baserat på 43.066 Hz/bin
    const int K_LOW = 3;
    const int K_HIGH = 7;
    const size_t N = FFT_SIZE;
    
    // 2. Nollställ alla frekvenser utanför bandet
    for (size_t k = 0; k < N; ++k) {
        // Nollställ för k < K_LOW (lågpass-del)
        // Nollställ för k > K_HIGH OCH k < N - K_HIGH (området mellan positiva och negativa spektrumet)
        // Nollställ för k > N - K_LOW (högpass-del i det negativa spektrumet)
        
        if ((k < K_LOW) || (k > K_HIGH && k < N - K_HIGH) || (k > N - K_LOW)) {
            // Nollställ om utanför bandet [K_LOW, K_HIGH] U [N - K_HIGH, N - K_LOW]
            spectrum[k] = 0.0;
        }
    }

    // 3. IFFT
    return ifft(spectrum);
}


// Beräknar TDOA i samples för ett par
double TdoaCalculator::calculateGccPhat(const RealVector& signal1, const RealVector& signal2) {
    
    // Steg 1: Filtrering 
    RealVector filtered1 = bandpassFilter(signal1);
    RealVector filtered2 = bandpassFilter(signal2);

    // Steg 2: FFT
    ComplexVector spectrum1 = fft(filtered1);
    ComplexVector spectrum2 = fft(filtered2);

    // Steg 3: Beräkna PHAT-vikten och korskorrelationsspektrumet (R_ij^PHAT)
    ComplexVector cross_spectrum_phat;
    cross_spectrum_phat.reserve(FFT_SIZE);

    for (size_t k = 0; k < FFT_SIZE; ++k) {
        // Y_i[k] * Y_j^*[k]
        std::complex<double> product = spectrum1[k] * std::conj(spectrum2[k]);
        
        double magnitude = std::abs(product);
        if (magnitude > 1e-10) {
            // R_ij^PHAT[k] = (Y_i[k] * Y_j^*[k]) / |Y_i[k] * Y_j^*[k]|
            cross_spectrum_phat.push_back(product / magnitude);
        } else {
            cross_spectrum_phat.push_back(0.0);
        }
    }

    // Steg 4: IFFT
    RealVector correlation = ifft(cross_spectrum_phat);

    // Steg 5: Toppdetektering (Hitta indexet med högst magnitud)
    double max_corr = 0.0;
    int max_index = 0;

    // Vi letar efter toppen i korrelationsvektorn
    for (size_t n = 0; n < FFT_SIZE; ++n) {
        if (correlation[n] > max_corr) {
            max_corr = correlation[n];
            max_index = n;
        }
    }

    // Hantera cirkulär skiftning (wrap-around) av IFFT-resultatet
    // Om indexet är större än N/2, betyder det en negativ fördröjning
    int shift_samples = max_index;
    if (max_index > FFT_SIZE / 2) {
        shift_samples = max_index - (int)FFT_SIZE;
    }

    // Konvertera samples till tid (sekunder)
    return (double)shift_samples / SAMPLE_RATE;
}

// Huvudmetod
TdoaVector TdoaCalculator::calculate(
    const std::array<RealVector, MIC_COUNT_TDOA>& audioBuffers
) {
    TdoaVector result;

    // Indexmappning för de 6 TDOA-paren: (i, j)
    std::array<std::pair<int, int>, MEASUREMENT_DIM_TDOA> mic_pairs = {{
        {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
    }};

    for (int k = 0; k < MEASUREMENT_DIM_TDOA; ++k) {
        int i = mic_pairs[k].first;
        int j = mic_pairs[k].second;

        if (audioBuffers[i].size() != FFT_SIZE || audioBuffers[j].size() != FFT_SIZE) {
            std::cerr << "Fel: Bufferstorlek (" << audioBuffers[i].size() << ") matchar inte FFT_SIZE (" << FFT_SIZE << ")." << std::endl;
            // Returnerar en nolla, men detta borde hanteras bättre i en riktig applikation
            result(k) = 0.0; 
            continue;
        }
        
        // Beräkna TDOA för paret i och j
        result(k) = calculateGccPhat(audioBuffers[i], audioBuffers[j]);
    }
    
    return result;
}