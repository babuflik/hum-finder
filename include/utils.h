#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <complex>

// Constants
const double PI = 3.14159265358979323846;

// Utility functions
std::vector<std::complex<double>> performDFT(const std::vector<double>& input);
void logMessage(const std::string& message);
double calculateMagnitude(const std::complex<double>& complexNumber);
std::vector<double> calculatePhase(const std::vector<std::complex<double>>& frequencyData);

#endif // UTILS_H