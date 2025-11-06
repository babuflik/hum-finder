#ifndef LOCALIZER_H
#define LOCALIZER_H

#include <vector>
#include <complex>

class Localizer {
public:
    Localizer(const std::vector<std::pair<double, double>>& microphonePositions);
    std::pair<double, double> triangulate(const std::vector<double>& magnitudes, const std::vector<double>& phases);

private:
    std::vector<std::pair<double, double>> microphonePositions;
    double calculateDistance(const std::pair<double, double>& mic1, const std::pair<double, double>& mic2);
};

#endif // LOCALIZER_H