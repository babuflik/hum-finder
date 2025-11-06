#ifndef LOCALIZER_H
#define LOCALIZER_H

#include "microphone.h"
#include <array>
#include <vector>

class Localizer {
public:
    explicit Localizer(const std::array<Microphone, 3>& microphones);
    Localizer() = delete;  // Prevent default construction
    
    // Returns the estimated source position [x, y]
    std::array<double, 2> locateSource(const std::vector<double>& magnitudes);

private:
    std::array<Microphone, 3> microphones_;
    
    // Helper method to calculate intersection points
    std::array<double, 2> calculateIntersection(
        const std::array<double, 2>& pos1,
        const std::array<double, 2>& pos2,
        double r1,
        double r2);
};

#endif // LOCALIZER_H