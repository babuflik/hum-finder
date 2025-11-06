#include "localizer.h"
#include <cmath>
#include <stdexcept>

Localizer::Localizer(const std::array<Microphone, 3>& microphones) 
    : microphones_(microphones) {
}

std::array<double, 2> Localizer::locateSource(const std::vector<double>& magnitudes) {
    if (magnitudes.size() != 3) {
        throw std::invalid_argument("Need exactly 3 magnitude readings");
    }

    // Calculate intersections between each pair of circles
    auto pos1 = calculateIntersection(
        microphones_[0].getPosition(),
        microphones_[1].getPosition(),
        magnitudes[0],
        magnitudes[1]
    );

    auto pos2 = calculateIntersection(
        microphones_[1].getPosition(),
        microphones_[2].getPosition(),
        magnitudes[1],
        magnitudes[2]
    );

    // Return average position as estimated source location
    return std::array<double, 2>{(pos1[0] + pos2[0]) / 2.0, (pos1[1] + pos2[1]) / 2.0};
}

std::array<double, 2> Localizer::calculateIntersection(
    const std::array<double, 2>& pos1,
    const std::array<double, 2>& pos2,
    double r1,
    double r2) {
    
    double dx = pos2[0] - pos1[0];
    double dy = pos2[1] - pos1[1];
    double d = std::sqrt(dx * dx + dy * dy);

    if (d > r1 + r2 || d < std::abs(r1 - r2)) {
        throw std::runtime_error("No intersection exists");
    }

    double a = (r1 * r1 - r2 * r2 + d * d) / (2 * d);
    double h = std::sqrt(r1 * r1 - a * a);

    double x = pos1[0] + (dx * a) / d;
    double y = pos1[1] + (dy * a) / d;

    return std::array<double, 2>{x, y};
}