#include "microphone.h"

Microphone::Microphone(double x, double y) 
    : position_{x, y} {
    samples_.reserve(BUFFER_SIZE);
}

std::array<double, 2> Microphone::getPosition() const {
    return position_;
}

void Microphone::addSample(double sample) {
    if (samples_.size() >= BUFFER_SIZE) {
        samples_.erase(samples_.begin());
    }
    samples_.push_back(sample);
}

std::vector<double> Microphone::getSamples() const {
    return samples_;
}