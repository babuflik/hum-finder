#ifndef MICROPHONE_H
#define MICROPHONE_H

#include <vector>
#include <array>
#include <cstddef>

class Microphone {
public:
    Microphone(double x, double y);
    
    // Position getter
    std::array<double, 2> getPosition() const;
    
    // Sample handling
    void addSample(double sample);
    std::vector<double> getSamples() const;

private:
    std::array<double, 2> position_;
    std::vector<double> samples_;
    static constexpr size_t BUFFER_SIZE = 1024;
};

#endif // MICROPHONE_H