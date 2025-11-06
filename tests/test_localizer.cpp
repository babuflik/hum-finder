#include "localizer.h"
#include <gtest/gtest.h>
#include <array>
#include <vector>

class LocalizerTest : public ::testing::Test {
protected:
    LocalizerTest() : 
        microphones{
            Microphone(0.0, 0.0),      // First microphone at origin
            Microphone(1.0, 0.0),      // Second microphone 1m to the right
            Microphone(0.5, 0.866)     // Third microphone to form triangle
        },
        localizer(microphones) {}

    void SetUp() override {}

    std::array<Microphone, 3> microphones;
    Localizer localizer;
};

TEST_F(LocalizerTest, TestTriangulateSingleSource) {
    // Simulate magnitudes that would come from a source at (0.5, 0.5)
    std::vector<double> magnitudes = {0.7071, 0.7071, 0.7071}; // Example values
    
    auto result = localizer.locateSource(magnitudes);
    
    // Check if the result is close to expected position
    EXPECT_NEAR(result[0], 0.5, 0.1); // x coordinate
    EXPECT_NEAR(result[1], 0.5, 0.1); // y coordinate
}