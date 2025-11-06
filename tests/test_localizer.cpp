#include "localizer.h"
#include <gtest/gtest.h>

class LocalizerTest : public ::testing::Test {
protected:
    Localizer localizer;

    void SetUp() override {
        // Initialize the localizer with microphone positions
        std::vector<std::pair<double, double>> microphone_positions = {
            {0.0, 0.0},   // Microphone 1
            {1.0, 0.0},   // Microphone 2
            {0.5, 1.0}    // Microphone 3
        };
        localizer.setMicrophonePositions(microphone_positions);
    }
};

TEST_F(LocalizerTest, TestTriangulateSingleSource) {
    // Simulated frequency data from microphones
    std::vector<double> magnitudes = {10.0, 15.0, 12.0};
    std::vector<double> phases = {0.0, 0.5, 0.3};

    auto result = localizer.triangulate(magnitudes, phases);

    // Check if the result is within expected bounds
    EXPECT_NEAR(result.first, 0.5, 0.1); // x-coordinate
    EXPECT_NEAR(result.second, 0.5, 0.1); // y-coordinate
}

TEST_F(LocalizerTest, TestTriangulateMultipleSources) {
    // Simulated frequency data for multiple sources
    std::vector<double> magnitudes = {5.0, 20.0, 10.0};
    std::vector<double> phases = {0.0, 1.0, 0.5};

    auto result = localizer.triangulate(magnitudes, phases);

    // Check if the result is within expected bounds
    EXPECT_NEAR(result.first, 0.75, 0.1); // x-coordinate
    EXPECT_NEAR(result.second, 0.25, 0.1); // y-coordinate
}

TEST_F(LocalizerTest, TestInvalidData) {
    // Simulated invalid frequency data
    std::vector<double> magnitudes = {0.0, 0.0, 0.0};
    std::vector<double> phases = {0.0, 0.0, 0.0};

    auto result = localizer.triangulate(magnitudes, phases);

    // Expect an invalid result or a specific behavior
    EXPECT_EQ(result.first, -1.0); // Indicating an error
    EXPECT_EQ(result.second, -1.0); // Indicating an error
}