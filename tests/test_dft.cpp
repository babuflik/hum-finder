#include <gtest/gtest.h>
#include "dft.h"
#include <vector>
#include <complex>

class DFTTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Set up test data
        input_size = 8;
        dft = new DFT(input_size);
    }

    void TearDown() override {
        delete dft;
    }

    DFT* dft;
    size_t input_size;
};

TEST_F(DFTTest, ComputeSineWave) {
    std::vector<double> input(input_size);
    // Generate simple sine wave
    for (size_t i = 0; i < input_size; i++) {
        input[i] = std::sin(2.0 * M_PI * i / input_size);
    }
    
    auto result = dft->compute(input);
    
    // For a pure sine wave, we expect peak at frequency 1
    std::vector<double> magnitudes = dft->getMagnitude(result);
    size_t peak_idx = 0;
    double peak_value = 0;
    
    for (size_t i = 0; i < magnitudes.size(); i++) {
        if (magnitudes[i] > peak_value) {
            peak_value = magnitudes[i];
            peak_idx = i;
        }
    }
    
    EXPECT_EQ(peak_idx, 1);
}