#include <gtest/gtest.h>
#include "../include/dft.h"

TEST(DFTTest, BasicFunctionality) {
    DFT dft;

    // Test input signal: a simple sine wave
    const int N = 8;
    double input[N] = {0.0, 0.707, 1.0, 0.707, 0.0, -0.707, -1.0, -0.707};
    std::vector<std::complex<double>> output = dft.compute(input, N);

    // Expected output (magnitude) for the sine wave
    double expected_magnitude[N] = {0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    for (int i = 0; i < N; ++i) {
        EXPECT_NEAR(std::abs(output[i]), expected_magnitude[i], 1e-5);
    }
}

TEST(DFTTest, ZeroInput) {
    DFT dft;

    const int N = 8;
    double input[N] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    std::vector<std::complex<double>> output = dft.compute(input, N);

    for (int i = 0; i < N; ++i) {
        EXPECT_NEAR(output[i].real(), 0.0, 1e-5);
        EXPECT_NEAR(output[i].imag(), 0.0, 1e-5);
    }
}

TEST(DFTTest, ConstantInput) {
    DFT dft;

    const int N = 8;
    double input[N] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    std::vector<std::complex<double>> output = dft.compute(input, N);

    EXPECT_NEAR(output[0].real(), 8.0, 1e-5);
    for (int i = 1; i < N; ++i) {
        EXPECT_NEAR(output[i].real(), 0.0, 1e-5);
        EXPECT_NEAR(output[i].imag(), 0.0, 1e-5);
    }
}