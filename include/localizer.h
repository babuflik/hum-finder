#ifndef LOCALIZER_H
#define LOCALIZER_H

#include "microphone.h"
#include <Eigen/Dense> 
#include <array>
#include <vector>
#include <cmath>

// Definiera tillståndsdimensionen (3D-position + 3D-hastighet)
static constexpr int STATE_DIM = 6; 
// Antal mikrofoner
static constexpr int MIC_COUNT = 4;
// Antal TDOA-mätningar (6 par: M1-M2, M1-M3, M1-M4, M2-M3, M2-M4, M3-M4)
static constexpr int MEASUREMENT_DIM = 6;

class Localizer {
public:
    explicit Localizer(
        const std::array<Microphone, MIC_COUNT>& microphones,
        double initial_dt,
        double sound_speed = 343.0 // Ljudhastighet (m/s)
    );
    Localizer() = delete; 

    // Huvudfunktion för att köra ett EKF-steg
    std::array<double, 3> locateSource(
        const std::array<std::vector<double>, MIC_COUNT>& audioBuffers
    );

private:
    // Mikrofon-data
    std::array<Microphone, MIC_COUNT> microphones_;
    std::array<Eigen::Vector3d, MIC_COUNT> mic_positions_;

    // Konstant Ljudhastighet
    const double C_; 
    // Tidssteg mellan uppdateringar (sampling interval)
    double dt_; 

    // --- EKF Tillstånd och Matriser ---

    // x (State Vector): [x, y, z, vx, vy, vz]
    Eigen::Matrix<double, STATE_DIM, 1> x_;

    // P (Covariance Matrix): Osäkerhet i tillståndet
    Eigen::Matrix<double, STATE_DIM, STATE_DIM> P_;

    // Q (Process Noise): Osäkerhet i rörelsemodellen (konstant hastighet)
    Eigen::Matrix<double, STATE_DIM, STATE_DIM> Q_;
    
    // R (Measurement Noise): Osäkerhet i TDOA-mätningarna
    Eigen::Matrix<double, MEASUREMENT_DIM, MEASUREMENT_DIM> R_;

    // F (State Transition Matrix): Linjäriserad modell för tillståndsövergång
    Eigen::Matrix<double, STATE_DIM, STATE_DIM> F_;

    // --- EKF Steg och Hjälpfunktioner ---

    // Beräknar TDOA-mätningarna (z) från råa ljuddata
    Eigen::Matrix<double, MEASUREMENT_DIM, 1> calculateTDOA(
        const std::array<std::vector<double>, MIC_COUNT>& audioBuffers
    );
    
    // Implementerar EKF Prediktionssteget: x_k|k-1 = F * x_k-1|k-1
    void predictionStage();

    // Beräknar den förväntade mätningen (h(x)) baserat på det predikterade tillståndet
    Eigen::Matrix<double, MEASUREMENT_DIM, 1> calculateExpectedTDOA(
        const Eigen::Matrix<double, STATE_DIM, 1>& x_pred
    );

    // Beräknar H (Jacobian av h(x))
    Eigen::Matrix<double, MEASUREMENT_DIM, STATE_DIM> calculateJacobianH();

    // Implementerar EKF Uppdateringssteget (Korrigerar med TDOA-mätningen)
    void updatingStage(const Eigen::Matrix<double, MEASUREMENT_DIM, 1>& z);
    
    // Enkel hjälpfunktion för att beräkna avstånd
    double distance(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2) const {
        return (p1 - p2).norm();
    }
};

#endif // LOCALIZER_H