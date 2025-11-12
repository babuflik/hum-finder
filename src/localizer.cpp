#include "Localizer.h"
#include <iostream>
// Inkludera en dummy-funktion för GCC-PHAT (MÅSTE BYTAS UT mot en riktig implementering)
#include "gcc_phat_stub.h" 

// Konstruktor
Localizer::Localizer(
    const std::array<Microphone, MIC_COUNT>& microphones,
    double initial_dt,
    double sound_speed
) : microphones_(microphones), C_(sound_speed), dt_(initial_dt) {

    // 1. Initiera Mikrofonpositioner
    for (int i = 0; i < MIC_COUNT; ++i) {
        auto pos = microphones[i].getPosition();
        mic_positions_[i] << pos[0], pos[1], pos[2];
    }

    // 2. Initiera EKF Tillstånd (x) och Kovarians (P)
    x_.setZero();
    // Startposition mitt i arrayen (1.5, 1.5, 1.5)
    x_(0) = 1.5; x_(1) = 1.5; x_(2) = 1.5;
    
    P_.setIdentity(); 
    // Position har hög osäkerhet, hastighet lägre
    P_.block<3, 3>(0, 0) *= 5.0; // Större osäkerhet i position
    P_.block<3, 3>(3, 3) *= 0.1; // Mindre osäkerhet i hastighet

    // 3. Initiera Process Noise (Q) och Measurement Noise (R)
    // Q: Vi antar konstant hastighet (vx, vy, vz) plus lite slumpmässigt brus
    Q_.setIdentity();
    Q_.block<3, 3>(0, 0) *= 0.01; // Litet brus på position
    Q_.block<3, 3>(3, 3) *= 0.05; // Mer brus på hastighet (för att tillåta långsam rörelse/drift)
    
    // R: Osäkerhet i TDOA-mätningen (bör vara baserad på samplingsfrekvens och brus)
    R_.setIdentity();
    R_ *= 1.0e-6; // Exempel: 1 microsekund^2 osäkerhet

    // 4. Initiera State Transition Matrix (F)
    // Detta är för en Constant Velocity (CV) modell: x(t+dt) = x(t) + v(t)*dt
    F_.setIdentity();
    F_.block<3, 3>(0, 3) = Eigen::Matrix3d::Identity() * dt_;
}

// Huvudfunktion för att köra ett EKF-steg
std::array<double, 3> Localizer::locateSource(
    const std::array<std::vector<double>, MIC_COUNT>& audioBuffers
) {
    // 1. Mätning: Beräkna TDOA-vektorn (z)
    Eigen::Matrix<double, MEASUREMENT_DIM, 1> z = calculateTDOA(audioBuffers);
    
    // 2. Prediktion
    predictionStage();
    
    // 3. Uppdatering
    updatingStage(z);
    
    // Returnera den filtrerade positionen [x, y, z]
    return {x_(0), x_(1), x_(2)};
}

// Implementerar EKF Prediktionssteget
void Localizer::predictionStage() {
    // 1. Prediktera Tillstånd: x_k|k-1 = F * x_k-1|k-1
    x_ = F_ * x_;

    // 2. Prediktera Kovarians: P_k|k-1 = F * P_k-1|k-1 * F^T + Q
    P_ = F_ * P_ * F_.transpose() + Q_;
}

// Beräknar den förväntade mätningen h(x)
// Detta är den icke-linjära modellen: TDOA = (Avstånd_till_Mi - Avstånd_till_Mj) / C
Eigen::Matrix<double, MEASUREMENT_DIM, 1> Localizer::calculateExpectedTDOA(
    const Eigen::Matrix<double, STATE_DIM, 1>& x_pred
) {
    // Positionen av ljudkällan baserat på det predikterade tillståndet
    Eigen::Vector3d source_pos = x_pred.block<3, 1>(0, 0); 
    Eigen::Matrix<double, MEASUREMENT_DIM, 1> h_x;

    // Funktion för att beräkna avstånd till en mikrofon Mi
    auto dist_to_mic = [&](int i) {
        return (mic_positions_[i] - source_pos).norm();
    };
    
    // Beräkna de 6 TDOA-värdena (i sekunder)
    h_x(0) = (dist_to_mic(0) - dist_to_mic(1)) / C_; // M1-M2
    h_x(1) = (dist_to_mic(0) - dist_to_mic(2)) / C_; // M1-M3
    h_x(2) = (dist_to_mic(0) - dist_to_mic(3)) / C_; // M1-M4
    h_x(3) = (dist_to_mic(1) - dist_to_mic(2)) / C_; // M2-M3
    h_x(4) = (dist_to_mic(1) - dist_to_mic(3)) / C_; // M2-M4
    h_x(5) = (dist_to_mic(2) - dist_to_mic(3)) / C_; // M3-M4

    return h_x;
}

// Beräknar H (Jacobian av h(x))
// H är derivatan av h(x) med avseende på tillståndsvektorn x: H = dh/dx
Eigen::Matrix<double, MEASUREMENT_DIM, STATE_DIM> Localizer::calculateJacobianH() {
    Eigen::Matrix<double, MEASUREMENT_DIM, STATE_DIM> H = Eigen::MatrixXd::Zero(MEASUREMENT_DIM, STATE_DIM);

    // Positionen av ljudkällan (x, y, z)
    Eigen::Vector3d source_pos = x_.block<3, 1>(0, 0); 

    // Hjälpfunktion för att beräkna derivatan av avståndet R_i
    // Derivatan av R_i med avseende på x, y, z är: (x_s - x_i) / R_i
    auto dR_dx = [&](int i) {
        // Positionen för mikrofon i
        const Eigen::Vector3d& m_i = mic_positions_[i];
        
        // Avstånd R_i
        double R_i = (source_pos - m_i).norm(); 
        
        // Derivatans delar dR/d(x,y,z)
        Eigen::RowVector3d derivative_vector;
        
        // dR_i / dx_s
        derivative_vector(0) = (source_pos.x() - m_i.x()) / R_i; 
        // dR_i / dy_s
        derivative_vector(1) = (source_pos.y() - m_i.y()) / R_i;
        // dR_i / dz_s
        derivative_vector(2) = (source_pos.z() - m_i.z()) / R_i;
        
        return derivative_vector;
    };

    // Indexmappning för de 6 TDOA-paren: (i, j)
    std::array<std::pair<int, int>, MEASUREMENT_DIM> mic_pairs = {{
        {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
    }};

    // Fyll Jacobianska matrisen H
    for (int k = 0; k < MEASUREMENT_DIM; ++k) {
        int i = mic_pairs[k].first;
        int j = mic_pairs[k].second;

        // h_k = (R_i - R_j) / C.
        // d(h_k)/d(x,y,z) = (1/C) * (dR_i/d(x,y,z) - dR_j/d(x,y,z))
        
        // H(k, 0:2) är derivatan med avseende på position (x, y, z)
        H.block<1, 3>(k, 0) = (1.0 / C_) * (dR_dx(i) - dR_dx(j));
        
        // H(k, 3:5) är derivatan med avseende på hastighet (vx, vy, vz), vilket är 0 i h(x)
        // Dessa är redan nollställda, så inget mer att göra här.
    }
    
    return H;
}

// Implementerar EKF Uppdateringssteget
void Localizer::updatingStage(const Eigen::Matrix<double, MEASUREMENT_DIM, 1>& z) {
    // 1. Beräkna den förväntade mätningen h(x)
    Eigen::Matrix<double, MEASUREMENT_DIM, 1> h_x = calculateExpectedTDOA(x_);
    
    // 2. Beräkna Jacobian H
    Eigen::Matrix<double, MEASUREMENT_DIM, STATE_DIM> H = calculateJacobianH();

    // 3. Beräkna Innovationskovarians S
    // S = H * P * H^T + R
    Eigen::Matrix<double, MEASUREMENT_DIM, MEASUREMENT_DIM> S = 
        H * P_ * H.transpose() + R_;

    // 4. Beräkna Kalman Gain K
    // K = P * H^T * S^-1
    Eigen::Matrix<double, STATE_DIM, MEASUREMENT_DIM> K = 
        P_ * H.transpose() * S.inverse();

    // 5. Uppdatera Tillstånd
    // x_k|k = x_k|k-1 + K * (z - h(x))
    Eigen::Matrix<double, MEASUREMENT_DIM, 1> innovation = z - h_x;
    x_ = x_ + K * innovation;

    // 6. Uppdatera Kovarians
    // P_k|k = (I - K * H) * P_k|k-1
    Eigen::Matrix<double, STATE_DIM, STATE_DIM> I = Eigen::MatrixXd::Identity(STATE_DIM, STATE_DIM);
    P_ = (I - K * H) * P_;
    
    // Säkerställ symmetri (p.g.a. flyttalsfel)
    P_ = 0.5 * (P_ + P_.transpose());
}

// -------------------------------------------------------------
// DUMMY IMPLEMENTATION AV TDOA - MÅSTE BYTAS UT
// -------------------------------------------------------------

// Detta är en platshållarfil som du MÅSTE BYTA UT mot en robust
// implementering av GCC-PHAT, t.ex. med ett FFT-bibliotek som FFTW.
#include <random>

Eigen::Matrix<double, MEASUREMENT_DIM, 1> Localizer::calculateTDOA(
    const std::array<std::vector<double>, MIC_COUNT>& audioBuffers
) {
    // Denna DUMMY returnerar bara nollor med lite brus.
    // Ljudlokaliseringsnoggrannheten är helt beroende av denna funktion!
    Eigen::Matrix<double, MEASUREMENT_DIM, 1> z;
    z.setZero();
    
    // Lägg till brus för att simulera en riktig mätning
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0.0, 1.0e-4); // Brus med std.avv 0.1 ms (mycket bra)

    for (int i = 0; i < MEASUREMENT_DIM; ++i) {
        z(i) += d(gen);
    }
    
    // Hårdkodad TDOA från en känd källa (för testning)
    // T.ex. Om källan är vid (1, 1, 1), skulle detta vara de sanna värdena
    // z(0) = TDOA(M1(0,0,0) - M2(3,0,0) vid (1,1,1))
    // Vi använder nollor som en neutral start.
    
    std::cout << "VARNING: Använder dummy TDOA. Byt ut mot GCC-PHAT-kod." << std::endl;

    return z;
}