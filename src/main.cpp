// src/main.cpp
#include <iostream>
#include <array>
#include <vector>
#include "microphone.h"
#include "localizer.h"

int main() {
    // OBS: MIC_COUNT är 4 i din header - skapa 4 mikrofoner
    std::array<Microphone, MIC_COUNT> microphones = {
        Microphone(0.0, 0.0, 0.0),      // M1
        Microphone(0.20, 0.0, 0.0),     // M2 (20 cm till höger)
        Microphone(0.10, 0.1732, 0.0),  // M3 (triangel)
        Microphone(-0.20, 0.0, 0.0)     // M4 (ytterligare mikrofon)
    };

    // Skapa Localizer (använd samma dt som du vill, t.ex. 0.01s)
    Localizer localizer(microphones, 0.01);

    // Bygg upp test-audio-buffers: varje buffer måste vara Microphone::BUFFER_SIZE långa
    std::array<std::vector<double>, MIC_COUNT> audioBuffers;
    for (size_t i = 0; i < MIC_COUNT; ++i) {
        audioBuffers[i].assign(Microphone::BUFFER_SIZE, 0.0); // fyll med nollor (dummy)
        // För testning kan du lägga in en impuls i mitten:
        audioBuffers[i][Microphone::BUFFER_SIZE/2] = (i==0 ? 1.0 : 0.9); // liten offset
    }

    // Kör lokaliseraren (använder din dummy-TDOA om du inte bytt ut den)
    auto sourceLocation = localizer.locateSource(audioBuffers);

    // Skriv ut x, y, z
    std::cout << "Source located at: ("
              << sourceLocation[0] << ", "
              << sourceLocation[1] << ", "
              << sourceLocation[2] << ")" << std::endl;

    return 0;
}
