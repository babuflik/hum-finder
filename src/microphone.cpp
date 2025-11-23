#include "microphone.h"
#include <algorithm> // För std::fill

// Konstruktor
Microphone::Microphone(double x, double y, double z) 
    : position_{x, y, z} {
    
    // Allokera minne för hela bufferten en gång vid start
    samples_.resize(BUFFER_SIZE);
    
    // Nollställ bufferten
    std::fill(samples_.begin(), samples_.end(), 0.0);
}

// Position getter
std::array<double, 3> Microphone::getPosition() const {
    return position_;
}

void Microphone::setPosition(const std::array<double,3>& pos) {
    position_ = pos;
}

// Lägger till ett nytt ljudprov och hanterar ringbufferten
void Microphone::addSample(double sample) {
    // Lägg in provet vid det aktuella skrivindexet
    samples_[write_index_] = sample;
    
    // Öka skrivindexet och rulla över till 0 om vi når slutet
    write_index_ = (write_index_ + 1) % BUFFER_SIZE;
    
    // Räkna hur många prover som har fyllts. Ökar tills den når maxstorleken.
    if (samples_filled_count_ < BUFFER_SIZE) {
        samples_filled_count_++;
    }
}

// Returnerar hela den aktuella bufferten
std::vector<double> Microphone::getSamples() const {
    
    // Om bufferten inte är full returnerar vi bara de insamlade proverna
    if (samples_filled_count_ < BUFFER_SIZE) {
        return samples_;
    }

    // --- Ringbuffert läslogik ---
    
    // Om bufferten är full: De äldsta proverna ligger vid write_index_.
    // Vi måste läsa ut bufferten i RÄTT TIDSORDNING: 
    // [Äldsta prover] följt av [Nyaste prover].
    
    std::vector<double> ordered_samples;
    ordered_samples.reserve(BUFFER_SIZE);
    
    // 1. Kopiera de äldsta proverna (från write_index_ till slutet)
    if (write_index_ < BUFFER_SIZE) {
        ordered_samples.insert(
            ordered_samples.end(),
            samples_.begin() + write_index_,
            samples_.end()
        );
    }
    
    // 2. Kopiera de nyaste proverna (från början till write_index_)
    if (write_index_ > 0) {
        ordered_samples.insert(
            ordered_samples.end(),
            samples_.begin(),
            samples_.begin() + write_index_
        );
    }

    // Observera: Detta returnerar en Kopia. För realtidsoptimering kan man överväga 
    // att returnera en referens eller en specialiserad iterator.
    return ordered_samples;
}