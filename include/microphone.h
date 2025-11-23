#ifndef MICROPHONE_H
#define MICROPHONE_H

#include <vector>
#include <array>
#include <cstddef>
#include <iostream>

class Microphone {
public:
    // Storlek på datablocket som mikrofonen ska buffra.
    // T.ex. 4096 prover. Detta bestämmer storleken på TDOA-fönstret.
    static constexpr size_t BUFFER_SIZE = 1024; 

    // Konstruktor: Position x, y, z
    Microphone(double x, double y, double z);

    // Position getter
    std::array<double, 3> getPosition() const;

    void setPosition(const std::array<double,3>& pos);
    
    // Lägger till ett nytt ljudprov till bufferten (lägger till i slutet och skriver över äldst om full)
    void addSample(double sample);
    
    // Returnerar hela den nuvarande bufferten (BUFFERT_SIZE prover).
    // OBS: Returnerar en KOPIA av bufferten för trådsäkerhet och renhet.
    std::vector<double> getSamples() const;
    
    // Returnerar om bufferten är full och redo att behandlas
    bool isBufferFull() const { return samples_filled_count_ == BUFFER_SIZE; }

private:
    // Fysisk position i 3D-rymden
    std::array<double, 3> position_; 

    // Ringbuffert för ljudprover. std::vector används, men logiken behandlar den som fix storlek.
    std::vector<double> samples_; 
    
    // Index där nästa prov ska läggas in
    size_t write_index_ = 0; 
    
    // Antal prover som faktiskt finns i bufferten (för att veta när den är "full")
    size_t samples_filled_count_ = 0; 
};

#endif // MICROPHONE_H