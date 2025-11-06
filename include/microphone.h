#ifndef MICROPHONE_H
#define MICROPHONE_H

#include <vector>

class Microphone {
public:
    Microphone(int id);
    void initialize();
    void startCapture();
    void stopCapture();
    std::vector<double> getAudioData() const;

private:
    int microphoneId;
    std::vector<double> audioData;
    bool isCapturing;
};

#endif // MICROPHONE_H