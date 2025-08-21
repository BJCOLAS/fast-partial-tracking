//
//  WavExporter.hpp
//  Real-Time-Fast-Partial-Tracking
//
//  Created by Colas Benjamin on 07/08/2025.
//

//
//  WavExporter.hpp
//  Real-Time-Fast-Partial-Tracking
//
//  Created by Colas Benjamin on 07/08/2025.
//

#ifndef WAVEXPORTER_HPP
#define WAVEXPORTER_HPP

#include <vector>
#include <string>
#include <fstream>
#include <cstdint>
#include <algorithm>

namespace PartialTracking {

inline void exportWAV(const std::vector<double>& signal, int sampleRate, const std::string& filename) {
    // === Normalisation du signal dans [-1.0, 1.0] ===
    std::vector<int16_t> pcmData(signal.size());
    double maxVal = 1e-9;
    for (double s : signal) maxVal = std::max(maxVal, std::abs(s));
    double gain = (maxVal > 1.0) ? (1.0 / maxVal) : 1.0;

    for (size_t i = 0; i < signal.size(); ++i) {
        double scaled = signal[i] * gain;
        scaled = std::clamp(scaled, -1.0, 1.0);
        pcmData[i] = static_cast<int16_t>(scaled * 32767.0);
    }

    // === Paramètres WAV ===
    int numChannels = 1;
    int bitsPerSample = 16;
    int byteRate = sampleRate * numChannels * bitsPerSample / 8;
    int blockAlign = numChannels * bitsPerSample / 8;
    int dataSize = static_cast<int>(pcmData.size()) * blockAlign;
    int chunkSize = 36 + dataSize;

    // === Écriture du fichier WAV ===
    std::ofstream out(filename, std::ios::binary);
    if (!out) {
        std::cerr << "Erreur : impossible d’ouvrir le fichier WAV pour écriture.\n";
        return;
    }

    // RIFF header
    out.write("RIFF", 4);
    out.write(reinterpret_cast<const char*>(&chunkSize), 4);
    out.write("WAVE", 4);

    // fmt chunk
    out.write("fmt ", 4);
    int subchunk1Size = 16;
    int audioFormat = 1; // PCM
    out.write(reinterpret_cast<const char*>(&subchunk1Size), 4);
    out.write(reinterpret_cast<const char*>(&audioFormat), 2);
    out.write(reinterpret_cast<const char*>(&numChannels), 2);
    out.write(reinterpret_cast<const char*>(&sampleRate), 4);
    out.write(reinterpret_cast<const char*>(&byteRate), 4);
    out.write(reinterpret_cast<const char*>(&blockAlign), 2);
    out.write(reinterpret_cast<const char*>(&bitsPerSample), 2);

    // data chunk
    out.write("data", 4);
    out.write(reinterpret_cast<const char*>(&dataSize), 4);
    out.write(reinterpret_cast<const char*>(pcmData.data()), dataSize);

    out.close();
}

} // namespace PartialTracking

#endif // WAVEXPORTER_HPP
