//
//  PartialsStorage.hpp
//  Real-Time-Fast-Partial-Tracking
//
//  Created by Colas Benjamin on 27/06/2025.
//

// PartialsStorage.hpp
#pragma once
#include <vector>
#include <array>
#include <stdexcept>

namespace PartialTracking {

template<typename T>
struct PerFramePartials {
    // PartialParameter = [omega , amp, phase, damp, trackID]
    std::vector<std::array<T,4>> buffer;
    std::size_t size = 0;

    PerFramePartials(std::size_t maxPartials)
      : buffer(maxPartials), size(0) {}

    // Redimensionnement
    void resize(std::size_t n) {
        if (n > buffer.size())
            throw std::runtime_error("Dépassement capacité partiels");
        size = n;
    }

    std::array<T,4>&       operator[](std::size_t i)       { return buffer[i]; }
    const std::array<T,4>& operator[](std::size_t i) const { return buffer[i]; }
};

template<typename T>
struct PartialsStorage {
    std::vector<PerFramePartials<T>> perFrame;
    PartialsStorage(std::size_t M, std::size_t maxPartials)
    {
        perFrame.reserve(M);
        for (std::size_t m = 0; m < M; ++m)
            perFrame.emplace_back(maxPartials);
    }

    PerFramePartials<T>& at(std::size_t m)       { return perFrame[m]; }
    const PerFramePartials<T>& at(std::size_t m) const { return perFrame[m]; }
};

} // namespace PartialTracking
