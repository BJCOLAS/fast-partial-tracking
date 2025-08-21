//
//  AdditiveSynthesis.hpp
//  Real-Time-Fast-Partial-Tracking
//
//  Created by Colas Benjamin on 06/08/2025.
//

#pragma once
#include <vector>
#include <cmath>

namespace PartialTracking {

template<typename T>
class AdditiveSynthesizer {
public:
  AdditiveSynthesizer(size_t hopSize)
    : hopSize_(hopSize)
  {}

  /**
   * Calcule la trame synthétisée à partir des partiels
   *
   * @param logAmps  : [numTracks][hopSize] log-amplitudes
   * @param phases   : [numTracks][hopSize] phases
   * @param outFrame : [hopSize] sortie temporelle synthétisée
   */
  void synthesizeFrame(
    const std::vector<std::vector<T>>& logAmps,
    const std::vector<std::vector<T>>& phases,
    std::vector<T>&                    outFrame
  ) const {
    size_t numTracks = logAmps.size();
    outFrame.assign(hopSize_, T(0));

    for (size_t i = 0; i < numTracks; ++i) {
      const auto& a = logAmps[i];
      const auto& p = phases [i];

      for (size_t n = 0; n < hopSize_; ++n) {
        // exp(a) * cos(phase)
        T partial = std::exp(a[n]) * std::cos(p[n]);
        outFrame[n] += partial;
      }
    }
  }

private:
  size_t hopSize_;
};

} // namespace PartialTracking
