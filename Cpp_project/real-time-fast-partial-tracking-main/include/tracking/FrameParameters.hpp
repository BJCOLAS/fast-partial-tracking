// FrameParameters.hpp
#pragma once

#include <vector>
#include <algorithm>  // std::copy_n
#include <cstddef>    // std::size_t

namespace PartialTracking {
  template<typename T>
  class ParameterEstimation;
}

namespace PartialTracking {

template<typename T>
struct FrameParameters
{
  std::vector<T> frequencies;
  std::vector<T> logAmplitudes;
  std::vector<T> phases;
  std::size_t    nTrue = 0;

  // Constructor réserve une fois pour toute la capacité max
  FrameParameters(std::size_t maxCapacity = 0)
  {
    frequencies.reserve(maxCapacity);
    logAmplitudes.reserve(maxCapacity);
    phases.reserve(maxCapacity);
  }

  // Copie les nTrue premières valeurs du tracker sans ré-allocation
  void fromTracker(const ParameterEstimation<T>& tracker)
  {
    nTrue = tracker.getNTrue();
    frequencies.resize(nTrue);
    logAmplitudes.resize(nTrue);
    phases.resize(nTrue);

    auto freqSpan  = tracker.getTrueFrequencies();
    auto logAmpSpan= tracker.getLogAmplitudes();
    auto phaseSpan = tracker.getTruePhases();

    std::copy_n(freqSpan.begin(),    nTrue, frequencies.begin());
    std::copy_n(logAmpSpan.begin(),  nTrue, logAmplitudes.begin());
    std::copy_n(phaseSpan.begin(),   nTrue, phases.begin());
  }
};

} // namespace PartialTracking

