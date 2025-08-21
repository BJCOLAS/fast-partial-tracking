//
//  PartialTracking.h
//  Real-Time-Fast-Partial-Tracking
//
//  Created by Colas Benjamin on 31/03/2025.
//

#pragma once
#include <vector>
#include "SnailFrame.hpp"
#include "Partial.hpp"
#include "TrackingWorkspace.hpp"
#include "utilities/TrackingUtils.hpp"
#include <span>
#include "FrameParameters.hpp"
#include "TrackingParams.hpp"


namespace PartialTracking {
template <typename T>
class ParameterEstimation
{
public:
    using PartialType = Partial<T>;

    ParameterEstimation() = default;

    void reset()
    {
        frameCount = 0;
        // reset autres états si besoin
    }
    
    void initialize(std::size_t Ndft, T fs)
    {
        workspace.initialize(Ndft, fs);
    }

    // On appelle cette méthode à chaque nouvelle trame reçue
    void processFrame(const SnailFrame<T>& newFrame, const TrackingParams<T>& params)
    {

        // Initalisation
        workspace.clearForFrame();

        // Determination des paramètres
        
        PartialTracking::find_f_true<T>(
            newFrame.deltaPhi,
            workspace.frequencies,
            workspace.true_f_,
            workspace.true_f_index,
            workspace.n_true
        );
        PartialTracking::interp_amp<T>(
            newFrame.amplitudes,
            workspace.true_f_,
            workspace.true_f_index,
            workspace.frequencies,
            workspace.amp_true_,
            workspace.n_true
        );
        
        PartialTracking::amplitude_thresholding<T>(
            params.Peak_dB,
            workspace.amp_true_,
            workspace.true_f_,
            workspace.true_f_index,
            workspace.n_true
        );
        PartialTracking::log_amplitude<T>(
            workspace.amp_true_,
            workspace.log_amp_true_,
            workspace.n_true
        );
        PartialTracking::interp_phase<T>(
            newFrame.phases,
            workspace.true_f_,
            workspace.true_f_index,
            workspace.frequencies,
            workspace.phase_true_,
            workspace.n_true
        );

        // Incrémente le compteur de trame
        ++frameCount;
    }

    
    std::span<const T> getTrueFrequencies() const {
        return std::span<const T>(workspace.true_f_.data(), workspace.n_true);
    }
    
    std::span<const T> getTrueAmplitudes() const {
        return std::span<const T>(workspace.amp_true_.data(), workspace.n_true);
    }
    
    std::span<const T> getLogAmplitudes() const {
        return std::span<const T>(workspace.log_amp_true_.data(), workspace.n_true);
    }
    
    std::span<const T> getTruePhases() const {
        return std::span<const T>(workspace.phase_true_.data(), workspace.n_true);
    }
        
    std::size_t getNTrue() const { return workspace.n_true; }

    std::size_t getFrameCount() const { return frameCount; }



private:
    std::size_t frameCount{0};
    
    TrackingWorkspace<T> workspace;
};

template<typename T>
  FrameParameters<T> extractFrameParameters(const ParameterEstimation<T>& tracker)
  {
    return {
      tracker.getTrueFrequencies(),
      tracker.getLogAmplitudes(),
      tracker.getTruePhases(),
      tracker.getNTrue()
    };
  }

}
