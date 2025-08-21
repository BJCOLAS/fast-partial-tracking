//
//  TrackingWorkspace.hpp
//  Real-Time-Fast-Partial-Tracking
//
//  Created by Colas Benjamin on 25/06/2025.
//

#pragma once

#include <vector>
#include <utility>
#include "utilities/trackingUtils.hpp"

template<typename T>
struct TrackingWorkspace
{
  std::vector<T>                 frequencies;
  std::vector<T>                 true_f_;        // taille = maxCandidates
  std::vector<std::pair<int,int>> true_f_index;  // idem
  std::vector<T>                 amp_true_;
  std::size_t                    n_true = 0;     // count logique
  std::vector<T>                 log_amp_true_;
  std::vector<T>                 phase_true_;
  void initialize(std::size_t Ndft, T fs)
  {
    PartialTracking::compute_frequencies_fourier(frequencies, Ndft, fs);
    std::size_t max_cand = Ndft/2;
    true_f_.assign(max_cand, T{});
    true_f_index.assign(max_cand, {0,0});
    amp_true_.assign(max_cand, T{});
    log_amp_true_.assign(max_cand, T{});
    phase_true_.assign(max_cand,T{});
  }

  void clearForFrame()
  {
    n_true = 0;
  }
};




