//
//  PartialTrackingRunner.hpp
//  Real-Time-Fast-Partial-Tracking
//
//  Created by Colas Benjamin on 26/06/2025.
//

// PartialTrackingRunner.hpp
// PartialTrackingRunner.hpp
// PartialTrackingRunner.hpp

#pragma once
#include "ParameterEstimation.hpp"   // PartialTracker<T>
#include "FrameParameters.hpp"   // FrameParameters<T>
#include "CostMatrix.hpp"        // CostMatrix<T>
#include <cstddef>
#include <optional>
#include <cmath>                 // std::exp, std::min
#include "TrackingParams.hpp"
#include "RectangularLsap.h"

/*
 Apariement des paramètres current et previous pour former les trajectoires
 1. Calcul de la matrice de coût
 2. Resolution du linear assignement problem
 
 */
namespace PartialTracking {

template<typename T>
class PartialTrackingRunner
{
public:
    PartialTrackingRunner(std::size_t Ndft, T fs)
      : maxCandidates_(Ndft/2 + 1),
        prevParams_(maxCandidates_),
        curParams_(maxCandidates_),
        costMatrix_(maxCandidates_),
        hasPrev_(false)
    {
        tracker_.initialize(Ndft, fs);
    }

    void processFrame(const SnailFrame<T>& frame,const TrackingParams<T>& params)
    {
        tracker_.processFrame(frame,params);

        if (hasPrev_)
            prevParams_ = curParams_;
        else
            hasPrev_ = true;

        curParams_.fromTracker(tracker_);
    }

    bool hasPair()  const { return hasPrev_; }
    auto& previous() const { return prevParams_; }
    auto& current()  const { return curParams_; }

    // Remplit costMatrix_ et renvoie sa référence
    CostMatrix<T>& computeCostMatrix(const TrackingParams<T>& p) {
      std::size_t R = prevParams_.nTrue;
      std::size_t C = curParams_.nTrue;
      costMatrix_.resize(R, C);
      for (std::size_t i = 0; i < R; ++i) {
        T f1 = prevParams_.frequencies[i], a1 = prevParams_.logAmplitudes[i];
        for (std::size_t j = 0; j < C; ++j) {
          T f2 = curParams_.frequencies[j], a2 = curParams_.logAmplitudes[j];
          T df = f1 - f2, da = a1 - a2;
          T Ause = 1 - std::exp(-(df*df)/p.varF - (da*da)/p.varA);
          T Bspa = 1 - (1 - p.delta)*Ause;
          T c    = std::min(Ause, Bspa);
          auto idx = i*C + j;              // **row-major**
          costMatrix_.costs [idx] = c;
          costMatrix_.labels[idx] = (Ause <= Bspa ? 'A' : 'B');
        }
      }
      return costMatrix_;
    }

    
    std::vector<std::pair<std::size_t,std::size_t>>
    assignFromCostMatrix(const TrackingParams<T>& params)
    {
      auto& cm = computeCostMatrix(params);
      std::size_t R = cm.rows, C = cm.cols;

      // arrays de sortie row->col et col->row
      std::vector<int64_t> row2col(R), col2row(C);

      int status = solve_rectangular_linear_sum_assignment(
        static_cast<intptr_t>(R),      // rows
        static_cast<intptr_t>(C),      // cols
        cm.costs.data(),               // row-major R×C
        false,                         // minimize
        row2col.data(),                // chaque i de [0..R) → j
        col2row.data()                 // optionnel
      );
      if (status != 0) {
        std::cerr << "LSAP failed (code " << status << ")\n";
        return {};
      }

        std::vector<std::pair<std::size_t, std::size_t>> assignments;
        assignments.reserve(std::min(R, C));
        int max_loop = std::min(R, C);
        for (std::size_t j = 0; j < max_loop; ++j) {
            int64_t i = col2row[j];
            int64_t l = row2col[j];
            if (i < 0 || i >= static_cast<int64_t>(R)) continue;

            std::size_t idx = l * C + i;
            if (cm.labels[idx] == 'A') {
                assignments.emplace_back(
                    static_cast<std::size_t>(i),
                    static_cast<std::size_t>(l)
                );
            }
        }

        // Inversion des indices pour obtenir résultat similaire au Python
        std::vector<std::pair<std::size_t, std::size_t>> reordered;
        reordered.reserve(assignments.size());
        for (const auto& [i, j] : assignments) {
            reordered.emplace_back(j, i);
        }

        // Tri croissant (fait implicitement en Python)
        std::sort(reordered.begin(), reordered.end(),
                  [](const auto& a, const auto& b) {
                      return a.first < b.first;
                  });
        
        return reordered;
    }
    
private:
    ParameterEstimation<T>   tracker_;
    std::size_t          maxCandidates_;
    FrameParameters<T>   prevParams_, curParams_;
    CostMatrix<T>        costMatrix_;
    bool                 hasPrev_;
};

} // namespace PartialTracking
