// TrackManager.hpp
#pragma once
#include <vector>
#include <tuple>
#include <unordered_set>
#include <stdexcept>
#include <cmath>                // pour M_PI
#include "PartialsStorage.hpp"  // PerFramePartials<T> doit être array<T,4>
#include "FrameParameters.hpp"  // FrameParameters<T>

namespace PartialTracking {

/// ---------- Ici on retire le template<typename T> ----------
/// match_sets : retourne (match, nomatchA, nomatchB)
static inline
std::tuple<
  std::vector<std::pair<std::size_t,std::size_t>>,
  std::vector<std::size_t>,
  std::vector<std::size_t>
>
match_sets(
  const std::vector<std::size_t>& A,            // indices i de la frame m-1
  const std::vector<std::size_t>& B_indices      // B = assignments[].first
) {
    std::size_t nA = A.size(), nB = B_indices.size();
    std::unordered_set<std::size_t> setA(A.begin(), A.end());
    if (setA.size() != nA)
        throw std::runtime_error("match_sets: A must be unique");

    std::vector<bool> validA(nA, true), validB(nB, true);
    std::vector<std::pair<std::size_t,std::size_t>> match;
    
    for (std::size_t a = 0; a < nA; ++a) {
        if (!validA[a]) continue;
        for (std::size_t b = 0; b < nB; ++b) {
            if (validB[b] && B_indices[b] == A[a]) {
                match.emplace_back(a, b);
                validA[a] = validB[b] = false;
                break;
            }
        }
    }

    std::vector<std::size_t> nomatchA, nomatchB;
    for (std::size_t a = 0; a < nA; ++a)
        if (validA[a]) nomatchA.push_back(a);
    for (std::size_t b = 0; b < nB; ++b)
        if (validB[b]) nomatchB.push_back(b);

    return { match, nomatchA, nomatchB };
}

/*
 Gestion de du stockage des partiels dans PartialsStorage 
 */
template<typename T>
class TrackManager
{
public:
    TrackManager()
      : numTracks_(0), numActiveLast_(0)
    {}

    void update(
      std::size_t                                        m,
      const std::vector<std::pair<std::size_t,std::size_t>>& assignments,
      const FrameParameters<T>&                         lastFP,
      const FrameParameters<T>&                         curFP,
      PartialsStorage<T>&                               storage
    ) {
        // 1) Préparer B = [i_prev for each assignment]
        std::vector<std::size_t> B;
        B.reserve(assignments.size());
        for (auto &pr : assignments)
            B.push_back(pr.first);

        // 2) Découper en continuations vs naissances
        auto [toContinue, nomatchA, nomatchB] = match_sets(pairsLast_, B);
        std::size_t numC      = toContinue.size();
        std::size_t numB      = nomatchB.size();
        std::size_t numActive = numC + numB;

        // 3) Reconstituer pairs_ et tracks_
        pairs_.resize(numActive);
        tracks_.resize(numActive);
        for (std::size_t idx = 0; idx < numC; ++idx) {
            auto [aIdx, bIdx] = toContinue[idx];
            pairs_[idx]  = assignments[bIdx].second; // j_cur
            tracks_[idx] = tracksLast_[aIdx];
        }
        for (std::size_t b = 0; b < numB; ++b) {
            std::size_t bIdx = nomatchB[b];
            pairs_[numC + b]  = assignments[bIdx].second;
            tracks_[numC + b] = numTracks_ + b;
        }
        
        // +etape : Initialisation du premier point de stockage à 0 :
        if (m==0){
            auto& bufFirst = storage.at(0);
            bufFirst.resize(1);
            bufFirst[0]={
                0,
                0,
                0,
                0
            };
        }
        // 4) Écrire les naissances dans m-1
        if (numB > 0 && m > 0) {
            auto& bufLast = storage.at(m-1);
            bufLast.resize(numActiveLast_ + numB);
            for (std::size_t b = 0; b < numB; ++b) {
                auto [iPrev, jCur] = assignments[nomatchB[b]];
                bufLast[numActiveLast_ + b] = {
                    T(2*M_PI)/48000 * lastFP.frequencies[iPrev],
                    lastFP.logAmplitudes[iPrev],
                    lastFP.phases[iPrev],
                    T(tracks_[numC + b])
                };
            }
        }

        // 5) Écrire tous les partiels actifs dans m
        {
            auto& bufCur = storage.at(m);
            bufCur.resize(numActive);
            for (std::size_t k = 0; k < numActive; ++k) {
                std::size_t assignIdx = (k < numC
                  ? toContinue[k].second
                  : nomatchB[k - numC]);
                auto [iPrev, jCur] = assignments[assignIdx];
                bufCur[k] = {
                    T(2*M_PI)/48000 * curFP.frequencies[jCur],
                    curFP.logAmplitudes[jCur],
                    curFP.phases[jCur],
                    T(tracks_[k])
                };
            }
            
        }

        // 6) Mise à jour d’état
        pairsLast_     = pairs_;
        tracksLast_    = tracks_;
        numTracks_    += numB;
        numActiveLast_ = numActive;
    }

private:
    std::size_t               numTracks_     = 0;
    std::size_t               numActiveLast_ = 0;
    std::vector<std::size_t>  pairsLast_;       // i de m-1
    std::vector<std::size_t>  tracksLast_;      // trackID m-1

    std::vector<std::size_t>  pairs_, tracks_;  // buffers temporaires
};

} // namespace PartialTracking

