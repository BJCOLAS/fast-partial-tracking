#pragma once
#include <vector>
#include <array>
#include <tuple>
#include <cmath>
#include "tracking/FrameParameters.hpp"
#include "tracking/TrackManager.hpp"   // TrackManager<T> est défini ici

namespace PartialTracking {

/// Paramètres à extraire pour chaque piste à la trame m
template<typename T>
struct SynthTrackParam {
    std::size_t     trackId;
    std::array<T,2> Omega_m;   // fréquence normalisée à m-1 et m
    std::array<T,2> Amp_m;     // log-amplitude à m-1 et m
    std::array<T,2> Phase_m;   // phase centrée à m-1 et m
};

/// Décodeur consolidé sans stockage additionnel
template<typename T>
class SynthesisDecoder {
public:
  explicit SynthesisDecoder(size_t windowSize, T sampleRate)
    : Tc_(T(windowSize) * T(0.5))
    , fs_(sampleRate)
    , numTracks_(0)
    , numActiveLast_(0)
  {}

  /**
   * Traite une paire de trames (m-1, m) et génère les SynthTrackParam
   *
   * @param m             : indice de la trame courante
   * @param assignments   : vector de pairs (i_prev, j_cur)
   * @param lastFP, curFP : FrameParameters à m-1 et m
   * @param hopSize       : saut en échantillons
   * @param outParams     : vecteur pré-réservé pour stocker le résultat
   */
    void processFrame(
      size_t m,
      const std::vector<std::pair<size_t,size_t>>& assignments,
      const FrameParameters<T>& lastFP,
      const FrameParameters<T>& curFP,
      size_t hopSize,

      // sorties
      std::vector<SynthTrackParam<T>>&     outParams,
      std::vector<std::vector<T>>&         outLogAmps,
      std::vector<std::vector<T>>&         outPhases,
      std::vector<std::vector<T>>&         outFreqs
    ){
    
        
    
    // 1) Prépare B = {i_prev pour chaque assignment}
    std::vector<std::size_t> B;
    B.reserve(assignments.size());
    for (auto &pr : assignments)
      B.push_back(pr.first);

    // 2) Matching et attribution d’ID via matchs set dans track manager
      auto [matches, nomatchA, nomatchB]
            = match_sets(pairsLast_, B); //
      
    std::size_t numC      = matches.size();
    std::size_t numB      = nomatchB.size();
    std::size_t numActive = numC + numB;

    // Reconstruction de pairs_ et tracks_
    pairs_.resize(numActive);
    tracks_.resize(numActive);

    
    // Continuations existantes
    for (std::size_t idx = 0; idx < numC; ++idx) {
      auto [ia, ib] = matches[idx];
      pairs_[idx]  = assignments[ib].second;     // j_cur
      tracks_[idx] = tracksLast_[ia];            // même trackId
    }
    // Nouvelles pistes (naissances)
    for (std::size_t b = 0; b < numB; ++b) {
      std::size_t bi = nomatchB[b];
      pairs_[numC + b]  = assignments[bi].second;
      tracks_[numC + b] = numTracks_ + b;        // nouveaux IDs
    }

    // 3) Extraction des SynthTrackParam
    outParams.clear();
    outParams.reserve(numActive);

    // Temps centrés
    T t0 = T((m-1)*hopSize) + Tc_;
    T t1 = T(m*hopSize)   + Tc_;
    T H  = std::abs(t1 - t0);

     
        // === Naissances ===
        for (std::size_t k = 0; k < nomatchB.size(); ++k) {
          auto bi = nomatchB[k];                // index dans assignments
          auto j  = assignments[bi].second;     // index dans curFP
          T f     = curFP.frequencies[j];
          T a     = curFP.logAmplitudes[j];
          T ph    = curFP.phases[j];

          T omega = (T(2) * M_PI * f) / fs_;
          T phi1  = ph + omega * Tc_;
          T phi0  = phi1 - omega * H;

          outParams.push_back({
            tracks_[numC + k],             // nouvel ID
            {omega, omega},
            {T(-1e10), a},
            {phi0, phi1}
          });
        }

        // === Continuations ===
        for (auto [ia, ib] : matches) {
          // ia → index dans pairsLast_  (cad index dans lastFP)
          // ib → offset dans assignments (cad index dans curFP)
          std::size_t prevIdx = pairsLast_[ia];
          std::size_t curIdx  = assignments[ib].second;

          T f0 = lastFP.frequencies[prevIdx];
          T f1 = curFP   .frequencies[curIdx];
          T a0 = lastFP.logAmplitudes[prevIdx];
          T a1 = curFP   .logAmplitudes[curIdx];

          T omega0 = (T(2) * M_PI * f0) / fs_;
          T omega1 = (T(2) * M_PI * f1) / fs_;    // **ATTENTION** : utilise f1 ici
          T ph0    = lastFP.phases[prevIdx] + omega0 * Tc_;
          T ph1    = curFP.phases[curIdx]     + omega1 * Tc_;

          outParams.push_back({
            tracks_[ib],                 // même trackId que continué
            {omega0, omega1},
            {a0, a1},
            {ph0, ph1}
          });
        }

        // === Morts ===
        for (std::size_t ia : nomatchA) {
          // ia → index dans pairsLast_  (cad index dans lastFP)
          std::size_t prevIdx = pairsLast_[ia];

          T fLast   = lastFP.frequencies[prevIdx];
          T aLast   = lastFP.logAmplitudes[prevIdx];
          T phiLast = lastFP.phases[prevIdx];

          T omega = (T(2) * M_PI * fLast) / fs_;
          T phi0  = phiLast       + omega * Tc_;
          T phi1  = phi0          + omega * H;

          outParams.push_back({
            tracksLast_[ia],          // trackId inchangé
            {omega, omega},           // fréquence constante
            {aLast, T(-1e10)},        // extinction en log-amp
            {phi0,  phi1}
          });
        }
    
    

    // Echantillonage :
        
    // Nombre de piste active / frame :
    size_t N = outParams.size();
        
    // Durée d’un hop en échantillons et en secondes
    T Hs     = T(hopSize);       // taille de saut en échantillons
    T Ts     = Hs / fs_;         // durée du saut en secondes
    T two_pi = T(2) * M_PI;      // constante 2π
        
    outLogAmps.resize(N);
    outPhases .resize(N);
    outFreqs  .resize(N);
        
    for (size_t i = 0; i < N; ++i) {
      // 1) Récupère les 4 “points” frontières
      auto &P  = outParams[i];
      T f0     = P.Omega_m[0];    // freq à m-1
      T f1     = P.Omega_m[1];    // freq à m
      T p0     = P.Phase_m [0];   // phase à m-1
      T p1     = P.Phase_m [1];   // phase à m
      T A0     = P.Amp_m   [0];   // log-amp à m-1
      T A1     = P.Amp_m   [1];   // log-amp à m
        
      T amp0    = std::exp(A0);
      T amp1    = std::exp(A1);

      // 2) “Unwrap” de la phase pour éviter les sauts 2π
      T Mstar  = std::floor((Hs*(f0+f1) + 2*(p0-p1)) / (4*M_PI) + T(0.5));
      T p1_adj = p1 + two_pi * Mstar;

      // 3) Calcul des deltas pour la spline cubique
      //    D1 = p1_adj - p0 - f0*Ts,   D2 = f1 - f0
      T D1 = p1_adj - p0 - f0 * Ts;
      T D2 = f1     - f0;

      // 4) Coefficients du polynôme φ(t) = c0 + c1 t + c2 t² + c3 t³
      T c0 = p0;
      T c1 = f0;

      T c2 = (6 * Mstar * M_PI) / (Hs * Hs)
               - (2 * f0) / Hs
               - (3 * p0) / (Hs * Hs)
               - f1 / Hs
               + (3 *p1) / (Hs * Hs);

      T c3 = (-4 * Mstar * M_PI) / (Hs * Hs * Hs)
               + f0 / (Hs * Hs)
               + (2 * p0) / (Hs * Hs * Hs)
               + f1 / (Hs * Hs)
               - (2 * p1) / (Hs * Hs * Hs);

      // 5) Pente linéaire pour la log-amplitude
      T ampSlope = (amp1 - amp0) / Hs;


      // 6) Échelle pour convertir φ′(t) en Hz
      T scale = fs_ / two_pi;
        

      outLogAmps[i].resize(hopSize);
      outPhases [i].resize(hopSize);
      outFreqs  [i].resize(hopSize);


      // 7) Échantillonnage en Horner (une seule boucle n)
      for (size_t n = 0; n < hopSize; ++n) {
        T t = T(n);

        // log-amplitude interpolée
        T amp_lin = amp0 + ampSlope * t;
        outLogAmps[i][n] = std::log(amp_lin);

        // phase interpolée via Horner
        // φ(t) = ((c3*t + c2)*t + c1)*t + c0
        outPhases[i][n]  = ((c3*t + c2)*t + c1)*t + c0;

        // fréquence instantanée = scale * φ′(t)
        // φ′(t) = c1 + 2·c2·t + 3·c3·t²
        outFreqs[i][n]   = scale * (c1 + T(2)*c2 * t + T(3)*c3 * t * t);
      }
    }
        
    

    // 4) Mise à jour de l’état interne
    pairsLast_     = pairs_;
    tracksLast_    = tracks_;
    numActiveLast_ = numActive;
    numTracks_    += numB;
  }

private:
  T                                 Tc_;
  T                                 fs_ ;
  std::size_t                        numTracks_;
  std::size_t                        numActiveLast_;
  std::vector<std::size_t>          pairsLast_, tracksLast_;
  std::vector<std::size_t>          pairs_,     tracks_;
};

} // namespace PartialTracking

