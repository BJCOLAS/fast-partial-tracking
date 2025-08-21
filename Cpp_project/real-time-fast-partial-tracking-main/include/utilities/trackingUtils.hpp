//
//  TrackingUtils.hpp
//  Real-Time-Fast-Partial-Tracking
//
//  Created by Colas Benjamin on 25/06/2025.
//

#pragma once

#include <vector>
#include <utility>  // std::pair

namespace PartialTracking
{
template <typename T>

void compute_frequencies_fourier(std::vector<T>& frequencies, int Ndft, T fs)
{
    std::size_t N = Ndft / 2 + 1;  // rfft size
    frequencies.resize(N);
    
    for (std::size_t k = 0; k < N; ++k)
    {
        frequencies[k] = static_cast<T>(k) * fs / Ndft;
    }
}

/*
 Fonction de determination des fréquences des partiels
 */
template<typename T>
void find_f_true(
                 const std::vector<T>& deltaPhi,
                 const std::vector<T>& freq,
                 std::vector<T>&       true_f,
                 std::vector<std::pair<int,int>>& true_f_idx,
                 std::size_t&          n_true)
{
    n_true = 0;
    auto N = deltaPhi.size();
    for(std::size_t k=0; k+1<N; ++k)
    {
        if (deltaPhi[k]>=0 && deltaPhi[k+1]<0)
        {
            T x0 = deltaPhi[k], x1 = deltaPhi[k+1];
            T y0 = freq[k],     y1 = freq[k+1];
            T f = y0 + (y1-y0)*(-x0)/(x1-x0);
            
            // on vérifie qu'on n'excède pas la capacité
            if (n_true < true_f.size())
            {
                true_f     [n_true] = f;
                true_f_idx [n_true] = {int(k), int(k+1)};
                ++n_true;
            }
        }
    }
}

/*
 Fonction  d'interpolation d'amplitude linéaire
 */
template<typename T>
void interp_amp(
                const std::vector<T>&         trame_amp,
                const std::vector<T>&         true_f,
                const std::vector<std::pair<int,int>>& true_f_idx,
                const std::vector<T>&         freq,
                std::vector<T>&               amp_out,
                std::size_t                   n_true)
{
    for(std::size_t i=0; i<n_true; ++i)
    {
        int   i0 = true_f_idx[i].first,
        i1 = true_f_idx[i].second;
        T x0 = freq[i0], x1 = freq[i1];
        T y0 = trame_amp[i0], y1 = trame_amp[i1];
        amp_out[i] = y0 + (y1-y0)*(true_f[i]-x0)/(x1-x0);
    }
}

/*
 Fonction de threshold d'amplitude
 */
template<typename T>
void amplitude_thresholding(
    T threshold,
    std::vector<T>& amp_true,
    std::vector<T>& f_true,
    std::vector<std::pair<int, int>>& true_f_index,
    std::size_t& n_true)
{
    std::size_t dst = 0;

    for (std::size_t k = 0; k < n_true; ++k)
    {
        if (amp_true[k] > threshold)
        {
            amp_true[dst]       = amp_true[k];
            f_true[dst]         = f_true[k];
            true_f_index[dst]   = true_f_index[k];
            ++dst;
        }
    }

    n_true = dst;
}

/*
 Fonction de conversion amplitude dB -> logAmplitude
 */
template<typename T>
void log_amplitude(
    const std::vector<T>& amp_dB,
    std::vector<T>& log_amp_out,
    std::size_t n_true)
{
    constexpr double kLog10Over20 = 0.11512925464970229;
    for (std::size_t i = 0; i < n_true; ++i)
    {
        log_amp_out[i] = amp_dB[i] * kLog10Over20;
    }
}

/*
 Fonction d'interpolation de phase
 */
template <typename T>
void interp_phase(
    const std::vector<T>& trame_phase,
    const std::vector<T>& true_f_,
    const std::vector<std::pair<int, int>>& true_f_index,
    const std::vector<T>& frequencies,
    std::vector<T>& phase_true_,
    std::size_t n_true)
{
    phase_true_.resize(n_true); // Tu peux aussi préallouer ça au max_candidates dans initialize()

    for (std::size_t k = 0; k < n_true; ++k)
    {
        int idx_low  = true_f_index[k].first;
        int idx_high = true_f_index[k].second;

        T phi0 = trame_phase[idx_low];
        T phi1 = trame_phase[idx_high];

        // Phase unwrap (équivalent de np.angle(exp(1j*(phi1 - phi0))))
        T delta = std::remainder(phi1 - phi0, T(2.0 * M_PI));
        T phi1_unwrapped = phi0 + delta;

        T f0 = frequencies[idx_low];
        T f1 = frequencies[idx_high];
        T f  = true_f_[k];

        T w = (f1 != f0) ? (f - f0) / (f1 - f0) : T(0.0);
        phase_true_[k] = (1.0 - w) * phi0 + w * phi1_unwrapped;
    }
}

}



