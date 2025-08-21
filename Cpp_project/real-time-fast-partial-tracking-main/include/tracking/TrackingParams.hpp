
//
//  TrackingParams.hpp
//  Real-Time-Fast-Partial-Tracking
//
//  Created by Colas Benjamin on 26/06/2025.
//

#pragma once

#include <cmath>     // std::log
#include <cstddef>   // std::size_t

/*
 Structure de données contenant les paramètres utiles au tracking des partiels
 */

namespace PartialTracking {

template<typename T>
struct TrackingParams
{
    // Paramètres initiaux du modèle
    std::size_t windowSizeMs = 50;
    std::size_t stepSizeMs   = 5;
    
    T           Peak_dB      = -50;

    T           delta        = 0.2;

    T           zetaF        = 50.0;
    T           zetaA        = 15.0;
    T           zetaA_log     = zetaA / T(20.0) * std::log(T(10.0));

    std::size_t oversample   = 2;

    // Variances calculées automatiquement
    T varF;
    T varA;

    // Constructeur
    TrackingParams()
    {
        // Évite division par zéro si delta mal défini
        T den = (delta - 2.0);
        T num = (delta - 1.0);
        T logTerm = std::log(num / den);

        varF = -zetaF * zetaF * logTerm;
        varA = -zetaA_log * zetaA_log * logTerm;
    }
};

} // namespace PartialTracking

