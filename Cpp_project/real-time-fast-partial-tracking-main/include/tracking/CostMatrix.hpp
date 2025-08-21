//
//  CostMatrix.hpp
//  Real-Time-Fast-Partial-Tracking
//
//  Created by Colas Benjamin on 26/06/2025.
//

#pragma once
#include <vector>
#include <cstddef>

/*
 Defintion  d'une structure CostMatrix pour stocker les matrices de coûts calculée à chaque itération
 */
namespace PartialTracking {

template<typename T>
struct CostMatrix {
    std::vector<T>    costs;   // valeurs de coût, rangées en row-major
    std::vector<char> labels;  // 'A' ou 'B'
    std::size_t       rows;    // nombre de lignes = prev.nTrue
    std::size_t       cols;    // nombre de colonnes = cur.nTrue

    CostMatrix(std::size_t maxN = 0)
      : costs(maxN*maxN),
        labels(maxN*maxN),
        rows(0),
        cols(0)
    {}

    // redimensionne la sous‐matrice active
    void resize(std::size_t r, std::size_t c) {
        rows = r;
        cols = c;
        costs.resize(r*c);
        labels.resize(r*c);
    }

};

} // namespace PartialTracking

