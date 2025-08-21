//
//  partial.hpp
//  Real-Time-Fast-Partial-Tracking
//
//  Created by Colas Benjamin on 13/06/2025.
//


#pragma once
#include <cstddef>



template <typename T>
struct Partial{
    T omega {0};
    T amplitude {0};
    T phase {0};
    std::size_t index {0};
    
    Partial() = default;
    
    Partial(T omega_,T amp_, T phase_, std::size_t idx): omega(omega_), amplitude(amp_), phase(phase_), index(idx) {}
};
