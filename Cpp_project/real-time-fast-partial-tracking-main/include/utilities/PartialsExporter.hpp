//
//  PartialsExporter.hpp
//  Real-Time-Fast-Partial-Tracking
//
//  Created by Colas Benjamin on 27/06/2025.
//


#pragma once
#include <fstream>
#include <iomanip>
#include <string>
#include "PartialsStorage.hpp"

namespace PartialTracking {

/// Exporte storage en CSV : frame,partial,omega,amp,phase,trackID
template<typename T>
class PartialsExporter {
public:
  static void toCSV(const PartialsStorage<T>& storage,
                    const std::string&       filename,
                    int                      precision = 8)
  {
    std::ofstream ofs(filename);
    if (!ofs) throw std::runtime_error("Impossible d'ouvrir " + filename);
    ofs << "frame,partial,omega,amp,phase,trackID\n";
    ofs << std::fixed << std::setprecision(precision); //

    for (size_t m = 0; m < storage.perFrame.size(); ++m) {
      const auto& buf = storage.at(m);
      for (size_t k = 0; k < buf.size; ++k) {
        auto& e = buf.buffer[k];  // {ω, amp, phase, trackID}
        ofs
         << m << ','
         << k << ','
         << e[0] << ','
         << e[1] << ','
         << e[2] << ','
         << e[3]
         << '\n';
      }
    }
  }
};

} // namespace PartialTracking


