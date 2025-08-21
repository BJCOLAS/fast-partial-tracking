//
//  SignalExporter.hpp
//  Real-Time-Fast-Partial-Tracking
//
//  Created by Colas Benjamin on 07/08/2025.
//

#include <fstream>
#include <string>
#include <vector>

inline void exportSignalToCSV(const std::vector<double>& signal, const std::string& filename) {
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Erreur : impossible d’ouvrir le fichier CSV pour écriture.\n";
        return;
    }

    for (size_t i = 0; i < signal.size(); ++i) {
        out << signal[i] << "\n";
    }

    out.close();
}
