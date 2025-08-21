//
//  test_csv.cpp
//  Real-Time-Fast-Partial-Tracking
//
//  Created by Colas Benjamin on 28/05/2025.
//

#include <iostream>
#include <iomanip>
#include "CSVReader.hpp"

// You can define your test function here
void testCSVReading(const std::string& path, int trame_index) {
    auto data = CSVReader::readCSV(path);
    if (data.empty()) {
        std::cerr << "Fichier vide ou introuvable\n";
        return;
    }
    if (trame_index >= data.size()) {
        std::cerr << "Index de trame invalide\n";
        return;
    }

    std::cout << "Trame " << trame_index << ": [";
    for (double val : data[trame_index]) {
        std::cout << std::fixed << std::setprecision(15) << val << " ";
    }
    std::cout << "]\n";
    std::cout << "Shape: (" << data.size() << ", " << data[0].size() << ")\n";
}




