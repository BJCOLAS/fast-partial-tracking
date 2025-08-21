//
//  CSVReader.cpp
//  Real-Time-Fast-Partial-Tracking
//
//  Created by Colas Benjamin on 28/05/2025.
//

#include "utilities/CSVReader.hpp"
#include <fstream>
#include <sstream>
#include <iostream>

std::vector<std::vector<double>> CSVReader::readCSV(const std::string& filepath){
    std::vector<std::vector<double>> data;
    std::ifstream file(filepath);
    std::string line;
    
    if (!file.is_open()) {
        std::cerr << "Erreur : impossible d'ouvrir le fichier " << filepath << std::endl;
        return data;
    }

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string cell;
        std::vector<double> row;

        while (std::getline(ss, cell, ',')) {
            try {
                row.push_back(std::stod(cell));
            } catch (...) {
                row.push_back(0.0);  // fallback
            }
        }

        data.push_back(row);
    }

    file.close();
    return data;
}
    

