//
//  CSVRead.hpp
//  Real-Time-Fast-Partial-Tracking
//
//  Created by Colas Benjamin on 28/05/2025.
//

#pragma once

#include <string>
#include <vector>


class CSVReader {
public :
    static std::vector<std::vector<double>> readCSV(const std::string& filepath);
};

