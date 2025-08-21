//
//  main.cpp
//  Real-Time-Fast-Partial-Tracking
//
//  Created by Colas Benjamin on 31/03/2025.
//

#include <iostream>
#include <iomanip>
#include "utilities/CSVReader.hpp"
#include "utilities/SnailFrame.hpp"
#include "tracking/ParameterEstimation.hpp"
#include "SnailFrame.hpp"
#include <algorithm>
#include "tracking/FrameParameters.hpp"
#include "tracking/PartialTrackingRunner.hpp"
#include "tracking/PartialsStorage.hpp"
#include "tracking/TrackManager.hpp"
#include "utilities/PartialsExporter.hpp"
#include <chrono>
#include "RectangularLsap.h"    
#include "synthesis/SynthesisDecoder.hpp"
#include "synthesis/AdditiveSynthesis.hpp"
#include "utilities/WavExporter.hpp"
#include "utilities/SignalExporter.hpp"


int main() {
        
    // Lectures des fichiers CSV
    
    std::string path_amp = "/Users/colas/Documents/Programmation/Python/Snail_estimator/oboe/oboe_impro_loudness.mat.csv";
    std::string path_phi = "/Users/colas/Documents/Programmation/Python/Snail_estimator/oboe/oboe_impro_phiDem.mat.csv";
    std::string path_dphi = "/Users/colas/Documents/Programmation/Python/Snail_estimator/oboe/oboe_impro_deltaPhi.mat.csv";
    
    //std::string path_amp = "/Users/colas/Documents/Programmation/Python/Snail_estimator/3sine/sine3_loudness.mat.csv";
    //std::string path_phi = "/Users/colas/Documents/Programmation/Python/Snail_estimator/3sine/sine3_phiDem.mat.csv";
    //std::string path_dphi = "/Users/colas/Documents/Programmation/Python/Snail_estimator/3sine/sine3_deltaPhi.mat.csv";
        
    std::vector<std::vector<double>> amp_dB = CSVReader::readCSV(path_amp);
    std::vector<std::vector<double>> phi_fourrier = CSVReader::readCSV(path_phi);
    std::vector<std::vector<double>> delta_phi = CSVReader::readCSV(path_dphi);
    
    //  Parameters :
    int Ndft = 4096;
    int Nw = 2400;
    int fs = 48000;
    int hopSize = 240;
    
    // Instanciation du runner
    const size_t M   = amp_dB.size(); //Frame number
    PartialTracking::TrackingParams<double> params;
    PartialTracking::PartialTrackingRunner<double> runner(Ndft, fs);
    
    // Instanciation du stockage des partiels
    const size_t maxPartials = Ndft/2 + 1;
    PartialTracking::PartialsStorage<double> storage(M, maxPartials);
    PartialTracking::TrackManager<double>    manager;
    
    
    // Instanciation du décodeur (plage 2 frames à chaque fois)
    PartialTracking::SynthesisDecoder<double>       decoder(Nw, fs);

    
    // Réserve les buffers de sortie
    std::vector<PartialTracking::SynthTrackParam<double>> rawParams;
    rawParams.reserve(maxPartials);
    std::vector<std::vector<double>> outLogAmps, outPhases, outFreqs;
    
    
    // Instanciation du synthétiseur :
    
    PartialTracking::AdditiveSynthesizer<double> synthesizer(hopSize);
    std::vector<double> globalSignal; // Signal final
    globalSignal.reserve(M * hopSize); // Allocation
    
    double totalTime = 0.0;
    int frameCount = 0;
    size_t totalPartials = 0;
    double maxTimeUs   = 0.0;
    
    // Boucle sur les trames
    
    for (size_t t = 0; t < M; ++t)
    {
        
        auto start = std::chrono::high_resolution_clock::now();
        
        //Instentation objet frame
        SnailFrame<double> frame;
        
        // Remplissage des valeurs
        frame.setFromVectors(amp_dB[t], phi_fourrier[t], delta_phi[t]);
        
        // Estimation des paramètres
        runner.processFrame(frame,params);
        
        // Verification que l'on a 2 trames pour lancer l'appariement
        if (! runner.hasPair()) {
            continue;
        }
        
        // Appariement des paramètres
        auto matches = runner.assignFromCostMatrix(params);
        
        totalPartials += matches.size();
        
        // Mise à jour des paramètres previous et current
        auto lastFP  = runner.previous();
        auto curFP   = runner.current();
        
        // Mise à jour du stockage des partiels dans le tableau storage
        //manager.update(t, matches, lastFP, curFP, storage);
        
        // Extraction synthétique en un seul appel
        
        decoder.processFrame(
                    /*m*/        t,
                    /*assigns*/  matches,
                    /*lastFP*/   lastFP,
                    /*curFP*/    curFP,
                    /*hopSize*/  hopSize,
                    /*outParams*/ rawParams,
                    /*outLogAmps*/ outLogAmps,
                    /*outPhases*/  outPhases,
                    /*outFreqs*/   outFreqs
                );
        
        
        // Synthèse
        
        // Synthèse additive pour cette trame
        std::vector<double> synthFrame;
        synthesizer.synthesizeFrame(outLogAmps, outPhases, synthFrame);

        // Ajout au signal global
        globalSignal.insert(globalSignal.end(), synthFrame.begin(), synthFrame.end());

        // Timer
        auto end = std::chrono::high_resolution_clock::now();
        
        std::chrono::duration<double, std::micro> elapsed_us = end - start;
        std::cout << "⏱ Trame " << t << " en "
        << elapsed_us.count() << " µs\n";
        

        if (t > 0) {
            totalTime += elapsed_us.count();
            frameCount++;
            if (elapsed_us.count() > maxTimeUs) {
                maxTimeUs = elapsed_us.count();
            }
        }
        
    }
    
    double avgTimePerFrame = totalTime / frameCount;
    double avgPartialsPerFrame = static_cast<double>(totalPartials) / M;

    std::cout << "\n Temps moyen par trame (hors première) : "
              << avgTimePerFrame << " µs\n";
    std::cout << " Temps max sur une trame : "
              << maxTimeUs << " µs\n";
    std::cout << "Nombre moyen de partiels par trame : "
              << avgPartialsPerFrame << "\n";
    
    
    // Export WAV final
    std::string outputPath = "/Users/colas/Documents/Programmation/Cpp/Real-Time-Fast-Partial-Tracking/data/wav/signal_synthetise_50dB.wav";
    PartialTracking::exportWAV(globalSignal, fs, outputPath);
    std::cout << "Signal synthétisé exporté vers : " << outputPath << "\n";

    // Export CSV
    PartialTracking::PartialsExporter<double>::toCSV(
      storage,
      "/Users/colas/Documents/Programmation/Cpp/Real-Time-Fast-Partial-Tracking/data/Partials/partials_test.csv",
      12 // 12 chiffres après la virgule
    );
    
    std::cout << "CSV généré : partials.csv\n";
    
    // Export Signal
    
    std::string csvPath = "/Users/colas/Documents/Programmation/Cpp/Real-Time-Fast-Partial-Tracking/data/signal_synthetise.csv";
    exportSignalToCSV(globalSignal, csvPath);
    std::cout << "Signal exporté vers : " << csvPath << "\n";


    return 0;
}

    


