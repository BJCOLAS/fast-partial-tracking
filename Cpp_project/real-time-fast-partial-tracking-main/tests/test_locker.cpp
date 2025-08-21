//
//  test_locker.cpp
//  Real-Time-Fast-Partial-Tracking
//
//  Created by Colas Benjamin on 26/06/2025.
//

/*
tracker.initialize(Ndft, fs);
tracker.processFrame(frame);

auto trueFrequencies = tracker.getTrueFrequencies();
auto trueAmplitudes = tracker.getTrueAmplitudes();
auto trueLogAmplitudes = tracker.getLogAmplitudes();
auto truePhases = tracker.getTruePhases();

std::cout << "\nFréquences interpolées f_true :\n";
for (size_t i = 0; i < std::min(trueFrequencies.size(), size_t(10)); ++i) {
    
    std::cout << std::setprecision(16) << trueFrequencies[i] << " Hz\n";
}

std::cout << "\nAmpltiudes interpolées amp_true :\n";
for (size_t i = 0; i < std::min(trueAmplitudes.size(), size_t(10)); ++i) {
    
    std::cout << std::setprecision(16) << trueAmplitudes[i] << " dB\n";
}

std::cout << "\nn_true : " << tracker.getNTrue() << "\n";

std::cout << "\nLogAmpltiudes interpolées Logamp_true :\n";
for (size_t i = 0; i < std::min(trueLogAmplitudes.size(), size_t(10)); ++i) {
    
    std::cout << std::setprecision(16) << trueLogAmplitudes[i] << " dB\n";
}

std::cout << "\nPhases interpolées Phases_true :\n";
for (size_t i = 0; i < std::min(truePhases.size(), size_t(10)); ++i) {
    
    std::cout << std::setprecision(16) << truePhases[i] << " dB\n";
}

std::cout << "\nn_true : " << tracker.getNTrue() << "\n";
 
 
 
 void estimateParameters(const SnailFrame<double>& frame){
     std::cout << "Estimation des paramètres de la trame : \n";
     size_t N = frame.amplitudes.size();
     
     for (size_t i = 0; i < std::min(N, size_t(5)); i++){
         std::cout << "Partial " << i << ": "
         << "amp = " << frame.amplitudes[i] << ", "
         << "phi = " << frame.phases[i] << ", "
         << "deltaPhi = " << frame.deltaPhi[i] << "\n";
     }
 }
 
 */

/*if (t == 0)
{
    if (cur.nTrue > 0)
    {
        std::cout << "=== Trame 0 ===\n";
        std::cout << " cur.nTrue = " << cur.nTrue
                  << ", freq[0] = " << cur.frequencies[0]
                  << ", logAmp[0] = " << cur.logAmplitudes[0]
                  << ", phase[0] = " << cur.phases[0] << "\n\n";
    }
    else
    {
        std::cout << "=== Trame 0 === (aucun peak détecté)\n\n";
    }
}
else if (runner.hasPair())
{
    const auto& prev = runner.previous();

    if (prev.nTrue > 0 && cur.nTrue > 0)
    {
        std::cout << "=== Trame " << (t - 1) << " vs " << t << " ===\n";
        std::cout << " prev.freq[0] = " << prev.frequencies[0]
                  << ", logAmp[0] = " << prev.logAmplitudes[0]
                  << ", phase[0] = " << prev.phases[0] << "\n";
        std::cout << " cur.freq[0]  = " << cur.frequencies[0]
                  << ", logAmp[0] = " << cur.logAmplitudes[0]
                  << ", phase[0]  = " << cur.phases[0]  << "\n\n";
    }
    else
    {
        std::cout << "=== Trame " << (t - 1) << " vs " << t << " === (aucun peak détecté)\n\n";
    }
}
// après avoir affiché prev vs cur
if (runner.hasPair()) {
    PartialTracking::TrackingParams<double> params;
    auto& C = runner.computeCostMatrix(params);
    std::cout << "CostMatrix (" << C.rows << "×" << C.cols << "):\n";
    for (std::size_t i = 0; i < C.rows; ++i) {
        for (std::size_t j = 0; j < C.cols; ++j) {
            auto  cost  = C.costs[i*C.cols + j];
            auto& label = C.labels[i*C.cols + j];
            std::cout << "(" << label << "," << cost << ") ";
        }
        std::cout << "\n";
    }
}*/
