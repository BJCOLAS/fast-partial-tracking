//
//  SnailFrame.hpp
//  Real-Time-Fast-Partial-Tracking
//
//  Created by Colas Benjamin on 28/05/2025.
//

#pragma once

#include <vector>
#include <stdexcept>
#include <cmath>  // pour M_PI
/*
 Classe definissant un objet Frame du Snail-Analyse contenant : Amplitudes, Phases et DeltaPhi.
 Conversion de la phase [0,1] -> [0,2pi]
 */
template <typename FloatType = float>
struct SnailFrameInterface
{
    virtual const std::vector<FloatType>& getAmplitudes() const = 0;
    virtual const std::vector<FloatType>& getPhases() const = 0;
    virtual const std::vector<FloatType>& getDeltaPhi() const = 0;

    virtual std::vector<FloatType>& getAmplitudes() = 0;
    virtual std::vector<FloatType>& getPhases() = 0;
    virtual std::vector<FloatType>& getDeltaPhi() = 0;

    virtual ~SnailFrameInterface() = default;
};

template <typename FloatType = float>
struct SnailFrame : public SnailFrameInterface<FloatType>
{
    std::vector<FloatType> amplitudes;
    std::vector<FloatType> phases;     // en radians
    std::vector<FloatType> deltaPhi;

    void resize(size_t size)
    {
        amplitudes.assign(size, FloatType(0));
        phases.assign(size, FloatType(0));
        deltaPhi.assign(size, FloatType(0));
    }

    void setFromVectors(const std::vector<FloatType>& amps,
                        const std::vector<FloatType>& phis_normalized,
                        const std::vector<FloatType>& dphis)
    {
        size_t size = amps.size();
        if (phis_normalized.size() != size || dphis.size() != size)
            throw std::runtime_error("setFromVectors: input vectors must have same size.");

        amplitudes = amps;
        deltaPhi = dphis;

        phases.resize(size);
        const FloatType twoPi = FloatType(2.0) * FloatType(M_PI);
        for (size_t i = 0; i < size; ++i)
            phases[i] = phis_normalized[i] * twoPi;
    }

    const std::vector<FloatType>& getAmplitudes() const override { return amplitudes; }
    const std::vector<FloatType>& getPhases() const override     { return phases; }
    const std::vector<FloatType>& getDeltaPhi() const override   { return deltaPhi; }

    std::vector<FloatType>& getAmplitudes() override { return amplitudes; }
    std::vector<FloatType>& getPhases() override     { return phases; }
    std::vector<FloatType>& getDeltaPhi() override   { return deltaPhi; }
};
