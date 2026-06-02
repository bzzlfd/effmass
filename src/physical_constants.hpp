#pragma once


constexpr double BOHR_RADIUS_ANGSTROM = 0.52917721067;  // 1 Bohr = 0.52917721067 Å
constexpr double HARTREE_TO_EV = 27.211386245988;        // CODATA 2018

// Energy unit conversions: 1 Hartree = 2 Rydberg
constexpr double HARTREE_TO_RYDBERG = 2.0;

enum class EnergyUnit { Hartree, Rydberg, eV };

// Factor to convert 1 unit of `from` into `to`:  x [from] * factor = x [to]
constexpr auto convertEnergy(EnergyUnit from, EnergyUnit to) -> double {
    auto perHartree = [](EnergyUnit u) -> double {
        using enum EnergyUnit;
        switch (u) {
            case Rydberg: return 1.0 / HARTREE_TO_RYDBERG;  // 1 Ry = 0.5 Ha
            case eV:      return 1.0 / HARTREE_TO_EV;
            default:      return 1.0;  // Hartree
        }
    };
    return perHartree(from) / perHartree(to);
}

constexpr auto convertEnergy(double value, EnergyUnit from, EnergyUnit to) -> double {
    return value * convertEnergy(from, to);
}
