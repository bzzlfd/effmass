module;
#include "pugixml.hpp"

export module pseudo.paw_upf;

import std;

// === PAW Header structure (placeholder) ===
export struct PAWHeader {
    std::string element;
    int z_valence = 0;
    int mesh_size = 0;
    int number_of_proj = 0;
    int number_of_wfc = 0;
    bool core_correction = false;
    std::string functional;
};

// === PAW class (not yet implemented) ===
export class PAW {
public:
    explicit PAW(const std::string& filename);

    auto header() const -> const PAWHeader& { return header_; }

private:
    PAWHeader header_;
};

PAW::PAW(const std::string& /*filename*/) {
    throw std::runtime_error("PAW reader not yet implemented");
}
