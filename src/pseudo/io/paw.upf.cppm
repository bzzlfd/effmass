module;
#include "pugixml.hpp"

export module pseudo.io.paw_upf;

import std;

// === PAWUPF Header structure (placeholder) ===
export struct PAWUPFHeader {
    std::string element;
    int z_valence = 0;
    int mesh_size = 0;
    int number_of_proj = 0;
    int number_of_wfc = 0;
    bool core_correction = false;
    std::string functional;
};

// === PAWUPF class (not yet implemented) ===
export class PAWUPF {
public:
    explicit PAWUPF(const std::string& filename);

    auto header() const -> const PAWUPFHeader& { return header_; }

private:
    PAWUPFHeader header_;
};

PAWUPF::PAWUPF(const std::string& /*filename*/) {
    throw std::runtime_error("PAWUPF reader not yet implemented");
}
