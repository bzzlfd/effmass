module;
#include "pugixml.hpp"

export module pseudo.uspp_upf;

import std;

// === USPP Header structure (placeholder) ===
export struct USPPHeader {
    std::string element;
    int z_valence = 0;
    int mesh_size = 0;
    int number_of_proj = 0;
    int number_of_wfc = 0;
    bool core_correction = false;
    std::string functional;
};

// === USPP class (not yet implemented) ===
export class USPP {
public:
    explicit USPP(const std::string& filename);

    auto header() const -> const USPPHeader& { return header_; }

private:
    USPPHeader header_;
};

USPP::USPP(const std::string& /*filename*/) {
    throw std::runtime_error("USPP reader not yet implemented");
}
