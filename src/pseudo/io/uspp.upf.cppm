module;
#include "pugixml.hpp"

export module pseudo.io.uspp_upf;

import std;


// === USPPUPF Header structure (placeholder) ===
export struct USPPUPFHeader {
    std::string element;
    int z_valence = 0;
    int mesh_size = 0;
    int number_of_proj = 0;
    int number_of_wfc = 0;
    bool core_correction = false;
    std::string functional;
};


// === USPPUPF class (not yet implemented) ===
export class USPPUPF {
public:
    explicit USPPUPF(const std::string& filename);

    auto header() const -> const USPPUPFHeader& { return header_; }

private:
    USPPUPFHeader header_;
};


USPPUPF::USPPUPF(const std::string& /*filename*/) {
    throw std::runtime_error("USPPUPF reader not yet implemented");
}
