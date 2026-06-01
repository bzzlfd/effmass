import pseudo;
import std;


auto main() -> int {
    try {
        // Load UPF data first
        UPF upf("test/data_io_upf/Ge-spd-high.PD04.PBE.UPF");

        // Construct NCPP potential operator from UPF data
        NCPP pot(upf);

        std::println("All NCPP potential tests passed!");
        return 0;
    } catch (const std::exception& e) {
        std::println(std::cerr, "Test failed: {}", e.what());
        return 1;
    }
}
