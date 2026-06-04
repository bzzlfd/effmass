import io;
import std;


auto testVR() -> void {
    RHO vr("test/data_io-local/OUT.VR");

    const auto& m = vr.meta;

    // Metadata checks
    if (m.n1 != 20) {
        throw std::runtime_error("VR n1 mismatch: expected 20, got " + std::to_string(m.n1));
    }
    if (m.n2 != 20) {
        throw std::runtime_error("VR n2 mismatch: expected 20, got " + std::to_string(m.n2));
    }
    if (m.n3 != 32) {
        throw std::runtime_error("VR n3 mismatch: expected 32, got " + std::to_string(m.n3));
    }
    if (m.nnode != 2) {
        throw std::runtime_error("VR nnodes mismatch: expected 2, got " + std::to_string(m.nnode));
    }
    if (m.nstate != 1) {
        throw std::runtime_error("VR nstate mismatch: expected 1, got " + std::to_string(m.nstate));
    }

    // Spot-check first value
    if (std::abs(vr[0, 0, 0] - 0.1545698540504107) > 1e-9) {
        throw std::runtime_error("VR[0,0,0] mismatch: expected ~0.1545698540504107, got " +
            std::to_string(vr[0, 0, 0]));
    }

    // Spot-check some other values
    if (std::abs(vr[0, 0, 1] - 0.1682704177691634) > 1e-9) {
        throw std::runtime_error("VR[0,0,1] mismatch: expected ~0.1682704177691634, got " +
            std::to_string(vr[0, 0, 1]));
    }
    if (std::abs(vr[0, 0, 2] - 0.2065417110366990) > 1e-9) {
        throw std::runtime_error("VR[0,0,2] mismatch: expected ~0.2065417110366990, got " +
            std::to_string(vr[0, 0, 2]));
    }

    // Access with explicit state=0 (same because nstate=1)
    if (std::abs(vr[0, 0, 0, 0] - vr[0, 0, 0]) > 1e-15) {
        throw std::runtime_error("VR[0,0,0,0] != VR[0,0,0]");
    }

    vr.print_info();
}


auto testRHO() -> void {
    RHO rho("test/data_io-local/OUT.RHO");

    const auto& m = rho.meta;

    if (m.n1 != 20) {
        throw std::runtime_error("RHO n1 mismatch: expected 20, got " + std::to_string(m.n1));
    }
    if (m.n2 != 20) {
        throw std::runtime_error("RHO n2 mismatch: expected 20, got " + std::to_string(m.n2));
    }
    if (m.n3 != 32) {
        throw std::runtime_error("RHO n3 mismatch: expected 32, got " + std::to_string(m.n3));
    }
    if (m.nnode != 2) {
        throw std::runtime_error("RHO nnodes mismatch: expected 2, got " + std::to_string(m.nnode));
    }
    if (m.nstate != 1) {
        throw std::runtime_error("RHO nstate mismatch: expected 1, got " + std::to_string(m.nstate));
    }

    // RHO has different values from VR
    if (std::abs(rho[0, 0, 0] - 0.01370603011048329) > 1e-12) {
        throw std::runtime_error("RHO[0,0,0] mismatch: expected ~0.01370603011048329, got " +
            std::to_string(rho[0, 0, 0]));
    }

    rho.print_info();
}


auto main() -> int {
    try {
        testVR();
        testRHO();
        std::println("All RHO/VR tests passed!");
        return 0;
    } catch (const std::exception& e) {
        std::println(std::cerr, "Test failed: {}", e.what());
        return 1;
    }
}
