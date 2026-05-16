module;

#include "pocketfft.h"

export module transform.fft3d;

import std;


namespace transform {

export auto fft1d(std::span<std::complex<double>> data, bool forward) -> void {
    int n = static_cast<int>(data.size());
    if (n <= 1) return;

    std::vector<double> buf(2 * n);
    for (int i = 0; i < n; ++i) {
        buf[2 * i]     = data[i].real();
        buf[2 * i + 1] = data[i].imag();
    }

    cfft_plan plan = make_cfft_plan(n);
    if (forward) {
        cfft_forward(plan, buf.data(), 1.0);
    } else {
        cfft_backward(plan, buf.data(), 1.0 / n);
    }
    destroy_cfft_plan(plan);

    for (int i = 0; i < n; ++i) {
        data[i] = std::complex<double>(buf[2 * i], buf[2 * i + 1]);
    }
}


auto fft1d_strided(
    std::complex<double>* data,
    int n,
    int stride,
    bool forward
) -> void {
    std::vector<std::complex<double>> buf(n);
    for (int i = 0; i < n; ++i) buf[i] = data[i * stride];
    fft1d(buf, forward);
    for (int i = 0; i < n; ++i) data[i * stride] = buf[i];
}


export auto fft3d(
    std::span<std::complex<double>> grid,
    int n1,
    int n2,
    int n3,
    bool forward
) -> void {
    int s23 = n2 * n3;
    int s3 = n3;

    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
            auto slice = grid.subspan((i * n2 + j) * n3, n3);
            fft1d(slice, forward);
        }
    }

    for (int i = 0; i < n1; ++i) {
        for (int k = 0; k < n3; ++k) {
            fft1d_strided(&grid[i * s23 + k], n2, s3, forward);
        }
    }

    for (int j = 0; j < n2; ++j) {
        for (int k = 0; k < n3; ++k) {
            fft1d_strided(&grid[j * s3 + k], n1, s23, forward);
        }
    }
}

} // namespace transform
