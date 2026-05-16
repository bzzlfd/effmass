module;

#ifdef EFFMASS_USE_POCKETFFT
#include "pocketfft.h"
#endif

export module transform.fft3d;

import std;

namespace transform {

namespace detail {

    auto bitReverse(int x, int log2n) -> int {
        int y = 0;
        for (int i = 0; i < log2n; ++i) {
            y = (y << 1) | (x & 1);
            x >>= 1;
        }
        return y;
    }

    auto isPowerOf2(int n) -> bool {
        return n > 0 && (n & (n - 1)) == 0;
    }

    auto nextPowerOf2(int n) -> int {
        int m = 1;
        while (m < n) m <<= 1;
        return m;
    }

    auto fft1d_radix2(std::span<std::complex<double>> data, bool forward) -> void {
        int n = static_cast<int>(data.size());
        int log2n = 0;
        while ((1 << log2n) < n) ++log2n;

        for (int i = 0; i < n; ++i) {
            int j = bitReverse(i, log2n);
            if (i < j) std::swap(data[i], data[j]);
        }

        for (int len = 2; len <= n; len <<= 1) {
            double angle = (forward ? -2.0 : 2.0) * std::numbers::pi / len;
            std::complex<double> wlen(std::cos(angle), std::sin(angle));
            for (int i = 0; i < n; i += len) {
                std::complex<double> w = 1.0;
                int half = len >> 1;
                for (int j = 0; j < half; ++j) {
                    auto u = data[i + j];
                    auto v = data[i + j + half] * w;
                    data[i + j] = u + v;
                    data[i + j + half] = u - v;
                    w *= wlen;
                }
            }
        }

        if (!forward) {
            double inv_n = 1.0 / n;
            for (int i = 0; i < n; ++i) data[i] *= inv_n;
        }
    }

    auto fft1d_direct(std::span<std::complex<double>> data, bool forward) -> void {
        int n = static_cast<int>(data.size());
        std::vector<std::complex<double>> out(n);
        double sign = forward ? -1.0 : 1.0;
        for (int k = 0; k < n; ++k) {
            std::complex<double> sum = 0.0;
            for (int j = 0; j < n; ++j) {
                double angle = sign * 2.0 * std::numbers::pi * j * k / n;
                sum += data[j] * std::complex<double>(std::cos(angle), std::sin(angle));
            }
            out[k] = sum;
        }
        if (!forward) {
            for (int i = 0; i < n; ++i) out[i] /= n;
        }
        for (int i = 0; i < n; ++i) data[i] = out[i];
    }

    auto fft1d_bluestein(std::span<std::complex<double>> data, bool forward) -> void {
        int n = static_cast<int>(data.size());
        int M = nextPowerOf2(2 * n - 1);

        std::vector<std::complex<double>> a(M, 0.0);
        std::vector<std::complex<double>> b(M, 0.0);

        for (int i = 0; i < n; ++i) {
            double angle = std::numbers::pi * static_cast<double>(i * i) / n;
            std::complex<double> w(std::cos(angle), std::sin(angle));
            a[i] = data[i] * w;
            b[i] = std::conj(w);
        }

        fft1d_radix2(a, false);
        fft1d_radix2(b, false);

        for (int i = 0; i < M; ++i) a[i] *= b[i];

        fft1d_radix2(a, true);

        for (int i = 0; i < n; ++i) {
            double angle = std::numbers::pi * static_cast<double>(i * i) / n;
            std::complex<double> w(std::cos(angle), -std::sin(angle));
            data[i] = a[i] * w;
        }
        if (!forward) {
            double inv_n = 1.0 / n;
            for (int i = 0; i < n; ++i) data[i] *= inv_n;
        }
    }

} // namespace detail


export auto fft1d_custom(std::span<std::complex<double>> data, bool forward) -> void {
    int n = static_cast<int>(data.size());
    if (n <= 1) return;

    if (detail::isPowerOf2(n)) {
        detail::fft1d_radix2(data, forward);
    } else if (n <= 64) {
        detail::fft1d_direct(data, forward);
    } else {
        detail::fft1d_bluestein(data, forward);
    }
}

#ifdef EFFMASS_USE_POCKETFFT

export auto fft1d_pocket(std::span<std::complex<double>> data, bool forward) -> void {
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

#endif


export auto fft1d(std::span<std::complex<double>> data, bool forward) -> void {
#ifdef EFFMASS_USE_POCKETFFT
    fft1d_pocket(data, forward);
#else
    fft1d_custom(data, forward);
#endif
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
