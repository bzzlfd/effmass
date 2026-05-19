module;

#include "pocketfft.h"

export module math.fft;

import std;


export constexpr bool R2G = true;
export constexpr bool G2R = false;

class FFT1D {
public:
    auto operator()(std::span<std::complex<double>> data, bool forward) const -> void {
        auto sz = data.size();
        if (sz > static_cast<std::size_t>(std::numeric_limits<int>::max()))
            throw std::runtime_error("FFT1D: span size exceeds int max");
        int n = static_cast<int>(sz);
        if (n <= 1) return;

        auto it = cache_.find(n);
        if (it == cache_.end()) {
            it = cache_.emplace(n,
                std::unique_ptr<cfft_plan_i, void(*)(cfft_plan_i*)>{
                    make_cfft_plan(n),
                    [](cfft_plan_i* p) { if (p) destroy_cfft_plan(p); }}).first;
        }
        auto& plan = it->second;

        auto* buf = reinterpret_cast<double*>(data.data());
        if (forward) {
            cfft_forward(plan.get(), buf, 1.0 / n);
        } else {
            cfft_backward(plan.get(), buf, 1.0);
        }
    }

private:
    inline static std::unordered_map<int, std::unique_ptr<cfft_plan_i, void(*)(cfft_plan_i*)>> cache_{};
};


export class FFT3D {
public:
    FFT3D(int n1, int n2, int n3)
        : scratch_(static_cast<std::size_t>(std::max({n1, n2, n3}))),
          n1_(n1), n2_(n2), n3_(n3) {}

    FFT3D(const FFT3D&) = default;
    auto operator=(const FFT3D&) -> FFT3D& = default;
    FFT3D(FFT3D&&) noexcept = default;
    auto operator=(FFT3D&&) noexcept -> FFT3D& = default;

    auto operator()(std::span<std::complex<double>> grid, bool forward) const -> void {
        int s23 = n2_ * n3_;
        int s3  = n3_;

        for (int i = 0; i < n1_; ++i) {
            for (int j = 0; j < n2_; ++j) {
                fft1d_(grid.subspan((i * n2_ + j) * n3_, n3_), forward);
            }
        }

        for (int i = 0; i < n1_; ++i) {
            for (int k = 0; k < n3_; ++k) {
                for (int j = 0; j < n2_; ++j) scratch_[j] = grid[i * s23 + j * s3 + k];
                fft1d_(std::span(scratch_.data(), n2_), forward);
                for (int j = 0; j < n2_; ++j) grid[i * s23 + j * s3 + k] = scratch_[j];
            }
        }

        for (int j = 0; j < n2_; ++j) {
            for (int k = 0; k < n3_; ++k) {
                for (int i = 0; i < n1_; ++i) scratch_[i] = grid[i * s23 + j * s3 + k];
                fft1d_(std::span(scratch_.data(), n1_), forward);
                for (int i = 0; i < n1_; ++i) grid[i * s23 + j * s3 + k] = scratch_[i];
            }
        }
    }

private:
    FFT1D fft1d_;
    mutable std::vector<std::complex<double>> scratch_;
    int n1_{0};
    int n2_{0};
    int n3_{0};
};

