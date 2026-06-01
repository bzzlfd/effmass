export module utils.array2d;

import std;


export template<typename T>
struct array2d {
    std::vector<T> data;
    int rows = 0;
    int cols = 0;

    array2d() = default;
    array2d(int r, int c) : data(static_cast<std::size_t>(r) * c), rows(r), cols(c) {}

    // Multi-parameter subscript: a[i, j]
    auto operator[](int i, int j) const -> const T& {
        return data[static_cast<std::size_t>(i) * cols + j];
    }
    auto operator[](int i, int j) -> T& {
        return data[static_cast<std::size_t>(i) * cols + j];
    }

    // Row-wise subscript: a[i][j]
    auto operator[](int i) -> std::span<T> {
        return std::span<T>(data.data() + static_cast<std::size_t>(i) * cols,
                            static_cast<std::size_t>(cols));
    }
    auto operator[](int i) const -> std::span<const T> {
        return std::span<const T>(data.data() + static_cast<std::size_t>(i) * cols,
                                  static_cast<std::size_t>(cols));
    }

    [[nodiscard]] auto size() const -> std::size_t { return data.size(); }
    [[nodiscard]] auto empty() const -> bool { return data.empty(); }
};
