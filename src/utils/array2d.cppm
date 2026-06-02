export module utils.array2d;

import std;


export template<typename T, int Rows = -1, int Cols = -1>
struct array2d;


// Dynamic-sized (heap) specialization: array2d<T>
export template<typename T>
struct array2d<T, -1, -1> {
    std::vector<T> data;
    int rows = 0;
    int cols = 0;

    array2d() = default;
    array2d(int r, int c) : data(static_cast<std::size_t>(r) * c), rows(r), cols(c) {}

    auto operator[](int i, int j) const -> const T& { return data[static_cast<std::size_t>(i) * cols + j]; }
    auto operator[](int i, int j) -> T& { return data[static_cast<std::size_t>(i) * cols + j]; }

    auto operator[](int i) -> std::span<T> {
        return {data.data() + static_cast<std::size_t>(i) * cols, static_cast<std::size_t>(cols)};
    }
    auto operator[](int i) const -> std::span<const T> {
        return {data.data() + static_cast<std::size_t>(i) * cols, static_cast<std::size_t>(cols)};
    }

    [[nodiscard]] auto size() const -> std::size_t { return data.size(); }
    [[nodiscard]] auto empty() const -> bool { return data.empty(); }
};


// Static-size (stack) specialization: array2d<T, Rows, Cols>  where Rows > 0 && Cols > 0
export template<typename T, int Rows, int Cols>
    requires (Rows > 0 && Cols > 0)
struct array2d<T, Rows, Cols> {
    static constexpr int rows = Rows;
    static constexpr int cols = Cols;
    std::array<T, static_cast<std::size_t>(Rows) * Cols> data{};

    array2d() = default;

    auto operator[](int i, int j) const -> const T& {
        return data[static_cast<std::size_t>(i) * cols + j];
    }
    auto operator[](int i, int j) -> T& {
        return data[static_cast<std::size_t>(i) * cols + j];
    }

    auto operator[](int i) -> std::span<T> {
        return {data.data() + static_cast<std::size_t>(i) * cols, static_cast<std::size_t>(cols)};
    }
    auto operator[](int i) const -> std::span<const T> {
        return {data.data() + static_cast<std::size_t>(i) * cols, static_cast<std::size_t>(cols)};
    }

    [[nodiscard]] auto size() const -> std::size_t { return data.size(); }
    [[nodiscard]] auto empty() const -> bool { return false; }
};
