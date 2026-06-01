export module utils.matrix;

import std;


export template<typename T>
struct DenseMatrix {
    std::vector<T> data;
    int rows = 0;
    int cols = 0;

    auto operator[](int i, int j) const -> const T& {
        return data[static_cast<std::size_t>(i) * cols + j];
    }
    auto operator[](int i, int j) -> T& {
        return data[static_cast<std::size_t>(i) * cols + j];
    }

    [[nodiscard]] auto size() const -> std::size_t { return data.size(); }
    [[nodiscard]] auto empty() const -> bool { return data.empty(); }
};
