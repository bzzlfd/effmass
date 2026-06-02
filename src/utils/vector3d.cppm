export module utils.vector3d;

import std;


export template<typename T>
struct vector3d {
    T x{}, y{}, z{};

    vector3d() = default;
    constexpr vector3d(T x_, T y_, T z_) : x(x_), y(y_), z(z_) {}

    [[nodiscard]] auto operator[](int i) const -> const T& {
        return i == 0 ? x : (i == 1 ? y : z);
    }
    [[nodiscard]] auto operator[](int i) -> T& {
        return i == 0 ? x : (i == 1 ? y : z);
    }

    [[nodiscard]] auto data() const -> const T* { return &x; }
    [[nodiscard]] auto data() -> T* { return &x; }

    [[nodiscard]] auto norm_squared() const -> T { return x * x + y * y + z * z; }
    [[nodiscard]] auto norm() const -> T { using std::sqrt; return sqrt(norm_squared()); }

    auto operator+=(const vector3d& rhs) -> vector3d& {
        x += rhs.x; y += rhs.y; z += rhs.z; return *this;
    }
    auto operator-=(const vector3d& rhs) -> vector3d& {
        x -= rhs.x; y -= rhs.y; z -= rhs.z; return *this;
    }
    auto operator*=(T scalar) -> vector3d& {
        x *= scalar; y *= scalar; z *= scalar; return *this;
    }
    auto operator/=(T scalar) -> vector3d& {
        x /= scalar; y /= scalar; z /= scalar; return *this;
    }
};


export template<typename T>
auto operator+(const vector3d<T>& a, const vector3d<T>& b) -> vector3d<T> {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

export template<typename T>
auto operator-(const vector3d<T>& a, const vector3d<T>& b) -> vector3d<T> {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

export template<typename T>
auto operator*(const vector3d<T>& v, T scalar) -> vector3d<T> {
    return {v.x * scalar, v.y * scalar, v.z * scalar};
}

export template<typename T>
auto operator*(T scalar, const vector3d<T>& v) -> vector3d<T> {
    return v * scalar;
}

export template<typename T>
auto operator/(const vector3d<T>& v, T scalar) -> vector3d<T> {
    return {v.x / scalar, v.y / scalar, v.z / scalar};
}

export template<typename T>
auto operator-(const vector3d<T>& v) -> vector3d<T> {
    return {-v.x, -v.y, -v.z};
}

export template<typename T>
auto dot(const vector3d<T>& a, const vector3d<T>& b) -> T {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

export template<typename T>
auto cross(const vector3d<T>& a, const vector3d<T>& b) -> vector3d<T> {
    return {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

export template<typename T>
auto operator==(const vector3d<T>& a, const vector3d<T>& b) -> bool {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}

export template<typename T>
auto operator!=(const vector3d<T>& a, const vector3d<T>& b) -> bool {
    return !(a == b);
}

export template<typename T>
auto operator<<(std::ostream& os, const vector3d<T>& v) -> std::ostream& {
    return os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
}
