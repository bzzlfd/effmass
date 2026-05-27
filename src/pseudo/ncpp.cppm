export module pseudo.ncpp;

import std;
import pseudo.io.ncpp_upf;

export {
    class NCPP;
}


class NCPP {
public:
    explicit NCPP(const NCPPUPF& upf);

    // Placeholder for future physical interfaces
    auto upfData() const -> const NCPPUPF& { return upf_; }

private:
    NCPPUPF upf_;
};


NCPP::NCPP(const NCPPUPF& upf) : upf_(upf) {}
