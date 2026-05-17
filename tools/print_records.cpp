#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

struct Record {
    std::int64_t length;
};

template<typename T>
auto full_pass(const std::string& filename,
               std::int64_t file_size,
               std::vector<Record>& records) -> bool
{
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open())
        return false;

    records.clear();

    while (true) {
        T lead{};
        file.read(reinterpret_cast<char*>(&lead), sizeof(lead));
        if (!file) {
            if (file.eof() && file.gcount() == 0)
                break;
            return false;
        }

        if (lead < 0)
            return false;

        auto pos = file.tellg();
        if (pos == std::streampos(-1))
            return false;

        auto bytes_left = file_size - static_cast<std::int64_t>(pos);
        if (bytes_left < static_cast<std::int64_t>(sizeof(T)))
            return false;
        if (lead > bytes_left - static_cast<std::int64_t>(sizeof(T)))
            return false;

        file.seekg(lead, std::ios::cur);
        if (!file)
            return false;

        T trail{};
        file.read(reinterpret_cast<char*>(&trail), sizeof(trail));
        if (!file)
            return false;

        if (trail != lead)
            return false;

        records.push_back({static_cast<std::int64_t>(lead)});
    }

    return true;
}


auto main(int argc, char* argv[]) -> int {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <fortran_unformatted_file>\n";
        return 1;
    }

    const std::string filename = argv[1];

    std::ifstream size_file(filename, std::ios::binary | std::ios::ate);
    if (!size_file.is_open()) {
        std::cerr << "Error: cannot open file '" << filename << "'\n";
        return 1;
    }
    const auto file_size = static_cast<std::int64_t>(size_file.tellg());
    size_file.close();

    if (file_size < 8) {
        std::cerr << "Error: file is too small (" << file_size << " bytes)\n";
        return 1;
    }

    std::vector<Record> records_int32, records_int64;
    bool int32_ok = full_pass<std::int32_t>(filename, file_size, records_int32);
    bool int64_ok = full_pass<std::int64_t>(filename, file_size, records_int64);

    int marker_size = 0;
    std::vector<Record>* records = nullptr;

    if (int32_ok && !int64_ok) {
        marker_size = 4;
        records = &records_int32;
    } else if (int64_ok && !int32_ok) {
        marker_size = 8;
        records = &records_int64;
    } else if (!int32_ok && !int64_ok) {
        std::cerr << "Error: cannot detect record marker size"
                     " (neither int32 nor int64 markers passed full verification)\n";
        return 1;
    } else {
        std::cerr << "Error: both int32 and int64 markers passed verification,"
                     " cannot determine which is correct\n";
        return 1;
    }

    std::cout << "File: " << filename << '\n';
    std::cout << "Marker size: " << marker_size << " bytes\n";

    std::int64_t total_data_bytes = 0;
    for (std::size_t i = 0; i < records->size(); ++i) {
        total_data_bytes += (*records)[i].length;
        std::cout << "record " << (i + 1) << ": " << (*records)[i].length
                  << " bytes\n";
    }

    auto nrec = static_cast<std::int64_t>(records->size());
    auto total_overhead = nrec * marker_size * 2;
    auto computed_file_size = total_data_bytes + total_overhead;

    std::cout << "\nsummary: " << nrec << " records\n";
    std::cout << "total data bytes: " << total_data_bytes << '\n';
    std::cout << "file size (actual): " << file_size << '\n';
    if (computed_file_size != file_size) {
        std::cout << "file size (from records): " << computed_file_size << '\n';
        std::cout << "difference: " << (file_size - computed_file_size)
                  << " bytes\n";
    }

    return 0;
}
