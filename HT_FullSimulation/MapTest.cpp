#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <tuple>
#include <cstdint>
#include <limits>
#include <chrono>
#include <cstdlib>

// =========================================================
// Dimensions (edit and recompile if needed)
// =========================================================
constexpr int DIM_X = 360;
constexpr int DIM_Y = 15;
constexpr int DIM_Z = 6;

// =========================================================
// Index compression utilities
// =========================================================
namespace IndexCodec {

    constexpr int DIM_YZ = DIM_Y * DIM_Z;

    // Safety check — call once at program startup
    inline void checkDimensions() {
        constexpr uint64_t total =
            static_cast<uint64_t>(DIM_X) *
            static_cast<uint64_t>(DIM_Y) *
            static_cast<uint64_t>(DIM_Z);

        if (total > std::numeric_limits<uint32_t>::max()) {
            std::cerr << "Error: dimensions (" << DIM_X << " x " << DIM_Y
                      << " x " << DIM_Z << ") = " << total
                      << " exceed uint32_t capacity ("
                      << std::numeric_limits<uint32_t>::max() << ")\n";
            std::exit(EXIT_FAILURE);
        }
    }

    inline uint32_t compress(int i, int j, int k) {
        return static_cast<uint32_t>(i * DIM_YZ + j * DIM_Z + k);
    }

    inline std::tuple<int, int, int> decompress(uint32_t id) {
        int i = id / DIM_YZ;
        int j = (id / DIM_Z) % DIM_Y;
        int k = id % DIM_Z;
        return {i, j, k};
    }

} // namespace IndexCodec

// =========================================================
// Map definition
// =========================================================
using Cell    = uint32_t;
using CellMap = std::map<unsigned, std::vector<Cell>>;

// =========================================================
// Helper to add a cell
// =========================================================
inline void addCell(CellMap& m, unsigned key, int i, int j, int k) {
    Cell c = IndexCodec::compress(i, j, k);
    m[key].push_back(c);
}

// =========================================================
// Binary I/O with timing
// =========================================================
void writeMapBinary(const CellMap& m, const std::string& filename) {
    auto start = std::chrono::high_resolution_clock::now();

    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs) throw std::runtime_error("Cannot open file for writing");

    uint64_t mapSize = m.size();
    ofs.write(reinterpret_cast<const char*>(&mapSize), sizeof(mapSize));

    for (const auto& [key, vec] : m) {
        ofs.write(reinterpret_cast<const char*>(&key), sizeof(key));
        uint64_t vecSize = vec.size();
        ofs.write(reinterpret_cast<const char*>(&vecSize), sizeof(vecSize));
        ofs.write(reinterpret_cast<const char*>(vec.data()), vecSize * sizeof(Cell));
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Wrote map to " << filename
              << " in " << elapsed.count() << " seconds\n";
}

CellMap readMapBinary(const std::string& filename) {
    auto start = std::chrono::high_resolution_clock::now();

    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs) throw std::runtime_error("Cannot open file for reading");

    CellMap m;

    uint64_t mapSize = 0;
    ifs.read(reinterpret_cast<char*>(&mapSize), sizeof(mapSize));

    for (uint64_t i = 0; i < mapSize; ++i) {
        unsigned key;
        uint64_t vecSize;
        ifs.read(reinterpret_cast<char*>(&key), sizeof(key));
        ifs.read(reinterpret_cast<char*>(&vecSize), sizeof(vecSize));

        std::vector<Cell> vec(vecSize);
        ifs.read(reinterpret_cast<char*>(vec.data()), vecSize * sizeof(Cell));
        m[key] = std::move(vec);
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Read map from " << filename
              << " in " << elapsed.count() << " seconds\n";

    return m;
}

// =========================================================
// Demo
// =========================================================
//#include <iostream>
//#include <string>
//#include "CellMap.cpp"   // or better: split into CellMap.h / CellMap.cpp

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <CellMap.dat> <N>\n";
        return 1;
    }

    std::string filename = argv[1];
    int N = std::stoi(argv[2]);

    IndexCodec::checkDimensions();

    // Load map
    CellMap m = readMapBinary(filename);

    std::cout << "Map contains " << m.size() << " keys\n";
    std::cout << "Printing first " << N << " keys:\n";

    int count = 0;
    for (const auto& [key, vec] : m) {
        std::cout << "Key " << key << " -> " << vec.size() << " cells\n";
        if (++count >= N) break;
    }

    return 0;
}
