
 #ifndef MDHT_CELLMAP_CPP
 #define MDHT_CELLMAP_CPP

#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <tuple>
#include <cstdint>
#include <limits>
#include <chrono>
#include <cstdlib>

// =========================================================
// Dimensions (edit and recompile if needed)
// =========================================================
constexpr int DIM_X = 15;
constexpr int DIM_Y = 360;
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
using CellMap = std::map<unsigned, std::set<Cell>>;

// =========================================================
// Helper to add a cell
// =========================================================
inline void addCell(CellMap& m, unsigned key, int i, int j, int k) {
    Cell c = IndexCodec::compress(i, j, k);
    m[key].insert(c);  // set ensures no duplicates
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

    for (const auto& [key, s] : m) {
        ofs.write(reinterpret_cast<const char*>(&key), sizeof(key));
        uint64_t setSize = s.size();
        ofs.write(reinterpret_cast<const char*>(&setSize), sizeof(setSize));

        for (Cell c : s) {
            ofs.write(reinterpret_cast<const char*>(&c), sizeof(Cell));
        }
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
        uint64_t setSize;
        ifs.read(reinterpret_cast<char*>(&key), sizeof(key));
        ifs.read(reinterpret_cast<char*>(&setSize), sizeof(setSize));

        std::set<Cell> s;
        for (uint64_t j = 0; j < setSize; ++j) {
            Cell c;
            ifs.read(reinterpret_cast<char*>(&c), sizeof(Cell));
            s.insert(c);
        }
        m[key] = std::move(s);
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Read map from " << filename
              << " in " << elapsed.count() << " seconds\n";

    return m;
}

// =========================================================
// Dump first N keys with sizes and cell coordinates
// =========================================================
void dumpFirstKeys(const CellMap& m, size_t N) {
    size_t count = 0;
    for (const auto& [key, s] : m) {
        std::cout << "Key " << key << " has " << s.size() << " cells\n";
        for (Cell c : s) {
            auto [i, j, k] = IndexCodec::decompress(c);
            std::cout << "   (" << i << "," << j << "," << k << ")\n";
        }
        if (++count >= N) break;
    }
}

// =========================================================
// Demo
// =========================================================
#ifdef CELL_MAP_DEMO
int main() {
    IndexCodec::checkDimensions();

    CellMap m;

    // Add some cells (duplicates are ignored)
    addCell(m, 42, 1, 2, 3);
    addCell(m, 42, 10, 5, 0);
    addCell(m, 42, 1, 2, 3);  // duplicate
    addCell(m, 99, 123, 7, 4);

    // Write map to disk
    writeMapBinary(m, "map.bin");

    // Read it back
    CellMap m2 = readMapBinary("map.bin");

    // Print only the first 2 keys
    dumpFirstKeys(m2, 2);
}
#endif


#endif //MDHT_CELLMAP_CPP