/*
======================================================================
 DumpMapStats.cpp — Utility for inspecting CellMap binary files
======================================================================

This tool reads a CellMap binary file (e.g. CellMap_30M.dat) and
provides basic statistics and inspection capabilities.

----------------------------------------------------------------------
 Options:
----------------------------------------------------------------------

  --range kMin kMax   
      Restrict processing to keys in the interval [kMin, kMax].
      If kMax = 0, it means “until the end of the map”.
      Example:
          ./DumpMapStats CellMap_30M.dat --range 1000 2000
      Output:
          Processing keys from 1000 to 2000

      Example with kMax=0:
          ./DumpMapStats CellMap_30M.dat --range 5000 0
      Output:
          Processing keys from 5000 to the end

  --dump N
      Print the first N keys (in the selected range) and all their
      associated (i,j,k) triples after decompression.
      Example:
          ./DumpMapStats CellMap_30M.dat --dump 3

  --key KEY
      Print the contents (triples) of a specific key.
      Example:
          ./DumpMapStats CellMap_30M.dat --key 67109121

  --hist file.root   [Requires compilation with -DUSE_ROOT]
      Create a ROOT histogram of set sizes (number of cells per key)
      and save it to file.root.
      Example:
          ./DumpMapStats CellMap_30M.dat --hist sizes.root

  --help
      Show a help message describing all options.

----------------------------------------------------------------------
 Example combined usage:
----------------------------------------------------------------------

  ./DumpMapStats CellMap_30M.dat --range 1000 2000 --dump 2 --hist out.root

    → Reads the map file
    → Processes keys 1000–2000
    → Prints the first 2 keys and their triples
    → Saves histogram of set sizes to out.root

----------------------------------------------------------------------
 Notes:
----------------------------------------------------------------------
 - The map is stored as std::map<unsigned, std::set<Cell>>.
 - Each Cell is a compressed index, decoded into (i,j,k) triples.
 - Average set size across all processed keys is reported at the end.

======================================================================
*/

#include "CellMap.cpp"
#include <iostream>
#include <string>
#include <vector>

#ifdef USE_ROOT
#include "TH1D.h"
#include "TFile.h"
#endif

// ------------------------------------------------------------------
// Print runtime help
// ------------------------------------------------------------------
void printHelp(const char* progName) {
    std::cout << "Usage: " << progName << " <mapfile> [options...]\n\n"
              << "Options:\n"
              << "  --range kMin kMax   Restrict processed keys. kMax=0 → until end.\n"
              << "  --dump N            Print first N keys in the range.\n"
              << "  --key KEY           Print a specific key and its triples.\n"
#ifdef USE_ROOT
              << "  --hist file.root    Save histogram of set sizes (ROOT required).\n"
#endif
              << "  --help              Show this message.\n\n"
              << "Examples:\n"
              << "  " << progName << " CellMap_30M.dat --range 1000 2000 --dump 2\n"
              << "  " << progName << " CellMap_30M.dat --key 67109121\n"
#ifdef USE_ROOT
              << "  " << progName << " CellMap_30M.dat --hist sizes.root\n"
#endif
              << std::endl;
}

// ------------------------------------------------------------------
// Main program
// ------------------------------------------------------------------
int main(int argc, char** argv) {
    if (argc < 2) {
        printHelp(argv[0]);
        return 1;
    }

    std::string filename;
    size_t dumpN = 0;
    unsigned keyToPrint = 0;
    unsigned kMin = 0, kMax = 0;
#ifdef USE_ROOT
    std::string histFile;
#endif

    filename = argv[1];
    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help") {
            printHelp(argv[0]);
            return 0;
        } else if (arg == "--dump" && i + 1 < argc) {
            dumpN = std::stoul(argv[++i]);
        } else if (arg == "--key" && i + 1 < argc) {
            keyToPrint = std::stoul(argv[++i]);
        } else if (arg == "--range" && i + 2 < argc) {
            kMin = std::stoul(argv[++i]);
            kMax = std::stoul(argv[++i]);
        }
#ifdef USE_ROOT
        else if (arg == "--hist" && i + 1 < argc) {
            histFile = argv[++i];
        }
#endif
        else {
            std::cerr << "Unknown option: " << arg << "\n";
            printHelp(argv[0]);
            return 1;
        }
    }

    std::cout << "Reading map from " << filename << "\n";
    CellMap m = readMapBinary(filename);

    std::cout << "Processing keys from " << kMin;
    if (kMax == 0)
        std::cout << " to the end\n";
    else
        std::cout << " to " << kMax << "\n";

    size_t totalSize = 0;
    size_t keyCount = 0;
    size_t printed = 0;

#ifdef USE_ROOT
    TH1D* hSizes = nullptr;
    if (!histFile.empty()) {
        hSizes = new TH1D("hSizes", "Set sizes;size;count", 100, 0, 1000);
    }
#endif

    for (const auto& [key, s] : m) {
        if (key < kMin) continue;
        if (kMax != 0 && key > kMax) break;

        totalSize += s.size();
        keyCount++;

#ifdef USE_ROOT
        if (hSizes) hSizes->Fill(s.size());
#endif

        if (dumpN > 0 && printed < dumpN) {
            std::cout << "Key " << key << " size=" << s.size() << ":\n";
            for (Cell c : s) {
                auto [i, j, k] = IndexCodec::decompress(c);
                std::cout << "   (" << i << "," << j << "," << k << ")\n";
            }
            printed++;
        }

        if (keyToPrint != 0 && key == keyToPrint) {
            std::cout << "Key " << key << " size=" << s.size() << ":\n";
            for (Cell c : s) {
                auto [i, j, k] = IndexCodec::decompress(c);
                std::cout << "   (" << i << "," << j << "," << k << ")\n";
            }
        }
    }

    if (keyCount > 0) {
        double avgSize = double(totalSize) / keyCount;
        std::cout << "Total keys processed: " << keyCount << "\n";
        std::cout << "Average set size: " << avgSize << "\n";
    } else {
        std::cout << "No keys processed.\n";
    }

#ifdef USE_ROOT
    if (hSizes) {
        TFile fout(histFile.c_str(), "RECREATE");
        hSizes->Write();
        fout.Close();
        std::cout << "Histogram written to " << histFile << "\n";
    }
#endif

    return 0;
}
