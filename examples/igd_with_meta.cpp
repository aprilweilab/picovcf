/* Example loading an IGD file along with numeric metadata.
 */
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "picovcf.hpp"

using namespace picovcf;

template <typename T> T lineToValue(const std::string& line);

template <> std::string lineToValue(const std::string& line) { return line; }

template <> uint64_t lineToValue(const std::string& line) {
    char* endPtr = nullptr;
    auto result = static_cast<uint64_t>(std::strtoull(line.c_str(), &endPtr, 10));
    if (endPtr == line.c_str()) {
        throw std::runtime_error("Could not parse integer");
    }
    return result;
}

template <typename T> std::vector<T> readMeta(const std::string& filename) {
    std::vector<T> result;
    std::ifstream metaTextFile(filename);
    std::string line;
    while (std::getline(metaTextFile, line)) {
        // Skip comments
        if (!line.empty() && line[0] == '#') {
            continue;
        }
        result.emplace_back(lineToValue<T>(line));
    }
    return std::move(result);
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: igd_with_meta <igd file> <numerical meta file>" << std::endl;
        return 1;
    }

    IGDData igd(argv[1]);

    // Assume integer-base metadata for this example. Change to readMeta<std::string>() to
    // load arbitrary strings.
    std::vector<uint64_t> metadata = readMeta<uint64_t>(argv[2]);

    for (size_t i = 0; i < igd.numVariants(); i++) {
        const uint64_t position = igd.getPosition(i);
        std::cout << "Variant=" << i << ", Position=" << position << ", Metadata=" << metadata.at(i) << std::endl;
    }

    return 0;
}
