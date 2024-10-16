/* Testing tool that just outputs the GZIPd data as uncompressed.
 *
 * Usage:
 *  gzcat <file>
 */
#include <iostream>

#include "picovcf.hpp"

using namespace picovcf;

int main(int argc, char *argv[]) {
  if (argc < 1) {
    std::cerr << "Please pass in an input file" << std::endl;
    return 1;
  }

  const std::string filename(argv[1]);
  ZBufferedReader reader(filename, 128 * 1024);
  std::string line;
  while (reader.readline(line) > 0) {
    std::cout << line << std::endl;
  }

  return 0;
}