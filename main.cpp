#include "truss.h"
#include <fstream>
#include <iostream>

int main(int argc, char* argv[]) {
    const char* inPath  = argc > 1 ? argv[1] : "truss.txt";
    const char* outPath = argc > 2 ? argv[2] : "result.txt";

    std::ifstream in(inPath);
    if (!in) { std::cerr << "cannot open " << inPath << '\n'; return 1; }

    Truss t;
    if (!readTruss(in, t)) { std::cerr << "bad input format\n"; return 1; }
    std::cout << "read " << t.nodes.size() << " nodes, " << t.members.size() << " members\n";

    Result r = solve(t);

    std::ofstream out(outPath);
    writeResult(out, r);
    writeResult(std::cout, r);
    return 0;
}
