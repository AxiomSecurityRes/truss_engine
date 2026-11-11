#include "truss.h"
#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
    const char* inPath  = argc > 1 ? argv[1] : "truss.txt";
    const char* outPath = argc > 2 ? argv[2] : "result.txt";

    std::ifstream in(inPath);
    if (!in) {
        std::cerr << "cannot open " << inPath << " (is it in the working directory?)\n";
        return 1;
    }

    Truss t;
    std::string err;
    if (!readTruss(in, t, &err)) {
        std::cerr << inPath << ": " << err << '\n';
        return 1;
    }
    std::cout << "read " << t.nodes.size() << " nodes, " << t.members.size() << " members\n";

    const Result r = solve(t);

    std::ofstream out(outPath);
    if (!out) {
        std::cerr << "cannot write " << outPath << " (is it open in another program?)\n";
        return 1;
    }
    const bool ok = writeResult(out, r);
    out.close();
    if (!out) {
        std::cerr << "failed while writing " << outPath << '\n';
        return 1;
    }

    writeResult(std::cout, r);
    return ok ? 0 : 1;
}
