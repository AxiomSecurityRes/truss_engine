#include "truss.h"
#include <string>

// 도면판 TXT 형식을 읽는다. 형식이 틀리면 false.
bool readTruss(std::istream& in, Truss& t) {
    std::string tag;
    int version = 0, j = 0, m = 0;
    if (!(in >> tag >> version >> j >> m)) return false;
    if (tag != "TRUSS" || version != 1 || j < 0 || m < 0) return false;

    t.nodes.assign(j, Node{});
    t.members.assign(m, Member{});

    for (Node& n : t.nodes) {
        if (!(in >> n.x >> n.y >> n.fixX >> n.fixY >> n.fx >> n.fy)) return false;
        if ((n.fixX != 0 && n.fixX != 1) || (n.fixY != 0 && n.fixY != 1)) return false;
    }
    for (Member& e : t.members) {
        if (!(in >> e.i >> e.j >> e.E >> e.A >> e.I >> e.rho >> e.sigmaT >> e.sigmaC >> e.K)) return false;
        if (e.i < 0 || e.i >= j || e.j < 0 || e.j >= j || e.i == e.j) return false;
    }
    return true;
}

// 결과를 RESULT 1 형식으로 쓴다. 17자리면 double이 손실 없이 넘어간다.
void writeResult(std::ostream& out, const Result& r) {
    const auto oldPrecision = out.precision(17);

    out << "RESULT 1\n";
    out << static_cast<int>(r.status) << ' ' << r.ux.size() << ' ' << r.members.size() << '\n';
    for (std::size_t k = 0; k < r.ux.size(); ++k)
        out << r.ux[k] << ' ' << r.uy[k] << ' ' << r.rx[k] << ' ' << r.ry[k] << '\n';
    for (const MemberResult& e : r.members)
        out << e.N << ' ' << e.stress << ' ' << e.cap << ' ' << e.dcr << ' ' << static_cast<int>(e.mode) << '\n';
    out << r.loadFactor << ' ' << r.Pfail << ' ' << r.critical << ' '
        << r.mass << ' ' << r.eta << ' ' << r.residual << '\n';

    out.precision(oldPrecision);
}
