#include "truss.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <utility>

namespace {

// 메모장 등이 파일 앞에 붙이는 UTF-8 BOM(EF BB BF)을 건너뛴다.
void skipUtf8Bom(std::istream& in) {
    if (in.peek() != 0xEF) return;
    char bom[3];
    in.read(bom, 3);
    if (in.gcount() != 3 || bom[1] != '\xBB' || bom[2] != '\xBF') in.setstate(std::ios::failbit);
}

bool fail(std::string* err, const std::string& msg) {
    if (err) *err = msg;
    return false;
}

bool positiveFinite(double v) { return std::isfinite(v) && v > 0; }

} // namespace

bool readTruss(std::istream& in, Truss& t, std::string* err) {
    t = Truss{};
    skipUtf8Bom(in);

    std::string tag;
    int version = 0;
    if (!(in >> tag >> version)) return fail(err, "header: expected 'TRUSS 1'");
    if (tag != "TRUSS" || version != 1)
        return fail(err, "header: expected 'TRUSS 1', got '" + tag + " " + std::to_string(version) + "'");

    long long j = 0, m = 0;
    if (!(in >> j >> m)) return fail(err, "line 2: expected node count j and member count m");
    if (j < 0 || m < 0) return fail(err, "line 2: j and m must be 0 or more");
    constexpr long long kMaxCount = 100000;   // 조밀 강성행렬로는 이 규모도 이미 감당하기 어렵다
    if (j > kMaxCount || m > kMaxCount)
        return fail(err, "line 2: j and m must be at most " + std::to_string(kMaxCount));

    t.nodes.resize(static_cast<std::size_t>(j));
    t.members.resize(static_cast<std::size_t>(m));

    for (long long k = 0; k < j; ++k) {
        Node& n = t.nodes[k];
        const std::string where = "node " + std::to_string(k);
        if (!(in >> n.x >> n.y >> n.fixX >> n.fixY >> n.fx >> n.fy))
            return fail(err, where + ": expected 6 values (x y fix_x fix_y fx fy)");
        if (!std::isfinite(n.x) || !std::isfinite(n.y) || !std::isfinite(n.fx) || !std::isfinite(n.fy))
            return fail(err, where + ": x, y, fx, fy must be finite numbers");
        if ((n.fixX != 0 && n.fixX != 1) || (n.fixY != 0 && n.fixY != 1))
            return fail(err, where + ": fix_x and fix_y must be 0 or 1");
    }

    std::set<std::pair<int, int>> seen;
    for (long long k = 0; k < m; ++k) {
        Member& e = t.members[k];
        const std::string where = "member " + std::to_string(k);
        if (!(in >> e.i >> e.j >> e.E >> e.A >> e.I >> e.rho >> e.sigmaT >> e.sigmaC >> e.K))
            return fail(err, where + ": expected 9 values (i j E A I rho sigma_t sigma_c K)");
        if (e.i < 0 || e.i >= j || e.j < 0 || e.j >= j)
            return fail(err, where + ": node index out of range 0.." + std::to_string(j - 1));
        if (e.i == e.j)
            return fail(err, where + ": both ends are node " + std::to_string(e.i));
        const std::pair<int, int> key = std::minmax(e.i, e.j);
        if (!seen.insert(key).second)
            return fail(err, where + ": duplicates an earlier member between nodes "
                             + std::to_string(key.first) + " and " + std::to_string(key.second));
        const Node& a = t.nodes[e.i];
        const Node& b = t.nodes[e.j];
        if (a.x == b.x && a.y == b.y)
            return fail(err, where + ": zero length (nodes " + std::to_string(e.i) + " and "
                             + std::to_string(e.j) + " are at the same point)");
        if (!positiveFinite(e.E) || !positiveFinite(e.A) || !positiveFinite(e.I) || !positiveFinite(e.rho)
            || !positiveFinite(e.sigmaT) || !positiveFinite(e.sigmaC) || !positiveFinite(e.K))
            return fail(err, where + ": E, A, I, rho, sigma_t, sigma_c, K must all be positive numbers");
    }

    in >> std::ws;
    if (!in.eof())
        return fail(err, "extra data after member " + std::to_string(m - 1) + " (check j and m on line 2)");
    return true;
}

bool writeResult(std::ostream& out, const Result& r) {
    const std::size_t j = r.ux.size();
    const bool consistent = r.uy.size() == j && r.rx.size() == j && r.ry.size() == j;

    const auto oldPrecision = out.precision(17);   // 17자리면 double이 손실 없이 넘어간다
    out << "RESULT 1\n";
    if (!consistent) {
        // solve가 배열 크기를 잘못 채운 경우: 읽을 수 있는 실패 결과를 남긴다.
        out << static_cast<int>(Status::InvalidInput) << " 0 0\n0 0 -1 0 0 0\n";
        out.precision(oldPrecision);
        std::cerr << "writeResult: ux, uy, rx, ry sizes differ (" << r.ux.size() << ", " << r.uy.size()
                  << ", " << r.rx.size() << ", " << r.ry.size() << ")\n";
        return false;
    }

    out << static_cast<int>(r.status) << ' ' << j << ' ' << r.members.size() << '\n';
    for (std::size_t k = 0; k < j; ++k)
        out << r.ux[k] << ' ' << r.uy[k] << ' ' << r.rx[k] << ' ' << r.ry[k] << '\n';
    for (const MemberResult& e : r.members)
        out << e.N << ' ' << e.stress << ' ' << e.cap << ' ' << e.dcr << ' ' << static_cast<int>(e.mode) << '\n';
    out << r.loadFactor << ' ' << r.Pfail << ' ' << r.critical << ' '
        << r.mass << ' ' << r.eta << ' ' << r.residual << '\n';

    out.precision(oldPrecision);
    return true;
}
