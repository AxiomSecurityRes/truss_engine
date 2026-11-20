// 대표 4종 트러스 생성. 도면판(HTML)의 generate() 와 절점·부재 번호가 같다.
//   절점: 하현 B0..Bn 먼저, 그다음 상현(T)·중간(M) 절점
//   지점: B0 핀, Bn 롤러.  하중: config 의 load (center: B(n/2) 에 P, panel: B1..B(n−1) 에 P/(n−1)씩)
#include "model.h"
#include <cctype>
#include <cmath>
#include <map>
#include <string>

const Kind ALL_KINDS[4] = {Kind::Warren, Kind::Pratt, Kind::Howe, Kind::KTruss};

const char* kindName(Kind k) {
    switch (k) {
        case Kind::Warren: return "Warren";
        case Kind::Pratt: return "Pratt";
        case Kind::Howe: return "Howe";
        case Kind::KTruss: return "K-truss";
    }
    return "?";
}

bool parseKind(const std::string& s, Kind& k) {
    std::string t;
    for (char c : s) t += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (t == "warren") { k = Kind::Warren; return true; }
    if (t == "pratt") { k = Kind::Pratt; return true; }
    if (t == "howe") { k = Kind::Howe; return true; }
    if (t == "k" || t == "ktruss" || t == "k-truss") { k = Kind::KTruss; return true; }
    return false;
}

Truss makePreset(Kind kind, int n, double L, double H, const Params& p) {
    Truss t;
    std::map<std::string, int> idx;
    auto X = [&](int k) { return k * L / n; };   // 곱을 먼저 해서 패널 점이 정확히 떨어지게
    auto add = [&](const std::string& name, double x, double y) {
        idx[name] = static_cast<int>(t.nodes.size());
        Node q; q.x = x; q.y = y;
        t.nodes.push_back(q);
    };
    auto e = [&](const std::string& a, const std::string& b) { t.members.push_back(makeMember(idx.at(a), idx.at(b), p)); };
    auto B = [](int k) { return "B" + std::to_string(k); };
    auto T = [](int k) { return "T" + std::to_string(k); };
    auto M = [](int k) { return "M" + std::to_string(k); };
    const int half = n / 2;

    for (int k = 0; k <= n; ++k) add(B(k), X(k), 0);
    if (kind == Kind::Warren) {
        for (int k = 0; k < n; ++k) add(T(k), (2 * k + 1) * L / (2 * n), H);
        for (int k = 0; k < n; ++k) e(B(k), B(k + 1));
        for (int k = 0; k < n - 1; ++k) e(T(k), T(k + 1));
        for (int k = 0; k < n; ++k) { e(T(k), B(k)); e(T(k), B(k + 1)); }
    } else if (kind == Kind::Pratt || kind == Kind::Howe) {
        for (int k = 1; k < n; ++k) add(T(k), X(k), H);
        for (int k = 0; k < n; ++k) e(B(k), B(k + 1));
        for (int k = 1; k < n - 1; ++k) e(T(k), T(k + 1));
        e(B(0), T(1)); e(B(n), T(n - 1));
        for (int k = 1; k < n; ++k) e(T(k), B(k));
        for (int k = 1; k < n - 1; ++k) {
            const bool left = k < half;
            if (kind == Kind::Pratt) { if (left) e(T(k), B(k + 1)); else e(T(k + 1), B(k)); }
            else                     { if (left) e(T(k + 1), B(k)); else e(T(k), B(k + 1)); }
        }
    } else {   // K-truss
        for (int k = 1; k < n; ++k) { add(T(k), X(k), H); add(M(k), X(k), H / 2); }
        for (int k = 0; k < n; ++k) e(B(k), B(k + 1));
        for (int k = 1; k < n - 1; ++k) e(T(k), T(k + 1));
        e(B(0), T(1)); e(B(n), T(n - 1));
        for (int k = 1; k < n; ++k) { e(T(k), M(k)); e(M(k), B(k)); }
        for (int k = 1; k < n - 1; ++k) {
            if (k < half) { e(T(k), M(k + 1)); e(B(k), M(k + 1)); }
            else          { e(T(k + 1), M(k)); e(B(k + 1), M(k)); }
        }
        e(B(0), M(1)); e(B(n), M(n - 1));
    }
    Node& n0 = t.nodes[idx.at(B(0))];
    Node& nn = t.nodes[idx.at(B(n))];
    n0.fixX = 1; n0.fixY = 1; nn.fixY = 1;
    applyLoad(t, n, p);   // 하현 절점은 0..n 번이다
    return t;
}
