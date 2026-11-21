// 제작 가능성 검사 (빨대 부재 + 이쑤시개 핀 접합)
//
//   1. 부재 길이 ≤ L_max      빨대 한 개로 만들 수 있어야 한다
//   2. 부재 길이 ≥ L_min      너무 짧으면 핀 두 개를 꽂을 수 없다
//   3. 절점 차수 ≤ max_deg    한 핀에 꽂을 수 있는 빨대 수
//   4. 절점 사이각 ≥ min_angle 빨대 끝이 서로 겹치지 않아야 한다
//   5. 부재 교차 없음          평면 모형에서 두 빨대가 절점 없이 겹칠 수 없다
//
//   위반량 = Σ(초과 길이 / L_max) + Σ(부족 길이 / L_min) + Σ(초과 차수)
//          + Σ(부족 각 / min_angle) + 교차 쌍 수        (0 이면 제작 가능)
#include "model.h"
#include <algorithm>
#include <cmath>
#include <vector>

namespace {

double cross(double ax, double ay, double bx, double by) { return ax * by - ay * bx; }

// 두 선분이 끝점이 아닌 곳에서 만나거나 겹치는가
bool segmentsClash(const Node& p1, const Node& p2, const Node& q1, const Node& q2, double scale) {
    const double eps = 1e-9 * scale * scale;
    const double d1 = cross(p2.x - p1.x, p2.y - p1.y, q1.x - p1.x, q1.y - p1.y);
    const double d2 = cross(p2.x - p1.x, p2.y - p1.y, q2.x - p1.x, q2.y - p1.y);
    const double d3 = cross(q2.x - q1.x, q2.y - q1.y, p1.x - q1.x, p1.y - q1.y);
    const double d4 = cross(q2.x - q1.x, q2.y - q1.y, p2.x - q1.x, p2.y - q1.y);
    if (((d1 > eps && d2 < -eps) || (d1 < -eps && d2 > eps)) && ((d3 > eps && d4 < -eps) || (d3 < -eps && d4 > eps)))
        return true;   // 내부에서 교차
    // 한 직선 위에 있으면 겹치는 길이가 있는지 본다
    if (std::fabs(d1) <= eps && std::fabs(d2) <= eps) {
        const double dx = p2.x - p1.x, dy = p2.y - p1.y, len2 = dx * dx + dy * dy;
        double a = ((q1.x - p1.x) * dx + (q1.y - p1.y) * dy) / len2;
        double b = ((q2.x - p1.x) * dx + (q2.y - p1.y) * dy) / len2;
        if (a > b) std::swap(a, b);
        return std::min(b, 1.0) - std::max(a, 0.0) > 1e-9;
    }
    // 한 부재의 끝점이 다른 부재의 중간에 닿는 경우 (절점이 없는 T자 접촉)
    auto onSeg = [&](const Node& r, const Node& s1, const Node& s2, double d) {
        if (std::fabs(d) > eps) return false;
        const double dx = s2.x - s1.x, dy = s2.y - s1.y;
        const double u = ((r.x - s1.x) * dx + (r.y - s1.y) * dy) / (dx * dx + dy * dy);
        return u > 1e-9 && u < 1 - 1e-9;
    };
    return onSeg(q1, p1, p2, d1) || onSeg(q2, p1, p2, d2) || onSeg(p1, q1, q2, d3) || onSeg(p2, q1, q2, d4);
}

} // namespace

FabReport checkFab(const Truss& t, const Params& p) {
    FabReport f;
    const int j = static_cast<int>(t.nodes.size());
    const int m = static_cast<int>(t.members.size());
    const double pi = std::acos(-1.0);
    double scale = 1;
    for (const Node& q : t.nodes) scale = std::max({scale, std::fabs(q.x), std::fabs(q.y)});

    f.minLen = m ? 1e300 : 0;
    std::vector<std::vector<double>> ang(j);
    for (const Member& e : t.members) {
        const Node& a = t.nodes[e.i];
        const Node& b = t.nodes[e.j];
        const double L = std::hypot(b.x - a.x, b.y - a.y);
        f.maxLen = std::max(f.maxLen, L);
        f.minLen = std::min(f.minLen, L);
        if (L > p.Lmax * (1 + 1e-12)) { ++f.tooLong; f.violation += (L - p.Lmax) / p.Lmax; }
        if (L < p.Lmin * (1 - 1e-12)) { ++f.tooShort; f.violation += (p.Lmin - L) / p.Lmin; }
        ang[e.i].push_back(std::atan2(b.y - a.y, b.x - a.x));
        ang[e.j].push_back(std::atan2(a.y - b.y, a.x - b.x));
    }
    for (int q = 0; q < j; ++q) {
        const int d = static_cast<int>(ang[q].size());
        f.maxDeg = std::max(f.maxDeg, d);
        if (d > p.maxDeg) { ++f.overDeg; f.violation += d - p.maxDeg; }
        if (d < 2) continue;
        std::vector<double>& a = ang[q];
        std::sort(a.begin(), a.end());
        for (int k = 0; k < d; ++k) {
            double gap = (k + 1 < d ? a[k + 1] : a[0] + 2 * pi) - a[k];
            const double deg = gap * 180 / pi;
            f.minAngle = std::min(f.minAngle, deg);
            if (deg < p.minAngle - 1e-9) { ++f.sharp; f.violation += (p.minAngle - deg) / p.minAngle; }
        }
    }
    if (p.noCross) {
        for (int a = 0; a < m; ++a)
            for (int b = a + 1; b < m; ++b) {
                const Member& e1 = t.members[a];
                const Member& e2 = t.members[b];
                const bool share = e1.i == e2.i || e1.i == e2.j || e1.j == e2.i || e1.j == e2.j;
                if (share) continue;   // 끝을 공유하면 겹침은 사이각 검사가 잡는다
                if (segmentsClash(t.nodes[e1.i], t.nodes[e1.j], t.nodes[e2.i], t.nodes[e2.j], scale)) { ++f.crossings; f.violation += 1; }
            }
    }
    return f;
}
