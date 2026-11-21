// 연쇄 파괴 (사건-사건 해석)
//
//   하중을 같은 비율로 키워 가며 (하중계수 λ)
//   1. 남은 부재로 푼다. DCR 은 하중에 비례하므로 부재 k 가 용량에 닿는 λ = 1 / DCR_k 이다.
//   2. 가장 먼저 닿는 부재를 제거한다 (좌굴·파단 후 힘을 잃는다고 본다).
//   3. 같은 λ 에서 다시 풀어, 용량을 넘는 부재가 있으면 하중을 올리지 않고 이어서 제거한다 (연쇄).
//   4. 강성행렬이 특이해지면 (메커니즘) 붕괴로 끝낸다.
//   P_ult = (도달한 최대 λ) × |하중 합력|
//
//   동시에 용량에 닿는 부재가 여럿이면 번호가 작은 부재를 먼저 제거하고 다시 판단한다.
#include "model.h"
#include <cmath>
#include <vector>

namespace {

// alive 인 부재만 남긴 트러스. 부재가 하나도 없고 하중도 없는 절점은 뺀다 (남겨 두면 강성이 0이라 메커니즘으로 잘못 판정된다).
// 하중이 걸린 절점이 고립되면 그 하중을 받을 길이 없으므로 false (메커니즘).
bool subTruss(const Truss& t, const std::vector<char>& alive, Truss& s, std::vector<int>& memberMap) {
    const int j = static_cast<int>(t.nodes.size());
    std::vector<int> deg(j, 0);
    for (size_t k = 0; k < t.members.size(); ++k)
        if (alive[k]) { ++deg[t.members[k].i]; ++deg[t.members[k].j]; }
    std::vector<int> nodeMap(j, -1);
    s = Truss{};
    for (int q = 0; q < j; ++q) {
        const Node& nd = t.nodes[q];
        if (deg[q] == 0) {
            if (nd.fx != 0 || nd.fy != 0) return false;
            continue;
        }
        nodeMap[q] = static_cast<int>(s.nodes.size());
        s.nodes.push_back(nd);
    }
    memberMap.clear();
    for (size_t k = 0; k < t.members.size(); ++k) {
        if (!alive[k]) continue;
        Member e = t.members[k];
        e.i = nodeMap[e.i]; e.j = nodeMap[e.j];
        s.members.push_back(e);
        memberMap.push_back(static_cast<int>(k));
    }
    return !s.members.empty();
}

} // namespace

CollapseResult progressiveCollapse(const Truss& t) {
    CollapseResult c;
    double fx = 0, fy = 0;
    for (const Node& q : t.nodes) { fx += q.fx; fy += q.fy; }
    const double Fmag = std::hypot(fx, fy);
    std::vector<char> alive(t.members.size(), 1);
    double lam = 0;
    const double tol = 1e-9;

    for (;;) {
        Truss s;
        std::vector<int> map;
        Result r;
        if (!subTruss(t, alive, s, map) || (r = solve(s)).status != Status::Ok) {
            c.end = c.events.empty() ? "unstable" : "mechanism";
            break;
        }
        // 지금 λ 에서 가장 많이 넘친 부재 (같으면 원래 번호가 작은 부재)
        int worst = -1;
        double wr = 0;
        for (size_t k = 0; k < r.members.size(); ++k) {
            const double d = r.members[k].dcr;
            if (d > wr * (1 + 1e-12) || (worst < 0 && d > 0)) { wr = d; worst = static_cast<int>(k); }
        }
        if (worst < 0) { c.end = "all-removed"; break; }   // 힘을 받는 부재가 없다

        CollapseEvent ev;
        if (lam > 0 && lam * wr > 1 + tol) {
            ev.cascade = true;                // 하중을 올리지 않고 이어진 파괴
        } else {
            lam = std::max(lam, 1 / wr);      // 다음 부재가 용량에 닿는 하중계수
        }
        ev.lam = lam;
        ev.P = lam * Fmag;
        ev.member = map[worst];
        ev.mode = r.members[worst].mode;
        c.events.push_back(ev);
        alive[ev.member] = 0;
        c.lamUlt = std::max(c.lamUlt, lam);
    }
    c.ok = !c.events.empty();
    if (c.ok) {
        c.lamFirst = c.events.front().lam;
        c.Pfirst = c.lamFirst * Fmag;
        c.Pult = c.lamUlt * Fmag;
    }
    return c;
}
