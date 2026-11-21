// 유전 알고리즘: 트러스의 위상(어느 절점·부재를 쓸지)과 형상(상현 절점 높이)을 함께 찾는다.
//
// [설계 공간 — 바탕 구조 (ground structure)]
//   하현: B0..Bn (간격 L/n), 항상 있음. B0 핀, Bn 롤러, B(n/2) 에 하중 P.
//   상현 후보 절점: 반 패널 간격의 슬롯 s = 1 .. 2n−1, x = s·L/(2n).
//     슬롯 s 는 하현 절점 중 수평 거리가 한 패널(L/n) 이하인 것과만 이을 수 있다.
//       s 홀수 (반 패널 점): 양옆 하현 절점 2개
//       s 짝수 (패널 점 k = s/2): B(k−1), B(k), B(k+1) 3개
//     쓰는 상현 절점끼리는 x 순서대로 상현재로 잇는다.
//   좌우 대칭: 하중·지점이 대칭이고, 좌우를 같게 만들면 제작이 쉬우므로
//     슬롯 s 와 2n−s 는 같은 유전자를 쓴다 (왼쪽 절반 + 중앙만 유전자로 둔다).
//   이 공간에는 Warren, Pratt, Howe 와 수직재 Warren, 높이가 다른 상현 (Parker 꼴) 등이 모두 들어 있다.
//   K-truss 는 중간 높이 절점이 있어 이 공간 밖이다 (따로 비교한다).
//
// [유전자]  act[s] (1 bit), h[s] (높이 단계), web[g] (1 bit)
//   n = 6: act 6 bit + h 6개 + web 14 bit
//
// [평가]  해독 → 제작 제약 검사(fab.cpp) → 직접강성법 해석(solve.cpp)
//   적합도 = η = P_fail / 질량  (fitness = etaUlt 이면 P_ult / 질량)
//   제약 처리: Deb 규칙 (벌점 계수 없음)
//     ① 제약을 지킨 해는 어긴 해보다 항상 낫다
//     ② 둘 다 지키면 적합도가 큰 쪽
//     ③ 둘 다 어기면 위반량이 작은 쪽  (메커니즘은 위반량 1 + 부족한 부재 수)
//
// [연산자]  토너먼트 선택, 슬롯 단위 교차 (slot_cross = 0 이면 유전자별 균등 교차), 비트 반전 / 높이 ±1 단계 또는 새로 뽑기,
//           엘리트 보존, 25세대마다 국소 탐색, 정체 시 재시작
// [난수]    splitmix64 를 직접 구현 — 컴파일러가 달라도 같은 seed 면 같은 결과가 나온다
#include "model.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>
#include <sstream>
#include <unordered_map>

namespace {

struct Rng {
    unsigned long long& s;
    unsigned long long next() {
        unsigned long long z = (s += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }
    double uni() { return static_cast<double>(next() >> 11) * (1.0 / 9007199254740992.0); }   // [0, 1)
    int below(int k) { return static_cast<int>(next() % static_cast<unsigned long long>(k)); }
};

std::string key(const Genome& g) {
    std::string s(g.act.begin(), g.act.end());
    for (int v : g.h) s += static_cast<char>(v & 0xFF), s += static_cast<char>((v >> 8) & 0xFF);
    s.append(g.web.begin(), g.web.end());
    return s;
}

int clampInt(int v, int lo, int hi) { return v < lo ? lo : v > hi ? hi : v; }   // std::clamp 는 C++17 이라 직접 둔다

// a 가 b 보다 나은가 (Deb 규칙). 탐색에는 부드러운 점수(score), 보고에는 참 적합도(fit) 를 쓴다
bool better(const Eval& a, const Eval& b) {
    if (a.feasible != b.feasible) return a.feasible;
    if (a.feasible) return a.score > b.score;
    return a.viol < b.viol;
}
bool betterTrue(const Eval& a, const Eval& b) {
    if (a.feasible != b.feasible) return a.feasible;
    if (a.feasible) return a.fit > b.fit;
    return a.viol < b.viol;
}

} // namespace

GaSpace::GaSpace(const Params& p) : p_(p), n_(p.n) {
    levels_ = static_cast<int>(std::floor((p.gaHmax - p.gaHmin) / p.gaHstep + 1e-9)) + 1;
    const int n = n_, h = n / 2;
    for (int s = 1; s <= n; ++s) {
        if (s == n) {   // 중앙 슬롯 (패널 점 n/2): 수직재 1개, 양쪽 사재 한 쌍
            genes_.push_back({s, {{s, h}}});
            genes_.push_back({s, {{s, h - 1}, {s, h + 1}}});
        } else if (s % 2) {
            for (int k : {(s - 1) / 2, (s + 1) / 2}) genes_.push_back({s, {{s, k}, {2 * n - s, n - k}}});
        } else {
            const int k = s / 2;
            for (int kk : {k - 1, k, k + 1}) genes_.push_back({s, {{s, kk}, {2 * n - s, n - kk}}});
        }
    }
}

double GaSpace::log10Size() const {
    return (n_ + webGenes()) * std::log10(2.0) + n_ * std::log10(static_cast<double>(levels_));
}

double GaSpace::height(int idx) const { return p_.gaHmin + idx * p_.gaHstep; }

Truss GaSpace::decode(const Genome& g) const {
    const int n = n_;
    // 슬롯이 실제로 쓰이려면 act 가 1 이고 웨브 부재가 하나 이상 있어야 한다 (상현재에만 매달린 절점은 쓸모없는 메커니즘)
    std::vector<char> hasWeb(n + 1, 0);
    for (size_t q = 0; q < genes_.size(); ++q)
        if (g.web[q]) hasWeb[genes_[q].slot] = 1;
    auto base = [&](int s) { return s <= n ? s : 2 * n - s; };
    auto active = [&](int s) { const int b = base(s); return g.act[b - 1] && hasWeb[b]; };

    Truss t;
    for (int k = 0; k <= n; ++k) { Node q; q.x = k * p_.L / n; t.nodes.push_back(q); }
    std::vector<int> slotNode(2 * n, -1);
    for (int s = 1; s <= 2 * n - 1; ++s) {
        if (!active(s)) continue;
        Node q;
        q.x = s * p_.L / (2 * n);
        q.y = height(g.h[base(s) - 1]);
        slotNode[s] = static_cast<int>(t.nodes.size());
        t.nodes.push_back(q);
    }
    t.nodes[0].fixX = 1; t.nodes[0].fixY = 1;
    t.nodes[n].fixY = 1;
    applyLoad(t, n, p_);

    for (int k = 0; k < n; ++k) t.members.push_back(makeMember(k, k + 1, p_));         // 하현재
    int prev = -1;
    for (int s = 1; s <= 2 * n - 1; ++s) {                                            // 상현재
        if (slotNode[s] < 0) continue;
        if (prev >= 0) t.members.push_back(makeMember(slotNode[prev], slotNode[s], p_));
        prev = s;
    }
    std::set<std::pair<int, int>> seen;
    for (size_t q = 0; q < genes_.size(); ++q) {                                     // 웨브
        if (!g.web[q] || !active(genes_[q].slot)) continue;
        for (const auto& link : genes_[q].links) {
            const int s = link.first, k = link.second;
            const int a = slotNode[s];
            if (a < 0 || !seen.insert({a, k}).second) continue;
            t.members.push_back(makeMember(k, a, p_));
        }
    }
    return t;
}

Eval GaSpace::evaluate(const Genome& g) const {
    Eval e;
    const Truss t = decode(g);
    e.j = static_cast<int>(t.nodes.size());
    e.m = static_cast<int>(t.members.size());
    e.deg = maxwellDegree(t);
    const FabReport f = checkFab(t, p_);
    const Result r = solve(t);
    e.status = static_cast<int>(r.status);
    e.viol = f.violation;
    if (r.status != Status::Ok || !(r.loadFactor > 0) || !std::isfinite(r.loadFactor)) {
        e.viol += 1 + std::max(0, -e.deg);
        return e;
    }
    e.eta = r.eta; e.Pfail = r.Pfail; e.mass = r.mass;
    e.feasible = e.viol == 0;
    e.fit = e.eta;
    if (e.feasible && p_.fitness == "etaUlt") {
        const CollapseResult c = progressiveCollapse(t);
        e.Pult = c.Pult;
        e.fit = c.Pult / r.mass;
    }
    // 탐색용 점수: η 는 DCR 최댓값 하나로만 정해져서, 거의 같은 힘을 받는 부재가 여럿이면 (활 모양 상현재처럼)
    // 그중 하나만 고쳐서는 η 가 오르지 않는다 (평평한 땅). 그래서 최댓값 대신 p-노름 (Σ DCRᵖ)^(1/p) 를 쓴
    // 부드러운 η 로 고른다. p → ∞ 이면 η 와 같다. 보고하는 값과 최종 결과는 참 η 다.
    e.score = e.fit;
    if (p_.softP > 0) {
        double mx = 0;
        for (const MemberResult& q : r.members) mx = std::max(mx, q.dcr);
        double sum = 0;
        for (const MemberResult& q : r.members) sum += std::pow(q.dcr / mx, p_.softP);
        e.score = e.fit / std::pow(sum, 1.0 / p_.softP);   // = fit × max / ‖dcr‖ₚ  (≤ fit)
    }
    return e;
}

Genome GaSpace::random(unsigned long long& state) const {
    Rng rng{state};
    Genome g;
    g.act.resize(n_); g.h.resize(n_); g.web.resize(genes_.size());
    // 개체마다 밀도(절점·부재를 쓸 확률)를 따로 뽑아 성긴 트러스부터 촘촘한 트러스까지 고루 시작한다.
    // 높이는 개체마다 "가장자리 높이 → 중앙 높이" 를 뽑아 그 사이를 포물선으로 잇고 ±2단계 잡음을 더한다
    // (평행 현·활 모양·거꾸로 된 활이 모두 나온다. 슬롯마다 완전히 따로 뽑으면 상현재가 톱니처럼 꺾여 거의 다 제작 불가가 된다).
    const double pAct = 0.2 + 0.8 * rng.uni(), pWeb = 0.2 + 0.8 * rng.uni();
    const int hEdge = rng.below(levels_), hMid = rng.below(levels_);
    for (int s = 0; s < n_; ++s) {
        g.act[s] = rng.uni() < pAct;
        const double u = static_cast<double>(s + 1) / n_;            // 가장자리 쪽 0 → 중앙 1
        const double base = hEdge + (hMid - hEdge) * (1 - (1 - u) * (1 - u));
        g.h[s] = clampInt(static_cast<int>(std::lround(base)) + rng.below(5) - 2, 0, levels_ - 1);
    }
    for (auto& w : g.web) w = rng.uni() < pWeb;
    return g;
}

bool GaSpace::encodePreset(Kind kind, double H, Genome& g) const {
    if (kind == Kind::KTruss) return false;
    const double fi = (H - p_.gaHmin) / p_.gaHstep;
    const int hi = static_cast<int>(std::lround(fi));
    if (std::fabs(fi - hi) > 1e-9 || hi < 0 || hi >= levels_) return false;
    const int n = n_, half = n / 2;
    g.act.assign(n, 0); g.h.assign(n, hi); g.web.assign(genes_.size(), 0);
    for (size_t q = 0; q < genes_.size(); ++q) {
        const int s = genes_[q].slot;
        const int k = genes_[q].links[0].second;
        bool on = false;
        if (kind == Kind::Warren) {
            on = s % 2 == 1;                                      // 반 패널 점, 양쪽 사재
        } else if (s % 2 == 0) {
            const int pk = s / 2;                                 // 패널 점 번호
            if (s == n) on = genes_[q].links.size() == 1 ? true : kind == Kind::Howe;   // 중앙: 수직재, Howe 는 사재 쌍도
            else if (k == pk) on = true;                          // 수직재
            else if (pk == 1 && k == 0) on = true;                // 끝 경사재
            else if (kind == Kind::Pratt) on = k == pk + 1;       // Pratt: T(k)–B(k+1) (중앙을 향해 내려감)
            else on = k == pk - 1 && pk >= 2;                     // Howe: T(k)–B(k−1)
        }
        g.web[q] = on;
        if (on) g.act[s - 1] = 1;
    }
    (void)half;
    return true;
}

std::vector<int> GaSpace::slotGenes(int s) const {
    std::vector<int> v;
    for (size_t q = 0; q < genes_.size(); ++q) if (genes_[q].slot == s) v.push_back(static_cast<int>(q));
    return v;
}

std::string GaSpace::describe(const Genome& g) const {
    const Truss t = decode(g);
    std::ostringstream o;
    o << "절점 " << t.nodes.size() << "개, 부재 " << t.members.size() << "개, m+r-2j = " << maxwellDegree(t) << "\n상현 절점 (x, y) mm:";
    for (size_t k = n_ + 1; k < t.nodes.size(); ++k) o << " (" << t.nodes[k].x << ", " << t.nodes[k].y << ")";
    return o.str();
}

GaRun runGa(const GaSpace& space, const Params& p, unsigned seed) {
    unsigned long long state = 0x5DEECE66DULL * (seed + 1) + 0xB;
    Rng rng{state};
    const double pm = p.pm > 0 ? p.pm : 1.0 / space.totalGenes();
    const int L = space.levels();

    std::unordered_map<std::string, Eval> cache;
    GaRun run;
    int gen = 0;
    bool haveBest = false;
    // 지금까지 평가한 모든 개체 중 참 적합도가 가장 높은 것을 따로 보관한다
    auto eval = [&](const Genome& g) -> Eval {
        const std::string k = key(g);
        auto it = cache.find(k);
        if (it != cache.end()) return it->second;
        ++run.evals;
        const Eval e = cache.emplace(k, space.evaluate(g)).first->second;
        if (!haveBest || betterTrue(e, run.bestEval)) { run.best = g; run.bestEval = e; run.foundGen = gen; haveBest = true; }
        return e;
    };
    // 국소 탐색: 아래 이웃을 모두 보고
    // 가장 좋은 이웃으로 옮기기를 더 나아지지 않을 때까지 반복한다. 옮긴 횟수를 돌려준다.
    auto localSearch = [&](Genome& g, Eval& e, bool (*cmp)(const Eval&, const Eval&)) {
        int steps = 0;
        for (;;) {
            Genome bestN = g;
            Eval bestE = e;
            auto tryN = [&](const Genome& c) { const Eval x = eval(c); if (cmp(x, bestE)) { bestE = x; bestN = c; } };
            // (1) 슬롯 다시 잇기: 한 슬롯의 사용 여부와 웨브 유전자 조합을 한꺼번에 바꾼다
            //     (절점 하나를 새로 넣으려면 절점과 부재를 같이 바꿔야 하므로 한 유전자씩으로는 못 넘는 골짜기를 넘는다)
            for (int s = 1; s <= space.slots(); ++s) {
                const std::vector<int> gs = space.slotGenes(s);
                for (int a = 0; a < 2; ++a)
                    for (int mask = 0; mask < (1 << gs.size()); ++mask) {
                        Genome c = g;
                        c.act[s - 1] = static_cast<unsigned char>(a);
                        for (size_t b = 0; b < gs.size(); ++b) c.web[gs[b]] = (mask >> b) & 1;
                        tryN(c);
                    }
            }
            // (2) 높이: 한 절점 ±1·±2 단계
            for (size_t q = 0; q < g.h.size(); ++q)
                for (int d : {-2, -1, 1, 2}) {
                    const int v = g.h[q] + d;
                    if (v < 0 || v >= L) continue;
                    Genome c = g; c.h[q] = v; tryN(c);
                }
            // (3) 연속한 여러 절점 (슬롯 a..b) 을 함께 ±1 단계: 활 모양 전체를 조금 올리고 내리는 이동
            for (size_t a = 0; a < g.h.size(); ++a)
                for (size_t b = a + 1; b < g.h.size(); ++b)
                    for (int d : {-1, 1}) {
                        Genome c = g;
                        bool ok = true;
                        for (size_t q = a; q <= b && ok; ++q) { c.h[q] += d; ok = c.h[q] >= 0 && c.h[q] < L; }
                        if (ok) tryN(c);
                    }
            if (!cmp(bestE, e)) return steps;
            g = bestN; e = bestE; ++steps;
        }
    };

    std::vector<std::vector<int>> slotGenesOf(space.slots() + 1);
    for (int s = 1; s <= space.slots(); ++s) slotGenesOf[s] = space.slotGenes(s);
    std::vector<Genome> pop(p.pop);
    for (auto& g : pop) g = space.random(state);
    int lastImprove = 0;
    double lastBest = -1;

    for (gen = 0; gen < p.gens; ++gen) {
        std::vector<Eval> ev(pop.size());
        int bi = 0;
        for (size_t i = 0; i < pop.size(); ++i) {
            ev[i] = eval(pop[i]);
            if (better(ev[i], ev[bi])) bi = static_cast<int>(i);
        }
        // 몇 세대마다 그 세대 최고 개체를 국소 탐색으로 다듬어 되돌려 넣는다 (메모틱 GA)
        if (p.lsEvery > 0 && gen % p.lsEvery == p.lsEvery - 1) {
            if (localSearch(pop[bi], ev[bi], better) > 0) ++run.lsSteps;
        }
        const double cur = run.bestEval.feasible ? run.bestEval.fit : -1;
        if (cur > lastBest) { lastBest = cur; lastImprove = gen; }
        double sum = 0;
        int feas = 0;
        for (const Eval& e : ev) if (e.feasible) { sum += e.fit; ++feas; }
        run.bestFit.push_back(run.bestEval.feasible ? run.bestEval.fit : 0);
        run.meanFit.push_back(feas ? sum / feas : 0);
        run.feasRate.push_back(static_cast<double>(feas) / pop.size());
        if (gen == p.gens - 1) break;

        std::vector<int> order(pop.size());
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return better(ev[a], ev[b]); });
        std::vector<Genome> next;
        next.reserve(pop.size());
        for (int e = 0; e < p.elite; ++e) next.push_back(pop[order[e]]);

        // 정체: restart 세대 동안 최고값이 그대로면 엘리트만 남기고 나머지를 새로 뿌린다 (조기 수렴 탈출)
        if (p.restart > 0 && gen - lastImprove >= p.restart) {
            while (next.size() < pop.size()) next.push_back(space.random(state));
            pop.swap(next);
            lastImprove = gen;
            ++run.restarts;
            continue;
        }

        auto pick = [&]() {
            int w = rng.below(static_cast<int>(pop.size()));
            for (int r = 1; r < p.tour; ++r) {
                const int c = rng.below(static_cast<int>(pop.size()));
                if (better(ev[c], ev[w])) w = c;
            }
            return w;
        };
        while (next.size() < pop.size()) {
            const Genome& a = pop[pick()];
            const Genome& b = pop[pick()];
            Genome c = a;
            if (rng.uni() < p.pc) {
                if (p.slotCross) {
                    // 슬롯 단위 교차: 한 슬롯의 사용 여부·높이·웨브 부재를 한 부모에게서 통째로 받는다
                    // (절점과 그 절점에 붙은 부재는 함께 있어야 의미가 있으므로 따로 섞으면 좋은 조합이 깨진다)
                    for (int s = 1; s <= space.slots(); ++s) {
                        if (rng.uni() >= 0.5) continue;
                        c.act[s - 1] = b.act[s - 1];
                        c.h[s - 1] = b.h[s - 1];
                        for (int q : slotGenesOf[s]) c.web[q] = b.web[q];
                    }
                } else {   // 유전자별 균등 교차
                    for (size_t q = 0; q < c.act.size(); ++q) if (rng.uni() < 0.5) c.act[q] = b.act[q];
                    for (size_t q = 0; q < c.h.size(); ++q) if (rng.uni() < 0.5) c.h[q] = b.h[q];
                    for (size_t q = 0; q < c.web.size(); ++q) if (rng.uni() < 0.5) c.web[q] = b.web[q];
                }
            }
            for (auto& v : c.act) if (rng.uni() < pm) v ^= 1;   // 돌연변이
            for (auto& v : c.web) if (rng.uni() < pm) v ^= 1;
            for (auto& v : c.h)
                if (rng.uni() < pm) {
                    if (rng.uni() < 0.7) v = clampInt(v + (rng.uni() < 0.5 ? -1 : 1), 0, L - 1);
                    else v = rng.below(L);
                }
            next.push_back(std::move(c));
        }
        pop.swap(next);
    }
    // 마지막으로 최고 개체를 국소 탐색으로 다듬는다
    run.gaFit = run.bestEval.feasible ? run.bestEval.fit : 0;
    gen = p.gens;                                    // 여기서 나아지면 "마지막 다듬기에서 찾음" 으로 기록된다
    Genome g = run.best;
    Eval e = run.bestEval;
    run.lsSteps += localSearch(g, e, betterTrue);   // 참 η 로 마지막 다듬기 (eval 이 보관 최고를 갱신한다)
    return run;
}
