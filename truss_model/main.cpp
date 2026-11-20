// truss_model — 논문 「트러스 형태에 따른 구조 효율성 분석 및 GA를 통한 최적의 트러스 형태 도출」 계산 프로그램
//
//   truss_model                      truss.txt 를 풀어 result.txt 를 쓴다 (원래 엔진과 같음)
//   truss_model solve  [in] [out]    한 트러스 해석
//   truss_model preset <형태> [n] [H] 대표 4종 하나를 만들어 풀고 .txt/.svg 로 쓴다 (warren|pratt|howe|k)
//   truss_model collapse [in]        연쇄 파괴 해석 → collapse.csv
//   truss_model efit [csv]           3절 E 회귀 → efit_result.txt
//   truss_model sweep                대표 4종 × 패널 수 × 높이 → sweep.csv, sweep_best.csv
//   truss_model compare              GA 와 같은 조건에서 대표 4종의 최적 높이 → compare.csv
//   truss_model ga                   유전 알고리즘 → ga_log.csv, ga_runs.csv, ga_best.txt/.svg, ga_summary.txt
//   truss_model verify               손계산·정답·불변식 검사
//   truss_model all                  verify, efit, sweep, compare, ga 를 차례로
//
//   설정: 같은 폴더의 config.txt (없으면 기본값). 다른 파일은 -c 파일이름
#include "model.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

namespace {

Params P;   // 설정

// CSV 칸 안에 쉼표가 있으면 따옴표로 감싼다
std::string csvSafe(const std::string& v) { return v.find(',') == std::string::npos ? v : '"' + v + '"'; }

bool loadTruss(const std::string& path, Truss& t) {
    std::ifstream in(path);
    if (!in) { std::cerr << "cannot open " << path << " (is it in the working directory?)\n"; return false; }
    std::string err;
    if (!readTruss(in, t, &err)) { std::cerr << path << ": " << err << '\n'; return false; }
    return true;
}

std::string fmt(double v, int prec = 6) {
    std::ostringstream o;
    o << std::setprecision(prec) << v;
    return o.str();
}

// 대칭 구조에서는 여러 부재가 같은 DCR 로 "동시에" 용량에 닿는다. 이론상 같은 값도 마지막 자리 반올림이
// 컴퓨터(컴파일러)마다 달라서, 그중 하나만 "최초 파괴 부재"로 고르면 컴퓨터마다 답이 달라진다.
// 그래서 최댓값과 상대오차 1e-9 안에 드는 부재를 모두 동시 파괴로 묶는다.
const double TIE = 1e-9;
std::vector<int> tiedCritical(const Result& r) {
    std::vector<int> v;
    double mx = 0;
    for (const MemberResult& m : r.members) mx = std::max(mx, m.dcr);
    for (size_t k = 0; k < r.members.size(); ++k)
        if (r.members[k].dcr >= mx * (1 - TIE)) v.push_back(static_cast<int>(k));
    return v;
}
// 예: "상현 좌굴", "사재 좌굴 8개 동시", "상현 압괴 1개 + 하현 인장 2개 동시"
std::string criticalText(const Truss& t, const Result& r) {
    if (r.status != Status::Ok || r.members.empty()) return "-";
    std::vector<std::pair<std::string, int>> g;
    const std::vector<int> tied = tiedCritical(r);
    for (int k : tied) {
        const std::string key = std::string(memberType(t, k)) + " " + modeName(r.members[k].mode);
        bool found = false;
        for (auto& x : g) if (x.first == key) { ++x.second; found = true; }
        if (!found) g.push_back({key, 1});
    }
    std::string out;
    for (const auto& x : g) out += (out.empty() ? "" : " + ") + x.first + (tied.size() > 1 ? " " + std::to_string(x.second) + "개" : "");
    return tied.size() > 1 ? out + " 동시" : out;
}
// 상대오차 1e-9 안의 값은 같은 값으로 본다 (순위가 컴퓨터마다 뒤집히지 않게)
bool tieEq(double a, double b) { return std::fabs(a - b) <= TIE * std::max(std::fabs(a), std::fabs(b)); }
// 큰 순서로 늘어놓고 같은 값은 "=" 로 잇는다. 같은 값끼리는 넣은 순서를 지킨다.
std::string rankText(const std::vector<std::pair<std::string, double>>& v) {
    std::vector<size_t> ord(v.size());
    for (size_t i = 0; i < v.size(); ++i) ord[i] = i;
    std::stable_sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return !tieEq(v[a].second, v[b].second) && v[a].second > v[b].second; });
    std::string out;
    for (size_t i = 0; i < ord.size(); ++i) {
        if (i) out += tieEq(v[ord[i - 1]].second, v[ord[i]].second) ? " = " : " > ";
        out += v[ord[i]].first;
    }
    return out;
}
// 순위 (같은 값은 같은 순위)
int rankOfValue(double x, const std::vector<double>& all) {
    int k = 1;
    for (double y : all) if (!tieEq(x, y) && y > x) ++k;
    return k;
}

std::vector<double> grid(double a, double b, double step) {
    std::vector<double> v;
    const int k = static_cast<int>(std::floor((b - a) / step + 1e-9));
    for (int i = 0; i <= k; ++i) v.push_back(a + i * step);   // 누적 덧셈 대신 곱으로 (오차 누적 방지)
    return v;
}

// 한 트러스의 요약을 화면에 쓴다
void summary(const Truss& t, const Result& r, std::ostream& o) {
    o << "  절점 j = " << t.nodes.size() << ", 부재 m = " << t.members.size() << ", m + r − 2j = " << maxwellDegree(t)
      << (maxwellDegree(t) == 0 ? " (정정)" : maxwellDegree(t) > 0 ? " (부정정)" : " (불안정)") << '\n';
    if (r.status != Status::Ok) { o << "  해석 불가: " << (r.status == Status::Mechanism ? "메커니즘 (강성행렬 특이)" : "입력 오류") << '\n'; return; }
    o << "  하중계수 λ = " << fmt(r.loadFactor) << ", 파괴하중 P_fail = " << fmt(r.Pfail) << " N, 질량 = " << fmt(r.mass) << " g\n"
      << "  구조 효율 η = P_fail / 질량 = " << fmt(r.eta) << " N/g\n";
    if (r.critical >= 0) {
        const std::vector<int> tied = tiedCritical(r);
        o << "  최초 파괴: " << criticalText(t, r) << " (|N| = " << fmt(std::fabs(r.members[r.critical].N)) << " N, 용량 " << fmt(r.members[r.critical].cap) << " N; 부재";
        for (int k : tied) o << " m" << k;
        o << ")\n";
    }
    o << "  평형 잔차 = " << r.residual << " N\n";
}

void fabLine(const FabReport& f, std::ostream& o) {
    o << "  제작 제약: " << (f.violation == 0 ? "만족" : "위반") << " (최장 " << fmt(f.maxLen, 4) << " mm, 최단 " << fmt(f.minLen, 4)
      << " mm, 최대 차수 " << f.maxDeg << ", 최소 사이각 " << fmt(f.minAngle, 4) << "°, 교차 " << f.crossings << ")\n";
}

// ---------------------------------------------------------------- solve
int cmdSolve(const std::string& inPath, const std::string& outPath) {
    Truss t;
    if (!loadTruss(inPath, t)) return 1;
    std::cout << "read " << t.nodes.size() << " nodes, " << t.members.size() << " members\n";
    const Result r = solve(t);
    std::ofstream out(outPath);
    if (!out) { std::cerr << "cannot write " << outPath << " (is it open in another program?)\n"; return 1; }
    const bool ok = writeResult(out, r);
    out.close();
    if (!out) { std::cerr << "failed while writing " << outPath << '\n'; return 1; }
    writeResult(std::cout, r);
    return ok ? 0 : 1;
}

// ---------------------------------------------------------------- preset
int cmdPreset(const std::vector<std::string>& a) {
    Kind k;
    if (a.empty() || !parseKind(a[0], k)) { std::cerr << "usage: preset <warren|pratt|howe|k> [n] [H]\n"; return 1; }
    const int n = a.size() > 1 ? std::stoi(a[1]) : P.n;
    const double H = a.size() > 2 ? std::stod(a[2]) : P.H;
    if (n < 2 || n % 2) { std::cerr << "n must be even\n"; return 1; }
    const Truss t = makePreset(k, n, P.L, H, P);
    const Result r = solve(t);
    const std::string base = std::string(kindName(k)) + "_n" + std::to_string(n) + "_H" + fmt(H);
    writeTruss(base + ".txt", t);
    writeSvg(base + ".svg", t, r, std::string(kindName(k)) + "  n = " + std::to_string(n) + ", H = " + fmt(H) + " mm, η = " + fmt(r.eta, 4) + " N/g");
    std::cout << kindName(k) << " (L = " << P.L << " mm, n = " << n << ", H = " << H << " mm, P = " << P.P << " N)\n";
    summary(t, r, std::cout);
    fabLine(checkFab(t, P), std::cout);
    const CollapseResult c = progressiveCollapse(t);
    if (c.ok) std::cout << "  연쇄 파괴: 최초 " << fmt(c.Pfirst) << " N → 붕괴 P_ult = " << fmt(c.Pult) << " N (사건 " << c.events.size() << "개)\n";
    std::cout << "  → " << base << ".txt (도면판에서 불러오기 가능), " << base << ".svg\n";
    return 0;
}

// ---------------------------------------------------------------- collapse
int cmdCollapse(const std::string& inPath) {
    Truss t;
    if (!loadTruss(inPath, t)) return 1;
    const Result r0 = solve(t);
    const CollapseResult c = progressiveCollapse(t);
    if (!c.ok) { std::cerr << "처음 트러스가 불안정해서 연쇄 파괴를 볼 수 없습니다\n"; return 1; }
    std::ofstream f;
    if (!openCsv(f, "collapse.csv")) { std::cerr << "cannot write collapse.csv\n"; return 1; }
    f << "순서,하중계수 λ,하중 P (N),파괴 부재,부재 종류,파괴 모드,연쇄\n";
    std::cout << "연쇄 파괴 (" << inPath << ")\n";
    for (size_t i = 0; i < c.events.size(); ++i) {
        const CollapseEvent& e = c.events[i];
        f << i + 1 << ',' << g17(e.lam) << ',' << g17(e.P) << ",m" << e.member << ',' << memberType(t, e.member) << ',' << modeName(e.mode)
          << ',' << (e.cascade ? "연쇄" : "") << '\n';
        std::cout << "  " << i + 1 << ". λ = " << fmt(e.lam) << ", P = " << fmt(e.P) << " N: m" << e.member << " " << memberType(t, e.member)
                  << " " << modeName(e.mode) << (e.cascade ? " (연쇄)" : "") << '\n';
    }
    std::cout << "  끝: " << (c.end == "mechanism" ? "메커니즘 (붕괴)" : c.end) << "\n  P_first = " << fmt(c.Pfirst) << " N, P_ult = " << fmt(c.Pult)
              << " N, P_ult / P_first = " << fmt(c.Pult / c.Pfirst) << '\n';
    if (r0.status == Status::Ok) std::cout << "  P_ult / 질량 = " << fmt(c.Pult / r0.mass) << " N/g\n";
    std::cout << "  → collapse.csv\n";
    return 0;
}

// ---------------------------------------------------------------- efit
int cmdEfit(const std::string& path) {
    std::ostringstream o;
    if (!runEfit(path, P, o)) return 1;
    std::cout << o.str();
    std::ofstream f("efit_result.txt");
    f << o.str();
    std::cout << "  → efit_result.txt\n";
    return 0;
}

// ---------------------------------------------------------------- sweep
int cmdSweep() {
    std::ofstream f, b;
    if (!openCsv(f, "sweep.csv") || !openCsv(b, "sweep_best.csv")) { std::cerr << "cannot write sweep csv\n"; return 1; }
    f << "형태,패널 수 n,높이 H (mm),H/L,절점 j,부재 m,m+r-2j,하중계수 λ,파괴하중 P_fail (N),질량 (g),η (N/g),최초 파괴 부재,부재 종류,파괴 모드,"
         "붕괴 하중 P_ult (N),P_ult/질량 (N/g),제작 가능,제작 위반량\n";
    b << "형태,패널 수 n,최적 H (mm),H/L,최고 η (N/g),최초 파괴 모드,제작 가능한 최적 H (mm),제작 가능한 최고 η (N/g)\n";
    const std::vector<double> Hs = grid(P.Hmin, P.Hmax, P.Hstep);
    int cases = 0;
    const auto t0 = std::chrono::steady_clock::now();
    std::cout << "높이 스윕: L = " << P.L << " mm, n = " << P.nMin << "~" << P.nMax << ", H = " << P.Hmin << "~" << P.Hmax << " mm (간격 " << P.Hstep << ")\n";
    for (Kind k : ALL_KINDS) {
        for (int n = std::max(2, P.nMin + P.nMin % 2); n <= P.nMax; n += 2) {
            double bestEta = -1, bestH = 0, fBestEta = -1, fBestH = 0;
            Mode bestMode = Mode::Tension;
            for (double H : Hs) {
                const Truss t = makePreset(k, n, P.L, H, P);
                const Result r = solve(t);
                const FabReport fr = checkFab(t, P);
                ++cases;
                f << kindName(k) << ',' << n << ',' << g17(H) << ',' << g17(H / P.L) << ',' << t.nodes.size() << ',' << t.members.size() << ',' << maxwellDegree(t) << ',';
                if (r.status != Status::Ok) { f << ",,,,,,,,," << (fr.violation == 0 ? 1 : 0) << ',' << g17(fr.violation) << '\n'; continue; }
                const CollapseResult c = progressiveCollapse(t);
                f << g17(r.loadFactor) << ',' << g17(r.Pfail) << ',' << g17(r.mass) << ',' << g17(r.eta) << ",m" << r.critical << ','
                  << memberType(t, r.critical) << ',' << modeName(r.members[r.critical].mode) << ',' << g17(c.Pult) << ',' << g17(c.Pult / r.mass) << ','
                  << (fr.violation == 0 ? 1 : 0) << ',' << g17(fr.violation) << '\n';
                if (r.eta > bestEta) { bestEta = r.eta; bestH = H; bestMode = r.members[r.critical].mode; }
                if (fr.violation == 0 && r.eta > fBestEta) { fBestEta = r.eta; fBestH = H; }
            }
            b << kindName(k) << ',' << n << ',' << g17(bestH) << ',' << g17(bestH / P.L) << ',' << g17(bestEta) << ',' << modeName(bestMode) << ',';
            if (fBestEta > 0) b << g17(fBestH) << ',' << g17(fBestEta) << '\n'; else b << ",\n";
            if (n == P.n)
                std::cout << "  " << std::left << std::setw(8) << kindName(k) << " n = " << n << ": 최고 η = " << fmt(bestEta, 5) << " N/g (H = " << bestH
                          << " mm, " << modeName(bestMode) << "), 제작 가능 최고 η = " << (fBestEta > 0 ? fmt(fBestEta, 5) : "-") << " (H = " << fBestH << " mm)\n";
        }
    }
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::cout << "  " << cases << "개 해석, " << fmt(sec, 3) << " s → sweep.csv, sweep_best.csv\n";
    return 0;
}

// ---------------------------------------------------------------- compare
struct RepBest { Kind k; double H = 0, eta = -1, Pfail = 0, mass = 0, Pult = 0; int crit = -1; Mode mode = Mode::Tension; int deg = 0; bool found = false; };

// GA 와 같은 조건 (n, 높이 격자, 제작 제약) 에서 대표 4종의 최적 높이를 찾는다
std::vector<RepBest> bestReps(bool fabOnly, int n) {
    std::vector<RepBest> out;
    for (Kind k : ALL_KINDS) {
        RepBest rb; rb.k = k;
        for (double H : grid(P.gaHmin, P.gaHmax, P.gaHstep)) {
            const Truss t = makePreset(k, n, P.L, H, P);
            const Result r = solve(t);
            if (r.status != Status::Ok) continue;
            if (fabOnly && checkFab(t, P).violation > 0) continue;
            const double fit = P.fitness == "etaUlt" ? progressiveCollapse(t).Pult / r.mass : r.eta;
            if (fit > rb.eta) {
                rb.H = H; rb.eta = fit; rb.Pfail = r.Pfail; rb.mass = r.mass; rb.crit = r.critical; rb.mode = r.members[r.critical].mode;
                rb.deg = maxwellDegree(t); rb.Pult = progressiveCollapse(t).Pult; rb.found = true;
            }
        }
        out.push_back(rb);
    }
    return out;
}

int cmdCompare(std::vector<RepBest>* keep = nullptr, std::vector<RepBest>* keep2 = nullptr) {
    // GA 의 상현 후보 절점은 반 패널 간격이므로, 같은 하현 패널 수 n 과 상현 절점 간격이 같은 2n 두 경우와 비교한다
    const std::vector<RepBest> fab = bestReps(true, P.n), any = bestReps(false, P.n), fab2 = bestReps(true, 2 * P.n);
    std::ofstream f;
    if (!openCsv(f, "compare.csv")) { std::cerr << "cannot write compare.csv\n"; return 1; }
    f << "형태,조건,최적 H (mm),H/L," << (P.fitness == "etaUlt" ? "P_ult/질량 (N/g)" : "η (N/g)") << ",P_fail (N),P_ult (N),질량 (g),최초 파괴 부재,파괴 모드,m+r-2j\n";
    std::cout << "대표 4종 비교 (n = " << P.n << ", L = " << P.L << " mm, H = " << P.gaHmin << "~" << P.gaHmax << " mm 간격 " << P.gaHstep << ")\n";
    auto row = [&](const RepBest& rb, const char* cond) {
        f << kindName(rb.k) << ',' << cond << ',';
        if (!rb.found) { f << ",,,,,,,,\n"; return; }
        f << g17(rb.H) << ',' << g17(rb.H / P.L) << ',' << g17(rb.eta) << ',' << g17(rb.Pfail) << ',' << g17(rb.Pult) << ',' << g17(rb.mass)
          << ",m" << rb.crit << ',' << modeName(rb.mode) << ',' << rb.deg << '\n';
    };
    for (size_t i = 0; i < 4; ++i) {
        row(fab[i], ("n = " + std::to_string(P.n) + ", 제작 제약 적용").c_str());
        row(any[i], ("n = " + std::to_string(P.n) + ", 제약 없음").c_str());
        row(fab2[i], ("n = " + std::to_string(2 * P.n) + ", 제작 제약 적용").c_str());
        std::cout << "  " << std::left << std::setw(8) << kindName(fab[i].k) << " 제작 가능 최적: "
                  << (fab[i].found ? "H = " + fmt(fab[i].H) + " mm, " + fmt(fab[i].eta, 5) + " N/g (" + modeName(fab[i].mode) + ")" : std::string("없음"))
                  << "   | 제약 없이: H = " << any[i].H << " mm, " << fmt(any[i].eta, 5) << " N/g"
                  << "   | n = " << 2 * P.n << ": " << (fab2[i].found ? "H = " + fmt(fab2[i].H) + " mm, " + fmt(fab2[i].eta, 5) + " N/g" : std::string("없음")) << '\n';
    }
    std::cout << "  → compare.csv\n";
    if (keep) *keep = fab;
    if (keep2) *keep2 = fab2;
    return 0;
}

// ---------------------------------------------------------------- loads
// 하중 조건에 따라 형태 순위가 바뀌는가 (목적 1: 형태별 효율 비교의 조건 의존성)
int cmdLoads() {
    std::ofstream f;
    if (!openCsv(f, "loads.csv")) { std::cerr << "cannot write loads.csv\n"; return 1; }
    f << "하중 조건,형태,H = " << P.H << " mm 에서 η (N/g),순위,최초 파괴 (동시에 용량에 닿는 부재),최적 H (mm),최적 η (N/g),최적 순위\n";
    const std::string keep = P.load;
    std::cout << "하중 조건별 대표 4종 순위 (n = " << P.n << ", L = " << P.L << " mm, 제작 제약 적용)\n";
    for (const char* ld : {"center", "panel"}) {
        P.load = ld;
        const std::vector<RepBest> opt = bestReps(true, P.n);
        struct Row { Kind k; double eta; std::string crit; double optH, optEta; };
        std::vector<Row> rows;
        for (size_t i = 0; i < 4; ++i) {
            const Truss t = makePreset(ALL_KINDS[i], P.n, P.L, P.H, P);
            const Result r = solve(t);
            rows.push_back({ALL_KINDS[i], r.eta, criticalText(t, r), opt[i].H, opt[i].eta});
        }
        std::vector<double> e1, e2;
        std::vector<std::pair<std::string, double>> v1, v2;
        for (const Row& r : rows) { e1.push_back(r.eta); e2.push_back(r.optEta); v1.push_back({kindName(r.k), r.eta}); v2.push_back({kindName(r.k), r.optEta}); }
        const std::string name = std::string(ld) == "center" ? "하현 중앙 1점" : "하현 패널 점 균등";
        for (const Row& r : rows)
            f << name << ',' << kindName(r.k) << ',' << g17(r.eta) << ',' << rankOfValue(r.eta, e1) << ',' << r.crit << ','
              << g17(r.optH) << ',' << g17(r.optEta) << ',' << rankOfValue(r.optEta, e2) << '\n';
        std::cout << "  [" << name << "]\n    H = " << P.H << " mm: " << rankText(v1) << "\n    높이 최적화: " << rankText(v2) << '\n';
        for (const Row& r : rows)
            std::cout << "      " << std::left << std::setw(8) << kindName(r.k) << " 최적 H = " << fmt(r.optH) << " mm, η = " << fmt(r.optEta, 4) << " N/g, 최초 파괴: " << r.crit << '\n';
    }
    P.load = keep;
    std::cout << "  → loads.csv\n";
    return 0;
}

// ---------------------------------------------------------------- sens
// 입력값이 흔들려도 결론(형태 순위, GA 최적이 가장 좋다)이 유지되는가 (목적 2: 실험값 오차가 결론에 주는 영향)
int cmdSens() {
    struct Case { std::string name; Params q; };
    std::vector<Case> cases;
    const Params base = P;
    cases.push_back({"기준", base});
    double E0, Elo, Ehi;
    if (efitRange("deflection.csv", base, E0, Elo, Ehi)) {
        Params a = base; a.E = Elo; cases.push_back({"E = " + fmt(Elo, 5) + " MPa (기울기 +1 표준오차)", a});
        Params b = base; b.E = Ehi; cases.push_back({"E = " + fmt(Ehi, 5) + " MPa (기울기 −1 표준오차)", b});
    } else {
        Params a = base; a.E *= 0.95; cases.push_back({"E −5 %", a});
        Params b = base; b.E *= 1.05; cases.push_back({"E +5 %", b});
    }
    // 3절 실험이 실제로 잰 것은 E 가 아니라 휨강성 EI 다 (기울기 s = L³/(48EI)).
    // 그래서 단면 치수가 틀렸다면 E 도 같은 기울기에서 다시 계산해야 맞다 → EI 는 그대로, A(질량) 와 σA 만 바뀐다.
    auto section = [&](const std::string& name, double d, double t) {
        Params a = base;
        double A0, I0, A1, I1;
        tubeSection(base.tubeD, base.tubeT, A0, I0);
        tubeSection(d, t, A1, I1);
        a.tubeD = d; a.tubeT = t; a.E = base.E * I0 / I1;
        cases.push_back({name + " (E 다시 계산, EI 일정)", a});
    };
    section("빨대 두께 −10 %", base.tubeD, base.tubeT * 0.9);
    section("빨대 두께 +10 %", base.tubeD, base.tubeT * 1.1);
    section("빨대 외경 −0.1 mm", base.tubeD - 0.1, base.tubeT);
    { Params a = base; a.sigmaT *= 0.8; a.sigmaC *= 0.8; cases.push_back({"강도 σt·σc −20 %", a}); }
    { Params a = base; a.K = 0.8; cases.push_back({"유효길이계수 K = 0.8 (핀이 조금 고정)", a}); }
    { Params a = base; a.load = base.load == "center" ? "panel" : "center"; cases.push_back({std::string("하중 조건 ") + (a.load == "panel" ? "패널 점 균등" : "중앙 1점"), a}); }

    Truss gaT;
    bool haveGa = false;
    {
        std::ifstream in("ga_best.txt");
        haveGa = in && readTruss(in, gaT, nullptr);
    }
    std::ofstream f;
    if (!openCsv(f, "sens.csv")) { std::cerr << "cannot write sens.csv\n"; return 1; }
    f << "조건,형태,η (N/g),기준 대비 (%),순위,최초 파괴 모드\n";
    std::cout << "민감도 (n = " << P.n << ", H = " << P.H << " mm" << (haveGa ? ", GA 최적 = ga_best.txt" : ", ga_best.txt 없음") << ")\n";
    std::vector<double> baseEta;
    std::string baseOrder;
    for (const Case& c : cases) {
        std::vector<std::pair<std::string, double>> v;
        std::vector<Mode> modes;
        for (Kind k : ALL_KINDS) {
            const Result r = solve(makePreset(k, c.q.n, c.q.L, c.q.H, c.q));
            v.push_back({kindName(k), r.eta});
            modes.push_back(r.members[r.critical].mode);
        }
        if (haveGa) {
            Truss t = gaT;   // 같은 형상에 바뀐 재료·하중을 넣는다
            for (Member& e : t.members) e = makeMember(e.i, e.j, c.q);
            int nb = 0;
            for (const Node& q : t.nodes) if (std::fabs(q.y) < 1e-9) ++nb;
            for (Node& q : t.nodes) q.fx = q.fy = 0;
            applyLoad(t, nb - 1, c.q);
            const Result r = solve(t);
            v.push_back({"GA 최적", r.eta});
            modes.push_back(r.members[r.critical].mode);
        }
        if (baseEta.empty()) for (auto& x : v) baseEta.push_back(x.second);
        const std::string order = rankText(v);
        if (baseOrder.empty()) baseOrder = order;
        std::vector<double> all;
        for (auto& x : v) all.push_back(x.second);
        for (size_t i = 0; i < v.size(); ++i) {
            const int rank = rankOfValue(v[i].second, all);
            f << csvSafe(c.name) << ',' << v[i].first << ',' << g17(v[i].second) << ',' << g17((v[i].second / baseEta[i] - 1) * 100) << ',' << rank << ',' << modeName(modes[i]) << '\n';
        }
        std::cout << "  " << c.name << "\n      " << order << (order == baseOrder ? "" : "   ← 순위 바뀜") << '\n';
    }
    std::cout << "  → sens.csv\n";
    return 0;
}

// ---------------------------------------------------------------- materials
// 같은 형상·같은 단면에서 재료만 바꿔 대표 4종과 GA 를 다시 한다 (빨대 결과가 실제 트러스 재료에도 통하는가)
int cmdMaterials() {
    std::ofstream f;
    if (!openCsv(f, "materials.csv")) { std::cerr << "cannot write materials.csv\n"; return 1; }
    f << "재료,E (MPa),ρ (g/cm³),σt (MPa),σc (MPa),E/ρ,경계 세장비 π√(E/σc),형태,최적 H (mm),η (N/g),최초 파괴 (동시에 용량에 닿는 부재),"
         "압괴가 좌굴보다 먼저인 압축재 비율,순위,GA 대비 (%),출처\n";
    double A, I;
    tubeSection(P.tubeD, P.tubeT, A, I);
    const double r = std::sqrt(I / A);
    const int runs = std::min(P.runs, 3);
    std::cout << "재료별 비교 (n = " << P.n << ", L = " << P.L << " mm, 같은 형상·같은 빨대 단면, 하중 " << (P.load == "panel" ? "패널 점 균등" : "중앙 1점")
              << ", GA 는 재료마다 " << runs << "회)\n"
              << "  빨대 단면 회전반경 r = " << fmt(r, 4) << " mm → 길이 50~180 mm 부재의 세장비 KL/r = " << fmt(P.K * 50 / r, 3) << "~" << fmt(P.K * 180 / r, 3) << '\n';
    const Params keep = P;
    for (int mi = 0; mi < MATERIAL_COUNT; ++mi) {
        const Material& m = MATERIALS[mi];
        Params q = keep;
        q.material = m.key; q.E = m.E; q.rho = m.rho; q.sigmaT = m.sigmaT; q.sigmaC = m.sigmaC;
        P = q;   // bestReps 가 전역 P 를 쓴다
        const std::vector<RepBest> reps = bestReps(true, q.n);
        const GaSpace space(q);
        GaRun best;
        bool have = false;
        for (int k = 0; k < runs; ++k) {
            const GaRun g = runGa(space, q, q.seed + k);
            if (!have || (g.bestEval.feasible && (!best.bestEval.feasible || g.bestEval.fit > best.bestEval.fit))) { best = g; have = true; }
        }
        P = keep;
        const double slim = std::acos(-1.0) * std::sqrt(q.E / q.sigmaC);
        // 형태별 행 (GA 포함) 을 η 순으로
        struct Row { std::string name; double H, eta; std::string crit; double crushFrac; };
        std::vector<Row> rows;
        auto describe = [&](const Truss& t, const std::string& name, double H) {
            const Result res = solve(t);
            int comp = 0, crush = 0;
            for (size_t k = 0; k < res.members.size(); ++k) {
                if (res.members[k].N >= 0) continue;
                ++comp;
                const Member& e = t.members[k];
                const double L = std::hypot(t.nodes[e.j].x - t.nodes[e.i].x, t.nodes[e.j].y - t.nodes[e.i].y);
                if (e.K * L / r < slim) ++crush;   // 경계보다 짧으면 압괴가 먼저
            }
            rows.push_back({name, H, res.eta, criticalText(t, res), comp ? static_cast<double>(crush) / comp : 0});
        };
        for (const RepBest& rb : reps) if (rb.found) describe(makePreset(rb.k, q.n, q.L, rb.H, q), kindName(rb.k), rb.H);
        Truss gt;
        if (best.bestEval.feasible) {
            gt = space.decode(best.best);
            describe(gt, "GA 최적", 0);
            const std::string base = std::string("ga_best_") + m.key;
            writeTruss(base + ".txt", gt);
            writeSvg(base + ".svg", gt, solve(gt), std::string("GA 최적 (") + m.name + ")  η = " + fmt(best.bestEval.fit, 4) + " N/g");
        }
        std::stable_sort(rows.begin(), rows.end(), [](const Row& a, const Row& b) { return !tieEq(a.eta, b.eta) && a.eta > b.eta; });
        const double gaEta = best.bestEval.feasible ? best.bestEval.fit : 0;
        std::vector<std::pair<std::string, double>> rv;
        std::vector<double> all;
        for (const Row& w : rows) { rv.push_back({w.name, w.eta}); all.push_back(w.eta); }
        const std::string order = rankText(rv);
        for (size_t i = 0; i < rows.size(); ++i) {
            const Row& w = rows[i];
            f << csvSafe(m.name) << ',' << g17(q.E) << ',' << g17(q.rho) << ',' << g17(q.sigmaT) << ',' << g17(q.sigmaC) << ',' << g17(q.E / q.rho) << ',' << g17(slim) << ','
              << w.name << ',' << (w.name == "GA 최적" ? std::string("") : g17(w.H)) << ',' << g17(w.eta) << ',' << w.crit << ','
              << g17(w.crushFrac) << ',' << rankOfValue(w.eta, all) << ',' << (gaEta > 0 ? g17((gaEta / w.eta - 1) * 100) : std::string("")) << ',' << csvSafe(m.source) << '\n';
        }
        std::cout << "\n  [" << m.name << "]  E/ρ = " << fmt(q.E / q.rho, 5) << ", 경계 세장비 π√(E/σc) = " << fmt(slim, 4) << '\n'
                  << "    순위: " << order << '\n';
        for (const Row& w : rows)
            std::cout << "      " << std::left << std::setw(10) << w.name << " η = " << std::setw(9) << fmt(w.eta, 5) << " N/g"
                      << (w.name == "GA 최적" ? std::string("              ") : "  (H = " + fmt(w.H) + " mm)") << "  최초 파괴: " << w.crit
                      << ",  압괴가 먼저인 압축재 " << fmt(100 * w.crushFrac, 3) << " %\n";
    }
    std::cout << "\n  → materials.csv, ga_best_<재료>.txt/.svg\n";
    return 0;
}

// ---------------------------------------------------------------- report
// 유효숫자 sig 자리로 (끝의 0 도 남긴다: 1.230)
std::string sf(double v, int sig = 4) {
    if (!std::isfinite(v)) return "-";
    if (v == 0) return "0";
    const int mag = static_cast<int>(std::floor(std::log10(std::fabs(v))));
    const int dec = std::max(0, sig - 1 - mag);
    char b[64];
    std::snprintf(b, sizeof b, "%.*f", dec, v);
    return b;
}
std::string pct(double v) { char b[32]; std::snprintf(b, sizeof b, "%+.1f %%", v); return b; }

// 앞의 명령들이 만든 파일(ga_best*.txt, ga_runs.csv)과 빠른 계산으로 논문용 요약을 단위와 함께 쓴다
int cmdReport() {
    std::ostringstream o;
    double A, I;
    tubeSection(P.tubeD, P.tubeT, A, I);
    const char* loadName = P.load == "panel" ? "하현 안쪽 절점에 균등 분배" : "하현 중앙 1점";
    const Material* mat = findMaterial(P.material);
    o << "논문용 요약 (truss_model report)\n"
      << "==========================================================\n"
      << "조건\n"
      << "  스팬 L = " << fmt(P.L) << " mm, 패널 수 n = " << P.n << ", 하중 P = " << fmt(P.P) << " N (" << loadName << "), 지점: 왼쪽 핀 · 오른쪽 롤러\n"
      << "  재료: " << (mat ? mat->name : "직접 입력") << " — E = " << sf(P.E) << " MPa, ρ = " << sf(P.rho, 3) << " g/cm³, σt = " << sf(P.sigmaT, 3)
      << " MPa, σc = " << sf(P.sigmaC, 3) << " MPa, K = " << P.K << '\n'
      << "  단면: 원형 관 외경 " << P.tubeD << " mm × 두께 " << P.tubeT << " mm → A = " << sf(A) << " mm², I = " << sf(I) << " mm⁴\n"
      << "  표기: 유효숫자 4자리. 상대오차 1e-9 안의 같은 값은 \"=\" (순위) 와 \"동시\" (파괴 부재) 로 적는다.\n\n";

    // [1] E
    EfitInfo ea, ek;
    if (efitInfo("deflection.csv", P, ea, ek)) {
        o << "[1] 탄성계수 E 회귀 (3절, deflection.csv, 단순보 스팬 " << fmt(P.beamL) << " mm)\n"
          << "  전체 " << ea.n << "점: 기울기 s = " << sf(ea.s) << " mm/N, R² = " << sf(ea.r2) << ", E = " << sf(ea.E) << " MPa\n"
          << "  판독 잔차 ≤ " << P.rmsMax << " px, " << ek.n << "점: s = " << sf(ek.s) << " ± " << sf(ek.se) << " mm/N (표준오차), R² = " << sf(ek.r2)
          << ", E = " << sf(ek.E) << " MPa (기울기 ±1 표준오차: " << sf(ek.lo) << " ~ " << sf(ek.hi) << " MPa)\n"
          << "  시뮬레이션에 쓴 E = " << sf(P.E) << " MPa (회귀값과 차이 " << pct((P.E / ek.E - 1) * 100) << ")\n\n";
    }

    // [2] H 고정
    o << "[2] 대표 4종, H = " << fmt(P.H) << " mm (n = " << P.n << ")\n"
      << "  형태 | 절점 j | 부재 m | m+r−2j | 하중계수 λ | P_fail (N) | 질량 (g) | η (N/g) | 최초 파괴\n";
    std::vector<std::pair<std::string, double>> rk;
    for (Kind k : ALL_KINDS) {
        const Truss t = makePreset(k, P.n, P.L, P.H, P);
        const Result r = solve(t);
        o << "  " << kindName(k) << " | " << t.nodes.size() << " | " << t.members.size() << " | " << maxwellDegree(t) << " | " << sf(r.loadFactor) << " | "
          << sf(r.Pfail) << " | " << sf(r.mass) << " | " << sf(r.eta) << " | " << criticalText(t, r) << '\n';
        rk.push_back({kindName(k), r.eta});
    }
    o << "  순위 (η): " << rankText(rk) << "\n\n";

    // [3] 높이 최적화
    const std::vector<RepBest> b1 = bestReps(true, P.n), b2 = bestReps(true, 2 * P.n);
    o << "[3] 높이 최적화 (H = " << fmt(P.gaHmin) << " ~ " << fmt(P.gaHmax) << " mm, " << fmt(P.gaHstep) << " mm 간격, 제작 제약 적용)\n"
      << "  형태 | n = " << P.n << ": 최적 H (mm), η (N/g), 최초 파괴 | n = " << 2 * P.n << ": 최적 H (mm), η (N/g)\n";
    std::vector<std::pair<std::string, double>> r1, r2;
    for (size_t i = 0; i < 4; ++i) {
        const Truss t = makePreset(b1[i].k, P.n, P.L, b1[i].H, P);
        o << "  " << kindName(b1[i].k) << " | " << fmt(b1[i].H) << ", " << sf(b1[i].eta) << ", " << criticalText(t, solve(t)) << " | "
          << fmt(b2[i].H) << ", " << sf(b2[i].eta) << '\n';
        r1.push_back({kindName(b1[i].k), b1[i].eta});
        r2.push_back({kindName(b2[i].k), b2[i].eta});
    }
    o << "  순위 n = " << P.n << ": " << rankText(r1) << "\n  순위 n = " << 2 * P.n << ": " << rankText(r2) << "\n\n";

    // [4] 하중 조건
    o << "[4] 하중 조건별 순위 (n = " << P.n << ", 높이 최적화)\n";
    {
        const std::string keep = P.load;
        for (const char* ld : {"center", "panel"}) {
            P.load = ld;
            const std::vector<RepBest> b = bestReps(true, P.n);
            std::vector<std::pair<std::string, double>> v;
            o << "  " << (std::string(ld) == "center" ? "하현 중앙 1점" : "하현 균등 분배") << ":";
            for (const RepBest& x : b) { v.push_back({kindName(x.k), x.eta}); o << "  " << kindName(x.k) << " " << sf(x.eta) << " N/g (H " << fmt(x.H) << " mm)"; }
            o << "\n    → " << rankText(v) << '\n';
        }
        P.load = keep;
    }
    o << '\n';

    // [5] 연쇄 파괴
    o << "[5] 연쇄 파괴 (H = " << fmt(P.H) << " mm)\n  형태 | P_first (N) | P_ult (N) | P_ult / P_first | 파괴 순서\n";
    for (Kind k : ALL_KINDS) {
        const Truss t = makePreset(k, P.n, P.L, P.H, P);
        const CollapseResult c = progressiveCollapse(t);
        o << "  " << kindName(k) << " | " << sf(c.Pfirst) << " | " << sf(c.Pult) << " | " << sf(c.Pult / c.Pfirst) << " |";
        const size_t tied = tiedCritical(solve(t)).size();
        for (size_t i = 0; i < c.events.size(); ++i) {
            const CollapseEvent& e = c.events[i];
            o << (i ? " → " : " ") << "m" << e.member << "(" << memberType(t, e.member) << " " << modeName(e.mode) << (e.cascade ? ", 연쇄" : "");
            if (i == 0 && tied > 1) o << ", 동시에 닿는 " << tied << "개 중 번호가 가장 작은 부재를 먼저 제거";
            o << ")";
        }
        o << " → " << (c.end == "mechanism" ? "메커니즘" : c.end) << '\n';
    }
    o << '\n';

    // [6] GA
    Truss gt;
    {
        std::ifstream in("ga_best.txt");
        if (in && readTruss(in, gt, nullptr)) {
            const Result r = solve(gt);
            const FabReport fr = checkFab(gt, P);
            o << "[6] 유전 알고리즘 최적 (ga_best.txt)\n"
              << "  절점 j = " << gt.nodes.size() << ", 부재 m = " << gt.members.size() << ", m+r−2j = " << maxwellDegree(gt) << "\n  상현 절점 (x, y) mm:";
            for (const Node& q : gt.nodes) if (q.y > 1e-9) o << " (" << fmt(q.x) << ", " << fmt(q.y) << ")";
            o << "\n  하중계수 λ = " << sf(r.loadFactor) << ", P_fail = " << sf(r.Pfail) << " N, 질량 = " << sf(r.mass) << " g, η = " << sf(r.eta) << " N/g\n"
              << "  최초 파괴: " << criticalText(gt, r) << "\n"
              << "  제작: 최장 " << sf(fr.maxLen) << " mm, 최단 " << sf(fr.minLen) << " mm, 최대 차수 " << fr.maxDeg << ", 최소 사이각 " << sf(fr.minAngle) << "°, 교차 "
              << fr.crossings << (fr.violation == 0 ? " → 제작 가능" : " → 제약 위반") << '\n';
            // 반복 결과
            std::ifstream rs("ga_runs.csv");
            std::string line;
            int runs = 0, same = 0;
            std::getline(rs, line);
            while (std::getline(rs, line)) {
                std::istringstream ls(line);
                std::string a, b, c;
                std::getline(ls, a, ','); std::getline(ls, b, ','); std::getline(ls, c, ',');
                if (c.empty()) continue;
                ++runs;
                if (tieEq(std::stod(c), r.eta)) ++same;
            }
            if (runs) o << "  반복: " << runs << "회 중 " << same << "회가 같은 최적에 도달\n";
            double t1 = 0, t2 = 0; std::string k1, k2;
            for (const RepBest& x : b1) if (x.eta > t1) { t1 = x.eta; k1 = kindName(x.k); }
            for (const RepBest& x : b2) if (x.eta > t2) { t2 = x.eta; k2 = kindName(x.k); }
            o << "  대표 최고 (n = " << P.n << "): " << k1 << " " << sf(t1) << " N/g → GA " << pct((r.eta / t1 - 1) * 100) << '\n'
              << "  대표 최고 (n = " << 2 * P.n << ", 상현 간격 같음): " << k2 << " " << sf(t2) << " N/g → GA " << pct((r.eta / t2 - 1) * 100) << "\n\n";
        } else {
            o << "[6] 유전 알고리즘: ga_best.txt 가 없습니다 (ga 를 먼저 실행)\n\n";
        }
    }

    // [7] 민감도
    {
        std::ostringstream so;
        auto* old = std::cout.rdbuf(so.rdbuf());
        cmdSens();
        std::cout.rdbuf(old);
        std::string sv = so.str();
        const auto cut = sv.find("  → sens.csv");
        if (cut != std::string::npos) sv.erase(cut);
        o << "[7] " << sv.substr(sv.find("민감도")) << '\n';
    }

    // [8] 재료별
    o << "[8] 재료별 (같은 형상·같은 빨대 단면, 높이 최적화, GA 는 ga_best_<재료>.txt)\n";
    {
        const Params keep = P;
        for (int mi = 0; mi < MATERIAL_COUNT; ++mi) {
            const Material& m = MATERIALS[mi];
            Params q = keep;
            q.material = m.key; q.E = m.E; q.rho = m.rho; q.sigmaT = m.sigmaT; q.sigmaC = m.sigmaC;
            P = q;
            const std::vector<RepBest> b = bestReps(true, q.n);
            P = keep;
            std::vector<std::pair<std::string, double>> v;
            double top = 0;
            for (const RepBest& x : b) { v.push_back({kindName(x.k), x.eta}); top = std::max(top, x.eta); }
            o << "  " << m.name << " (E = " << sf(m.E) << " MPa, ρ = " << sf(m.rho, 3) << " g/cm³, 경계 세장비 π√(E/σc) = " << sf(std::acos(-1.0) * std::sqrt(m.E / m.sigmaC), 3) << ")\n";
            std::ifstream in(std::string("ga_best_") + m.key + ".txt");
            Truss t;
            if (in && readTruss(in, t, nullptr)) {
                for (Member& e : t.members) e = makeMember(e.i, e.j, q);
                const Result r = solve(t);
                v.insert(v.begin(), {"GA", r.eta});
                o << "    GA η = " << sf(r.eta) << " N/g (" << criticalText(t, r) << "), 대표 최고 " << sf(top) << " N/g → " << pct((r.eta / top - 1) * 100) << '\n';
            }
            for (const RepBest& x : b) {
                const Truss tt = makePreset(x.k, q.n, q.L, x.H, q);
                o << "    " << kindName(x.k) << " η = " << sf(x.eta) << " N/g (H = " << fmt(x.H) << " mm, " << criticalText(tt, solve(tt)) << ")\n";
            }
            o << "    순위: " << rankText(v) << '\n';
        }
    }

    std::cout << o.str();
    std::ofstream f("summary.txt", std::ios::binary);
    f << "\xEF\xBB\xBF" << o.str();
    std::cout << "\n  → summary.txt (메모장으로 열 수 있음)\n";
    return 0;
}

// ---------------------------------------------------------------- ga
int cmdGa() {
    const GaSpace space(P);
    std::cout << "유전 알고리즘\n"
              << "  설계 공간: 상현 슬롯 " << space.slots() << "개 (대칭) × 높이 " << space.levels() << "단계 (" << P.gaHmin << "~" << P.gaHmax << " mm), 웨브 유전자 "
              << space.webGenes() << "개 → 유전자 " << space.totalGenes() << "개, 경우의 수 약 10^" << fmt(space.log10Size(), 3) << '\n'
              << "  개체 " << P.pop << ", 세대 " << P.gens << ", 토너먼트 " << P.tour << ", 엘리트 " << P.elite << ", 교차 " << P.pc
              << ", 돌연변이 " << (P.pm > 0 ? P.pm : 1.0 / space.totalGenes()) << ", 국소 탐색 " << P.lsEvery << "세대마다, 정체 " << P.restart << "세대면 재시작, 반복 " << P.runs << "회 (seed " << P.seed << "~" << P.seed + P.runs - 1 << ")\n"
              << "  적합도: " << (P.fitness == "etaUlt" ? "P_ult / 질량" : "η = P_fail / 질량") << ", 제작 제약: 길이 " << P.Lmin << "~" << P.Lmax << " mm, 차수 ≤ " << P.maxDeg
              << ", 사이각 ≥ " << P.minAngle << "°, 교차 " << (P.noCross ? "금지" : "허용") << '\n';

    std::ofstream log, runs;
    if (!openCsv(log, "ga_log.csv") || !openCsv(runs, "ga_runs.csv")) { std::cerr << "cannot write ga csv\n"; return 1; }
    log << "반복,seed,세대,지금까지 최고 적합도,제작 가능 개체 평균 적합도,제작 가능 개체 비율\n";
    runs << "반복,seed,최고 적합도,η (N/g),P_fail (N),질량 (g),절점 j,부재 m,m+r-2j,찾은 세대,서로 다른 평가 수,상현 절점 (x:y)\n";

    GaRun best;
    bool have = false;
    std::vector<double> fits;
    const auto t0 = std::chrono::steady_clock::now();
    for (int r = 0; r < P.runs; ++r) {
        const unsigned seed = P.seed + r;
        const GaRun run = runGa(space, P, seed);
        for (size_t g = 0; g < run.bestFit.size(); ++g)
            log << r + 1 << ',' << seed << ',' << g << ',' << g17(run.bestFit[g]) << ',' << g17(run.meanFit[g]) << ',' << g17(run.feasRate[g]) << '\n';
        const Truss t = space.decode(run.best);
        runs << r + 1 << ',' << seed << ',' << g17(run.bestEval.fit) << ',' << g17(run.bestEval.eta) << ',' << g17(run.bestEval.Pfail) << ',' << g17(run.bestEval.mass)
             << ',' << run.bestEval.j << ',' << run.bestEval.m << ',' << run.bestEval.deg << ',' << run.foundGen << ',' << run.evals << ',';
        for (size_t k = P.n + 1; k < t.nodes.size(); ++k) runs << (k > static_cast<size_t>(P.n) + 1 ? " " : "") << g17(t.nodes[k].x) << ':' << g17(t.nodes[k].y);
        runs << '\n';
        std::cout << "  반복 " << r + 1 << " (seed " << seed << "): " << (run.bestEval.feasible ? fmt(run.bestEval.fit, 6) + " N/g" : std::string("제작 가능한 해 없음"))
                  << " (GA " << fmt(run.gaFit, 6) << " → 마지막 국소 탐색 후), " << (run.foundGen >= P.gens ? std::string("마지막 다듬기에서") : std::to_string(run.foundGen) + "세대에")
                  << " 찾음, 재시작 " << run.restarts << "회, 서로 다른 평가 " << run.evals << "개\n";
        fits.push_back(run.bestEval.feasible ? run.bestEval.fit : 0);
        if (!have || (run.bestEval.feasible && (!best.bestEval.feasible || run.bestEval.fit > best.bestEval.fit))) { best = run; have = true; }
    }
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (!best.bestEval.feasible) { std::cerr << "제작 가능한 해를 찾지 못했습니다. 제약이나 높이 범위를 확인하세요.\n"; return 1; }

    const Truss t = space.decode(best.best);
    const Result r = solve(t);
    const FabReport fr = checkFab(t, P);
    const CollapseResult c = progressiveCollapse(t);
    writeTruss("ga_best.txt", t);
    writeSvg("ga_best.svg", t, r, "GA 최적  n = " + std::to_string(P.n) + ", η = " + fmt(r.eta, 4) + " N/g, 질량 " + fmt(r.mass, 4) + " g");

    std::vector<RepBest> reps, reps2;
    std::ostringstream cmpOut;
    {
        auto* old = std::cout.rdbuf(cmpOut.rdbuf());
        cmdCompare(&reps, &reps2);
        std::cout.rdbuf(old);
    }
    const RepBest* top = nullptr;
    for (const RepBest& rb : reps) if (rb.found && (!top || rb.eta > top->eta)) top = &rb;
    const RepBest* top2 = nullptr;
    for (const RepBest& rb : reps2) if (rb.found && (!top2 || rb.eta > top2->eta)) top2 = &rb;

    int same = 0, near = 0;
    for (double v : fits) {
        if (std::fabs(v - best.bestEval.fit) <= 1e-9 * best.bestEval.fit) ++same;
        if (v >= 0.97 * best.bestEval.fit) ++near;
    }
    double mean = 0, sd = 0;
    for (double v : fits) mean += v;
    mean /= fits.size();
    for (double v : fits) sd += (v - mean) * (v - mean);
    sd = fits.size() > 1 ? std::sqrt(sd / (fits.size() - 1)) : 0;

    std::ostringstream o;
    o << "GA 결과 요약\n"
      << "  조건: L = " << P.L << " mm, n = " << P.n << ", P = " << P.P << " N, E = " << P.E << " MPa, 빨대 외경 " << P.tubeD << " mm × 두께 " << P.tubeT << " mm\n"
      << "  설계 공간: 유전자 " << space.totalGenes() << "개, 경우의 수 약 10^" << fmt(space.log10Size(), 3) << '\n'
      << "  반복 " << P.runs << "회: 최고 " << fmt(best.bestEval.fit, 6) << " N/g, 평균 " << fmt(mean, 6) << ", 표준편차 " << fmt(sd, 3)
      << "\n  최고값에 도달한 반복 " << same << "/" << P.runs << "회, 최고값의 97% 이상 " << near << "/" << P.runs << "회, 걸린 시간 " << fmt(sec, 3) << " s\n\n"
      << "최적 트러스 (ga_best.txt, ga_best.svg)\n"
      << "  " << space.describe(best.best) << '\n';
    summary(t, r, o);
    fabLine(fr, o);
    if (c.ok) o << "  연쇄 파괴: 최초 " << fmt(c.Pfirst) << " N → P_ult = " << fmt(c.Pult) << " N\n";
    o << "\n대표 4종 (같은 n, 같은 높이 격자, 같은 제작 제약에서 높이를 최적화)\n" << cmpOut.str();
    if (top) {
        const double gain = (best.bestEval.fit / top->eta - 1) * 100;
        o << "\n  가장 좋은 대표 형태: " << kindName(top->k) << " (H = " << top->H << " mm, " << fmt(top->eta, 6) << " N/g)\n"
          << "  GA 최적 / 대표 최고 = " << fmt(best.bestEval.fit / top->eta, 6) << "  →  효율 " << fmt(gain, 4) << " % " << (gain >= 0 ? "향상" : "감소") << '\n';
    }
    if (top2) {
        const double gain = (best.bestEval.fit / top2->eta - 1) * 100;
        o << "  (상현 절점 간격이 같은 n = " << 2 * P.n << " 대표 최고: " << kindName(top2->k) << ", H = " << top2->H << " mm, " << fmt(top2->eta, 6)
          << " N/g → GA 가 " << fmt(gain, 4) << " % " << (gain >= 0 ? "높음" : "낮음") << ")\n";
    }
    // 대표 형태가 GA 공간 안에 있음을 확인 (같은 트러스를 유전자로 만들어 같은 값이 나오는지)
    o << "\n검증: 대표 형태를 GA 유전자로 표현해 같은 값이 나오는가\n";
    for (Kind k : {Kind::Warren, Kind::Pratt, Kind::Howe}) {
        Genome g;
        const double H = P.H;
        if (!space.encodePreset(k, H, g)) { o << "  " << kindName(k) << ": H = " << H << " mm 이 GA 높이 격자에 없음\n"; continue; }
        const Result a = solve(makePreset(k, P.n, P.L, H, P));
        const Eval e = space.evaluate(g);
        o << "  " << std::left << std::setw(7) << kindName(k) << " H = " << H << ": 직접 생성 η = " << g17(a.eta) << ", 유전자 해독 η = " << g17(e.eta)
          << (std::fabs(a.eta - e.eta) <= 1e-12 * a.eta ? "  일치" : "  불일치") << '\n';
    }
    o << "  K-truss 는 중간 높이 절점이 있어 GA 공간 밖이며, 위 비교에서 따로 다룬다.\n";
    std::cout << '\n' << o.str();
    std::ofstream s("ga_summary.txt");
    s << o.str();
    std::cout << "  → ga_log.csv, ga_runs.csv, ga_best.txt, ga_best.svg, ga_summary.txt, compare.csv\n";
    return 0;
}

// ---------------------------------------------------------------- verify
int fails = 0, passes = 0;
void check(bool ok, const std::string& what) {
    (ok ? passes : fails)++;
    std::cout << (ok ? "  [통과] " : "  [실패] ") << what << '\n';
}
bool rel(double a, double b, double tol) { return std::fabs(a - b) <= tol * std::max(1.0, std::fabs(b)); }

int cmdVerify() {
    std::cout << "검증\n";
    // 1. 기본 삼각형: README 의 정답과 비교
    {
        Truss t;
        std::istringstream in("TRUSS 1\n3 3\n0 0 1 1 0 0\n200 0 0 1 0 0\n100 141.42135623730951 0 0 0 -10\n"
                              "0 1 1500 2.7570351643589585 11.795788893216003 0.905 30 40 1\n"
                              "0 2 1500 2.7570351643589585 11.795788893216003 0.905 30 40 1\n"
                              "1 2 1500 2.7570351643589585 11.795788893216003 0.905 30 40 1\n");
        std::ifstream file("truss.txt");
        const bool useFile = static_cast<bool>(file) && readTruss(file, t, nullptr);
        if (!useFile) readTruss(in, t, nullptr);
        const Result r = solve(t);
        check(r.status == Status::Ok && rel(r.members[0].N, 3.5355339059327369, 1e-12) && rel(r.members[1].N, -6.1237243569579443, 1e-12)
                  && rel(r.loadFactor, 0.95095134548908589, 1e-12) && rel(r.eta, 6.9757989126382771, 1e-12),
              "기본 삼각형: 부재력·λ·η 가 README 정답과 상대오차 1e-12 이내");
        const double th = std::atan(std::sqrt(2.0));
        check(rel(r.members[0].N, 5 / std::tan(th), 1e-12) && rel(r.members[1].N, -5 / std::sin(th), 1e-12), "기본 삼각형: 손계산 (P/2)cotθ, −P/(2 sinθ) 와 일치");
    }
    // 2. 부재 수 일반식 (보조정리 1.1.1) 과 정정 여부
    {
        bool ok = true;
        std::string bad;
        for (int n = 2; n <= 20; n += 2)
            for (Kind k : ALL_KINDS) {
                const Truss t = makePreset(k, n, 600, 100, P);
                int j = 0, m = 0, d = 0;
                switch (k) {
                    case Kind::Warren: j = 2 * n + 1; m = 4 * n - 1; d = 0; break;
                    case Kind::Pratt: case Kind::Howe: j = 2 * n; m = 4 * n - 3; d = 0; break;
                    case Kind::KTruss: j = 3 * n - 1; m = 6 * n - 4; d = 1; break;
                }
                if (n == 2 && k == Kind::KTruss) continue;   // n = 2 K-truss 는 정의상 의미가 없다
                if (static_cast<int>(t.nodes.size()) != j || static_cast<int>(t.members.size()) != m || maxwellDegree(t) != d) {
                    ok = false; bad += std::string(" ") + kindName(k) + " n=" + std::to_string(n);
                }
            }
        check(ok, "부재 수 일반식: Warren j=2n+1, m=4n−1 / Pratt·Howe j=2n, m=4n−3 / K j=3n−1, m=6n−4, 1차 부정정 (n = 4~20)" + bad);
    }
    // 3. 평형: 반력 합 = 하중, 잔차
    {
        bool ok = true;
        double worst = 0;
        for (Kind k : ALL_KINDS)
            for (int n : {4, 6, 12})
                for (double H : {50.0, 100.0, 200.0}) {
                    const Truss t = makePreset(k, n, 600, H, P);
                    const Result r = solve(t);
                    double ry = 0;
                    for (double v : r.ry) ry += v;
                    worst = std::max(worst, r.residual);
                    if (r.status != Status::Ok || !rel(ry, P.P, 1e-12) || r.residual > 1e-9) ok = false;
                }
        check(ok, "대표 4종 36경우: 반력 합 = P, 평형 잔차 최대 " + fmt(worst, 3) + " N");
    }
    // 4. 중앙 상현재 힘 = PL/(4H)  (Pratt, 단면법)
    {
        const double H = 100;
        const Truss t = makePreset(Kind::Pratt, 6, 600, H, P);
        const Result r = solve(t);
        double mx = 0;
        for (size_t k = 0; k < t.members.size(); ++k) if (std::string(memberType(t, static_cast<int>(k))) == "상현") mx = std::max(mx, -r.members[k].N);
        check(rel(mx, P.P * 600 / (4 * H), 1e-12), "Pratt n=6, H=100: 상현재 최대 압축력 = PL/(4H) = " + fmt(P.P * 600 / (4 * H)) + " N (단면법)");
        const double sh = H / std::hypot(100.0, H);
        double endPost = 0;
        for (size_t k = 0; k < t.members.size(); ++k) if (std::string(memberType(t, static_cast<int>(k))) == "사재") endPost = std::max(endPost, -r.members[k].N);
        check(rel(endPost, P.P / 2 / sh, 1e-12), "Pratt n=6, H=100: 끝 경사재 압축력 = (P/2)/sinθ = " + fmt(P.P / 2 / sh) + " N");
    }
    // 5. 연쇄 파괴: 첫 사건 = 선형 해석의 λ
    {
        bool ok = true;
        for (Kind k : ALL_KINDS) {
            const Truss t = makePreset(k, 6, 600, 100, P);
            const Result r = solve(t);
            const CollapseResult c = progressiveCollapse(t);
            // 대칭 구조는 같은 DCR 부재가 둘 (마지막 자리만 다름) → 번호 대신 DCR 이 최댓값과 같은지 본다
            if (!c.ok || !rel(c.lamFirst, r.loadFactor, 1e-12) || !rel(r.members[c.events.front().member].dcr, r.members[r.critical].dcr, 1e-12) || c.lamUlt < c.lamFirst) ok = false;
        }
        check(ok, "연쇄 파괴: 첫 사건의 λ와 부재 DCR 이 선형 해석의 최댓값과 같고, P_ult ≥ P_first (4종)");
        Params q = P; q.P = 6.6;
        const CollapseResult c = progressiveCollapse(makePreset(Kind::KTruss, 6, 600, 100, q));
        check(c.ok && c.events.size() >= 2 && c.events[1].cascade && c.end == "mechanism",
              "K-truss n=6, H=100: 첫 좌굴 뒤 연쇄 파괴로 메커니즘 (도면판과 같은 결과, P_first = " + fmt(c.Pfirst) + " N)");
    }
    // 6. GA 공간이 대표 3종을 포함하고, 해독 결과가 직접 생성과 같다
    {
        Params q = P; q.gaHmin = 30; q.gaHmax = 200; q.gaHstep = 5;
        const GaSpace sp(q);
        bool ok = true;
        for (Kind k : {Kind::Warren, Kind::Pratt, Kind::Howe})
            for (double H : {50.0, 100.0, 150.0}) {
                Genome g;
                if (!sp.encodePreset(k, H, g)) { ok = false; continue; }
                const Truss a = makePreset(k, q.n, q.L, H, q), b = sp.decode(g);
                const Result ra = solve(a), rb = solve(b);
                if (a.members.size() != b.members.size() || !rel(ra.eta, rb.eta, 1e-12) || !rel(ra.Pfail, rb.Pfail, 1e-12)) ok = false;
            }
        check(ok, "GA 해독: Warren·Pratt·Howe 를 유전자로 표현하면 직접 생성한 트러스와 부재 수·η·P_fail 이 같다");
        // 대칭: 해독한 트러스는 좌우 대칭
        unsigned long long st = 7;
        bool sym = true;
        for (int i = 0; i < 200; ++i) {
            const Truss t = sp.decode(sp.random(st));
            for (const Node& a : t.nodes) {
                bool found = false;
                for (const Node& b : t.nodes) if (rel(b.x, q.L - a.x, 1e-12) && rel(b.y, a.y, 1e-12)) found = true;
                if (!found) sym = false;
            }
        }
        check(sym, "GA 해독: 무작위 유전자 200개의 트러스가 모두 좌우 대칭");
    }
    // 7. 제작 제약 검사기
    {
        Params q = P;
        const FabReport ok = checkFab(makePreset(Kind::Pratt, 6, 600, 100, q), q);
        Truss x;   // X 자로 교차하는 두 부재
        const double pts[4][2] = {{0, 0}, {100, 0}, {0, 100}, {100, 100}};
        for (const auto& pt : pts) { Node nd; nd.x = pt[0]; nd.y = pt[1]; x.nodes.push_back(nd); }
        x.members = {makeMember(0, 3, q), makeMember(1, 2, q)};
        const FabReport cr = checkFab(x, q);
        q.Lmax = 120;
        const FabReport lg = checkFab(makePreset(Kind::Pratt, 6, 600, 100, q), q);   // 사재 141 mm > 120
        check(ok.violation == 0 && cr.crossings == 1 && lg.tooLong > 0, "제작 제약: Pratt(H=100) 통과, X 교차 1건 검출, L_max = 120 mm 이면 긴 사재 검출");
    }
    // 8. 하중 조건: 패널 점 균등 하중도 반력 합 = P, 좌우 반력이 같다
    {
        Params q = P; q.load = "panel";
        bool ok = true;
        for (Kind k : ALL_KINDS) {
            const Result r = solve(makePreset(k, 6, 600, 100, q));
            const double ry0 = r.ry[0], ryn = r.ry[6];
            if (r.status != Status::Ok || !rel(ry0 + ryn, q.P, 1e-12) || !rel(ry0, ryn, 1e-9)) ok = false;
        }
        check(ok, "패널 점 균등 하중: 반력 합 = P, 좌우 반력 = P/2 (4종)");
    }
    std::cout << "  " << passes << "개 통과, " << fails << "개 실패\n";
    return fails ? 1 : 0;
}

} // namespace

int main(int argc, char* argv[]) {
    std::vector<std::string> args;
    std::string cfg = "config.txt";
    bool cfgGiven = false;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "-c" && i + 1 < argc) { cfg = argv[++i]; cfgGiven = true; }
        else args.push_back(a);
    }
    {
        std::ifstream probe(cfg);
        if (probe || cfgGiven) {
            std::string err;
            if (!readConfig(cfg, P, &err)) { std::cerr << err << '\n'; return 1; }
        }
    }
    const std::string cmd = args.empty() ? "solve" : args[0];
    const std::vector<std::string> rest(args.size() > 1 ? args.begin() + 1 : args.end(), args.end());
    auto arg = [&](size_t i, const char* def) { return rest.size() > i ? rest[i] : std::string(def); };

    if (cmd == "solve") return cmdSolve(arg(0, "truss.txt"), arg(1, "result.txt"));
    if (cmd == "preset") return cmdPreset(rest);
    if (cmd == "collapse") return cmdCollapse(arg(0, "truss.txt"));
    if (cmd == "efit") return cmdEfit(arg(0, "deflection.csv"));
    if (cmd == "sweep") return cmdSweep();
    if (cmd == "compare") return cmdCompare();
    if (cmd == "loads") return cmdLoads();
    if (cmd == "sens") return cmdSens();
    if (cmd == "materials") return cmdMaterials();
    if (cmd == "report") return cmdReport();
    if (cmd == "ga") return cmdGa();
    if (cmd == "verify") return cmdVerify();
    if (cmd == "all") {
        int rc = cmdVerify();
        std::cout << '\n'; rc |= cmdEfit("deflection.csv");
        std::cout << '\n'; rc |= cmdSweep();
        std::cout << '\n'; rc |= cmdLoads();
        std::cout << '\n'; rc |= cmdGa();
        std::cout << '\n'; rc |= cmdSens();
        std::cout << '\n'; rc |= cmdMaterials();
        std::cout << '\n'; rc |= cmdReport();
        return rc;
    }
    std::cerr << "unknown command '" << cmd << "'. commands: solve preset collapse efit sweep compare loads ga sens materials report verify all\n";
    return 1;
}
