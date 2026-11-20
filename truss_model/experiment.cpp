// 3절 실험: 하중–처짐 직선의 기울기로 탄성계수 E 를 구한다.
//
//   단순보 중앙 집중하중:  δ = P L³ / (48 E I)
//   δ = a + s·P 로 최소제곱 직선을 맞추면  E = L³ / (48 I s)
//   절편 a 에는 예비하중(빈 걸이 12 g)과 처음 휨이 들어가고, 기울기 s 에는 들어가지 않는다.
//   하중 P = m·g / 1000  [N]  (m: 물의 질량 [g], 1 mL = 1 g)
#include "model.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct Point { double mass, delta, rms; };

namespace {

struct Fit { int n = 0; double s = 0, a = 0, r2 = 0, E = 0, se = 0; };

// 정규방정식으로 δ = a + s P 를 맞춘다
Fit lineFit(const std::vector<double>& P, const std::vector<double>& d, double L, double I) {
    Fit f;
    f.n = static_cast<int>(P.size());
    if (f.n < 2) return f;
    double mp = 0, md = 0;
    for (int k = 0; k < f.n; ++k) { mp += P[k]; md += d[k]; }
    mp /= f.n; md /= f.n;
    double sxx = 0, sxy = 0, syy = 0;
    for (int k = 0; k < f.n; ++k) {
        sxx += (P[k] - mp) * (P[k] - mp);
        sxy += (P[k] - mp) * (d[k] - md);
        syy += (d[k] - md) * (d[k] - md);
    }
    f.s = sxy / sxx;
    f.a = md - f.s * mp;
    double sse = 0;
    for (int k = 0; k < f.n; ++k) { const double e = d[k] - (f.a + f.s * P[k]); sse += e * e; }
    f.r2 = 1 - sse / syy;
    f.se = f.n > 2 ? std::sqrt(sse / (f.n - 2) / sxx) : 0;   // 기울기의 표준오차
    f.E = L * L * L / (48 * I * f.s);
    return f;
}

} // namespace

bool efitInfo(const std::string& path, const Params& p, EfitInfo& all, EfitInfo& kept) {
    std::ifstream in(path);
    if (!in) return false;
    std::vector<double> Pa, da, Pk, dk;
    std::string line;
    bool header = false;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.size() >= 3 && line.compare(0, 3, "\xEF\xBB\xBF") == 0) line.erase(0, 3);
        if (line.empty() || line[0] == '#') continue;
        if (!header) { header = true; continue; }
        std::istringstream ls(line);
        Point q; char c1, c2;
        if (!(ls >> q.mass >> c1 >> q.delta >> c2 >> q.rms)) return false;
        const double P = q.mass * p.g / 1000;
        Pa.push_back(P); da.push_back(q.delta);
        if (q.rms <= p.rmsMax) { Pk.push_back(P); dk.push_back(q.delta); }
    }
    double A, I;
    tubeSection(p.tubeD, p.tubeT, A, I);
    const double L3 = p.beamL * p.beamL * p.beamL;
    auto fill = [&](const std::vector<double>& P, const std::vector<double>& d, EfitInfo& o) {
        const Fit f = lineFit(P, d, p.beamL, I);
        o.n = f.n; o.s = f.s; o.se = f.se; o.r2 = f.r2; o.E = f.E;
        o.lo = L3 / (48 * I * (f.s + f.se));
        o.hi = f.s > f.se ? L3 / (48 * I * (f.s - f.se)) : 0;
    };
    fill(Pa, da, all);
    fill(Pk, dk, kept);
    return kept.n >= 3;
}

bool efitRange(const std::string& path, const Params& p, double& E, double& lo, double& hi) {
    EfitInfo a, k;
    if (!efitInfo(path, p, a, k)) return false;
    E = k.E; lo = k.lo; hi = k.hi;
    return true;
}

bool runEfit(const std::string& path, const Params& p, std::ostream& out) {
    std::ifstream in(path);
    if (!in) { std::cerr << "cannot open " << path << '\n'; return false; }
    std::vector<Point> pts;
    std::string line;
    bool header = false;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.size() >= 3 && line.compare(0, 3, "\xEF\xBB\xBF") == 0) line.erase(0, 3);
        if (line.empty() || line[0] == '#') continue;
        if (!header) { header = true; continue; }   // mass_g,delta_mm,rms_px
        std::istringstream ls(line);
        Point q; char c1, c2;
        if (!(ls >> q.mass >> c1 >> q.delta >> c2 >> q.rms) || c1 != ',' || c2 != ',') {
            std::cerr << path << ": bad line '" << line << "'\n";
            return false;
        }
        pts.push_back(q);
    }
    double A, I;
    tubeSection(p.tubeD, p.tubeT, A, I);
    std::vector<double> Pa, da, Pk, dk;
    for (const Point& q : pts) {
        const double P = q.mass * p.g / 1000;
        Pa.push_back(P); da.push_back(q.delta);
        if (q.rms <= p.rmsMax) { Pk.push_back(P); dk.push_back(q.delta); }
    }
    const Fit all = lineFit(Pa, da, p.beamL, I);
    const Fit kept = lineFit(Pk, dk, p.beamL, I);

    out << "E 회귀 (3절)\n"
        << "  단면: 외경 " << p.tubeD << " mm, 두께 " << p.tubeT << " mm → A = " << A << " mm², I = " << I << " mm⁴\n"
        << "  스팬 L = " << p.beamL << " mm, g = " << p.g << " m/s²\n"
        << "  식: E = L³ / (48 I s),  s = Δδ/ΔP\n\n";
    auto show = [&](const char* name, const Fit& f) {
        const double Elo = p.beamL * p.beamL * p.beamL / (48 * I * (f.s + f.se));
        const double Ehi = p.beamL * p.beamL * p.beamL / (48 * I * (f.s - f.se));
        out << "  [" << name << "] 점 " << f.n << "개\n"
            << "    기울기 s = " << f.s << " mm/N (표준오차 " << f.se << "),  절편 = " << f.a << " mm,  R² = " << f.r2 << '\n'
            << "    E = " << f.E << " MPa  (기울기 ±1 표준오차 → " << Elo << " ~ " << Ehi << " MPa)\n";
    };
    show("전체", all);
    show(("판독 잔차 rms ≤ " + std::to_string(p.rmsMax).substr(0, 4) + " px").c_str(), kept);
    out << "  제외한 점:";
    for (const Point& q : pts) if (q.rms > p.rmsMax) out << ' ' << q.mass << " g (rms " << q.rms << " px)";
    out << '\n';
    return true;
}
