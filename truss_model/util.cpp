// 설정 읽기, 단면, 파일 쓰기 같은 공용 도구
#include "model.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>

namespace {

std::string trim(const std::string& s) {
    const auto a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    const auto b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

} // namespace

// "키 = 값" 줄을 읽는다. # 뒤는 주석. 모르는 키는 오류로 알린다 (오타 방지).
bool readConfig(const std::string& path, Params& p, std::string* err) {
    std::ifstream in(path);
    if (!in) { if (err) *err = "cannot open " + path; return false; }
    std::map<std::string, double*> dbl = {
        {"E", &p.E}, {"rho", &p.rho}, {"sigma_t", &p.sigmaT}, {"sigma_c", &p.sigmaC}, {"K", &p.K},
        {"tube_d", &p.tubeD}, {"tube_t", &p.tubeT},
        {"L", &p.L}, {"H", &p.H}, {"P", &p.P},
        {"H_min", &p.Hmin}, {"H_max", &p.Hmax}, {"H_step", &p.Hstep},
        {"L_max", &p.Lmax}, {"L_min", &p.Lmin}, {"min_angle", &p.minAngle},
        {"ga_H_min", &p.gaHmin}, {"ga_H_max", &p.gaHmax}, {"ga_H_step", &p.gaHstep},
        {"pc", &p.pc}, {"soft_p", &p.softP}, {"pm", &p.pm},
        {"beam_L", &p.beamL}, {"g", &p.g}, {"rms_max", &p.rmsMax},
    };
    std::map<std::string, int*> itg = {
        {"n", &p.n}, {"n_min", &p.nMin}, {"n_max", &p.nMax},
        {"max_deg", &p.maxDeg}, {"no_cross", &p.noCross},
        {"slot_cross", &p.slotCross}, {"ls_every", &p.lsEvery}, {"restart", &p.restart}, {"pop", &p.pop}, {"gens", &p.gens}, {"tour", &p.tour}, {"elite", &p.elite}, {"runs", &p.runs},
    };
    std::string line;
    int no = 0;
    while (std::getline(in, line)) {
        ++no;
        if (no == 1 && line.size() >= 3 && line.compare(0, 3, "\xEF\xBB\xBF") == 0) line.erase(0, 3);
        const auto hash = line.find('#');
        if (hash != std::string::npos) line.erase(hash);
        line = trim(line);
        if (line.empty()) continue;
        const auto eq = line.find('=');
        if (eq == std::string::npos) { if (err) *err = path + ":" + std::to_string(no) + ": expected 'key = value'"; return false; }
        const std::string key = trim(line.substr(0, eq)), val = trim(line.substr(eq + 1));
        std::istringstream vs(val);
        if (key == "seed") { if (!(vs >> p.seed)) goto bad; continue; }
        if (key == "material") {
            if (val != "custom" && !findMaterial(val)) { if (err) *err = path + ":" + std::to_string(no) + ": unknown material '" + val + "' (straw, steel, aluminum, douglas_fir, balsa, custom)"; return false; }
            p.material = val; continue;
        }
        if (key == "load") {
            if (val != "center" && val != "panel") { if (err) *err = path + ":" + std::to_string(no) + ": load must be center or panel"; return false; }
            p.load = val; continue;
        }
        if (key == "fitness") {
            if (val != "eta" && val != "etaUlt") { if (err) *err = path + ":" + std::to_string(no) + ": fitness must be eta or etaUlt"; return false; }
            p.fitness = val; continue;
        }
        if (dbl.count(key)) { if (!(vs >> *dbl[key])) goto bad; continue; }
        if (itg.count(key)) { if (!(vs >> *itg[key])) goto bad; continue; }
        if (err) *err = path + ":" + std::to_string(no) + ": unknown key '" + key + "'";
        return false;
    bad:
        if (err) *err = path + ":" + std::to_string(no) + ": bad value for '" + key + "'";
        return false;
    }
    // material 이 custom 이 아니면 목록의 값이 E·rho·sigma_t·sigma_c 줄보다 우선한다
    if (p.material != "custom") {
        const Material* m = findMaterial(p.material);
        p.E = m->E; p.rho = m->rho; p.sigmaT = m->sigmaT; p.sigmaC = m->sigmaC;
    }
    if (p.n < 2 || p.n % 2) { if (err) *err = "n must be an even number >= 2 (load sits on the middle bottom node)"; return false; }
    if (!(p.tubeT > 0 && 2 * p.tubeT <= p.tubeD)) { if (err) *err = "tube_t must be in (0, tube_d/2]"; return false; }
    if (!(p.gaHstep > 0 && p.gaHmax >= p.gaHmin && p.gaHmin > 0)) { if (err) *err = "need 0 < ga_H_min <= ga_H_max and ga_H_step > 0"; return false; }
    if (!(p.Hstep > 0 && p.Hmax >= p.Hmin && p.Hmin > 0)) { if (err) *err = "need 0 < H_min <= H_max and H_step > 0"; return false; }
    if (p.pop < 4 || p.gens < 1 || p.tour < 1 || p.elite < 0 || p.elite >= p.pop || p.runs < 1) { if (err) *err = "bad GA sizes (pop >= 4, gens >= 1, tour >= 1, 0 <= elite < pop, runs >= 1)"; return false; }
    return true;
}

// 원형 관: A = π t (d − t),  I = A (d² + dᵢ²) / 16,  dᵢ = d − 2t
// (π/4 (d² − dᵢ²) 를 인수분해한 꼴이라 얇은 벽에서도 뺄셈 오차가 없다)
void tubeSection(double d, double t, double& A, double& I) {
    const double pi = std::acos(-1.0);
    const double di = d - 2 * t;
    A = pi * t * (d - t);
    I = A * (d * d + di * di) / 16;
}

Member makeMember(int i, int j, const Params& p) {
    Member e;
    e.i = i; e.j = j;
    e.E = p.E; e.rho = p.rho; e.sigmaT = p.sigmaT; e.sigmaC = p.sigmaC; e.K = p.K;
    tubeSection(p.tubeD, p.tubeT, e.A, e.I);
    return e;
}

void applyLoad(Truss& t, int n, const Params& p) {
    const double P = std::fabs(p.P);
    if (p.load == "panel") {
        for (int k = 1; k < n; ++k) t.nodes[k].fy = -P / (n - 1);   // 바닥판 하중이 패널 점으로 고르게 전해지는 경우
    } else {
        t.nodes[n / 2].fy = -P;                                     // 하현 중앙 1점 (3절 실험과 같은 조건)
    }
}

std::string g17(double v) {
    char b[40];
    std::snprintf(b, sizeof b, "%.17g", v);
    return b;
}

bool writeTruss(const std::string& path, const Truss& t) {
    std::ofstream f(path);
    if (!f) return false;
    f << "TRUSS 1\n" << t.nodes.size() << ' ' << t.members.size() << '\n';
    for (const Node& q : t.nodes)
        f << g17(q.x) << ' ' << g17(q.y) << ' ' << q.fixX << ' ' << q.fixY << ' ' << g17(q.fx) << ' ' << g17(q.fy) << '\n';
    for (const Member& e : t.members)
        f << e.i << ' ' << e.j << ' ' << g17(e.E) << ' ' << g17(e.A) << ' ' << g17(e.I) << ' ' << g17(e.rho) << ' '
          << g17(e.sigmaT) << ' ' << g17(e.sigmaC) << ' ' << g17(e.K) << '\n';
    return static_cast<bool>(f);
}

bool openCsv(std::ofstream& f, const std::string& path) {
    f.open(path, std::ios::binary);
    if (!f) return false;
    f << "\xEF\xBB\xBF";
    return true;
}

const char* modeName(Mode m) {
    switch (m) {
        case Mode::Tension: return "인장";
        case Mode::Buckling: return "좌굴";
        case Mode::Crushing: return "압괴";
    }
    return "?";
}

const char* memberType(const Truss& t, int k) {
    const Node& a = t.nodes[t.members[k].i];
    const Node& b = t.nodes[t.members[k].j];
    const double tol = 1e-9 * (1 + std::fabs(a.x) + std::fabs(b.x));
    if (std::fabs(a.y - b.y) < tol) return std::fabs(a.y) < tol ? "하현" : "상현";
    if (std::fabs(a.x - b.x) < tol) return "수직";
    // 높이가 다른 상현 절점끼리 이은 경사 부재도 상현으로 본다.
    // 상현 절점 = 같은 x 에서 가장 높은 절점 (K-truss 의 중간 절점은 위에 상현 절점이 있어 제외된다)
    auto isTop = [&](const Node& q) {
        if (q.y <= tol) return false;
        for (const Node& o : t.nodes)
            if (std::fabs(o.x - q.x) < tol && o.y > q.y + tol) return false;
        return true;
    };
    if (isTop(a) && isTop(b)) return "상현";
    return "사재";
}

int maxwellDegree(const Truss& t) {
    int r = 0;
    for (const Node& q : t.nodes) r += q.fixX + q.fixY;
    return static_cast<int>(t.members.size()) + r - 2 * static_cast<int>(t.nodes.size());
}

// 부재력 그림: 인장 파랑, 압축 빨강, 굵기 ∝ |N|, 최초 파괴 부재는 굵은 테두리
bool writeSvg(const std::string& path, const Truss& t, const Result& r, const std::string& title) {
    if (t.nodes.empty()) return false;
    double x0 = 1e300, x1 = -1e300, y0 = 1e300, y1 = -1e300;
    for (const Node& q : t.nodes) { x0 = std::min(x0, q.x); x1 = std::max(x1, q.x); y0 = std::min(y0, q.y); y1 = std::max(y1, q.y); }
    const double W = 1000, pad = 60;
    const double sc = (W - 2 * pad) / std::max(x1 - x0, 1e-9);
    const double Hh = (y1 - y0) * sc + 2 * pad + 110;   // 아래쪽에 지점·하중 화살표·범례 자리
    auto X = [&](double x) { return pad + (x - x0) * sc; };
    auto Y = [&](double y) { return Hh - 90 - (y - y0) * sc; };
    double nmax = 0;
    const bool ok = r.status == Status::Ok && r.members.size() == t.members.size();
    if (ok) for (const auto& m : r.members) nmax = std::max(nmax, std::fabs(m.N));
    std::ofstream f(path);
    if (!f) return false;
    f << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << W << "\" height=\"" << Hh << "\" viewBox=\"0 0 " << W << ' ' << Hh
      << "\" font-family=\"sans-serif\"><rect width=\"100%\" height=\"100%\" fill=\"#fff\"/>\n";
    f << "<text x=\"" << pad << "\" y=\"28\" font-size=\"18\" fill=\"#1b2328\">" << title << "</text>\n";
    for (size_t k = 0; k < t.members.size(); ++k) {
        const Node& a = t.nodes[t.members[k].i];
        const Node& b = t.nodes[t.members[k].j];
        std::string col = "#8a9298";
        double w = 2;
        if (ok) {
            const double N = r.members[k].N;
            if (std::fabs(N) > 1e-9 * std::max(nmax, 1e-300)) col = N > 0 ? "#1f5fa8" : "#c0392b";
            w = 1.5 + 6 * std::fabs(N) / std::max(nmax, 1e-300);
        }
        if (ok && static_cast<int>(k) == r.critical)
            f << "<line x1=\"" << X(a.x) << "\" y1=\"" << Y(a.y) << "\" x2=\"" << X(b.x) << "\" y2=\"" << Y(b.y)
              << "\" stroke=\"#f2c14e\" stroke-width=\"" << w + 8 << "\" stroke-linecap=\"round\"/>\n";
        f << "<line x1=\"" << X(a.x) << "\" y1=\"" << Y(a.y) << "\" x2=\"" << X(b.x) << "\" y2=\"" << Y(b.y)
          << "\" stroke=\"" << col << "\" stroke-width=\"" << w << "\" stroke-linecap=\"round\"/>\n";
        f << "<text x=\"" << (X(a.x) + X(b.x)) / 2 + 4 << "\" y=\"" << (Y(a.y) + Y(b.y)) / 2 - 4
          << "\" font-size=\"11\" fill=\"#5b6870\">m" << k << "</text>\n";
    }
    for (size_t k = 0; k < t.nodes.size(); ++k) {
        const Node& q = t.nodes[k];
        if (q.fixX || q.fixY) {
            const double cx = X(q.x), cy = Y(q.y);
            f << "<path d=\"M" << cx << ' ' << cy << " L" << cx - 10 << ' ' << cy + 16 << " L" << cx + 10 << ' ' << cy + 16
              << " Z\" fill=\"" << (q.fixX ? "#1b2328" : "#fff") << "\" stroke=\"#1b2328\" stroke-width=\"1.5\"/>\n";
        }
        if (q.fx != 0 || q.fy != 0) {
            const double cx = X(q.x), cy = Y(q.y);
            f << "<line x1=\"" << cx << "\" y1=\"" << cy + 48 << "\" x2=\"" << cx << "\" y2=\"" << cy + 8
              << "\" stroke=\"#1b2328\" stroke-width=\"2\"/><path d=\"M" << cx << ' ' << cy + 18 << " L" << cx - 5 << ' ' << cy + 26
              << " L" << cx + 5 << ' ' << cy + 26 << " Z\" fill=\"#1b2328\" transform=\"rotate(180 " << cx << ' ' << cy + 22 << ")\"/>\n";
        }
        f << "<circle cx=\"" << X(q.x) << "\" cy=\"" << Y(q.y) << "\" r=\"4\" fill=\"#fff\" stroke=\"#1b2328\" stroke-width=\"1.5\"/>\n";
    }
    f << "<text x=\"" << pad << "\" y=\"" << Hh - 12 << "\" font-size=\"13\" fill=\"#5b6870\">파랑 = 인장, 빨강 = 압축, 굵기 ∝ |N|, 노란 테두리 = 최초 파괴 부재</text>\n";
    f << "</svg>\n";
    return static_cast<bool>(f);
}
