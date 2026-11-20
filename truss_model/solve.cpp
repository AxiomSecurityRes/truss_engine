// =====================================================================
//  solve.cpp — 평면 트러스 해석 (직접강성법)
// ---------------------------------------------------------------------
//  함수마다 [하는 일] [입력] [반환] [식] [확인값] 을 적어 두었다.
//  확인값은 모두 truss.txt (기본 삼각형) 기준이며, 전부 일치함을 확인했다.
//
//      절점 0 (0, 0)        핀      x, y 구속
//      절점 1 (200, 0)      롤러    y 만 구속
//      절점 2 (100, 141.42) 자유    fy = -10 N
//      부재 0: 0-1 (밑변)   부재 1: 0-2 (왼쪽 빗변)   부재 2: 1-2 (오른쪽 빗변)
//
//  자유도 번호 규칙 (이 파일 전체에서 사용):
//      절점 k 의 x 방향 = 2k,   y 방향 = 2k + 1
//      → 절점이 j 개면 자유도는 2j 개 (기본 삼각형: 0~5)
// =====================================================================

#include "truss.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace {

using Matrix = std::vector<std::vector<double>>;   // 조밀 행렬. K[행][열]
using Vector = std::vector<double>;

// =====================================================================
//  함수 1. geometry — 부재의 길이와 방향
// ---------------------------------------------------------------------
//  [하는 일] 부재 e 의 시작 절점 i 에서 끝 절점 j 로 가는 방향을 구한다.
//  [입력]    t: 트러스 전체 (절점 좌표를 꺼내려고),  e: 부재 하나
//  [반환]    Geo { L, c, s }
//              L = 부재 길이 [mm]
//              c = cosθ = dx / L   (x 방향 성분)
//              s = sinθ = dy / L   (y 방향 성분)
//  [식]      dx = x_j - x_i,  dy = y_j - y_i,  L = hypot(dx, dy)
//  [확인값]  부재 0: L = 200,     c =  1,        s = 0
//            부재 1: L = 173.205, c =  0.57735,  s = 0.816497
//            부재 2: L = 173.205, c = -0.57735,  s = 0.816497
// =====================================================================
struct Geo {
    double L;   // 길이 [mm]
    double c;   // cos
    double s;   // sin
};

Geo geometry(const Truss& t, const Member& e) {
    const Node& a = t.nodes[e.i];
    const Node& b = t.nodes[e.j];
    const double dx = b.x - a.x;
    const double dy = b.y - a.y;
    const double L  = std::hypot(dx, dy);   // io.cpp 가 길이 0 부재를 이미 거부했다
    return {L, dx / L, dy / L};
}

// =====================================================================
//  함수 2. memberDofs — 부재 양 끝의 자유도 번호 4개
// ---------------------------------------------------------------------
//  [하는 일] 부재의 국소 번호 0~3 이 전역 자유도 몇 번인지 알려 준다.
//  [입력]    e: 부재 하나
//  [반환]    { 2i, 2i+1, 2j, 2j+1 }   (i, j 는 e.i, e.j)
//              국소 0 = 절점 i 의 x,  국소 1 = 절점 i 의 y
//              국소 2 = 절점 j 의 x,  국소 3 = 절점 j 의 y
//  [확인값]  부재 2 (1-2): {2, 3, 4, 5}
// =====================================================================
std::array<int, 4> memberDofs(const Member& e) {
    return {2 * e.i, 2 * e.i + 1, 2 * e.j, 2 * e.j + 1};
}

// =====================================================================
//  함수 3. memberStiffness — 부재 하나의 4x4 강성행렬
// ---------------------------------------------------------------------
//  [하는 일] 용수철 하나(부재)가 양 끝 절점 변위 4개에 대해 내는 힘의 관계.
//  [입력]    e: 부재 (E, A 를 씀),  g: 함수 1 의 결과 (L, c, s 를 씀)
//  [반환]    4x4 행렬 ke  (ke[a][b] = 국소 b 번 변위를 1 줬을 때 국소 a 번에 생기는 힘)
//  [식]      k  = E·A / L                         (용수철 상수)
//            gv = { -c, -s, c, s }                (축 방향 벡터)
//            ke[a][b] = k · gv[a] · gv[b]         (a, b = 0..3)
//            → 풀어 쓰면 보조정리 1.2 의 [c², cs, -c², -cs; ...] 와 같다
//  [확인값]  부재 0: k = 20.6756,  gv = {-1, 0, 1, 0}
//            ke = 20.6756 × [ 1 0 -1 0 ;  0 0 0 0 ; -1 0 1 0 ;  0 0 0 0 ]
// =====================================================================
using Mat4 = std::array<std::array<double, 4>, 4>;

Mat4 memberStiffness(const Member& e, const Geo& g) {
    Mat4 ke{};
    const double k = e.E * e.A / g.L;
    const double gv[4] = {-g.c, -g.s, g.c, g.s};
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
            ke[a][b] = k * gv[a] * gv[b];
    return ke;
}

// =====================================================================
//  함수 4. assemble — 전역 강성행렬 K (2j x 2j)
// ---------------------------------------------------------------------
//  [하는 일] 모든 부재의 4x4 를 큰 행렬 K 의 제자리에 "더한다".
//  [입력]    t: 트러스 전체
//  [반환]    K (크기 2j x 2j)
//  [방법]    부재마다:
//              geo = geometry(t, e)
//              ke  = memberStiffness(e, geo)
//              d   = memberDofs(e)
//              K[d[a]][d[b]] += ke[a][b]     (a, b = 0..3)
//            ※ "=" 이 아니라 "+=" — 한 절점에 부재가 여러 개 붙기 때문
//  [확인값]  K[0][0] = 28.6337,  K[0][1] = 11.2544,  K[0][2] = -20.6756,
//            K[5][5] = 31.8322,  K[4][5] = 0
//            대칭: K[a][b] == K[b][a]
//            통째 이동: u = {1,0,1,0,1,0} 이면 K·u 의 모든 칸 = 0
// =====================================================================
Matrix assemble(const Truss& t) {
    const int n = 2 * static_cast<int>(t.nodes.size());
    Matrix K(n, Vector(n, 0.0));
    for (const Member& e : t.members) {
        const Geo g = geometry(t, e);
        const Mat4 ke = memberStiffness(e, g);
        const std::array<int, 4> d = memberDofs(e);
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                K[d[a]][d[b]] += ke[a][b];
    }
    return K;
}

// =====================================================================
//  함수 5. loadVector — 하중 벡터 F (크기 2j)
// ---------------------------------------------------------------------
//  [하는 일] 절점 하중을 자유도 번호 순서로 늘어놓는다.
//  [입력]    t: 트러스 전체
//  [반환]    F,  F[2k] = 절점 k 의 fx,  F[2k+1] = 절점 k 의 fy
//  [확인값]  F = {0, 0, 0, 0, 0, -10}
// =====================================================================
Vector loadVector(const Truss& t) {
    Vector F(2 * t.nodes.size(), 0.0);
    for (std::size_t k = 0; k < t.nodes.size(); ++k) {
        F[2 * k]     = t.nodes[k].fx;
        F[2 * k + 1] = t.nodes[k].fy;
    }
    return F;
}

// =====================================================================
//  함수 6. freeDofs — 구속되지 않은 자유도 번호 목록
// ---------------------------------------------------------------------
//  [하는 일] 지점에 묶이지 않아 "모르는 값"인 변위의 번호만 모은다.
//  [입력]    t: 트러스 전체
//  [반환]    오름차순 번호 목록
//              절점 k 의 fixX == 0 이면 2k 를,  fixY == 0 이면 2k+1 을 넣는다
//  [확인값]  {2, 4, 5}
//            (절점 0 은 둘 다 구속 → 없음, 절점 1 은 x 만 자유 → 2, 절점 2 → 4, 5)
// =====================================================================
std::vector<int> freeDofs(const Truss& t) {
    std::vector<int> f;
    for (int k = 0; k < static_cast<int>(t.nodes.size()); ++k) {
        if (t.nodes[k].fixX == 0) f.push_back(2 * k);
        if (t.nodes[k].fixY == 0) f.push_back(2 * k + 1);
    }
    return f;
}

// =====================================================================
//  함수 7. subMatrix / subVector — 자유 자유도만 뽑기
// ---------------------------------------------------------------------
//  [하는 일] K 와 F 에서 free 목록에 있는 행·열만 남긴다.
//  [입력]    K (또는 F),  free: 함수 6 의 결과 (크기 m)
//  [반환]    Kff (m x m):  Kff[p][q] = K[free[p]][free[q]]
//            Ff  (m)    :  Ff[p]     = F[free[p]]
//  [확인값]  Kff = [ 28.6337  -7.9580  11.2544 ;
//                    -7.9580  15.9161   0      ;
//                    11.2544   0       31.8322 ]
//            Ff  = {0, 0, -10}
// =====================================================================
Matrix subMatrix(const Matrix& K, const std::vector<int>& free) {
    const int m = static_cast<int>(free.size());
    Matrix Kff(m, Vector(m, 0.0));
    for (int p = 0; p < m; ++p)
        for (int q = 0; q < m; ++q)
            Kff[p][q] = K[free[p]][free[q]];
    return Kff;
}

Vector subVector(const Vector& F, const std::vector<int>& free) {
    Vector Ff(free.size(), 0.0);
    for (std::size_t p = 0; p < free.size(); ++p)
        Ff[p] = F[free[p]];
    return Ff;
}

// =====================================================================
//  함수 8. gaussSolve — A·x = b 를 푼다 (부분 피벗 가우스 소거)
// ---------------------------------------------------------------------
//  [하는 일] 연립방정식을 풀고, 풀 수 없으면(특이 = 메커니즘) 알려 준다.
//  [입력]    A (n x n), b (n)  ← 값으로 받으므로 안에서 마음대로 고쳐도 된다
//  [출력]    x (n)  ← 참조로 받은 변수에 답을 채운다
//  [반환]    true  = 풀었다
//            false = 피벗이 0 에 가까워 풀 수 없다 (구조물이 메커니즘)
//  [방법]
//    준비)  scale = max |A[i][i]|   (소거 전에 한 번)
//    전진 소거) col = 0 .. n-1 마다
//        ① 피벗 고르기: r = col..n-1 중 |A[r][col]| 이 가장 큰 행
//        ② 판정: 그 값 < 1e-12 * scale 이면 return false
//        ③ 행 교환: std::swap(A[col], A[r]);  std::swap(b[col], b[r]);
//        ④ 아래 행 지우기: r = col+1 .. n-1 마다
//               m = A[r][col] / A[col][col]
//               A[r][k] -= m * A[col][k]   (k = col .. n-1)
//               b[r]    -= m * b[col]
//    후진 대입) i = n-1 .. 0 (거꾸로)
//        x[i] = ( b[i] - Σ_{k>i} A[i][k]·x[k] ) / A[i][i]
//  [확인값]  A = {{2,1,-1},{-3,-1,2},{-2,1,2}}, b = {8,-11,-3} → true,  x = {2, 3, -1}
//            A = {{1,2},{2,4}},                 b = {3,6}      → false
// =====================================================================
bool gaussSolve(Matrix A, Vector b, Vector& x) {
    const int n = static_cast<int>(b.size());
    x.assign(n, 0.0);
    // 준비: 피벗 판정용 기준 크기
    double scale = 0.0;
    for (int i = 0; i < n; ++i) scale = std::max(scale, std::abs(A[i][i]));
    if (n > 0 && scale == 0.0) return false;          // 대각이 전부 0 → 강성이 전혀 없음
    const double tol = 1e-12 * scale;

    // 전진 소거
    for (int col = 0; col < n; ++col) {
        // ① 피벗 고르기
        int piv = col;
        for (int r = col + 1; r < n; ++r)
            if (std::abs(A[r][col]) > std::abs(A[piv][col])) piv = r;
        // ② 판정
        if (std::abs(A[piv][col]) < tol) return false;  // 메커니즘
        // ③ 행 교환
        std::swap(A[col], A[piv]);
        std::swap(b[col], b[piv]);
        // ④ 아래 행 지우기
        for (int r = col + 1; r < n; ++r) {
            const double m = A[r][col] / A[col][col];
            if (m == 0.0) continue;
            for (int k = col; k < n; ++k) A[r][k] -= m * A[col][k];
            b[r] -= m * b[col];
        }
    }

    // 후진 대입
    for (int i = n - 1; i >= 0; --i) {
        double sum = b[i];
        for (int k = i + 1; k < n; ++k) sum -= A[i][k] * x[k];
        x[i] = sum / A[i][i];
    }
    return true;
}

// =====================================================================
//  함수 9. multiply — 행렬 x 벡터
// ---------------------------------------------------------------------
//  [하는 일] y = K·u  (반력과 잔차 계산에 쓴다)
//  [입력]    K (n x n), u (n)
//  [반환]    y (n),  y[a] = Σ_b K[a][b]·u[b]
//  [확인값]  K·{1,0,1,0,1,0} = {0,0,0,0,0,0}   (함수 4 의 통째 이동 검사와 같음)
// =====================================================================
Vector multiply(const Matrix& K, const Vector& u) {
    Vector y(u.size(), 0.0);
    for (std::size_t a = 0; a < u.size(); ++a)
        for (std::size_t b = 0; b < u.size(); ++b)
            y[a] += K[a][b] * u[b];
    return y;
}

// =====================================================================
//  함수 10. axialForce — 부재 축력 (보조정리 1.3)
// ---------------------------------------------------------------------
//  [하는 일] 절점 변위로 부재가 얼마나 늘었는지 보고 축력을 구한다.
//  [입력]    e: 부재,  g: 함수 1 의 결과,  u: 전체 변위 벡터 (크기 2j)
//  [반환]    N [N]   (인장 +, 압축 -)
//  [식]      늘어난 길이 dL = c·(u[2j] - u[2i]) + s·(u[2j+1] - u[2i+1])
//            N = (E·A / L) · dL
//  [확인값]  부재 0:  3.5355339  (= 5/√2,  인장)
//            부재 1: -6.1237244  (압축)
//            부재 2: -6.1237244  (압축)
// =====================================================================
double axialForce(const Member& e, const Geo& g, const Vector& u) {
    const double dL = g.c * (u[2 * e.j]     - u[2 * e.i])
                    + g.s * (u[2 * e.j + 1] - u[2 * e.i + 1]);
    return e.E * e.A / g.L * dL;
}

// =====================================================================
//  함수 11. capacity — 부재가 버틸 수 있는 축력과 파괴 모드 (보조정리 1.4)
// ---------------------------------------------------------------------
//  [하는 일] 지금 축력의 부호(인장/압축)에 따라 용량을 정한다.
//  [입력]    e: 부재,  L: 길이,  N: 함수 10 의 축력
//  [반환]    Capacity { cap, mode }
//              cap  = 버틸 수 있는 축력의 크기 [N] (양수)
//              mode = 그 용량을 정한 파괴 방식
//  [식]      N >= 0 (인장):  cap = σt·A,                     mode = Tension
//            N <  0 (압축):  Pcr = π²·E·I / (K·L)²
//                            cap = min(Pcr, σc·A)
//                            Pcr <= σc·A 이면 mode = Buckling (좌굴)
//                            아니면           mode = Crushing (압괴)
//            π 는 std::acos(-1.0)   (3.14 를 쓰면 좌굴하중이 0.1% 작게 나온다)
//  [확인값]  부재 0: cap = 82.702427, Tension
//            부재 1, 2: Pcr = 5.8233639 < σc·A = 110.27 → cap = 5.8233639, Buckling
// =====================================================================
struct Capacity {
    double cap;
    Mode   mode;
};

Capacity capacity(const Member& e, double L, double N) {
    if (N >= 0.0)
        return {e.sigmaT * e.A, Mode::Tension};

    const double pi  = std::acos(-1.0);
    const double KL  = e.K * L;
    const double Pcr = pi * pi * e.E * e.I / (KL * KL);
    const double Pc  = e.sigmaC * e.A;
    if (Pcr <= Pc) return {Pcr, Mode::Buckling};
    return {Pc, Mode::Crushing};
}

} // namespace

// =====================================================================
//  solve — 위 함수들을 순서대로 불러 Result 를 채운다
// ---------------------------------------------------------------------
//  [입력]  t: readTruss 가 읽은 트러스 (입력 검사는 io.cpp 가 이미 했다)
//  [반환]  Result  (truss.h 참고. 각 칸의 뜻과 단위가 적혀 있다)
//
//  [최종 확인값]  README 의 정답 8 줄과 상대오차 1e-12 이내
//      RESULT 1
//      0 3 3
//      0 0 0 5
//      0.1710002499823563 0 0 5
//      0.08550012499117815 -0.37460523693241532 0 0
//      3.5355339059327369 1.2825018748676722 82.702426605751299 0.042750062495589075 0
//      -6.1237243569579443 -2.2213584080731512 5.8233639166534443 1.0515785110811193 1
//      -6.1237243569579443 -2.2213584080731512 5.8233639166534443 1.0515785110811193 1
//      0.95095134548908589 9.5095134548908593 1 1.3632149627567633 6.9757989126382771 0
// =====================================================================
Result solve(const Truss& t) {
    Result r;
    const int j = static_cast<int>(t.nodes.size());
    const int n = 2 * j;

    // ---- 단계 A. 행렬과 하중 만들기 ---------------------------------
    const Matrix K = assemble(t);
    const Vector F = loadVector(t);

    // ---- 단계 B. 자유 자유도만 남긴 작은 연립방정식 ------------------
    const std::vector<int> free = freeDofs(t);
    const Matrix Kff = subMatrix(K, free);
    const Vector Ff  = subVector(F, free);

    // ---- 단계 C. 풀기 -------------------------------------------------
    Vector uf;
    if (!gaussSolve(Kff, Ff, uf)) {
        r.status = Status::Mechanism;   // 불안정 구조: 나머지 칸은 비운 채로 돌려준다
        return r;
    }

    // ---- 단계 D. 전체 변위 복원 → r.ux, r.uy -------------------------
    //      확인: ux = {0, 0.1710003, 0.0855001},  uy = {0, 0, -0.3746052}
    Vector u(n, 0.0);
    for (std::size_t p = 0; p < free.size(); ++p) u[free[p]] = uf[p];
    r.ux.assign(j, 0.0);
    r.uy.assign(j, 0.0);
    for (int k = 0; k < j; ++k) {
        r.ux[k] = u[2 * k];
        r.uy[k] = u[2 * k + 1];
    }

    // ---- 단계 E. 반력 R = K·u − F (구속된 방향만) → r.rx, r.ry --------
    //      확인: rx = {0, 0, 0},  ry = {5, 5, 0}
    const Vector Ku = multiply(K, u);
    r.rx.assign(j, 0.0);
    r.ry.assign(j, 0.0);
    for (int k = 0; k < j; ++k) {
        if (t.nodes[k].fixX) r.rx[k] = Ku[2 * k]     - F[2 * k];
        if (t.nodes[k].fixY) r.ry[k] = Ku[2 * k + 1] - F[2 * k + 1];
    }

    // ---- 단계 F. 잔차: 자유 자유도에서 max |K·u − F| ------------------
    //      확인: 1e-12 이하
    r.residual = 0.0;
    for (int a : free) r.residual = std::max(r.residual, std::abs(Ku[a] - F[a]));

    // ---- 단계 G. 부재마다 축력, 응력, 용량, DCR → r.members ------------
    r.members.reserve(t.members.size());
    for (const Member& e : t.members) {
        const Geo g = geometry(t, e);
        MemberResult mr;
        mr.N      = axialForce(e, g, u);
        mr.stress = mr.N / e.A;
        const Capacity c = capacity(e, g.L, mr.N);
        mr.cap  = c.cap;
        mr.mode = c.mode;
        mr.dcr  = std::abs(mr.N) / mr.cap;
        r.members.push_back(mr);
    }

    // ---- 단계 H. 하중계수, 파괴하중, 위험 부재 -------------------------
    //      확인: loadFactor = 0.9509513, critical = 1, Pfail = 9.5095135
    double maxDcr = 0.0;
    for (std::size_t i = 0; i < r.members.size(); ++i) {
        if (r.members[i].dcr > maxDcr) {          // '>' 이므로 동점이면 작은 번호 유지
            maxDcr = r.members[i].dcr;
            r.critical = static_cast<int>(i);
        }
    }
    if (maxDcr > 0.0) {
        double sumFx = 0.0, sumFy = 0.0;
        for (const Node& nd : t.nodes) { sumFx += nd.fx; sumFy += nd.fy; }
        r.loadFactor = 1.0 / maxDcr;
        r.Pfail      = r.loadFactor * std::hypot(sumFx, sumFy);
    }
    // maxDcr == 0 (하중 없음): loadFactor = Pfail = 0, critical = -1 그대로

    // ---- 단계 I. 질량과 효율 ------------------------------------------
    //      확인: mass = 1.3632150 g,  eta = 6.9757989 N/g
    r.mass = 0.0;
    for (const Member& e : t.members)
        r.mass += e.rho * e.A * geometry(t, e).L / 1000.0;   // g/cm³ × mm³ / 1000 = g
    if (r.mass > 0.0) r.eta = r.Pfail / r.mass;

    // ---- 단계 J. 끝 ---------------------------------------------------
    r.status = Status::Ok;
    return r;
}
