#ifndef MODEL_H_INCLUDED
#define MODEL_H_INCLUDED
// 논문용 확장: 설정, 대표 4종 생성, 연쇄 파괴, 제작 제약, 유전 알고리즘, E 회귀
// 해석 본체는 solve.cpp (truss.h) 를 그대로 쓴다.
#include "truss.h"
#include <string>
#include <vector>

// ---------------------------------------------------------------------
//  설정 (config.txt 의 "키 = 값" 으로 바꿀 수 있다)
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
//  재료 목록 (출처는 materials.cpp). 모든 재료에 같은 단면(빨대 관)을 써서 형상과 세장비를 같게 둔다.
// ---------------------------------------------------------------------
struct Material {
    const char* key;      // config 의 material = 값
    const char* name;
    double E, rho, sigmaT, sigmaC;   // MPa, g/cm³, MPa, MPa
    const char* source;
};
extern const Material MATERIALS[];
extern const int MATERIAL_COUNT;
const Material* findMaterial(const std::string& key);

struct Params {
    // 재료·단면: 3절 실험의 PP 빨대
    double E = 1630;          // 탄성계수 [MPa] (3절 실측)
    double rho = 0.905;       // 밀도 [g/cm³]
    double sigmaT = 35.5;     // 인장강도 [MPa] (5,150 psi)
    double sigmaC = 39.6;     // 압축강도 [MPa] (5,750 psi)
    double K = 1;             // 유효길이계수
    std::string material = "straw";   // 재료 이름 (custom 이면 위 값을 그대로)
    double tubeD = 5;         // 빨대 외경 [mm]
    double tubeT = 0.1366;    // 빨대 벽두께 [mm]

    // 다리 (하현 중앙 1점 하중, 왼쪽 핀·오른쪽 롤러)
    double L = 600;           // 스팬 [mm]
    double H = 100;           // 높이 [mm]
    int    n = 6;             // 패널 수 (짝수)
    double P = 10;            // 하중 [N] (합력)
    std::string load = "center"; // center: 하현 중앙 1점, panel: 하현 안쪽 절점 n−1개에 P/(n−1)씩

    // 높이 스윕
    double Hmin = 30, Hmax = 300, Hstep = 5;
    int    nMin = 4, nMax = 20;

    // 제작 제약 (빨대 + 이쑤시개 핀)
    double Lmax = 180;        // 부재 최대 길이 [mm] (빨대 한 개 길이)
    double Lmin = 20;         // 부재 최소 길이 [mm]
    double minAngle = 15;     // 한 절점에서 두 부재 사이 최소 각 [deg]
    int    maxDeg = 6;        // 한 절점에 꽂을 수 있는 부재 수
    int    noCross = 1;       // 1: 부재끼리 교차 금지

    // 유전 알고리즘
    double gaHmin = 30, gaHmax = 200, gaHstep = 5;   // 상현 절점 높이 범위와 간격 [mm]
    int    pop = 200, gens = 1000, tour = 3, elite = 2, runs = 10;
    unsigned seed = 1;
    double pc = 0.9;          // 교차 확률
    double pm = 0.08;         // 유전자당 돌연변이 확률 (-1 이면 1 / 유전자 수)
    int    lsEvery = 25;      // 몇 세대마다 최고 개체를 국소 탐색할지 (0 = 끔)
    int    slotCross = 1;     // 1: 슬롯 단위 교차, 0: 유전자별 균등 교차
    double softP = 8;         // 탐색용 점수의 p-노름 지수 (0 = 끔, 참 η 로 고름)
    int    restart = 150;     // 이 세대 수 동안 나아지지 않으면 엘리트만 남기고 다시 뿌림 (0 = 끔)
    std::string fitness = "eta";   // eta = P_fail/질량, etaUlt = P_ult/질량

    // E 회귀 (3절)
    double beamL = 160;       // 단순보 스팬 [mm]
    double g = 9.80665;       // 중력가속도 [m/s²]
    double rmsMax = 3;        // 판독 잔차가 이보다 크면 제외 [px]
};

bool readConfig(const std::string& path, Params& p, std::string* err);
void tubeSection(double d, double t, double& A, double& I);   // A [mm²], I [mm⁴]
Member makeMember(int i, int j, const Params& p);
// 하현 절점 0..n (지점 0, n) 에 하중 조건대로 하중을 건다
void applyLoad(Truss& t, int n, const Params& p);

// TRUSS 1 형식으로 쓴다 (도면판에서 불러올 수 있다)
bool writeTruss(const std::string& path, const Truss& t);
// 부재력을 색과 굵기로 그린 SVG
bool writeSvg(const std::string& path, const Truss& t, const Result& r, const std::string& title);
// 엑셀에서 한글이 깨지지 않게 UTF-8 BOM 을 붙여 연다
bool openCsv(std::ofstream& f, const std::string& path);
std::string g17(double v);              // 유효숫자 17자리 (double 손실 없음)
const char* modeName(Mode m);           // 인장 / 좌굴 / 압괴
const char* memberType(const Truss& t, int k);   // 하현 / 상현 / 수직 / 사재
int  maxwellDegree(const Truss& t);     // m + r − 2j

// ---------------------------------------------------------------------
//  대표 4종 (도면판과 절점·부재 번호가 같다)
// ---------------------------------------------------------------------
enum class Kind { Warren, Pratt, Howe, KTruss };
extern const Kind ALL_KINDS[4];
const char* kindName(Kind k);
bool parseKind(const std::string& s, Kind& k);
Truss makePreset(Kind kind, int n, double L, double H, const Params& p);

// ---------------------------------------------------------------------
//  연쇄 파괴 (부재가 용량에 닿으면 제거하고 다시 푼다)
// ---------------------------------------------------------------------
struct CollapseEvent {
    double lam = 0, P = 0;     // 하중계수, 그때 하중 [N]
    int    member = -1;        // 원래 트러스의 부재 번호
    Mode   mode = Mode::Tension;
    bool   cascade = false;    // 하중을 더 올리지 않고 이어진 파괴
};
struct CollapseResult {
    bool ok = false;           // 처음 트러스가 안정해서 해석이 되었는가
    std::vector<CollapseEvent> events;
    double lamFirst = 0, lamUlt = 0, Pfirst = 0, Pult = 0;
    std::string end;           // mechanism | all-removed | unstable
};
CollapseResult progressiveCollapse(const Truss& t);

// ---------------------------------------------------------------------
//  제작 가능성 검사
// ---------------------------------------------------------------------
struct FabReport {
    double violation = 0;      // 0 이면 제작 가능
    int tooLong = 0, tooShort = 0, overDeg = 0, sharp = 0, crossings = 0;
    double maxLen = 0, minLen = 0, minAngle = 180;
    int maxDeg = 0;
};
FabReport checkFab(const Truss& t, const Params& p);

// ---------------------------------------------------------------------
//  유전 알고리즘 (위상 + 형상 최적화)
// ---------------------------------------------------------------------
struct Genome {
    std::vector<unsigned char> act;   // 상현 후보 절점 사용 여부 (왼쪽 절반 + 중앙)
    std::vector<int>           h;     // 상현 후보 절점 높이 단계
    std::vector<unsigned char> web;   // 웨브 부재 사용 여부
};
struct Eval {
    bool   feasible = false;
    double fit = 0;            // 적합도 (η 또는 P_ult/질량) — 보고하는 값
    double score = 0;          // 선택에 쓰는 부드러운 적합도 (soft_p 참고)
    double viol = 0;           // 제약 위반량 (작을수록 좋음)
    double eta = 0, Pfail = 0, mass = 0, Pult = 0;
    int    status = 2, deg = 0, j = 0, m = 0;
};
class GaSpace {
public:
    explicit GaSpace(const Params& p);
    int    slots() const { return n_; }
    int    levels() const { return levels_; }
    int    webGenes() const { return static_cast<int>(genes_.size()); }
    int    totalGenes() const { return 2 * n_ + webGenes(); }
    double log10Size() const;
    double height(int idx) const;
    Truss  decode(const Genome& g) const;
    Eval   evaluate(const Genome& g) const;
    Genome random(unsigned long long& state) const;
    // 대표 형태를 유전자로 나타낸다 (검증용). Warren, Pratt, Howe 만 가능 (K는 중간 절점이 없어 탐색 공간 밖)
    bool   encodePreset(Kind k, double H, Genome& g) const;
    std::string describe(const Genome& g) const;
    std::vector<int> slotGenes(int s) const;   // 슬롯 s (1..n) 에 딸린 웨브 유전자 번호
private:
    struct WebGene { int slot; std::vector<std::pair<int, int>> links; };   // (상현 슬롯, 하현 절점)
    Params p_;
    int n_ = 0, levels_ = 0;
    std::vector<WebGene> genes_;
};
struct GaRun {
    Genome best;
    Eval   bestEval;
    int    foundGen = 0;
    double gaFit = 0;          // 국소 탐색 전 GA 최고 적합도
    int    lsSteps = 0;        // 국소 탐색으로 나아진 횟수
    int    restarts = 0;       // 정체로 다시 뿌린 횟수
    long long evals = 0;
    std::vector<double> bestFit, meanFit, feasRate;   // 세대별
};
GaRun runGa(const GaSpace& space, const Params& p, unsigned seed);

// ---------------------------------------------------------------------
//  3절 E 회귀: deflection.csv (mass_g, delta_mm, rms_px) → 결과를 out 에 쓴다
// ---------------------------------------------------------------------
bool runEfit(const std::string& path, const Params& p, std::ostream& out);
struct EfitInfo { int n = 0; double s = 0, se = 0, r2 = 0, E = 0, lo = 0, hi = 0; };   // 점 수, 기울기·표준오차 [mm/N], R², E·범위 [MPa]
// 전체 점과 rms_max 이하 점 두 가지 회귀
bool efitInfo(const std::string& path, const Params& p, EfitInfo& all, EfitInfo& kept);
// 같은 회귀에서 E 와 기울기 ±1 표준오차에 해당하는 E 범위
bool efitRange(const std::string& path, const Params& p, double& E, double& lo, double& hi);

#endif // MODEL_H_INCLUDED
