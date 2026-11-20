#ifndef TRUSS_H_INCLUDED
#define TRUSS_H_INCLUDED
#include <istream>
#include <ostream>
#include <string>
#include <vector>

// 단위: 길이 mm, 힘 N, E·강도 MPa, A mm², I mm⁴, 밀도 g/cm³
// 좌표: x 오른쪽 +, y 위쪽 +. 지점 1 = 그 방향 구속

struct Node {
    double x = 0, y = 0;
    int    fixX = 0, fixY = 0;
    double fx = 0, fy = 0;
};

struct Member {
    int    i = 0, j = 0;       // 양 끝 절점 번호 (i -> j 방향)
    double E = 0, A = 0, I = 0, rho = 0;
    double sigmaT = 0, sigmaC = 0, K = 1;
};

struct Truss {
    std::vector<Node>   nodes;
    std::vector<Member> members;
};

enum class Status { Ok = 0, Mechanism = 1, InvalidInput = 2 };
enum class Mode   { Tension = 0, Buckling = 1, Crushing = 2 };

struct MemberResult {
    double N = 0;        // 축력 [N], 인장 +
    double stress = 0;   // N / A [MPa]
    double cap = 0;      // 용량 [N], 양수
    double dcr = 0;      // |N| / cap
    Mode   mode = Mode::Tension;
};

struct Result {
    Status status = Status::InvalidInput;
    std::vector<double> ux, uy;          // 절점 변위 [mm], 크기 j
    std::vector<double> rx, ry;          // 반력 [N], 크기 j (구속 안 된 방향은 0)
    std::vector<MemberResult> members;   // 크기 m, 입력과 같은 순서
    double loadFactor = 0;   // 1 / max dcr
    double Pfail = 0;        // loadFactor × |하중 합력| [N]
    int    critical = -1;    // dcr 최대 부재 (동점이면 작은 번호)
    double mass = 0;         // Σ ρAL / 1000 [g]
    double eta = 0;          // Pfail / mass [N/g]
    double residual = 0;     // 자유 자유도에서 max |K·u − F| [N]
};

// 도면판 TXT(TRUSS 1)를 읽는다. 실패하면 false, err가 있으면 이유를 적는다.
bool   readTruss  (std::istream& in, Truss& t, std::string* err = nullptr);
Result solve      (const Truss& t);
// RESULT 1 형식으로 쓴다. Result 안의 배열 크기가 서로 맞지 않으면
// 읽을 수 있는 InvalidInput 결과를 대신 쓰고 false를 돌려준다.
bool   writeResult(std::ostream& out, const Result& r);

#endif // TRUSS_H_INCLUDED
