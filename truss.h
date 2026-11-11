#ifndef TRUSS_H_INCLUDED
#define TRUSS_H_INCLUDED
#include <istream>
#include <ostream>
#include <vector>

struct Node   { double x, y; int fixX, fixY; double fx, fy; };
struct Member { int i, j; double E, A, I, rho, sigmaT, sigmaC, K; };
struct Truss  { std::vector<Node> nodes; std::vector<Member> members; };

enum class Status { Ok = 0, Mechanism = 1, InvalidInput = 2 };
enum class Mode   { Tension = 0, Buckling = 1, Crushing = 2 };

struct MemberResult { double N, stress, cap, dcr; Mode mode; };

struct Result {
    Status status = Status::InvalidInput;
    std::vector<double> ux, uy, rx, ry;
    std::vector<MemberResult> members;
    double loadFactor = 0, Pfail = 0;
    int    critical = -1;
    double mass = 0, eta = 0, residual = 0;
};

bool   readTruss  (std::istream& in, Truss& t);
Result solve      (const Truss& t);
void   writeResult(std::ostream& out, const Result& r);
#endif // TRUSS_H_INCLUDED
