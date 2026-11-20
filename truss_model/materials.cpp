// 재료 목록과 출처
//
// 모든 재료에 같은 단면(외경 5 mm × 두께 0.1366 mm 관)과 같은 형상을 쓴다.
// 그래서 부재의 세장비 KL/r 이 재료와 상관없이 같고, 재료에 따라 달라지는 것은
//   좌굴하중 ∝ E,  압괴·인장 용량 ∝ σ,  질량 ∝ ρ  뿐이다.
// 좌굴과 압괴가 바뀌는 경계 세장비 π√(E/σc) 가 재료마다 달라서, 같은 형상이라도 파괴 모드가 달라질 수 있다.
//
// 금속은 σt = σc = 항복강도로 둔다 (설계에서 흔히 쓰는 가정).
// 목재는 결 방향 값이고, 밀도 = 비중 × (1 + 함수율 0.12) 로 바꿨다 (Wood Handbook 의 비중은 전건 무게 기준).
#include "model.h"
#include <cstring>

const Material MATERIALS[] = {
    {"straw", "PP 빨대 (실측)", 1630, 0.905, 35.5, 39.6,
     "E: 3절 실측 / ρ: INEOS, Typical Engineering Properties of Polypropylene (2014) / σ: Plastics International, Polypropylene Homopolymer (5,150·5,750 psi)"},
    {"steel", "구조용 강 S275", 210000, 7.85, 275, 275,
     "E·항복강도: SteelConstruction.info, Steel material properties (EN 1993-1-1, t ≤ 16 mm) / ρ: 강재 일반값 7850 kg/m³"},
    {"aluminum", "알루미늄 6061-T6", 68900, 2.70, 276, 276,
     "Aerospace Metals, Aluminum 6061-T6 자료 (Aluminum Association 대표값): E 68.9 GPa, 항복 276 MPa, 2.70 g/cm³. σc 는 σt 와 같다고 가정"},
    {"douglas_fir", "더글러스퍼 (목재, 결 방향)", 13400, 0.538, 107.6, 49.9,
     "USDA Wood Handbook 표 5-3a (12 % 함수율): 비중 0.48, E 13.4 GPa, 압축 49.9 MPa / 인장: 표 5-7 (interior north, 생재) 107.6 MPa"},
    {"balsa", "발사 (목재, 결 방향)", 3400, 0.179, 14.9, 14.9,
     "USDA Wood Handbook 표 5-5a (12 % 함수율): 비중 0.16, E 3.4 GPa, 압축 14.9 MPa / 인장 자료 없음 → 압축과 같다고 가정 (보수적)"},
};
const int MATERIAL_COUNT = sizeof(MATERIALS) / sizeof(MATERIALS[0]);

const Material* findMaterial(const std::string& key) {
    for (int i = 0; i < MATERIAL_COUNT; ++i)
        if (key == MATERIALS[i].key) return &MATERIALS[i];
    return nullptr;
}
