# truss_engine

평면 트러스 직접강성법(DSM) 해석 엔진. 트러스 도면판에서 내보낸 TXT를 읽어 해석 결과를 씁니다.

## 파일

| 파일 | 내용 |
|---|---|
| `truss.h` | 구조체(`Node`, `Member`, `Truss`, `Result`)와 함수 선언 |
| `io.cpp` | `readTruss`, `writeResult` |
| `main.cpp` | `truss.txt` 읽기 → `solve` → `result.txt` 쓰기 |
| `solve.cpp` | `Result solve(const Truss&)` — 해석 본체 (작성 중) |
| `truss.txt` | 검증용 입력: 기본 삼각형 (밑변 200 mm, 높이 100√2 mm, 꼭짓점 하중 10 N) |

## 빌드

```
g++ -std=c++17 -Wall -O2 *.cpp -o truss_engine
./truss_engine [입력=truss.txt] [출력=result.txt]
```

## 입력 형식 (`TRUSS 1`)

```
TRUSS 1
j m
x y fix_x fix_y fx fy               ← 절점 j줄 (0부터 번호)
i j E A I rho sigma_t sigma_c K     ← 부재 m줄
```

단위: 길이 mm, 힘 N, E·강도 MPa, A mm², I mm⁴, 밀도 g/cm³. 좌표 x 오른쪽 +, y 위쪽 +. 지점 1 = 구속.

## 출력 형식 (`RESULT 1`)

```
RESULT 1
status j m                           ← status: 0 Ok, 1 Mechanism, 2 InvalidInput
ux uy rx ry                          ← 절점 j줄
N stress cap dcr mode                ← 부재 m줄 (mode: 0 인장, 1 좌굴, 2 압괴)
loadFactor Pfail critical mass eta residual
```
