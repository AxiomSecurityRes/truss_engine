# truss_engine

평면 트러스 직접강성법(DSM) 해석 엔진. 트러스 도면판에서 내보낸 TXT를 읽어 해석 결과를 씁니다.

## 파일

| 파일 | 내용 |
|---|---|
| `truss.h` | 구조체(`Node`, `Member`, `Truss`, `Result`)와 함수 선언 |
| `io.cpp` | `readTruss`, `writeResult` |
| `main.cpp` | `truss.txt` 읽기 → `solve` → `result.txt` 쓰기 |
| `solve.cpp` | `Result solve(const Truss&)` — 해석 본체 (작성 중, 아직 저장소에 없음) |
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

`readTruss`는 다음 입력을 거부하고, 어느 절점·부재에서 무엇이 틀렸는지 알려 줍니다.

- 값이 모자라거나 숫자가 아닌 칸, `j`·`m`보다 많은 데이터
- `fix_x`, `fix_y`가 0·1이 아닌 경우
- 범위 밖 절점 번호, 양 끝이 같은 부재, 같은 두 절점을 잇는 중복 부재, 길이 0인 부재
- `E`, `A`, `I`, `rho`, `sigma_t`, `sigma_c`, `K` 중 0 이하이거나 유한하지 않은 값

메모장이 붙이는 UTF-8 BOM과 Windows 줄바꿈(CRLF)은 그대로 읽힙니다.

## 출력 형식 (`RESULT 1`)

```
RESULT 1
status j m                           ← status: 0 Ok, 1 Mechanism, 2 InvalidInput
ux uy rx ry                          ← 절점 j줄
N stress cap dcr mode                ← 부재 m줄 (mode: 0 인장, 1 좌굴, 2 압괴)
loadFactor Pfail critical mass eta residual
```

실수는 유효숫자 17자리로 써서 double이 손실 없이 넘어갑니다. `solve`가 `ux`, `uy`, `rx`, `ry`를 서로 다른 크기로 채우면 `writeResult`는 `status 2`인 결과를 대신 쓰고 오류를 알립니다.

## 검증 — `truss.txt`의 정답

`solve`를 완성하면 아래와 같아야 합니다. 마지막 한두 자리는 연산 순서에 따라 달라질 수 있으니 상대오차 1e-12 이내면 맞습니다. `residual`은 1e-12 N 이하면 정상입니다.

```
RESULT 1
0 3 3
0 0 0 5
0.1710002499823563 0 0 5
0.08550012499117815 -0.37460523693241532 0 0
3.5355339059327369 1.2825018748676722 82.702426605751299 0.042750062495589075 0
-6.1237243569579443 -2.2213584080731512 5.8233639166534443 1.0515785110811193 1
-6.1237243569579443 -2.2213584080731512 5.8233639166534443 1.0515785110811193 1
0.95095134548908589 9.5095134548908593 1 1.3632149627567633 6.9757989126382771 0
```

손계산과의 대응: 밑변 축력 $\frac{P}{2}\cot\theta = 5/\sqrt{2}$, 빗변 축력 $-\frac{P}{2\sin\theta}$, 밑각 $\theta = \arctan\sqrt{2}$.
