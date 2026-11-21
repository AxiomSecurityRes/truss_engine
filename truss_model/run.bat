@echo off
chcp 65001 >nul
setlocal EnableExtensions EnableDelayedExpansion
cd /d "%~dp0"
set "ROOT=%~dp0"
title truss_model

rem ===== g++ 찾기: PATH -> Code::Blocks 에 딸린 MinGW -> 흔한 설치 위치 =====
set "GPP="
where g++ >nul 2>nul && set "GPP=g++"
if not defined GPP (
    for %%D in ("C:\Program Files\CodeBlocks\MinGW\bin" "C:\Program Files (x86)\CodeBlocks\MinGW\bin" "C:\MinGW\bin" "C:\msys64\ucrt64\bin" "C:\msys64\mingw64\bin") do (
        if not defined GPP if exist "%%~D\g++.exe" (
            set "GPP=%%~D\g++.exe"
            set "PATH=%%~D;!PATH!"
        )
    )
)

rem ===== 소스(.cpp/.h)가 exe 보다 새로우면 자동으로 다시 빌드 (예전 exe 로 새 config 를 읽는 오류 방지) =====
set "NEEDBUILD=0"
if not exist truss_model.exe set "NEEDBUILD=1"
if "%NEEDBUILD%"=="0" (
    powershell -NoProfile -Command "$e=(Get-Item 'truss_model.exe').LastWriteTime; if ((Get-ChildItem *.cpp,*.h | Where-Object { $_.LastWriteTime -gt $e }).Count -gt 0) { exit 1 } else { exit 0 }" >nul 2>nul
    if errorlevel 1 set "NEEDBUILD=1"
)
if "%NEEDBUILD%"=="1" (
    if exist truss_model.exe echo 소스가 바뀌어 다시 빌드합니다.
    call :build
)
if not exist truss_model.exe (
    echo.
    echo 프로그램을 만들지 못했습니다. 위 오류를 확인하세요.
    pause
    exit /b 1
)

:menu
cls
echo ================================================
echo   truss_model  -  트러스 구조 효율 / GA
echo ================================================
echo   결과는 results 폴더에 저장됩니다. 설정은 config.txt
echo.
echo   1. 전체 실행 (검증 - E - 스윕 - 하중 - GA - 민감도 - 재료 - 요약, 약 2분)
echo   2. 검증 (verify)
echo   3. E 회귀 (efit)
echo   4. 높이 스윕 (sweep)
echo   5. 대표 4종 비교 (compare)
echo   6. 하중 조건 비교 (loads)
echo   7. 유전 알고리즘 (ga, 약 40초)
echo   8. 민감도 (sens, GA 를 먼저 돌려야 GA 최적도 포함)
echo   9. 대표 형태 하나 만들기 (preset)
echo   M. 재료별 비교 (빨대, 강, 알루미늄, 목재 2종 / 약 1분)
echo   S. 논문용 요약 (summary.txt, 단위 포함. 1 또는 7, M 을 먼저 실행)
echo.
echo   C. config.txt 열기      F. results 폴더 열기
echo   R. 다시 빌드            0. 끝
echo.
set "SEL="
set /p "SEL=번호를 고르세요: "
if "%SEL%"=="1" call :run all & goto done
if "%SEL%"=="2" call :run verify & goto done
if "%SEL%"=="3" call :run efit & goto done
if "%SEL%"=="4" call :run sweep & goto done
if "%SEL%"=="5" call :run compare & goto done
if "%SEL%"=="6" call :run loads & goto done
if "%SEL%"=="7" call :run ga & goto done
if "%SEL%"=="8" call :run sens & goto done
if "%SEL%"=="9" goto preset
if /i "%SEL%"=="M" call :run materials & goto done
if /i "%SEL%"=="S" call :run report & (if exist results\summary.txt start "" notepad results\summary.txt) & goto done
if /i "%SEL%"=="C" start "" notepad config.txt & goto menu
if /i "%SEL%"=="F" (if not exist results mkdir results) & start "" explorer results & goto menu
if /i "%SEL%"=="R" call :build & goto done
if "%SEL%"=="0" exit /b 0
goto menu

:preset
set "KIND=warren"
set "PN="
set "PH="
set /p "KIND=형태 (warren / pratt / howe / k) [warren]: "
set /p "PN=패널 수 n, 짝수 [config 값]: "
set /p "PH=높이 H mm [config 값]: "
call :run preset !KIND! !PN! !PH!
goto done

:done
echo.
pause
goto menu

rem ===== 실행: results 폴더에서 돌린다 (설정과 입력 파일을 복사해 결과와 함께 남김) =====
:run
if not exist results mkdir results
copy /y config.txt results\ >nul
copy /y deflection.csv results\ >nul
copy /y truss.txt results\ >nul
pushd results
..\truss_model.exe %*
set "RC=!errorlevel!"
popd
echo.
if "!RC!"=="0" (echo 완료. 결과: !ROOT!results) else (echo 오류 코드 !RC!)
exit /b 0

rem ===== 빌드 =====
:build
if not defined GPP (
    echo g++ 를 찾지 못했습니다.
    echo Code::Blocks 를 MinGW 포함판으로 설치했거나 g++ 이 PATH 에 있어야 합니다.
    echo 아니면 Code::Blocks 에서 truss_model.cbp 를 열어 빌드하세요.
    exit /b 1
)
echo 빌드 중: "%GPP%"
"%GPP%" -std=c++14 -O2 -static *.cpp -o truss_model.exe
if errorlevel 1 (
    echo 빌드 실패
    exit /b 1
)
echo 빌드 완료: truss_model.exe
exit /b 0
