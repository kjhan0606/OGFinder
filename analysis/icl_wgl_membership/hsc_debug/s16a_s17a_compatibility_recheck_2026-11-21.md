# S16A–S17A 배경 선택과 렌즈별 적색편이 응답 재검증

## 결론

S17A Wide-COSMOS는 S16A weak-lensing shape 표본의 기본 색–등급 분포에 맞도록 SOM 가중치를 제공한다. 동일 이름의 Ephor_AB point-z와 유한한 risk를 요구하고 `z_best > z_lens + gap`을 양쪽에 적용할 수 있지만, 이 검증만으로 S16A와 S17A의 전체 선택 함수나 photo-z 버전의 동등성은 입증되지 않았다. 특히 COSMOS 결합표에는 로컬 Table 4 전체 선택과 accepted-pixel/별 마스크 및 photo-z 결측 선택을 재현할 입력이 충분하지 않다. risk 수치 절단은 적용하지 않았다.

[HSC S17A Wide-COSMOS 공식 문서](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s17a-wide-cosmos/)는 SOM과 shape 가중치를 함께 사용하고, 추가 적색편이 선택을 보정 표본에도 적용하도록 안내한다. [S16A shape 카탈로그 문서](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s16a-shape-catalog-pdr2/)는 Sirius mask와 별도의 selection bias 및 적색편이에 따른 responsivity 보정을 명시한다. 기본 색–등급 분포 일치가 본 연구의 조건부 선택까지 동일함을 보증하지는 않는다.

기존 COSMOS 결합은 기준 18,852개 중 18,768개를 연결했고, 양의 shape×SOM 가중치의 99.496%를 포함한다. 84개 미연결 객체는 원장에 남았다. 이는 전체 선택의 동등성을 확인한 수치가 아니다.

## 렌즈별 재계산

기존 Astrafek 실행기를 LagEunha CPU에서 기존 저장 입력만으로 재실행했다. 11개 렌즈, gap 0/0.1/0.2, COSMOS 전체 및 두 tract에 대한 99개 결과 행은 기존 산출물과 바이트 단위로 일치한다. 두 파일 모두 기존 manifest에 기재된 SHA-256 `7a15ce87b7f9e55fdd6a1a91794a7d11c4d5dbedaba6eaa0ad6ac5e6266ac898`와 일치한다. 응답과 그 역수의 곱은 최대 `1.11e-16`까지 1에서 벗어난다. 이는 같은 구현의 재현성과 산술 검사이며 독립적인 물리적 검증을 뜻하지 않는다.

gap 0.2와 전체 COSMOS 표본에서 탐색적 ΔΣ 곱셈계수는 다음과 같다. 관측 S16A와 COSMOS의 유지 가중치는 각각 같은 자료에서 gap 0과 비교한 값이다.

| 렌즈 | 렌즈 z | 계수 | S16A 유지율 | COSMOS 유지율 |
|---|---:|---:|---:|---:|
| S23B_3012 | 0.5420 | 1.02994 | 96.86% | 96.39% |
| S23B_9131 | 0.3920 | 1.01874 | 95.85% | 97.19% |
| S23B_9013 | 0.5413 | 1.02982 | 97.33% | 96.40% |
| S23B_8627 | 0.5252 | 1.02847 | 96.79% | 96.40% |
| S23B_6290 | 0.5186 | 1.02793 | 96.85% | 96.40% |
| S23B_10600 | 0.4140 | 1.02004 | 97.17% | 97.34% |
| S23B_9828 | 0.5360 | 1.02937 | 94.88% | 96.37% |
| S23B_3620 | 0.4020 | 1.01935 | 96.58% | 97.25% |
| S23B_8310 | 0.4720 | 1.02413 | 94.50% | 96.73% |
| S23B_8206 | 0.3454 | 1.01772 | 96.63% | 96.81% |
| S23B_9307 | 0.3720 | 1.01822 | 96.74% | 97.05% |

기존 반경별 민감도에서는 0.1–3 Mpc의 annulus별 재가중과 전체 반경 렌즈별 계수를 비교하면 스택의 유효 계수 차이가 −1.090%~+0.694%였다. 이는 고정된 진단 가중치와 선택된 표본에 대한 민감도이며, 측정된 전단이나 총 오차가 아니다. 별도 기존 거리 적분 검사에서 Astropy 직접 거리 대비 최대 차이는 0.0000863%였다.

기준 적색편이는 참값의 대용치이며, COSMOS 두 tract는 같은 하늘 영역으로 독립 우주분산 오차가 아니다. 소스 재사용, 은하단 구성원 오염, 선택 편향, 참조 적색편이 및 공분산이 남아 있다. 결과를 실제 전단이나 질량에 적용하지 않았다. `science_measurement=false`, `calibration_applied=false`.

## 재현과 원장

저장소 루트에서 `OPENBLAS_NUM_THREADS=1 PYTHONPATH=. python3 -B analysis/icl_wgl_membership/hsc_debug/audit_s16a_cosmos_lens_response.py --output NEW_DIRECTORY`를 실행한다. 입력과 산출물 해시는 생성되는 `summary.json`에 기록된다. 이번 재실행은 `/home/kjhan/BACKUP/Astrafek/s16a_s17a_recheck_20261121/`, 이전 결과는 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s16a_cosmos_lens_response_2026-11-19_01/`에 있다. 추가 HSC 조회나 데이터 다운로드는 없었다.
