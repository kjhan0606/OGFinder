# HSC Y3 tomographic weak-lensing consistency — 2026-10-20

## 실행과 입력

Astrafek `audit_wl_tomographic_consistency.py`는 이미 승인된 800개 CAMIRA 렌즈의
3,069,544 lens–source 행과 로컬 HSC Y3 SACC `n(z)`를 읽어 HSC source-bin 1–4를
따로 쌓았다. source bin 0의 814,898행은 공개 선택 규약대로 제외했고, bin 1–4에
속한 2,254,646 후보 중 품질·B-mode·거리·반경 선택을 통과한 쌍은 2,153,004개다.
각 bin 통과 수는 579,865, 688,189, 569,863, 315,087이다.

Slurm job **428830**은 12분 41초 만에 성공했고 최고 RSS는 1.83 GiB였다. 5도
equal-area sky cell 52개를 표적/원천 좌표에 공통 적용하고, 매 삭제에서 네 source
bin 모두의 shape·selection 보정을 다시 계산했다. 저장된 공분산은 56×56이며 대칭이다.
최소 고유값은 최대 고유값의 −1.6×10⁻¹⁷로, 부동소수점 수준의 작은 음수다.
처음 제출한 job 428828은 입력을 읽기 전에 `ogfkit` 경로 오류로 2초 만에 실패했다.
저장소 루트 경로를 명시하고 검사를 다시 통과한 뒤 428830을 실행했다.

기계 판독 결과는
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_tomographic_consistency_2026-10-20_01/summary.json`
에 있다. SHA-256은
`39b967e91d27ba3ebe1f05f24779245afff2fa197077469c8227d4d5efcfa643`이다.
그림 `tomographic_profiles.png` SHA-256은
`1bf67603e4631ed22f0d53199214277e110b828cfd3d3efff126fdf594ebd147`이다.
그림 생성기도 Astrafek의 `plot_wl_tomographic_consistency.py`로 저장했고,
합성 입력 smoke test를 포함한 관련 테스트 4개가 통과했다. 분석기·SACC·jackknife
코드 해시는 JSON 입력 항목에 기록되어 있다.

## 프로파일과 source-bin 일관성

표의 `DeltaSigma`와 jackknife 오차 단위는 10¹² M☉ Mpc⁻²다. 네 source-bin 추정치는
같은 lens population의 질량 프로파일을 보정해 비교하는 진단이다.

| R (Mpc) | bin 1 | bin 2 | bin 3 | bin 4 | 4-bin 일관성 χ²/dof | 7개 구간 Holm p |
|---|---:|---:|---:|---:|---:|---:|
| 0.1–0.2 | −3.9 ± 36.5 | 169.5 ± 19.6 | 179.0 ± 20.0 | 203.2 ± 28.7 | 25.80/3 | 6.29×10⁻⁵ |
| 0.2–0.4 | −10.3 ± 19.1 | 95.5 ± 9.7 | 146.5 ± 10.0 | 132.1 ± 13.3 | 45.15/3 | 6.02×10⁻⁹ |
| 0.4–0.7 | 35.2 ± 13.3 | 62.6 ± 5.7 | 73.5 ± 6.6 | 90.7 ± 5.7 | 21.22/3 | 4.74×10⁻⁴ |
| 0.7–1.0 | 57.8 ± 13.7 | 49.3 ± 5.4 | 58.4 ± 4.0 | 54.2 ± 6.3 | 2.18/3 | 1.000 |
| 1.0–1.5 | 38.3 ± 9.2 | 38.2 ± 4.0 | 39.8 ± 2.2 | 36.5 ± 4.4 | 0.51/3 | 1.000 |
| 1.5–2.0 | 29.6 ± 8.1 | 25.5 ± 3.2 | 25.6 ± 2.2 | 34.6 ± 3.7 | 6.26/3 | 0.398 |
| 2.0–3.0 | 9.7 ± 5.7 | 14.6 ± 1.7 | 15.1 ± 2.0 | 17.1 ± 1.9 | 3.03/3 | 1.000 |

공통 `DeltaSigma` GLS 검정은 네 bin과 공분산을 사용하며 일곱 반경 구간의 p-value에
Holm 보정을 적용했다. 0.1–0.7 Mpc의 세 구간에서 source-bin 간 차이가 보정 뒤에도
남는다. 나머지 바깥 구간은 이 진단에서 일관성과 맞는다. bin 1은 source redshift가
렌즈와 가까워 `n(z)` 평균 적색이동이 0.452이고, 렌즈 z=0.4에서 평균
`<D_ls/D_s>`가 약 0.129로 bin 2–4의 0.407, 0.541, 0.604보다 작다. 따라서 bin 1은
앞배경 희석과 `n(z)` 교정에 특히 민감하지만, bin 2–4만 남긴 사후 민감도에서도
0.2–0.7 Mpc 구간의 차이가 남아 이를 bin 1만의 문제로 돌릴 수 없다(raw
`p=0.00382, 0.00202`; 7개 반경에 대한 Holm `p=0.02294, 0.01414`). 이 사후
분석은 source bin 1의 불일치를 본 뒤 선택한 집합이며, 그 선택 효과를 포함해
보정한 검정이 아니므로 확증적 p-value로 해석하지 않는다.

전체 28개 cross 성분 공분산 검정은 χ²=74.75, rank=28, p=3.93×10⁻⁶이었다. 각 source
bin 내부의 7개 반경 cross 검정은 개별적으로 유의하지 않았으며, 네 bin에 대한 Holm
보정 뒤에도 유의한 bin은 없었다. 이 차이는 source-bin 사이의 공간 공분산을 포함한
joint test가 반경별 집계와 다르게 작동했음을 보여 준다. 더 정밀한 관측 영역 기반
공분산과 random 차감 전에는 체계오차 검토 신호로만 다룬다.

## 해석과 남은 관문

이 진단은 four-bin 보정 후에도 안쪽 `DeltaSigma`가 source redshift bin에 따라
일치하지 않음을 보인다. HSC SACC `n(z)`의 표본/선택 매칭, source-bin별 photo-z
dilution과 shear-selection 보정, cluster 주변 source incompleteness를 우선 점검해야
한다. 이 결과를 질량 추정이나 렌즈 검출 유의도로 보고하지 않는다.

이는 cluster-only stack이다. random-lens profile 차감, SACC posterior/systematic
불확실성 전파, reduced-shear 보정, 연결된 HSC 관측영역 기준 공분산은 포함하지 않았다.
raw B-mode 값 0을 제외했고 기존 좌표-영역 검증을 따른다. 최종 weak-lensing 결과로
승격하기 전에 matched random source 행을 확보하고 동일 jackknife 표본에서
`cluster − random` 차이의 공분산을 계산해야 한다.

모든 결과에는 `science_measurement=false`가 유지된다.

## 후속: cross-shear null 반경·source-bin 분해

전체 28성분 cross-shear null(χ²=74.75, dof=28, p=3.93×10⁻⁶)을 저장 covariance와 profile에서 독립 재구성하고, 반경 shell별 4-source-bin joint test와 source-bin별 7-shell joint test를 계산했다. 가장 낮은 shell별 원시 p는 1.5–2.0 Mpc의 0.0145였으나, 7개 반경 가족 Holm 보정 후 0.102였다. source-bin별 최소 p는 bin 2의 0.0775이며 4-bin 보정 후 0.310이다. 28개 개별 cross component의 최대 표준화 잔차는 2.11이고 28개 Holm 최소 p는 0.966이다.

따라서 이 진단에서는 단일 반경·source bin·profile 성분이 전체 결합 null 불일치를 보정 후 주도한다고 확인되지 않았다. 전체 joint signal은 유지되지만 원인은 특정하지 못한다. random-lens 차감, source 선택/shape 보정과 관측 footprint 기반 공분산을 확보하기 전에는 물리 신호로 해석하지 않는다. 이미 저장된 tomographic summary만 사용했고 추가 CAS 조회·자료 반출은 없었다.

Astrafek 재현기 [`audit_wl_cross_shear_null_decomposition.py`](audit_wl_cross_shear_null_decomposition.py), 테스트 3개 통과. 출력은 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_cross_shear_null_decomposition_2026-10-20_01/cross_shear_null_decomposition.json`; SHA-256 `94eb784e4a48e989c5ae0590a6049b2be289ae7544db9f2490c95a1ffc975e90`. 입력 summary SHA-256 `39b967e91d27ba3ebe1f05f24779245afff2fa197077469c8227d4d5efcfa643`, 분석기 SHA-256 `a8e424813afe90782a8224cc53baabd7993f6ae51d2d33e715762b62abc7dc46`다. `science_measurement=false` 유지.
