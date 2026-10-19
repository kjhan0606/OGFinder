# 위성 외곽광–ICL 연관성의 외곽 기준 민감도 — 2026-10-20

## 목적과 산출물

기존 800개 HSC coadd 표본에서 위성의 50–75 kpc 외곽광과 ICL의 200–300 kpc
차감광 사이 연관성을 확인하고, 각 target/control 쌍의 400–500 kpc 차감광을
같은 ICL 측정에서 빼면 연관성이 달라지는지 비교했다. 원시 ICL과 외곽 기준을 뺀
ICL은 동일한 cluster/filter/control 조인 표본에서 계산했다. 필터는 분리했다.

재현 분석기는 [`analyze_satellite_outer_icl_baseline_association.py`](analyze_satellite_outer_icl_baseline_association.py),
회귀검사는 [`test_analyze_satellite_outer_icl_baseline_association.py`](test_analyze_satellite_outer_icl_baseline_association.py)다.
실행은 5,000회 연결 하늘 블록 부트스트랩을 썼고, i 77개, i2 183개의 cluster/filter
쌍을 포함했다. 이는 외곽 기준 분석의 5% 유효면적 지지와 독립적으로 선택된 위성
50–75 kpc 프로필이 모두 존재하는 표본이다.

## 결과

| 필터 | 쌍/연결 블록 | ICL 항 | Spearman ρ [95% 블록 구간] | z·richness 조정 r [95% 블록 구간] |
|---|---:|---:|---:|---:|
| i | 77 / 21 | 원시 200–300 kpc | 0.001 [−0.262, 0.242] | 0.055 [−0.221, 0.300] |
| i | 77 / 21 | 400–500 kpc 차감 후 | 0.075 [−0.145, 0.291] | 0.113 [−0.071, 0.263] |
| i2 | 183 / 44 | 원시 200–300 kpc | 0.051 [−0.079, 0.178] | 0.040 [−0.076, 0.141] |
| i2 | 183 / 44 | 400–500 kpc 차감 후 | 0.086 [−0.063, 0.246] | 0.085 [−0.049, 0.224] |

차감 전후 상관계수 변화(차감 후 minus 원시) 역시 0과 양립한다.
i Spearman 변화는 `0.074 [−0.138, 0.295]`, 조정 상관 변화는
`0.058 [−0.149, 0.284]`다. i2는 각각 `0.034 [−0.119, 0.208]`와
`0.045 [−0.100, 0.209]`다. 현 표본에서는 위성 외곽광–ICL 연관성의 증거가 없으며,
외곽 기준 차감이 그 관계를 유의하게 바꾼다는 증거도 없다.

위성 표본 수와 ICL 면적 지지가 상관을 설명하는지 확인하기 위해 redshift·richness에
`log(1+n_satellites)`, 200–300 kpc 내부 지지율, 400–500 kpc 외곽 지지율을 더한
부분상관도 계산했다.

| 필터 | ICL 항 | 확장 조정 r [95% 블록 구간] |
|---|---|---:|
| i | 원시 | 0.029 [−0.252, 0.252] |
| i | 400–500 kpc 차감 후 | 0.077 [−0.141, 0.250] |
| i | 차감 후−원시 변화 | 0.048 [−0.199, 0.333] |
| i2 | 원시 | 0.033 [−0.082, 0.142] |
| i2 | 400–500 kpc 차감 후 | 0.083 [−0.051, 0.237] |
| i2 | 차감 후−원시 변화 | 0.050 [−0.097, 0.221] |

이 추가 조정에서도 모든 구간은 0을 포함한다. 이 결과는 위성 수나 면적 지지만으로
관계가 가려졌다는 근거를 주지 않는다. 다만 위성 수와 지원율은 측정 완전도와 함께
변할 수 있어 이 선형 조정 자체를 선택편향 보정으로 해석하지 않는다.

단, 400–500 kpc 공통 유효면적의 중앙값은 i 17.2%, i2 17.3%이고, 20% 미만인 쌍이
각각 47/77, 111/183이다. 외곽 기준광은 i에서 40%, i2에서 56%의 쌍에서 음수였다.
따라서 외곽 기준을 빼면 평균 ICL이 i `3.90 → 2.59`, i2 `2.16 → 3.62`
nJy arcsec⁻²로 달라진다. 이 차감은 검증된 sky 보정이 아니며, 실제 확장 은하단
빛과 외곽 면적·마스크의 잔차를 함께 포함할 수 있다.

## 해석 한계

위성 광 프로필은 50–75 kpc 구간에서 독립적으로 선택됐고 CAMIRA weight나 위성 질량을
적용하지 않았다. 자동 중심·마스크, 확장 PSF 날개, HSC 상류 sky 전달, 소광,
passband response, 전체 coadd 픽셀 공분산은 검증되지 않았다. 블록 구간은 선택한
5도 target/control 연결 블록만 반영한다. 두 필터의 분석은 탐색적이며
`science_measurement=false`를 유지한다. 이 결과는 외곽 기준 민감도 평가이지,
확정 ICL 프로필이나 weak-lensing 질량 측정이 아니다.

## 재현과 계보

산출물 디렉터리:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_outer_icl_baseline_assoc_2026-10-20_05/`.
요약 JSON SHA-256은 `3583d146a659843eea0ec0f0984e630b496a8694d32c150e6c8208d8ed40396d`,
쌍별 표 SHA-256은 `899f209a129c951aacfc4f9878ef6c0706b917f6bbf0574acd5798ea985eac48`이며,
분석기 SHA-256은 `3c5cd2229486c69437c4accfc07117969f9457cad8e1c5463ed1affe2172e41b`이다.
입력 자료의 경로와 SHA-256도 `summary.json`에 저장했다.
명령:

```bash
python3 analysis/icl_wgl_membership/hsc_debug/analyze_satellite_outer_icl_baseline_association.py \
  --work /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_outer_icl_baseline_assoc_2026-10-20_05 \
  --replicas 5000
```

분석기는 기존 FITS·catalog·PSF 산출물만 읽고 외부 질의를 수행하지 않는다. 새 HSC 자료를
받지 않았으며 다운로드 승인은 사용하지 않았다.
