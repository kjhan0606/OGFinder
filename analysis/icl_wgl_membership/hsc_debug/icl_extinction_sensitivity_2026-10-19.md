# ICL–위성 연관성의 소광 보정 민감도 (2026-10-19)

## 목적과 입력

기존 HSC-I 표면광 프로필에 은하수 소광 보정을 적용했을 때 200–300 kpc ICL 평균과 50–75 kpc 위성 외곽광 연관성이 얼마나 달라지는지 확인했다. HSC CAS를 호출하거나 영상을 새로 받지 않았다. 800개 대상의 SFD 중심값과 798개 대조장의 SFD 중심값을 따로 사용했고, HSC PDR3 `forced.a_i`는 대상 보정의 민감도 경로로만 사용했다.

분석 스크립트는 `analyze_extinction_sensitivity_icl_satellite.py`다. 완료된 PDR3 쿼리 결과, 입력 요약 해시, 기존 target/control 테이블 해시, 3초각 내 forced-source 최근접 대응, 원래 target-minus-control 값을 확인하고 실행한다. 산출물은 다음 경로에 저장됐다.

`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_extinction_sensitivity_2026-10-19_03/`

`annulus_extinction_corrected.tsv`에는 200–300 kpc 고정 절단 224쌍(i 71, i2 153)이 있고, `icl_satellite_extinction_corrected.tsv`에는 기존 위성 연관표와 정확히 겹치는 74쌍(i 21, i2 53)이 있다. 두 표의 기록 해시는 재검산했고 필요한 수치 열은 모두 유한값이다. 재현 실행은 lageunha에서 다음과 같다.

```bash
python3 analysis/icl_wgl_membership/hsc_debug/analyze_extinction_sensitivity_icl_satellite.py \
  --work /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_extinction_sensitivity_2026-10-19_03 \
  --replicas 5000 --seed 20261019
```

## 정량 결과

| 필터 | 고정 절단 쌍 | 원시 평균 ΔSB | SFD 중심 보정 평균 ΔSB | PDR3 `a_i` 대상 보정 평균 ΔSB |
|---|---:|---:|---:|---:|
| i | 71 | 3.226 | 3.380 | 3.379 |
| i2 | 153 | 3.787 | 4.001 | 4.002 |

단위는 nJy arcsec⁻²다. 공간 연결 5도 블록을 5,000회 재표집한 SFD 보정 평균의 95% 구간은 i `[2.196, 4.463]`, i2 `[2.824, 5.076]`이다. 보정 전 평균 대비 변화는 i `+0.154`, i2 `+0.214` nJy arcsec⁻²였다.

ICL 200–300 kpc와 위성 50–75 kpc 외곽광의 적색이동·richness 보정 부분 Pearson 상관은 다음과 같다. 구간은 같은 연결 5도 블록을 재표집한 95% 구간이다.

| 필터 | 쌍 수 | 원시 r (95% 구간) | SFD 보정 r (95% 구간) | PDR3 대상 보정 r (95% 구간) |
|---|---:|---:|---:|---:|
| i | 21 | 0.120 `[-0.238, 0.530]` | 0.112 `[-0.241, 0.522]` | 0.112 `[-0.241, 0.522]` |
| i2 | 53 | 0.234 `[0.021, 0.412]` | 0.239 `[0.023, 0.417]` | 0.239 `[0.023, 0.417]` |

가장 가까운 HSC PDR3 forced-source `a_i`와 필터별 SFD `A` 예측의 절대 차 중앙값은 i 표본 309개에서 0.000297 mag, i2 표본 491개에서 0.000580 mag였다. 95백분위는 각각 0.000635 mag와 0.001500 mag였고, 최대 차는 각각 0.001134 mag와 0.002730 mag였다.

## 해석 및 한계

소광을 보정해도 평균 프로필과 상관 추정은 거의 변하지 않는다. 따라서 이 보정은 현재 i2 양의 상관을 만들거나 없애는 주된 요인이 아닌 것으로 보인다. 하지만 i2의 표본은 53쌍이고 선택·절단 민감도 분석이며, 기존 familywise 검정에서 다중비교 보정 후 유의한 그룹은 없었다. 위 결과는 물리적 검출로 해석하지 않는다.

보정은 각 대상·대조장 중심의 소광값을 해당 프로필 전체에 적용했다. 픽셀별 SFD 변화, 광도 프로필 반경에 따른 소광 변화, 지도 계통오차는 전파하지 않았다. PDR3 `a_i`는 대상 보정의 민감도 비교이며 모든 반경에서 독립적으로 측정된 값이 아니다. 하늘 잔류량, PSF 외곽 날개, 마스크 누출, 주입·회수 편향도 해결되지 않았다. 따라서 산출물은 `science_measurement=false`인 탐색 결과다.

weak-lensing 쪽은 행 반출을 완료했지만 calibrated HSC Y3 `n(z)` SACC 입력과 B-mode 저장값 극성 확인이 남아 있어 stack을 확정하지 않았다.
