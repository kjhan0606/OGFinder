# ICL–위성 94쌍 지원표본 확장 민감도 (2026-11-05)

## 목적

기존 상관 진단 74쌍이 이미 보유한 전체 결합 지원표본 94쌍을 대표하는지 확인했다.
94쌍은 저장된 ICL 200–300 kpc 및 위성 50–75 kpc 자료의 지원 교집합이다(i 26쌍,
i2 68쌍). 이 계산은 보유한 요약 행과 좌표 manifest만 사용했다. 영상 픽셀 재독,
HSC 조회, 추가 자료 반출은 없었다.

## 방법

Astrafek의 [`analyze_icl_satellite_supported_94.py`](analyze_icl_satellite_supported_94.py)는
필터별로 raw target-minus-control ICL 200–300 kpc 표면밝기와 raw 위성 50–75 kpc
평균 대비의 부분 Pearson 상관을 계산했다. 공변량은 redshift와 CAMIRA richness다.
불확실도는 표적 또는 대조 endpoint 중 하나라도 같은 equal-area 5° 셀을 공유하는
쌍을 연결요소로 묶은 percentile block bootstrap 5,000회로 추정했다.

74쌍은 94쌍 안에 포함된 기존 pixel-extinction association 행으로 정했다. 두 표의
중첩 74행에서 ICL 값의 최대 절대 차이는 `1.36e-7 nJy arcsec⁻²`, 위성 값 차이는
0이었다. 표본 범위 변화의 paired interval은 94쌍의 연결블록 재표집을 74쌍과 94쌍에
같이 적용했다. 표본별 단독 bootstrap interval은 각 표본 안에서 블록을 다시 계산했다.

## 결과

| 필터 | 범위 | 쌍 | 연결 5° 블록 | 부분 r (95% block interval) |
|---|---:|---:|---:|---:|
| i | 기존 진단 | 21 | 13 | 0.120 `[-0.223, 0.515]` |
| i | 전체 지원 | 26 | 13 | 0.079 `[-0.183, 0.477]` |
| i2 | 기존 진단 | 53 | 28 | 0.234 `[0.012, 0.418]` |
| i2 | 전체 지원 | 68 | 33 | 0.293 `[0.039, 0.500]` |

동일한 94쌍 블록을 사용한 paired 범위 차이(94 minus 74)는 i에서 `−0.041`
(`[−0.245, 0.092]`), i2에서 `+0.059` (`[−0.188, 0.343]`)였다. 두 차이 구간 모두
0을 포함한다. i2의 점추정치는 추가 15쌍을 포함한 전체 지원 표본에서 커졌지만,
이 민감도만으로 표본 범위에 따른 효과 변화를 구별할 수 없다. i 추정치는 작고
불확실성이 크다.

## 해석과 다음 관문

확장해도 i2의 양의 탐색 상관은 관측되지만, 그 크기는 선택된 지원 표본과 공간 블록
불확실도에 민감하다. 지원행은 parent 800개에서 무작위 추출된 표본이 아니며, 이
분석은 기존 selection audit의 불균형을 보정하지 않는다. 보정된 extinction, coadd
공분산과 Poisson 입력, PSF wing, sky transfer, 위성 완전도 및 마스크 영향도 아직
과학 측정을 제한한다. 따라서 두 필터 모두 결과를 검출이나 인과 관계로 해석하지
않고 `science_measurement=false`를 유지한다.

현재 800개 target에서 충분한 확증 표본을 만들기보다, 다음 핵심 단계는 HSC frame
contributor count-only 관문을 승인 범위 안에서 완료하고 필요한 frame-level gain 및
exposure metadata의 별도 최소 반출 계획을 세우는 것이다. 기존 count-only 제출은
아직 승인 대기이며 본 분석은 그 대기와 무관하게 진행했다.

## 재현과 provenance

```bash
python3 analysis/icl_wgl_membership/hsc_debug/analyze_icl_satellite_supported_94.py \
  --output /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_supported_94_association_2026-11-05_01 \
  --replicas 5000
```

산출물은 `summary.json`과 `supported_pairs.tsv`이며, 요약 상태는
`complete_diagnostic`, `science_measurement=false`, `hsc_network_used=false`다.
요약 파일 SHA-256은 `291ffff417d4c1da33412ac1cf1aa3f02a3a3e5c2dc2c2677973fb2885e42ad6`,
쌍 표 SHA-256은 `79f7c4c6b18b29ad837e91da4a4359433b09a564cad62efaf60b4cce655685da`다.
재현기 SHA-256은 `dc8b14bcebd02566bc8375ce57575dce4986420ddbe4fd494131006d448e1e7e`다.

주요 입력 해시는 summary.json에 저장된다. 입력 경로는 기존 94쌍
`satellite_outer_icl_association_perbin_2026-10-18_03/matched_cluster_values.tsv`,
74쌍 `icl_pixel_extinction_2026-10-19_03/pixel_extinction_satellite_association.tsv`,
`sample_manifest.json`, 그리고 matched-control `download_manifest.json`이다.
