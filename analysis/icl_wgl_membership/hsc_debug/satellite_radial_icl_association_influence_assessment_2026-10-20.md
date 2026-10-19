# 위성 반경별 ICL 연관성 영향도 — 2026-10-20

## 결과

앞선 14개 common-support 분석군(i/i2, 위성 7개 반경)에서 은하단 한 개 또는 연결된
5도 target/control 블록 하나를 차례로 제외했다. 총 1,820개 은하단 삭제와 455개
연결 블록 삭제를 기록했다. 아래 수치는 400–500 kpc 고리 차감 뒤의 Spearman ρ와,
블록 하나를 뺀 뒤의 ρ 범위다. 삭제 범위는 신뢰구간이 아니다.

| 필터 | 위성 반경 (kpc) | 전체표본 ρ | 한 블록 제외 범위 | 원시 ICL ρ 부호가 바뀐 삭제 수 (블록/은하단) |
|---|---:|---:|---:|---:|
| i | 5–10 | 0.212 | 0.165–0.259 | 0/0 |
| i | 10–15 | 0.202 | 0.156–0.243 | 5/4 |
| i | 15–20 | 0.144 | 0.101–0.172 | 6/12 |
| i | 20–30 | 0.145 | 0.103–0.186 | 6/9 |
| i | 30–40 | 0.125 | 0.067–0.166 | 3/0 |
| i | 40–50 | 0.015 | −0.031–0.063 | 4/1 |
| i | 50–75 | 0.075 | 0.032–0.124 | 11/36 |
| i2 | 5–10 | 0.011 | −0.011–0.035 | 0/0 |
| i2 | 10–15 | 0.026 | 0.002–0.052 | 0/0 |
| i2 | 15–20 | 0.041 | 0.007–0.077 | 0/0 |
| i2 | 20–30 | 0.026 | −0.011–0.060 | 0/0 |
| i2 | 30–40 | 0.074 | 0.056–0.115 | 0/0 |
| i2 | 40–50 | 0.067 | 0.035–0.114 | 0/0 |
| i2 | 50–75 | 0.086 | 0.056–0.130 | 0/0 |

i-band raw ICL 상관은 10–50 kpc 구간에서 일부 삭제에 따라 부호가 바뀌지만,
외곽 차감 상관은 5–40 kpc에서 삭제 한 건씩에 대해서는 양수를 유지했다. 이는 선택한
표본 안에서 순위의 영향도가 제한적이라는 뜻일 뿐, 표본·광도·마스크·배경 오차가
작다는 뜻이 아니다. 앞 분석의 공간 블록 bootstrap 95% 구간은 모든 raw 상관 및
paired raw-to-corrected 변화에서 0을 포함했다. 이어서 7개 반경을 공동 재표집한
max-T family-wise 95% 구간도 14개 상관과 14개 paired 변화 모두 0을 포함했다.
따라서 반경별 ICL–위성 연관성을 검출했다고 판단하지 않는다. 상세 구간과 산출물 해시는
[family-wise 평가](satellite_radial_icl_association_familywise_assessment_2026-10-20.md)에
기록했다.

## 방법과 제한

Astrafek은 앞 단계의 고정 77개 i·183개 i2 cluster/filter 공통표본과 같은 ICL raw 및
외곽 차감값을 사용했다. 각 반경에서 satellite mean, 위성 수, redshift·richness, ICL
안·바깥 면적 지지율을 보존해 Spearman과 조정 상관을 재계산했다. 블록 삭제는 기존과
같은 연결 5도 target/control 규칙을 썼다. 출력은 영향도 진단이며 p값이나 확률적
신뢰구간을 추정하지 않는다.

원본 ICL 400–500 kpc 고리는 검증된 blank sky가 아니며, 위성 구성은 같은 cluster 안에서도
반경에 따라 달라진다. PSF 날개, 상류 sky 전달, source completeness, 소광 및 전체 coadd
픽셀 공분산은 계속 미해결이다. `science_measurement=false`다.

## 재현과 계보

분석기 [`audit_satellite_radial_icl_association_influence.py`](audit_satellite_radial_icl_association_influence.py),
시험 [`test_audit_satellite_radial_icl_association_influence.py`](test_audit_satellite_radial_icl_association_influence.py)
2개가 통과했다. 산출물은
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_radial_icl_association_influence_2026-10-20_01/`이다.

요약 JSON SHA-256은 `d899ba9bb42f71eda3e77693144fbeef3e826e6b5c44872f28d8d9b69388251b`,
블록 삭제 표 SHA-256은 `183db214d26203beb09ca8b21e55940b5aa9a401b430b793b2bf2e8eecea34f8`,
은하단 삭제 표 SHA-256은 `e4f2e94cf22542f5657f19419b4f453b2049ae2b942fdb2adccfc43a06c835a7`,
분석기 SHA-256은 `374d33b833d6699a7f31bc3848b64ee7c7825bab09cb693db11df7d50e9116f2`다.
분석은 기존 자료만 사용했고 HSC에 새 질의를 하지 않았다.
