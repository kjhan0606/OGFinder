# S23B ICL–위성–weak-lensing 고정표본 교집합 점검 (2026-11-11)

## 범위와 결과

기존 S23B CAMIRA member-position/ICL profile과 고정 52-lens weak-lensing per-lens profile을 해시 검증해 결합 가능한 표본을 확인했다. 원시 source CSV와 영상 픽셀은 읽지 않았고, CAS/네트워크 요청 및 데이터 다운로드도 없었다. 공동 관계 적합은 하지 않았으며 `science_measurement=false`를 유지한다.

| ICL 필터 | 고정 WL 표본 | ICL–위성 공통 지원 렌즈 | 안쪽·바깥쪽 모두 ≥10% 공통 픽셀 | 모두 ≥20% | 중앙 안쪽/바깥쪽 공통 픽셀 비율 |
|---|---:|---:|---:|---:|---:|
| i | 52 | 25 | 18 | 5 | 21.3% / 15.4% |
| i2 | 52 | 27 | 22 | 7 | 23.5% / 18.3% |

모든 52개 weak-lensing 렌즈가 저장된 ICL–위성 profile과 연결되지만, 각 렌즈는 i 또는 i2 한쪽에만 나타난다. 양쪽 필터가 모두 있는 렌즈는 없다. ICL contrast의 표본 중앙값은 i에서 7.43, i2에서 2.95 nJy arcsec⁻²다. 각 필터의 CAMIRA 후보 위성 밀도비 `log10(inner/outer)` 중앙값은 0.629와 0.586이다. 이 값들은 필터별 서로 다른 은하단 표본에서 나온 기술 통계이므로 두 필터 사이의 물리적 차이로 해석하지 않는다.

## 연구상 의미

52개 교집합은 저차원의 탐색적 ICL–위성–lensing 관계를 검토할 수 있는 표본이다. 다만 필터별 표본은 25개와 27개로 작고, 두 ICL annulus 모두 20% 공통 픽셀을 요구하면 각각 5개와 7개만 남는다. 공통 픽셀 선택 민감도와 약중력렌즈 표본이 높은 원시 형상원 가용성 및 다소 낮은 CAMIRA 풍부도 쪽으로 선택된 점을 함께 고려해야 한다.

따라서 이 교집합으로는 확인적 검출이나 전체 S23B 집단으로의 일반화를 주장하지 않는다. 다음 관계 계산은 약한 차원의 탐색 분석과 공간 공분산·영상 지지 민감도에 한정해야 한다. 넓은 표본으로 일반화하려면 추가 HSC shape-row 반출 계획과 별도 사용자 승인이 필요하다.

## 재현 정보

- 실행기: [audit_s23b_icl_wl_fixed_support_overlap.py](audit_s23b_icl_wl_fixed_support_overlap.py)
- 승인 row-plan SHA-256: `7bfb998149f30a50b5ee7a537bb703ed302b7ea269a3296c35204a836e75aabc`
- 고정 52 렌즈 ID SHA-256: `18312f27139b3750abdbf71543f22d9b7b18b2d3b24cbd29ed686a377524bb7d`
- ICL–위성 pair table SHA-256: `98ed9e0022bfa68d5921ead62dce55be7df1cfa169e6a777d66d1e9760f43994`
- WL per-lens profile SHA-256: `5bc599d97123264058011313b0ad35c3def49b08bba98a016cc5edf6ff03906c`
- 결과 디렉터리: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_wl_fixed_support_overlap_2026-11-11_01`
- 결과 `summary.json` SHA-256: `8145f969b60fe93075d86bd952c9cd0d0e38f9be71be2ff521252048754cf69e`
