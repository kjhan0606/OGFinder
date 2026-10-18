# Joint ICL–위성 주입 응답행렬 안정성 감사

**판정:** 현재 고정 마스크 주입 결과로 관측된 ICL–위성 상관을 역행렬 보정하지 않는다. 응답행렬에 음의 행렬식과 큰 조건수가 있고, 주입 표본도 기존 연관 표본 일부만 덮는다. `science_measurement=false`를 유지한다.

## 감사 범위

63개 공통 지원 클러스터-필터 쌍(i 17, i2 46)의 프레임별 단일 주입 결과로 아래 응답행렬을 구성했다.

\[
H = \begin{pmatrix}
\partial I_{200-300}/\partial A_{\rm diffuse} & \partial I_{200-300}/\partial A_{\rm satellite} \\
\partial S_{50-75}/\partial A_{\rm diffuse} & \partial S_{50-75}/\partial A_{\rm satellite}
\end{pmatrix}
\]

각 프레임 결과 JSON의 SHA-256과 프레임/필터 ID, 주입 코드 해시, 미승격 science gate를 확인했다. 결과의 2-노름 조건수와 행렬식을 계산했으며, 관측량을 역행렬로 보정하지 않았다.

## 결과

| 필터 | 주입 쌍 | 음의 행렬식 | 조건수 >10 | 조건수 >100 | 조건수 >1000 | 조건수 중앙값 / 90백분위 / 최댓값 |
|---|---:|---:|---:|---:|---:|---:|
| i | 17 | 7 | 6 | 0 | 0 | 5.09 / 25.41 / 84.87 |
| i2 | 46 | 8 | 14 | 3 | 1 | 7.38 / 53.03 / 5,282.92 |

ICL 행에서 위성 주입 응답 대 확산 주입 응답의 비율 중앙값은 i 1.99, i2 2.59였다. 이는 지정한 단위 주입 모양에서 위성 외곽광이 ICL 측정량에 크게 섞일 수 있음을 나타낸다. 음의 행렬식은 일부 프레임에서 응답이 입력 성분의 방향을 뒤집는다는 뜻이므로 단순한 양의 성분 기여로 해석하기 어렵다. 조건수가 큰 경우에는 작은 측정오차·모델오차가 역산된 성분을 크게 바꾼다. 표의 조건수 10/100/1000 구분은 규모를 보기 위한 기술통계이며 사전 정의한 통과 기준은 아니다.

기존 CAMIRA 가중 ICL–위성 표본은 74쌍(i 21, i2 53)이고, 주입 행렬이 있는 것은 52쌍(i 14, i2 38)뿐이다. 두 표본은 같은 크기와 구성이 아니므로 63개 주입값을 74개 연관 표본 전체의 보정으로 일반화할 수 없다.

## 해석 및 다음 조치

현재 결과로 단순한 2×2 역행렬을 적용하면 수치 불안정성과 잘못 지정된 신호 모양을 실제 ICL–위성 상관으로 오인할 위험이 있다. 따라서 앞서 계산한 가중/비가중 상관은 탐색적 관측 통계로만 남기고, 주입행렬을 이용해 수정하지 않는다.

과학적 분리를 진행하기 전에 같은 cluster/filter/coverage 표본에서 위성 수와 반경 분포를 실제 관측 분포에 맞춘 신호 주입, 관측 PSF convolution, sky 처리 전체 경로, source 재검출 및 마스크 갱신, 영상 잡음/공분산을 포함한 end-to-end 전달 검증이 필요하다. 그 검증 후에도 행렬이 음의 행렬식 또는 큰 조건수를 보이면 단일 선형 분해 대신 더 식별 가능한 관측량이나 사전 정의된 계층 모델을 선택해야 한다.

## 재현

```bash
python3 analysis/icl_wgl_membership/hsc_debug/audit_joint_transfer_matrix.py \\
  --batch-summary /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/joint_icl_satellite_batch_2026-10-18_01/results/summary.json \\
  --spatial-audit /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/joint_icl_satellite_batch_2026-10-18_01/results/spatial_audit.json \\
  --association-pairs /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/weighted_satellite_icl_association_2026-10-18_01/weighted_satellite_icl_pairs.tsv \\
  --output /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/joint_transfer_matrix_audit_2026-10-18_01/summary.json
```

JSON 산출물은 입력 해시와 63개 행렬별 행렬식·조건수를 보존한다. 이 계산은 작은 집계 JSON/TSV만 사용했으며 HSC 자료를 새로 조회하거나 다운로드하지 않았다.
