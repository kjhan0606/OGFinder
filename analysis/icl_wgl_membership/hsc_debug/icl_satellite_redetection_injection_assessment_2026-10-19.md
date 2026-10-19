# ICL–위성 재검출 주입/회수 감사 — 2026-10-19

## 범위

기존 63개 공동 지원 cluster/filter 쌍에서 ICL 확산 성분과 위성 외곽광을
메모리 영상에 주입하고, 각 조건(baseline, ICL only, satellite only, both)마다
SEP source detection과 Kron 2.5 source mask를 다시 계산했다. ICL은 raw
200–300 kpc와 600–900 kpc 평면 보정 프로필을, 위성은 50–75 kpc target-minus-
azimuthal-control 프로필을 측정했다. 원본 FITS는 수정하지 않았고 HSC 질의나
다운로드는 없었다.

실행은 Slurm job 428211, `a10`, GPU 1개 요청, 16 GiB 메모리로 3분 56초에 끝났고
최대 RSS는 825,000 KiB였다. 기존 고정 마스크 결과 폴더는 보존했으며 새 산출물은
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/joint_icl_satellite_redetection_2026-10-19_01/`
에 저장했다. 독립 산출물 감사는 20개 frame JSON의 SHA-256, 시나리오 수와
미승격 science 상태를 확인했다. 감사 상태는 `verified_partial`이고 요약 해시는
`6e4e84ad6ec030dd7ccd00cd06da1ab9a68cb1abfaa9393fd8dcc8e3c9141c72`다.

## 처리 결과

| 결과 | i | i2 | 합계 |
|---|---:|---:|---:|
| 성공 프레임 | 5 | 15 | 20 |
| 선택 프레임 | 17 | 46 | 63 |

43개 실패 중 41개는 재검출 후 중심에서 유일한 BCG 후보를 찾지 못해 차단됐다.
나머지 2개는 재검출 마스크 적용 후 50–75 kpc 위성 프로필 지원이 남지 않았다.
실패 프레임을 결과에 대입하거나 성공 표본을 전체 표본으로 일반화하지 않았다.

20개 성공 프레임을 원래 고정 마스크 결과와 같은 cluster/filter ID로 짝지어 비교한
기술 통계는 다음과 같다. 각 값은 필터별 성공 하위표본 중앙값이며, 차이는 각
프레임에서 `redetection − fixed mask`를 먼저 계산한 뒤의 중앙값이다.

| 응답량 | 필터 | n | 고정 마스크 중앙값 | 재검출 중앙값 | 짝지은 차이 중앙값 |
|---|---|---:|---:|---:|---:|
| ICL←확산 주입 | i | 5 | 0.2270 | 0.2264 | −0.0025 |
| ICL←확산 주입 | i2 | 15 | 0.2261 | 0.2211 | −0.0008 |
| ICL←위성 외곽광 | i | 5 | 0.4421 | 0.3758 | −0.0842 |
| ICL←위성 외곽광 | i2 | 15 | 0.4757 | 0.2732 | −0.0683 |
| 위성←확산광 | i | 5 | 0.0016 | −0.0371 | −0.0490 |
| 위성←확산광 | i2 | 15 | 0.0008 | 0.0374 | 0.0330 |
| 위성←위성 주입 | i | 5 | 0.2523 | 0.4490 | 0.1619 |
| 위성←위성 주입 | i2 | 15 | 0.3950 | 0.3913 | −0.0188 |

이는 재검출이 일부 교차 응답을 바꾸는지 보여주는 성공 프레임 한정 비교다.
짝지은 차이 범위는 위성←확산 응답에서 i `−0.252…+4.326`, i2
`−2.495…+0.320`, 위성←위성 응답에서 i2 `−7.330…+2.411`로 넓다.
고정 분산에서 마스크/분할 변화만 다시 계산한 것이므로 응답의 표본 중앙값을
보정계수나 실제 ICL–위성 공분산으로 사용할 수 없다.

## 판단과 다음 관문

현재 주입 측정기는 63쌍 모집단에서 유일 BCG 재검출과 위성 annulus 지원을 동시에
유지하지 못했다. 이는 주입행렬을 관측 ICL–위성 상관에 적용할 수 없다는 추가
근거다. `science_measurement=false`, `observational_validation_passed=false`,
`flux_covariance_available=false`를 유지한다.

이번 계산은 주입원 photon noise, 관측 PSF convolution, upstream HSC sky 처리,
상관된 coadd 분산을 포함하지 않는다. 다음 단계는 BCG 검출을 억지로 통과시키는
것이 아니라, 실패 원인을 중심 검출 후보 수·입력 BCG/중심 오프셋·Kron 마스크
coverage로 분리해 기록하고, 사전에 정한 마스크 정책으로 전체 표본의 주입 회수율을
다시 검증하는 것이다. 그 뒤 실제 PSF/잡음 및 source 재검출을 포함한 공분산 검증이
필요하다.

재현 도구는 [`joint_icl_satellite_injection.py`](joint_icl_satellite_injection.py),
[`run_joint_icl_satellite_batch.py`](run_joint_icl_satellite_batch.py),
[`audit_joint_icl_satellite_redetection.py`](audit_joint_icl_satellite_redetection.py)다.
Slurm 설정은 [`run_joint_icl_satellite_redetection.sbatch`](run_joint_icl_satellite_redetection.sbatch)에 있다.
