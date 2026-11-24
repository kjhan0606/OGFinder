# S16A photo-z 반경 지지와 Y1 보정 입력 점검

기존 S16A 파일럿에서 Ephor_AB best와 risk가 둘 다 유한한 소스를 집계했다. Table 4와 수치 조건을 통과한 169,518개 소스–중심 쌍에 대해 동일한 66개 중심과 8개 proper Mpc 반경 구간을 사용했다. 유한한 입력이 있는 것은 161,036쌍과 고유 소스 109,307개였다. 이 단계는 배경 소스 조건을 정하기 전에 photo-z 입력의 공간적 결측을 파악하기 위한 것이다.

| Proper Mpc | Table 4 및 수치 조건 이후 공통 그룹 | Ephor_AB 입력이 있는 공통 그룹 |
|---|---:|---:|
| 0–0.1 | 2 | 2 |
| 0.1–0.2 | 3 | 3 |
| 0.2–0.4 | 3 | 2 |
| 0.4–0.7 | 3 | 3 |
| 0.7–1 | 4 | 4 |
| 1–1.5 | 4 | 4 |
| 1.5–2 | 4 | 4 |
| 2–3 | 4 | 4 |

공통 그룹의 기준은 표적과 짝지은 대조 중심 5곳 모두에서 각 반경 구간에 소스가 한 개 이상 있는 것이다. 0.1–3 Mpc의 모든 구간을 만족하는 그룹은 S23B_3012, S23B_3620, S23B_9131 세 개에서 S23B_3012와 S23B_9131 두 개로 줄었다. S23B_3620의 대조 중심 `RC-cbf934142d614f978176`에서 0.2–0.4 Mpc의 65쌍 모두에 photo-z 입력이 없다. 같은 중심 0.4–0.7 Mpc에는 161쌍 중 6쌍만 있어 기존 shape weight의 3.98%만 남는다. S23B_3620의 다른 대조 중심 `RC-e4b4b443d83177619483`에서도 중앙 구간의 10쌍 모두에 입력이 없다. 기존 공통 지지 요건은 중앙 구간에서 이미 S23B_3620을 제외했다.

가중치 합 기준 완비율은 표적 95.11%, 대조 94.99%다. 두 평균은 가까워도 해당 반경과 중심에서의 결측을 가린다. 회색 셀은 기존 선택을 통과한 소스가 없는 구간이다. 보라색 최솟값은 선택을 통과한 소스는 있으나 유한한 photo-z 입력이 없는 구간이다. 이 파일럿의 공통 그룹 수로 전체 은하단 표본의 유지율을 추정할 수는 없다.

![Ephor_AB 입력의 반경 및 중심별 shape weight 완비율](/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s16a_photoz_radial_support_2026-11-18_01/photoz_radial_coverage.png)

실행 스크립트는 수신 계획, 파일 해시, 원본 pair 원장과 선택 원장을 확인하고 기존 object QA의 반경별 쌍 수를 다시 대조한다. 직교 단위 벡터의 chord 거리와 독립 Simpson 적분을 사용한 검산도 528개 구간의 쌍 수와 shape weight 합에서 일치했다. 결과와 파일 해시는 [summary.json](/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s16a_photoz_radial_support_2026-11-18_01/summary.json)에 있다. 재현할 때 `--output`은 아직 존재하지 않는 작업 디렉터리를 지정한다.

```bash
cd /home/kjhan/BACKUP/Astrafek/OGFinder
/home/kjhan/miniconda3/bin/python3 -B analysis/icl_wgl_membership/hsc_debug/audit_s16a_photoz_radial_support.py --output /path/to/new-output
```

## Y1 보정에 필요한 입력

| 항목 | 현재 확인 범위와 계산 전 필요한 입력 |
|---|---|
| 사용 가능 영역 | Shape 표의 멤버십으로 Y1 선별을 확인할 수 있다. 반경별 사용 면적에는 Wide2의 Y1 깊이·visit·coadd PSF 검사를 통과한 픽셀, 연결 영역과 Sirius bright-object mask의 위치·반경이 필요하다. 공개 random 표의 좌표와 visit/flag 열은 이 검증 전에는 면적의 대용 지표로만 쓸 수 있다. |
| Point-z 효율 보정 | S16A에 재가중한 S17A Wide-COSMOS의 참조 redshift, SOM weight, shape weight와 동일 방법의 Ephor_AB 추정량이 필요하다. 선택을 확정한 뒤 실자료와 보정 표본에 동일하게 적용한다. Proper Mpc 기준의 역임계표면밀도 단위와 lens/source 가중치도 일치시킨다. |
| Y1 shape 보정 | 현재 shape weight와 (e_{\rm rms}, m, c)가 있다. 반경별 최종 pair weight로 ensemble responsivity와 평균 (m)을 정하고 (e,c)를 같은 각도로 회전한다. Additive (c)는 shear 단위이므로 (2\mathcal R)로 나누지 않는다. 최종 선택에 대한 (R_2=0.3) 경계의 가중 밀도와 redshift 의존 responsivity를 따로 점검한다. |
| 공분산 | 같은 소스가 여러 대조 중심에 쓰인다. 사용 영역과 표본 선택을 확정한 뒤 공간 블록과 source 재사용을 함께 반영한다. |

Y1 공개 논문의 전체 면적 136.9 deg²는 개별 annulus의 면적을 제공하지 않는다. 시리우스 마스크는 Tycho-2/NOMAD를 사용한 조기 마스크다. 다른 공개 릴리스의 bright-star mask를 쓸 경우 그 결과를 Y1 면적과 등치할 근거가 필요하다. [Y1 Shape Catalog 논문](https://arxiv.org/html/1705.06745v3), [HSC S16A 공개 안내](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s16a-shape-catalog-pdr2/), [S16A random 공개 schema](https://hsc-release.mtk.nao.ac.jp/schema/schemata/pdr1.s16a_wide.random.json).

Point-z 효율 보정에서는 보정 표본의 참조 효율과 추정 효율의 비를 비교한다. (\beta=\Sigma_{\rm crit}^{-1})로 두면 보정 평균의 가중치는 SOM weight와 (w_{\rm shape}\beta_{\rm est}^2)의 곱이다. 후자는 실제 pair 선택·스택 정의와 맞춰야 한다. Y1의 resolution 선택 바이어스는 최종 가중 분포의 경계 밀도에 비례하며, 논문의 전체 표본 보정값을 이 파일럿에 고정값으로 적용할 수 없다. [Nakajima et al. §7.1](https://academic.oup.com/mnras/article/420/4/3240/973217), [Y1 shear calibration §5.3 및 §5.6](https://arxiv.org/html/1710.00885).

## 다음 자료 확인 범위

공개 안내와 dsigma 구현 예제가 공통으로 가리키는 최소 point-z 대응 파일 세 개의 주소·크기·ETag를 HSC의 인증된 `HEAD` 응답으로 확인했다. 앞 두 파일은 이름에 `pdf`가 있으나 예제에서 읽는 것은 point catalog다. S17A ID로 참조 redshift와 SOM/shape weight를 결합하는 데 사용한다. 현재 S16A CAS의 best/risk 버전과의 대응 및 동일 선택을 재현하는 데 필요한 추가 열은 FITS 수신 후 검증한다. 두 S17A 측광 파일의 필요 여부는 열 대조에 따라 결정한다. [HSC S17A Wide-COSMOS](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s17a-wide-cosmos/), [dsigma의 S16A 구현 예제](https://dsigma.readthedocs.io/en/v0.6.0/application_hsc.html).

| 파일 | 바이트 |
|---|---:|
| `pdf-s17a_wide-9812.cat.fits` | 5,008,320 |
| `pdf-s17a_wide-9813.cat.fits` | 63,452,160 |
| `Afterburner_reweighted_COSMOS_photoz_FDFC.fits` | 838,080 |
| 합계 및 수신 상한 | 69,298,560 (66.09 MiB) |

세 파일 본문 수신에는 사용자의 별도 승인이 필요하다. 준비 단계에서 받은 FITS 본문은 0 bytes다. 파일명과 주소, 버전, 상한은 [수신 계획](/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s16a_y1_ephor_ab_cosmos_plan_2026-11-18_01/plan.json)에 있다. 계획 SHA-256은 `73a28a83e9a2579b290ab68347702e0f2d28adfa9de162bf37e97df360727975`다. S16A Wide2의 Y1 영역과 Sirius 자산은 별도 검증 항목이다. `science_measurement=false`다. [HSC S16A 데이터 경로 안내](https://hsc-release.mtk.nao.ac.jp/doc/index.php/tools-2/).
