# HSC S23B 표본의 source-shape–PSF 정렬 진단

## 범위와 결론

승인된 HSC S19A shape-row export를 사용해 S23B 렌즈와 matched control 표본의 source ellipticity가 같은 위치의 PSF ellipticity와 선형 상관을 보이는지 점검했다. 새 CAS 질의, 자료 다운로드, 영상 픽셀 열람은 없었다. 계산은 `lageunha` CPU에서 수행했다.

프로필 분석과 같은 4/5 control-support 정책은 렌즈 세트 74개와 양의 source 지지를 가진 대조 중심 367개를 선택했다. HSC source bin 2–4, 원 catalog B-mode mask에서 통과한 source, `|e|<2`, `0≤sigma_e≤0.4`, `mag_A10≤25.5`, `R2≥0.3`, 양의 n(z) lensing efficiency, proper radius 0.1–3 Mpc를 적용했다. 공식 HSC 문서는 `e1/e2`를 하늘 좌표 ellipticity로 정의하고 PSF moments를 선택편향 보정에 사용한다고 설명한다. HSC calibration utility의 PSF ellipticity 정의도 검사 코드에 적용했다. [HSC S19A Shape Catalog 문서](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s19a-shape-catalog-pdr3/)와 [공식 문서에서 연결된 HSC Y3 calibration utility](https://github.com/mr-superonion/hsc-y3-calib/blob/main/utils.py)를 참조했다.

고유 source 단위 가중 회귀에서 렌즈와 대조의 PSF 정렬 계수 차이는 모두 약 1.17 공간 잭나이프 오차 이내였다. 2.5°와 5° 셀에서 계수 점추정치는 같았고 오차 크기도 비슷했다. 본 검사는 표본별 선형 정렬 null이므로 PSF 누출이 없음을 증명하지 않는다. `science_measurement=false`를 유지한다.

## 방법

각 렌즈 또는 대조 표본에서 `object_id`가 반복된 pair 행을 고유 source 한 개로 축약했다. 렌즈·대조 표본에 함께 나타난 source 3,824개는 역할 간 비교의 독립성을 유지하도록 양쪽 표본에서 제외했다. 한 `object_id`의 좌표, ellipticity, shape weight, `m/c`, 선택 변수, PSF moments, tomographic bin 값이 모든 pair 행에서 일치하는지도 검사했다.

각 source의 PSF ellipticity는 `e1_PSF=(PSF11−PSF22)/(PSF11+PSF22)` 및 `e2_PSF=2 PSF12/(PSF11+PSF22)`로 계산했다. raw `e1/e2`와 per-object `m/c` 보정을 적용한 `(e/(2R_z)−c)/(1+m)`를 각각 PSF ellipticity에 회귀했다. `R_z`는 선택된 pair 행에서 HSC shear weight `shape_weight×(D_l beta)^2`로 계산한 source-bin별 평균 responsivity이며, `R_z=1−sum(w rms_e^2)/sum(w)`로 정의했다. 모델은 HSC source bin별 절편과 두 PSF 성분을 포함했다. 두 회귀는 비교를 위해 같은 source-level 가중치 `shape_weight×[2R_z(1+m)]^2`를 사용했다. HSC stack의 radial-bin edge-selection 보정은 회귀식에 적용하지 않았다.

공간 오차는 source 좌표로 나눈 equal-area RA/sin(Dec) 셀을 하나씩 제외해 계산했다. 5° 셀은 현재 S23B profile의 주 공분산 분할과 맞췄고, 2.5° 셀은 분할 폭 민감도를 확인했다. 두 분할 모두 실제 HSC 연속 관측 영역을 검증한 것이 아니므로 p-value를 계산하지 않았다.

## 결과

선택 컷 뒤 렌즈 pair 행 174,064개는 고유 source 172,993개에 해당했다. 대조 pair 행 840,486개는 고유 source 460,289개에 해당했다. 교차 역할 source를 제외한 회귀 표본은 렌즈 169,169개와 대조 456,465개였다. 대조 pair에서 source가 여러 중심에 반복되는 빈도가 더 높아 고유 source 축약이 필수였다.

아래 표는 `e_corr = b_zbin + A e_PSF`의 계수 행렬 `A`를 보여준다. 행은 galaxy ellipticity 성분이고 열은 PSF ellipticity 성분이다. 오차는 같은 공간 셀 잭나이프에서 얻었다.

| 공간 셀 | 표본 | e1–PSF e1 | e1–PSF e2 | e2–PSF e1 | e2–PSF e2 |
|---|---|---:|---:|---:|---:|
| 5° | 렌즈 | −0.0225 ± 0.0249 | −0.0053 ± 0.0216 | 0.0014 ± 0.0200 | 0.0088 ± 0.0204 |
| 5° | 대조 | −0.0011 ± 0.0192 | 0.0064 ± 0.0128 | −0.0039 ± 0.0122 | −0.0133 ± 0.0114 |
| 5° | 렌즈−대조 | −0.0214 ± 0.0314 | −0.0117 ± 0.0261 | 0.0054 ± 0.0193 | 0.0221 ± 0.0199 |
| 2.5° | 렌즈−대조 | −0.0214 ± 0.0288 | −0.0117 ± 0.0242 | 0.0054 ± 0.0200 | 0.0221 ± 0.0190 |

렌즈−대조 계수의 최대 절대비는 2.5° 셀에서 잭나이프 오차의 약 1.16배였다. 5° 셀에서는 약 1.11배였다. 점추정치는 셀 폭에 거의 변하지 않았지만, 오차는 검증되지 않은 equal-area 분할에 조건부로 계산됐다.

## 해석의 한계

HSC PDR3 문서는 survey-level PSF model 정확도와 보정 뒤 B-mode·random-point null 검사를 설명한다. 본 감사는 현재 S23B 표본에서 galaxy shape와 catalog PSF moments 사이의 선형 관계만 측정했다. HSC star null test나 image-level PSF 재구성은 재현하지 않았다.

현재 결과는 표본에 남은 선형 PSF 정렬이 profile 오차보다 뚜렷하게 크다는 증거를 제공하지 않는다. 작은 계수만으로 PSF 영향을 배제할 수는 없다. PSF 오차와 공간 공분산, 대조 위치의 깊이·mask·sky 매칭, selection-response 불확실성, reduced shear, member dilution 및 magnification은 계속 독립 검증이 필요하다. 원래 profile과 마찬가지로 과학 측정 또는 질량 검출로 해석하지 않는다.

## 재현 정보

- 분석기 [audit_s23b_wl_psf_shape_alignment.py](audit_s23b_wl_psf_shape_alignment.py)
- aggregate 결과 [wl_s23b_psf_shape_alignment_2026-11-11.json](wl_s23b_psf_shape_alignment_2026-11-11.json)
- 입력 row-plan SHA-256 `7bfb998149f30a50b5ee7a537bb703ed302b7ea269a3296c35204a836e75aabc`
- row-status SHA-256 `fe2611d149e577c3ca3a6b44d0375848988eef63eb08f18702531925ccfc8bf4`
- HSC Y3 n(z) SACC SHA-256 `a28f9e2e088e92d96d2d87083a6eeeac4c0c84c6b48e033bd3d1dc1bf0d8958f`
- aggregate JSON SHA-256 `d735305cbe127bd83dc1ad50f120ce5634c0a9955eb436fe8a49c975e690e699`
- 분석기 SHA-256 `d468a9005489ae551a4e3010b4790b546586eff306a40f0ea4097c049a5d9761`
