# ICL–위성 측정의 공통 픽셀 잡음 공분산

**평가일:** 2026-10-18 (Asia/Seoul)<br>
**상태:** 대각 분산 기반 공분산 진단; `science_measurement=false`,
`coadd_covariance_available=false`

## 질문과 계산

실제 ICL 200–300 kpc 고리와 위성 50–75 kpc 고리에 공통 유효 픽셀이 많더라도, 위성
측정은 방위각 대조 고리 평균을 빼므로 두 추정기의 잡음 공분산은 단순 중첩 면적으로
판정할 수 없다. 저장된 target과 선택된 control 위치를 복원해 cluster-filter별 위성
target-minus-control 평균의 픽셀 가중치를 만들었다. target 가중치는 양수, 각 수락된
control 가중치는 음수다. 위성은 cluster 안에서 동일 가중 평균했다.

공통 coverage 표본은 63개 cluster-filter 쌍(17 i, 46 i2)과 162개 위성 프로필이다.
50–75 kpc 반경별 coverage 60% 조건만 적용해 저장된 per-bin 측정 선택을 재현했다.
모든 target 픽셀 수, control 개수 및 control 픽셀 수가 저장 표와 일치했다.

각 프레임의 HSC 분산 배열을 독립 픽셀 분산으로 보고 signed satellite weight와 ICL target
200–300 kpc 균일 가중치의 분산 가중 내적을 계산했다. 이는 공통 픽셀과 대조 고리를
포함하는 **명목 대각 잡음 공분산**이다.

## 결과

| 필터 | cluster-filter 쌍 | 음의 공분산 쌍 | 명목 픽셀 잡음 상관 중앙값 [16, 84 백분위] |
|---|---:|---:|---:|
| i | 17 | 11 | −0.0112 [−0.0209, +0.0080] |
| i2 | 46 | 25 | −0.0023 [−0.0243, +0.0058] |
| 전체 | 63 | 36 | 표본 합성 상관은 산출하지 않음 |

프레임별 명목 픽셀 잡음 상관 범위는 −0.0602에서 +0.0363이었다. 이전 실제 픽셀 감사에서
공통 유효 support가 확인됐지만 signed control 가중치를 반영하면 분산맵 기반 상관은 작고
대체로 약한 음수였다. 대조 위치가 공통 ICL annulus 픽셀의 일부를 빼는 효과와 일치하는
결과다.

## 해석과 한계

이 계산은 위성 측정의 target과 control을 모두 포함해 단순 중첩 진단보다 완전하지만,
HSC coadd의 픽셀 간 상관잡음은 무시한다. 또한 ICL의 random-control 항은 포함하지 않아
ICL **target-side** 추정량과의 공분산만 계산했다. sky 전달, PSF 날개, 보정오차, 잔류
천체광, 주입 신호의 천체물리 변동도 포함하지 않는다. 따라서 작은 명목 상관은 실제 측정
공분산이 0이라는 증거가 아니며, 관측된 ICL–위성 연관성을 보정하지 않는다.

다음 단계는 빈 하늘 위치를 이용해 coadd 상관잡음을 실측하고, ICL target/control과 위성
target/control을 함께 블록 재표집한 공분산을 추정하는 것이다. 그 뒤 PSF 및 sky 처리,
검출·마스크 재생성을 포함한 독립 주입·회수에서 이 결과를 검증해야 한다.

## 재현과 계보

```bash
python3 analysis/icl_wgl_membership/hsc_debug/audit_satellite_icl_shared_pixel_covariance.py \
  --output-dir /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_icl_shared_pixel_covariance_2026-10-18_01
```

프레임별 결과는 `pair_pixel_covariance.tsv`, 집계와 SHA-256은 같은 디렉터리의
`summary.json`에 있다. CPU 실행은 lageunha에서 완료했다.
