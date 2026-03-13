# AI-Assisted Source Merge Decision Plan

## 1. Problem Statement

### 1.1 Background
SEP/SExtractor의 deblending 알고리즘은 경험적 파라미터(`deblend_nthresh`, `deblend_mincont`)에 의존하며, 다음 두 가지 오류를 피할 수 없다:

| 오류 유형 | 설명 | 예시 |
|---|---|---|
| **Over-deblending** | 하나의 천체가 여러 소스로 분리됨 | 나선 팔, HII 영역, 조석 꼬리가 별도 소스로 검출 |
| **Under-deblending** | 물리적으로 다른 천체가 하나로 합쳐짐 | 시선 방향 겹침(chance alignment), 상호작용 은하쌍 |

### 1.2 목표
겹쳐 보이는 소스 쌍/그룹에 대해 **"merge 해야 하는가?"**를 AI 모델로 자동 판단한다.

### 1.3 출력
- 각 후보 그룹에 대해: `merge` / `separate` / `uncertain` 분류
- 신뢰도 점수 (0.0 ~ 1.0)
- DS9 패널에서 사용자 확인/거부 워크플로우

---

## 2. Merge Candidate 정의

### 2.1 후보 선별 기준
추출된 소스 중 다음 조건을 만족하는 **쌍(pair)** 또는 **그룹(group)**을 merge 후보로 선정:

```
Candidate Pair (i, j):
  d_ij < α × (ISO_RADIUS_i + ISO_RADIUS_j)    # 거리 조건
  where d_ij = sqrt((x_i - x_j)² + (y_i - y_j)²)
  α = 1.5 (tunable overlap factor)
```

그룹 확장: 쌍을 그래프로 구성하고, connected components를 그룹으로 묶는다.

### 2.2 후보 통계 예상
- 일반적인 deep field (HST/JWST): 전체 소스의 5~15%가 overlap 후보
- Crowded field (은하단 중심부): 20~40%까지 증가
- 각 그룹 크기: 대부분 2~3개, 드물게 5개 이상

---

## 3. Feature Engineering

### 3.1 Image Features (CNN 입력)

각 후보 그룹에 대해 **cutout 이미지**를 생성:

| Feature | 설명 | 크기 |
|---|---|---|
| **Raw cutout** | 그룹 bounding box + padding (2×ISO_RADIUS) | N×N pixels |
| **Background-subtracted cutout** | SEP background 제거 후 | N×N pixels |
| **Segmentation map** | SEP segmap에서 해당 소스들의 footprint | N×N pixels |
| **Residual map** | 원본 - 각 소스의 모델(elliptical Gaussian) 합 | N×N pixels |

Cutout 정규화:
- Asinh stretch: `I_norm = arcsinh(I / σ_bkg) / arcsinh(I_peak / σ_bkg)`
- 크기 표준화: 64×64 pixels로 리샘플링 (bilinear interpolation)
- 멀티채널: [raw, segmap, residual] → 3채널 입력

### 3.2 Catalog Features (MLP 입력)

각 후보 쌍 (i, j)에 대해 추출하는 수치 특성:

**공간적 관계 (Spatial)**
| Feature | 수식 | 의미 |
|---|---|---|
| `sep_ratio` | d_ij / (R_i + R_j) | 정규화된 분리도 |
| `overlap_frac` | A_overlap / min(A_i, A_j) | isophote 겹침 비율 |
| `angle_diff` | |θ_i - θ_j| mod 180° | 위치각 차이 |
| `pa_alignment` | cos(2 × (PA_ij - θ_i)) | 연결 방향과 장축 정렬도 |

**측광적 관계 (Photometric)**
| Feature | 수식 | 의미 |
|---|---|---|
| `mag_diff` | |m_i - m_j| | 등급 차이 |
| `flux_ratio` | F_min / F_max | 플럭스 비율 |
| `color_diff` | (m_i - m_j)_band1 - (m_i - m_j)_band2 | 색깔 차이 (다중 밴드 시) |
| `sb_bridge` | μ_bridge / μ_bkg | 소스 사이 surface brightness |

**형태학적 관계 (Morphological)**
| Feature | 수식 | 의미 |
|---|---|---|
| `ellip_i`, `ellip_j` | 1 - b/a | 각 소스 타원율 |
| `size_ratio` | R_min / R_max | 크기 비율 |
| `class_star_i`, `class_star_j` | CLASS_STAR | 점원/확산원 분류 |
| `kron_ratio` | R_kron / R_iso | Kron-to-iso 비율 (확산도) |

**Pixel-level (Segmap 기반)**
| Feature | 수식 | 의미 |
|---|---|---|
| `saddle_depth` | (I_saddle - I_bkg) / (I_peak - I_bkg) | 안장점 깊이 비율 |
| `npix_ratio` | N_overlap / N_total | 겹침 픽셀 비율 |
| `gradient_at_boundary` | |∇I| at segmap boundary | 경계 기울기 |

### 3.3 Context Features (Optional)

| Feature | 의미 |
|---|---|
| `local_density` | 반경 R 내 이웃 소스 수 |
| `psf_fwhm` | 해당 위치의 PSF FWHM |
| `snr_i`, `snr_j` | 각 소스의 Signal-to-Noise Ratio |

---

## 4. Model Architecture

### 4.1 Option A: Hybrid CNN + MLP (Recommended)

```
┌─────────────┐    ┌──────────────┐
│ Cutout Image │    │ Catalog      │
│ 64×64×3      │    │ Features     │
│              │    │ (~20 dims)   │
└──────┬───────┘    └──────┬───────┘
       │                   │
  ┌────▼────┐         ┌───▼───┐
  │  CNN    │         │  MLP  │
  │ Encoder │         │ 64→32 │
  │ →128dim │         │       │
  └────┬────┘         └───┬───┘
       │                   │
       └────────┬──────────┘
                │
          ┌─────▼─────┐
          │ Concat     │
          │ 128+32=160 │
          └─────┬──────┘
                │
          ┌─────▼─────┐
          │ FC 160→64  │
          │ FC 64→3    │  ← merge / separate / uncertain
          │ Softmax    │
          └────────────┘
```

**CNN Encoder:**
```
Conv2d(3, 32, 3, padding=1) → BN → ReLU → MaxPool(2)    # 64→32
Conv2d(32, 64, 3, padding=1) → BN → ReLU → MaxPool(2)   # 32→16
Conv2d(64, 128, 3, padding=1) → BN → ReLU → MaxPool(2)  # 16→8
Conv2d(128, 128, 3, padding=1) → BN → ReLU → AdaptiveAvgPool(1)  # 8→1
Flatten → 128
```

**MLP Branch:**
```
Linear(N_feat, 64) → BN → ReLU → Dropout(0.3)
Linear(64, 32) → BN → ReLU
```

**파라미터 수:** ~500K (경량 모델, CPU 추론 가능)

### 4.2 Option B: Pure MLP (Lightweight)

Cutout 없이 catalog features만 사용. 빠르지만 정확도 낮음.

```
Input(~20) → 128 → 64 → 32 → 3 (softmax)
```

파라미터 수: ~15K. 추론 속도: <1ms per pair.

### 4.3 Option C: Vision Transformer (High Accuracy)

ViT-Tiny를 cutout에 적용. 높은 정확도, 그러나 큰 모델.

```
Patch embedding (8×8 patches → 64 tokens) → 4 Transformer layers → CLS token → MLP head
```

파라미터 수: ~5M. GPU 권장.

### 4.4 권장 전략

1. **Phase 1**: Option B (MLP only)로 빠르게 프로토타입
2. **Phase 2**: Option A (Hybrid)로 정확도 향상
3. **Phase 3**: 필요 시 Option C로 업그레이드

---

## 5. Training Data

### 5.1 Data Sources

#### 5.1.1 Simulated Data (Primary)

**GalSim으로 합성 이미지 생성:**

| 시나리오 | Label | 생성 방법 |
|---|---|---|
| 단일 은하 (over-deblend) | `merge` | Sérsic 프로필 1개 → SEP로 추출 시 2개 이상 검출되는 경우 |
| 은하 + substructure | `merge` | 은하 + 밝은 knot/HII 영역 추가 → 별도 검출 |
| 물리적 은하쌍 | `separate` | 서로 다른 Sérsic 프로필 2개, 다른 색깔/크기 |
| 시선 방향 겹침 | `separate` | 서로 다른 z에서 우연히 겹침 (크기/SB 차이 큼) |
| 상호작용 은하 | `uncertain` | 조석 브릿지로 연결된 은하쌍 |

**시뮬레이션 파라미터 범위:**
```python
sersic_n   = [0.5, 1, 2, 4]          # Sérsic index
half_light = [2, 5, 10, 20] pixels    # Half-light radius
mag        = [20, 22, 24, 26, 28]     # Magnitude
ellip      = [0.0, 0.3, 0.5, 0.7]    # Ellipticity
separation = [0.5, 1.0, 1.5, 2.0] × R_sum  # 분리 거리 / 반경 합
snr        = [5, 10, 20, 50, 100]     # Signal-to-noise ratio
psf_fwhm   = [1.5, 3.0, 5.0] pixels  # PSF FWHM
```

**목표 데이터셋 크기:**
- Training: 50,000 pairs (simulated)
- Validation: 10,000 pairs (simulated)
- Test: 5,000 pairs (simulated) + 1,000 pairs (real, human-labeled)

#### 5.1.2 Real Data (Fine-tuning & Evaluation)

| Dataset | Field | Depth | Resolution | 용도 |
|---|---|---|---|---|
| HUDF (HST ACS/WFC3) | XDF | ~30 mag | 0.03"/pix | Primary test set |
| CANDELS (HST) | GOODS-S/N | ~28 mag | 0.06"/pix | Training augmentation |
| HSC-SSP (Subaru) | Wide/Deep/UD | ~26-28 mag | 0.17"/pix | Ground-based test |
| JWST JADES/CEERS | Deep fields | ~30 mag | 0.03"/pix | Next-gen test |

**Human labeling:**
- Galaxy Zoo 스타일 웹 인터페이스 또는 DS9 내 라벨링 도구
- 전문가 3인 이상의 majority vote
- 라벨: `merge` / `separate` / `uncertain`
- 목표: 최소 1,000 쌍 (real data)

#### 5.1.3 Data Augmentation

| Augmentation | 설명 |
|---|---|
| Rotation | 0°, 90°, 180°, 270° |
| Flip | Horizontal, Vertical |
| Noise injection | Gaussian noise (σ = 0.5~2.0 × σ_bkg) |
| PSF convolution | 다양한 FWHM으로 재 convolution |
| Brightness scaling | ×0.5 ~ ×2.0 |
| Crop jitter | Center ± 5 pixels |

---

## 6. Training Pipeline

### 6.1 Data Preparation

```
                    ┌──────────────┐
                    │ FITS Image   │
                    └──────┬───────┘
                           │
                    ┌──────▼───────┐
                    │ SEP Extract  │
                    │ + Segmap     │
                    └──────┬───────┘
                           │
                    ┌──────▼───────┐
                    │ Candidate    │
                    │ Pair Finder  │  ← α × (R_i + R_j) 거리 조건
                    └──────┬───────┘
                           │
              ┌────────────┼────────────┐
              │            │            │
       ┌──────▼──────┐ ┌──▼───┐ ┌──────▼──────┐
       │ Cutout Gen  │ │ Feat │ │ Label       │
       │ 64×64×3     │ │ Ext  │ │ (sim/human) │
       └──────┬──────┘ └──┬───┘ └──────┬──────┘
              │            │            │
              └────────────┼────────────┘
                           │
                    ┌──────▼───────┐
                    │  HDF5 / NPZ  │
                    │  Dataset     │
                    └──────────────┘
```

### 6.2 Training Configuration

```yaml
model: hybrid_cnn_mlp
optimizer: AdamW
learning_rate: 1e-3
weight_decay: 1e-4
scheduler: CosineAnnealingLR
epochs: 100
batch_size: 128
early_stopping_patience: 15

loss: CrossEntropyLoss
class_weights:              # 클래스 불균형 보정
  merge: 1.0
  separate: 1.0
  uncertain: 2.0            # uncertain 클래스에 높은 가중치

evaluation_metrics:
  - accuracy
  - precision (per class)
  - recall (per class)
  - F1-score
  - confusion_matrix
  - ROC-AUC (one-vs-rest)
```

### 6.3 Cross-validation

- 5-fold cross-validation on simulated data
- Hold-out real data for final evaluation
- Stratified split (merge/separate/uncertain 비율 유지)

---

## 7. Inference & Integration

### 7.1 Model Deployment

```
ds9_sextract (C binary)
    │
    ├── Extract sources (existing)
    ├── Find candidate pairs (NEW)
    ├── Generate cutouts + features (NEW)
    └── Run AI inference (NEW)
         │
         ├── Option A: Embedded (ONNX Runtime in C)
         │     - 모델을 ONNX 포맷으로 export
         │     - ds9_sextract에 ONNX Runtime 링크
         │     - 추가 의존성: libonnxruntime
         │
         ├── Option B: Python subprocess
         │     - ds9_merge_predict.py 별도 스크립트
         │     - exec로 호출, JSON 입출력
         │     - 추가 의존성: torch, numpy
         │
         └── Option C: Shared library (.so/.dylib/.dll)
               - PyTorch → TorchScript → C++ load
               - ds9_sextract에서 dlopen
```

**권장: Option B** (초기 프로토타입, 유연성 높음)
**장기: Option A** (배포 용이, 의존성 최소)

### 7.2 출력 포맷

ds9_sextract의 기존 TSV 출력에 AI merge suggestion을 추가:

```
# 기존 출력
NUMBER  X_IMAGE  Y_IMAGE  ...  MAG_AUTO  ...

# 추가 출력 (별도 섹션 또는 별도 파일)
# MERGE_SUGGESTIONS
GROUP_ID  MEMBERS       DECISION    CONFIDENCE  REASON
1         5641,8480     merge       0.94        substructure_of_single_galaxy
2         1234,5678     separate    0.87        distinct_colors_and_sizes
3         9012,3456     uncertain   0.52        possible_interaction
```

### 7.3 DS9 UI Integration

```
+-----------------------------------------------------------+
| Source Extractor  [Extract] [Settings...] [Mark All] ...   |
+-----------------------------------------------------------+
| [AI Merge]  ← 새 버튼                                      |
+-----------------------------------------------------------+
| Merge Suggestions (3 groups found):                        |
|                                                            |
| ○ Group 1: #5641 + #8480                                  |
|   Decision: MERGE (94%)                                    |
|   Reason: Substructure of single galaxy                    |
|   [Accept] [Reject] [Show]                                 |
|                                                            |
| ○ Group 2: #1234 + #5678                                  |
|   Decision: SEPARATE (87%)                                 |
|   [Accept] [Reject] [Show]                                 |
|                                                            |
| ○ Group 3: #9012 + #3456                                  |
|   Decision: UNCERTAIN (52%)                                |
|   [Merge] [Keep Separate] [Show]                           |
+-----------------------------------------------------------+
```

**워크플로우:**
1. Extract → 소스 검출
2. **AI Merge** 버튼 클릭 → 후보 쌍 탐색 + AI 추론
3. Suggestion 목록 표시 (신뢰도 순 정렬)
4. **Show**: 해당 소스 쌍을 이미지에서 하이라이트 (파란 박스)
5. **Accept/Merge**: merge 실행 (기존 `CatalogPanelMergeSources` 활용)
6. **Reject/Keep Separate**: skip, 다음 후보로
7. 모든 결정 후 카탈로그 갱신

### 7.4 Batch Mode (Non-interactive)

```bash
ds9_sextract --ai-merge --merge-threshold 0.8 image.fits > catalog.tsv
```

- `--ai-merge`: AI merge 자동 실행
- `--merge-threshold 0.8`: confidence ≥ 0.8인 쌍만 자동 merge
- threshold 미만은 그대로 유지 (사용자가 DS9에서 수동 확인)

---

## 8. Evaluation Metrics

### 8.1 Classification Metrics

| Metric | 목표 | 설명 |
|---|---|---|
| Accuracy | ≥ 90% | 전체 정확도 |
| Precision (merge) | ≥ 85% | merge로 판단 시 실제 merge 비율 |
| Recall (merge) | ≥ 90% | 실제 merge 중 검출 비율 |
| Precision (separate) | ≥ 90% | separate 판단 정확도 |
| F1 (macro) | ≥ 88% | 클래스 균형 F1 |

### 8.2 Astronomical Metrics

| Metric | 설명 |
|---|---|
| Photometric accuracy | merge 후 MAG_AUTO vs 실제 등급 차이 |
| Astrometric accuracy | merge 후 위치 vs 실제 중심 위치 오차 |
| Completeness | 전체 실제 merge 대상 중 AI가 찾은 비율 |
| Contamination | AI가 merge한 것 중 잘못된 merge 비율 |
| Number count impact | merge 전/후 source count 차이의 물리적 타당성 |

### 8.3 Benchmark Comparison

| Method | 비교 대상 |
|---|---|
| Default SEP deblending | 기존 알고리즘 baseline |
| SExtractor v2 | 고전적 SExtractor 결과 |
| ProFound (R) | 대안적 소스 추출기 |
| Human expert | 전문가 라벨 (gold standard) |

---

## 9. Active Learning (Iterative Improvement)

### 9.1 개념

초기 모델의 `uncertain` 판단 사례를 사용자에게 라벨링 요청 → 재학습 → 점진적 정확도 향상.

```
┌──────────┐     ┌──────────┐     ┌──────────┐
│ AI Model │────▶│ Uncertain│────▶│ User     │
│ v1       │     │ Cases    │     │ Labeling │
└──────────┘     └──────────┘     └────┬─────┘
                                       │
┌──────────┐     ┌──────────┐          │
│ AI Model │◀────│ Retrain  │◀─────────┘
│ v2       │     │ Pipeline │
└──────────┘     └──────────┘
```

### 9.2 DS9 내 라벨링 도구

- AI Merge 결과에서 `uncertain` 항목에 **Merge** / **Separate** 버튼 제공
- 사용자 결정을 `~/.ds9/merge_labels.tsv`에 누적 저장:
  ```
  IMAGE_FILE  SRC1  SRC2  LABEL  TIMESTAMP
  xdf.fits    5641  8480  merge  2026-03-13T14:30:00
  ```
- 누적 1000건 이상 → 재학습 트리거

---

## 10. Implementation Phases

### Phase 1: Data & Baseline (2~3주)

| Task | 내용 | 산출물 |
|---|---|---|
| 1.1 | GalSim 시뮬레이션 파이프라인 구축 | `simulate_pairs.py` |
| 1.2 | Candidate pair finder 구현 | `find_candidates.py` |
| 1.3 | Feature extractor 구현 | `extract_features.py` |
| 1.4 | Cutout generator 구현 | `generate_cutouts.py` |
| 1.5 | MLP baseline 학습 (Option B) | `model_mlp.pth` |
| 1.6 | 시뮬레이션 데이터 평가 | Accuracy, F1 리포트 |

### Phase 2: CNN Model & Real Data (3~4주)

| Task | 내용 | 산출물 |
|---|---|---|
| 2.1 | Hybrid CNN+MLP 모델 구현 | `model_hybrid.py` |
| 2.2 | HUDF 실제 데이터 라벨링 (1000쌍) | `labels_hudf.tsv` |
| 2.3 | Hybrid 모델 학습 + 평가 | `model_hybrid.pth` |
| 2.4 | Sim → Real transfer 성능 분석 | Fine-tuning 보고서 |
| 2.5 | ONNX export | `merge_model.onnx` |

### Phase 3: DS9 Integration (2~3주)

| Task | 내용 | 산출물 |
|---|---|---|
| 3.1 | `ds9_merge_predict.py` 추론 스크립트 | Python 스크립트 |
| 3.2 | `CatalogPanelAIMerge` Tcl 프로시저 | layout.tcl 수정 |
| 3.3 | AI Merge suggestion UI | DS9 패널 |
| 3.4 | Accept/Reject 워크플로우 | 사용자 상호작용 |
| 3.5 | Batch mode (`--ai-merge`) | CLI 지원 |

### Phase 4: Refinement & Deployment (2~3주)

| Task | 내용 | 산출물 |
|---|---|---|
| 4.1 | Active learning 파이프라인 | 라벨링 도구 + 재학습 |
| 4.2 | ONNX Runtime C 통합 (Option A) | ds9_sextract 수정 |
| 4.3 | 크로스플랫폼 테스트 | Linux/macOS/Windows |
| 4.4 | 문서화 (매뉴얼 추가) | main.tex 업데이트 |
| 4.5 | 벤치마크 논문 초안 | benchmark_results/ |

---

## 11. Directory Structure (Proposed)

```
SAOImageDS9/
├── ds9/library/
│   └── layout.tcl                    # AI Merge UI 추가
├── ai_merge/
│   ├── README.md
│   ├── requirements.txt              # torch, numpy, sep, galsim, onnx
│   ├── data/
│   │   ├── simulate_pairs.py         # GalSim 시뮬레이션
│   │   ├── find_candidates.py        # 후보 쌍 탐색
│   │   ├── extract_features.py       # 수치 특성 추출
│   │   └── generate_cutouts.py       # 이미지 cutout 생성
│   ├── models/
│   │   ├── mlp.py                    # MLP baseline
│   │   ├── hybrid_cnn_mlp.py         # Hybrid CNN+MLP
│   │   └── vit.py                    # Vision Transformer (optional)
│   ├── training/
│   │   ├── train.py                  # 학습 스크립트
│   │   ├── evaluate.py               # 평가 스크립트
│   │   └── config.yaml               # 학습 설정
│   ├── inference/
│   │   ├── ds9_merge_predict.py      # DS9 연동 추론 스크립트
│   │   ├── export_onnx.py            # ONNX export
│   │   └── merge_model.onnx          # 배포용 모델
│   └── tests/
│       ├── test_candidates.py
│       ├── test_features.py
│       └── test_inference.py
├── bin/
│   ├── ds9_sextract                  # 기존 바이너리
│   └── ds9_merge_predict             # AI 추론 (Python wrapper)
└── docs/
    └── manual/main.tex               # AI Merge 챕터 추가
```

---

## 12. Risk & Mitigation

| Risk | Impact | Mitigation |
|---|---|---|
| 시뮬레이션-현실 갭 (sim-to-real gap) | 실 데이터에서 낮은 정확도 | Fine-tuning + domain adaptation |
| 학습 데이터 부족 (real labels) | 과적합 | Data augmentation + active learning |
| 추론 속도 느림 (GPU 없는 환경) | 사용자 경험 저하 | ONNX Runtime CPU 최적화, batch processing |
| 모델 크기 (배포) | 설치 복잡성 | Quantization (INT8), ~2MB 모델 |
| Edge case: 상호작용 은하 | 명확한 정답 없음 | `uncertain` 클래스 + 사용자 확인 |
| 다중 밴드 지원 | 복잡도 증가 | Phase 1은 단일 밴드, 추후 확장 |

---

## 13. References

- Bertin & Arnouts 1996, SExtractor (A&AS 117, 393)
- Barbary 2016, SEP: Source Extraction and Photometry (JOSS 1, 58)
- Rowe et al. 2015, GalSim (Astronomy & Computing 10, 121)
- Hausen & Robertson 2020, Morpheus: Deep Learning for Galaxy Morphology (ApJS 249, 25)
- Bretonnière et al. 2022, SourceXtractor++ (A&A 657, A90)
- Melchior et al. 2018, SCARLET: Source Separation (Astronomy & Computing 24, 129)
