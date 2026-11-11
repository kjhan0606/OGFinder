# OGFinder / Astrafex Desktop 인수인계 (HANDOVER)

작성: 2026-10-08 (Asia/Seoul), 소유자 Juhan Kim <kjhan0606@gmail.com>.
추가: 2026-10-12 (Asia/Seoul). 6장. 1–5장은 2026-10-08 기록이다. 전체 서술은 ds10-web `HANDOVER.md` 11장이다.
추가: 2026-10-13 (Asia/Seoul). 7장. ICL–위성은하 외곽구조 장기 연구주제와 도구 설계.
추가: 2026-10-13 (Asia/Seoul). 8장. 분석 목적 기반 범용 파이프라인 요구사항.

**전체 인수인계 문서는 Astrafex Web 저장소의 [`HANDOVER.md`](https://github.com/kjhan0606/ds10-web/blob/master/HANDOVER.md)** (`kjhan0606/ds10-web`, private)이다.
프로젝트 목적, 이름(Astrafex), 아키텍처, 기능별 검증 결과, 작업 규칙, 남은 일은 모두 거기에 있다. 이 문서는 OGFinder에만 해당하는 부분만 다룬다.
두 저장소는 같은 상위 폴더에 나란히 clone한다 (`<WS>/OGFinder`, `<WS>/ds10-web`; 원래 box는 `<WS>=/workspace`).

## 1. 이 저장소는 무엇인가

* SAOImageDS9 8.7 소스 트리 + OGFinder 분석 기능 = **Astrafex Desktop** ("based on the OGFinder tool set"). 라이선스 GPL-3.0 (소유자의 잠정 결정, `README.md` License 절, `docs/licensing.md`, `docs/license_audit.md`). 번들된 third-party 구성요소는 각자 라이선스 유지.
* 저장소/폴더/라이선스 표기/ds9 쪽 GUI 문자열의 "OGFinder" 이름은 상표 검토가 끝날 때까지 그대로 둔다 (ds10-web `docs/naming.md`).
* 주요 디렉터리:

| 경로 | 내용 |
|---|---|
| `ds9/`, `tksao/`, `tcl8.6/`, `tk8.6/`, `tcllib/`, `ast/`, ... | ds9 8.7과 vendored Tcl/Tk/C 라이브러리 (빌드 대상) |
| `ds9/library/*.tcl` | ds9 GUI + OGFinder catalog panel, session recorder, layout (`docs/architecture.md`) |
| `plugins/<name>/` | 분석 plugin (manifest `plugin.json` + CLI `*.py` + Tcl GUI hook + `tests/`): extract, photometry, isophote, multifit, psfex, daophot, completeness, lsbg, icl, moving, trails, stacking, photoz_sed, sedcodes, lensmodel, spectra, xmatch, lightcurves, cluster, batch, repro, report, ai_services, ... (`docs/plugins.md`) |
| `plugins/ds10core/` | **Astrafex Core의 Desktop shell**. `ds10.py`가 ds10-web의 `server/ds10core`를 찾아 **별도 process**로 `python -m ds10core ...` 실행 (import하지 않음 = 라이선스 경계). 코어 탐색: `--core-dir`, `DS10_CORE`, `~/.ds9/ogf_params/ds10core.json`, `../ds10-web/server`, `/workspace/ds10-web/server`. GUI: Measure 탭 chip "Astrafex Core"(More 메뉴). 단계 목록과 검사: `docs/ds10core.md` |
| `ogfkit/` | 공용 Python 라이브러리 (cliexpand, trails, photerr, ...) |
| `moving/`, `icl/`, `lsbg/`, `photo_z/`, `ai_bridge/`, `ai_merge/`, ... | 분석 패키지 (이동천체, ICL, LSBG, photo-z, AI 서비스 연결) |
| `regression/` | 공개 데이터 회귀 검증 세트 (`docs/testing.md` 끝) |
| `scripts/` | `run_all_checks.sh`(모든 검사), `verify_*.tcl/.sh/.py`(GUI 검사), `astrafex`(Core CLI), 재현 스크립트 |
| `docs/` | 기능 문서, `testing.md`, `progress_log.md`(R1–R25 진행 기록), `phase2_web_design.md`(웹 설계), `shots/`(스크린샷) |

## 2. 빌드와 환경

```bash
# Debian/Ubuntu 빌드 의존성 (BUILD.txt)
sudo apt-get install -y automake autoconf libx11-dev zlib1g-dev libxml2-dev libxslt1-dev libxft-dev zip tcl8.6-dev tk8.6-dev \
     libcfitsio-dev xvfb x11-utils xdotool tcl
cd OGFinder
unix/configure && make            # -> bin/ds9   (원래 box: 8 core로 수십 분)
./build_sextract.sh               # -> bin/ds9_sextract (SEP 기반 추출기; batch/repro/report 테스트에 필요)
python3 -m venv --system-site-packages ~/ogf_venv
~/ogf_venv/bin/pip install --extra-index-url https://download.pytorch.org/whl/cpu -r scripts/ogf_venv_freeze.txt
export OGFINDER_PYTHON=~/ogf_venv/bin/python3     # run_all_checks 기본값은 /workspace/ogf_venv/bin/python3
export OGF_TEST_FITS=/path/to/fits                # m51.fits 가 있는 폴더 (기본 /workspace/fits); HUDF 검사는 hudf_f160w.fits 도 사용
```

* `scripts/ogf_venv_freeze.txt` = 원래 box에서 검사가 돌던 venv의 package 목록 (제품 dependency 목록이 아니라 기록). GPL(`rebound`, `assist`), LGPL(`sep`) 주의.
* **빌드가 tracked 파일을 바꾸거나 지운다** (`ast/config.h.in`, `tcllib/modules/*/build/*`, `tkimg/compat/libtiff/build/*` 등). 이런 변경은 **절대 commit하지 않는다**: 항상 `git add <지정 파일>`로만 commit. 특히 `ast/config.h.in`은 commit 금지.
* `m51.fits`는 600x600 2D 영상이다 (출처는 저장소에 기록되어 있지 않음). GUI 검사 일부(`cat_behavior` golden 등)는 원래 box의 그 파일 기준이라, 다른 영상이면 golden 비교가 다를 수 있다.
* `ds10core_*` 검사와 Astrafex Core 단계는 sibling `../ds10-web` checkout이 필요하다 (없으면 SKIP).

## 3. 검사 (`scripts/run_all_checks.sh`)

```bash
flock /tmp/astrafex_heavy.lock scripts/run_all_checks.sh            # 전체 (long 검사 포함)
flock /tmp/astrafex_heavy.lock scripts/run_all_checks.sh --quick    # long/network 검사(session_replay, link_bench, moving_session, regression_data) 제외
scripts/run_all_checks.sh --list                                     # 검사 이름
scripts/run_all_checks.sh --only geometry,mouse                      # 일부만
```

* GUI 검사는 `--display`(기본 `:77`)에 X가 없으면 Xvfb(1400x1000x24)를 직접 띄우고, 끝날 때(정상/에러/INT/TERM/HUP) PID로 정리한다. Xvfb는 fd>2를 닫고 띄워서 `flock`의 lock fd를 물려받지 않는다.
* 로그: `$OGF_CHECK_OUT` (기본 `/tmp/ogf_checks.<pid>/NAME.log`). 끝나면 지운다.
* 검사 표와 각 검사의 필요조건: `docs/testing.md`. live AI agent 테스트는 `OGF_LIVE=1`일 때만.
* `link_bench`는 `/workspace/work/inj*.json.pkl` 주입 세트, `moving_session`은 캐시된 BB89 노출(`~/.ds9/moving_cache`)이나 네트워크, `regression_data`는 `$OGF_DATA_CACHE`(기본 `~/.cache/ogfinder_regression`)에 공개 데이터를 받는다. 없으면 SKIP (실패 아님).
* 실행 중에는 저장소 파일을 수정하지 않는다.

### 3.1 최신 결과

인수인계 시점의 `run_all_checks.sh` 전체 표는 이 저장소에 없다. 검사 로그는 끝나면 지우고, `docs/progress_log.md` R24도 "full-run result is in the final report of this session (log kept at scripts run output, not in the repo)"라고 적는다.

저장소에 숫자가 있는 마지막 전체 실행은 `docs/progress_log.md`의 "Final full run"이다.

* HEAD `5d3271b7e` (2026-10-02): 0 failures, 35 PASS, 1 SKIP (`agent_real`, agent CLI 없음).
* 시간: newplugins 708 s, session_replay 183 s, link_bench 191 s, moving_session 629 s.
* 출력 파일은 원래 box의 `/workspace/work/run_all_checks_full.txt`이다. 저장소 밖이라 여기에는 없다.

그 뒤에 `ds10core`, `aperphot_gui`, `icl_tests`, `nightlink_tests`, trails stack이 들어왔으므로, 35 PASS는 인수인계 시점의 검사 목록과 같지 않다. 통과로 적힌 `geometry`는 R1 (HEAD `ac4fff4be`, table 181/769/154)까지다.

### 3.2 알려진 문제

* **geometry 검사의 Tcl `ALPHA_J2000` 오류**: `scripts/verify_geometry.tcl`이 추출 결과 행에서 `dict get $r ALPHA_J2000`을 하는데, 환경(추출 결과/WCS, 테스트 영상)에 따라 그 key가 없으면 Tcl 오류로 실패한다. 코드 회귀가 아니라 환경 문제로 분류해 왔다. 실패 로그 자체는 저장소에 없다. 스크립트는 `CatalogPanelExtract` 뒤 카탈로그의 3·4·5·7번째 행에 `ALPHA_J2000`과 `DELTA_J2000`이 있다고 가정한다 (`scripts/verify_geometry.tcl` 43–46행). 추출 결과에 그 열이 없으면 `dict get`이 Tcl 오류로 끝난다. 사용 영상은 `$OGF_TEST_FITS/m51.fits`(기본 `/workspace/fits/m51.fits`)이다.
* **`scripts/verify_ai_gui.py`, `scripts/verify_agent_gui.py`가 자기가 띄운 Xvfb를 남긴다**: 해당 display에 Xvfb가 없으면 `subprocess.Popen(['Xvfb', display, ...])`으로 띄우고 끝날 때 정리하지 않는다 (fd도 닫지 않아 lock fd를 물려받을 수 있음). `run_all_checks.sh` 안에서는 먼저 Xvfb를 띄우므로 문제 없고, **단독 실행할 때** 남는다. 끝난 뒤 `ps -ef | grep Xvfb`로 PID 확인 후 `kill <PID>`. 고칠 일: `run_all_checks.sh`처럼 띄운 Xvfb의 PID를 기억해 종료 시 kill, fd>2 닫기.
* trails: SDSS 후보 29개의 육안 확인은 2026-10-08에 끝났다 (spike 16, 별 5, 연속 streak 없음 8, 위성 trail로 확인된 것 없음). 표는 `plugins/trails/validation/trails_validation.md`. 실제 satellite-trail truth sample은 2026-10-08에 ACS/WFC flc 8장, 라벨 9개로 기록했다 (recall 6/9, precision 6/12). `jc8m32j5q`에서 이전에 적힌 3개는 그 노출에서 검출기가 낸 선분이다. 이동천체 linker의 검색은 two-body이고, 연결된 그룹의 N-body 보정은 별도 프로세스다 (`--nl-nbody`, `docs/moving_objects.md`). 이 머신에서는 ASSIST 없이 CODES로 입자 하나만 측정했다. `--nl-ranging`은 그 그룹의 궤도요소 분포다. 1시간 호는 `iod.ranging`(a 16/50/84% = 0.78/1.77/5.18 AU), 21일 호는 국소 Laplace 샘플(a 2.6941/2.6961/2.6985 AU, 스케일된 공분산으로 다시 잰 값)이다. OpenOrb 실행 파일은 쓰지 않았다. (1862) Apollo 관측 12줄의 예제는 `examples/apollo_1862_ranging/`이다. 자세한 남은 일은 ds10-web `HANDOVER.md` 8장.

## 4. 최근 commit (이번 인수인계 시점)

아래 목록은 2026-10-08, master `ad8761063` 기록이다. 2026-10-12의 HEAD `4dee078b467d6f1dac1aaef8ef2f80a55f577222`와 그 사이 24개 commit은 6장이다.

master `ad8761063` (2026-10-08, 이 문서) 기준 최근 commit. 한 줄은 `git log` 제목이다. 전체 이력은 약 5573 commit이고, 기능 단위 기록은 `docs/progress_log.md`다.

* `16798dd79` (2026-10-08) trails stack을 실제 HST dither (jc8m32010)에서 돌리며 WCS 왜곡, CR 제거, 스케일을 고침.
* `fa0c7cfe1` (2026-10-07) `run_all_checks`와 moving-session 검사가 띄운 Xvfb를 모든 종료 시그널에서 PID로 정리. lock fd를 물려받지 않게 fd>2를 닫고 띄움.
* `aaf0114a3`, `1fa6218ba`, `c6c5dc819` (2026-10-07) ICL: 공유 마스크와 밴드별 영점, hot+cold 마스크, μ가 정의된 픽셀만으로 f_ICL. Abell S1063 숫자를 설명서에 반영.
* `8f9e9bd28`, `885162e78` (2026-10-07) Astrafex Core 별 구경 측광 단계와 `aperphot_gui` 검사, 데스크톱 스크린샷.
* `2932ee850`, `9f6f71461`, `0b821cfa4` (2026-10-07) 하늘 카탈로그, MAST, 사용자 플러그인 단계와 `verify_ds10core` 검사.
* `f54dda631` (2026-10-06) PSF 측광 / CMD 단계.
* `60373f050`, `0c9d0fb8a`, `d28f5b6a3` (2026-10-05–06) LSBG finder와 `--classify`를 데스크톱에 연결.

## 5. 규칙 요약 (전체는 ds10-web `HANDOVER.md` 7장)

파일 지정 commit만; `ast/config.h.in` commit 금지; 무거운 작업은 `flock /tmp/astrafex_heavy.lock`으로 하나씩 (RAM 15 GB, swap 없음); process는 PID로 kill (`pkill -f` 금지); push는 사용자가 지시할 때만; git identity `Juhan Kim <kjhan0606@gmail.com>`; 기능 추가 시 `docs/`와 ds10-web 사용자 설명서도 갱신.

## 6. 2026-10-12 추가

이 머신의 작업 폴더는 `/home/kjhan/BACKUP/Astrafek`이다. HEAD는 `origin/master`와 같은 `4dee078b467d6f1dac1aaef8ef2f80a55f577222`이다. 스케일 ICL, 멤버 맞춤, HSC 디버그가 멈춘 지점, 커밋에서 빼 둘 파일은 ds10-web `HANDOVER.md` 11장에 있다. 여기에는 이 저장소에서 이어서 볼 것만 적는다.

### 6.1 `ad8761063` 다음 commit

`git log --oneline ad8761063..HEAD` 24개. 한 줄은 제목이다.

* `4dee078b4` Add the CAMIRA member match and the masked kpc stack step.
* `893b6e3af` Keep the full spectroscopic slice for the HSC stack list.
* `7956b45b6` Save the HSC stack debug procedure as scripts for the paper.
* `7e57dcd68` Wire the scaled BCG+diffuse stack into the ICL menu.
* `c8d3e44a0` Count lens sources at the named-data gate before any observed sentence.
* `0b81eb9bc` Draft the methods note for the scaled light stack and the lensing mass.
* `dd26ab3f0` Add a tangential postage-stamp diagnostic beside the shear catalogue.
* `5c72e470f` Average tangential shear profiles before the next NFW fit.
* `dd4a87e9c` Name the scaled stack BCG+diffuse and keep a kpc comparison beside it.
* `e441eb7f5` Record the synthetic ICL and weak-lensing self-consistency run.
* `e7cfadc2a` Label satellites and background galaxies for ICL masks and weak lensing.
* `1bbecb041` Scale cluster images by apparent size and M_200c before an ICL stack.
* `3bbaa6dd2` Add cluster weak lensing from metacalibration through an NFW mass.
* `9a4adcb8f` Add learned spectrum steps for a temperature letter, galaxy neighbours, and a user template match.
* `01e6610a3` Add a supernova type from one host redshift and the line depths.
* `8b2325e68` Add a synthetic spectrum example for preprocessing and the science step.
* `0282c122c` Classify one spectrum from the lines that actually fit.
* `92ef9a1f0` Draw the default galaxy stamp with an in-tree Sérsic, shear, and Gaussian PSF.
* `b7c00b6e0` Use in-tree models for the standalone measurements, orbits, and SED fits.
* `7b3b757cb` Record orbit ranging for linked groups and stop the GUI checks leaving Xvfb up.
* `e2f307324` Refine accepted nightlink groups in an N-body child process.
* `1d3edf8b6` Record the ACS satellite-trail truth sample (recall 6/9, precision 6/12)
* `2984e9971` Record the visual inspection of the 29 SDSS trail candidates
* `8984e74f7` HANDOVER.md: fill check-result, geometry, and git-log placeholders from the repo

### 6.2 이 저장소에서 이어서 볼 경로

* `analysis/icl_wgl_membership/`: 합성 자기일관성 (`run_analysis.py`, `history.json`), `implementation_plan.md`, `methods_draft.md`. `history.json`의 `parent_commit`은 `9a4adcb8f421838115acbff9c8612671d4fff723`이다. 그 파일은 고치지 않고, 스크립트는 다시 돌리지 않는다.
* `analysis/icl_wgl_membership/hsc_debug/`: 논문에 함께 제공하는 디버그 스크립트 01–08과 `run_debug.sh`. 실행 데이터는 `/gpfs/kjhan/Astrafek/data`에만 있다. 카탈로그, FITS, 로그는 커밋하지 않는다.
* 메뉴에 연결된 스택은 `plugins/icl/icl_stack.py`의 `Stack BCG+diffuse`이다. 산출물 키 `bcg+diffuse`는 BCG를 남긴 평균이다. 가린 kpc 평균의 이름은 `bcg-and-members-masked`이고 `hsc_debug/08_mask_stack.py`가 쓴다. 2026-10-12 실행은 컷아웃 HTTP 401에서 멈추었고, `cutouts/frames.tsv`가 없어 평균 파일은 쓰이지 않았다. `science_measurement`는 거짓이다.
* 접선 우표 진단은 `ogfkit/wlcluster.py`의 `stack_tangential_stamps`이다. 메뉴에는 연결되어 있지 않다.

### 6.3 작업 트리에 남은 파일

2026-10-12에 커밋하지 않고 남겨 둔 파일이다. 인수인계 커밋에 넣지 않는다.

* 수정: `ds9/library/ds9_psf_deconv.py`, `plugins/galaxy_model/galaxy_model.tcl`, `plugins/star_psf/plugin.json`, `plugins/star_psf/star_psf.tcl`, `psf_deconv/psf/sim_psf.py`
* 미추적: `plugins/galaxy_model/galaxy_fit.py`, `plugins/galaxy_model/tests/`, `psf_deconv/tests/test_diffraction_psf.py`

나선 은하는 exponential disc(preset `exp`)로 둔다. bulge+disk는 연결하지 않는다. `ast/config.h.in`과 빌드 산출물은 커밋하지 않는다.

### 6.4 2026-10-12 승인된 HSC kpc 시범과 인증 재확인

앞 절의 401은 최초 실행 기록이다. 사용자가 새 `hsc_debug/probe_login.py`를
터미널에서 실행한 뒤 HTTP 200을 받았고, 이어 `--global-sky`로 공식 매뉴얼과
`coadd/bg` 한 장의 모든 HDU 헤더를 받았다. 매뉴얼은 이를 `{rerun}/deepCoadd`
의 coadd 및 공격적 배경 차감 전 중간 자료로 설명한다. 헤더는 HSC-I, inherited
primary `FLUXMAG0=63095734448.01944` (AB 영점 27), HDU 1 IMAGE, 2 MASK, 3 VARIANCE를
확인한다. 이에 이 한 프레임의 제품 출처와 영점 매핑은 통과했다. 판정은 실행
데이터의 `logs/cutout_global_review.json`에 남겼다. 영상 픽셀은 probe가 폐기했으므로
대량 다운로드와 `cutouts/frames.tsv` 작성, 관측 스택은 아직 실행하지 않았다.
`MP_*` 헤더 값은 bit 위치이고 실제 제외 비트 값은 `2**n`이다. 검출 마스크에서
`DETECTED`만으로 BCG 픽셀을 제거하지 않는다.

승인된 첫 목표는 질량을 쓰지 않는 HSC-I BCG+확산광 kpc 프로파일이다.
`analysis/icl_wgl_membership/pilot_README.md`에 입력 스키마와 명령이 있고,
`hsc_pilot.py`의 prepare/inspect/measure, `icl/pilot.py`, `pilot_job.sh`가 새 경로다.
BCG를 남기고 모든 다른 검출을 가린다. 개별 원본 고리 프로파일의 동일 은하단
가중 평균이 주 추정량이며, 면적 보존 영상은 별도 대조다. matched random과
영역별 jackknife 경로가 있으나 실제 random·PSF·sky 전달·주입/회수 검증은
미완료이다. 모든 시범 산출물의 `science_measurement`는 false이다.

실행 데이터 아래 `pilot_2026-10-12/pilot_manifest.json`에 richness 순위로 12개를
준비했다. 영상·교정·품질 입력은 비어 있으며 `inspection_initial/inspection.json`
검사에서 12개 모두 멈춘다. 관측 스택은 실행하지 않았다. Slurm 작업 420654의
회귀 검사는 35 passed (8.94 s); 로그는 `logs/hsc-icl-pilot-420654.out`이다.
`mass_interpretation`은 입력 공분산을 기록할 뿐이며 `covariance_recorded`와
`covariance_propagated=false`로 구별한다. 자세한 전체 갱신은 Web 인수인계 11.9다.
이 추가 작업은 아직 커밋·푸시하지 않았다.

## 7. 2026-10-13 장기 연구주제

인증 조회 후 확인: catalogue_status는 161개 중 patch_qa 1개만 완료했고,
56,994행의 파일 해시는 일치했다. extinction_000(CAS ID 1059397)의 `.part`는
CAS CSV 첫 열의 `# ` 접두사 때문에 검사에서 중단됐다. 파서를 수정해 실제 파일
10개 대상/11행을 읽고 회귀 10건을 통과했다. 기존 `.part`는 완료로 가정하지 않고
동일 명령 재실행 시 기존 CAS ID에서 다시 받는다. patch_qa는 건너뛴다.
아직 소광 전체·random 목록은 완료되지 않았다. 상세 ancillary_README.md.

다음 준비에서 `validation_sample_2026-10-13_01`에 면적 조건 충족 40개와 마스크
stress 10개의 층화 예비 표본을 만들었다. 원본 800개는 유지한다. 활성 소광/대조장은
여전히 미입력이며 관측 측정을 실행하지 않았다. 표본·경로 독립 검증 및 선택 회귀
2건 통과. 사용자 지시로 숨김 암호를 실행당 1회 입력하는
`analysis/icl_wgl_membership/hsc_debug/fetch_ancillary_catalogues.py`를 추가했다.
161개 순차 CAS 쿼리의 결과/작업 ID/해시를 저장·재개하며 기존 HSC 영상 다운로더는
변경하지 않았다. 모의 CAS 회귀 5건 통과, 실제 인증 조회는 사용자 실행 대기다.
실행 명령과 결과 경로는 hsc_debug/ancillary_README.md에 있다.

후속 자료 상태: HSC-I 800개 다운로드 완료, 실제 파일 누락·크기 불일치 0,
42.24 GiB. `sample_800_2026-10-13/quality_2026-10-13`의 전체 CPU QA에서
구조·WCS 오류 0이다. 기존 마스크/고리별 면적 50% 기준으로 0–500 kpc 모든 고리를
통과하는 동일 표본 165개, 고리별 360–433개다. 800개를 모두 전 반경에 쓸 수 있다고
가정하지 않는다. `ancillary_2026-10-13_01`에는 SFD/i·i2 소광 후보 800행과 출처,
독립 조회 검증, 공식 random/patch_qa/a_i 조회 준비 파일이 있다. 소광 후보는 아직
활성 보정에 적용하지 않았다. 대조장 좌표·영상은 0이며 현재 환경에 HSC 인증이
남아 있지 않아 공식 CAS 조회는 미실행이다. 암호 재요청 없음.
`analysis/icl_wgl_membership/hsc_debug/ancillary_README.md` 및 실행 계획의
"후속: 다운로드 완료와 소광·대조장 준비"를 현재 기준으로 읽는다.

사용자는 **ICL과 위성은하 외곽구조의 관계를 공동 스태킹으로 연구**하고 장기
주제로 유지하도록 지시했다. 기준 문서는
[과학 연구계획](analysis/icl_wgl_membership/science_research_plan.md)이며,
[분석 도구 설계](analysis/icl_wgl_membership/research_software_design.md)와
[실행 계획](analysis/icl_wgl_membership/implementation_plan.md)을 함께 읽는다.

은하단 중심의 BCG+확산광과 위성 중심의 외곽광을 같은 은하단 표본에서 측정한다.
ICL의 조작적 정의, 공동 주입/회수, 광량 분배가 만드는 가짜 상관, 은하단 단위
공분산을 검증한 뒤 관계를 해석한다. 약중력은 후속 확장이다. 기존 213개 영상의
검사와 총 800개 후보 준비는 실행 계획에 있다. 실제 전송 상태는 실행 데이터의
상태 파일로 확인하며, 위 6장의 최초 401·헤더 시범 기록을 현재 상태로 쓰지 않는다.

다음 우선 작업은 추가 측광·PSF 입력 감사와 30–50개 예비 표본 설계다.
CPU 분석은 `lageunha`, 대용량 데이터는 `/gpfs/kjhan/Astrafek/data`에 둔다.
장기 연구주제는 등록됐으나 공동 분석 코드는 설계 단계이고 관측 과학 결과는
미완료다. 인증정보를 문서·설정에 남기거나 기존 다운로드 방식을 새로 작성하지 않는다.

## 8. 2026-10-13 목적 기반 파이프라인 도구

사용자는 “구상성단 별 측광”, “LSBG” 등 분석 목적에 따라 파이프라인을 자동
구성하고, 그 흐름을 따라 최종 결과를 얻는 도구를 요청했다. 요구사항·구조·
예시·개발 순서는 [science_workflow_designer.md](docs/science_workflow_designer.md)에
등록했다. `astrispex`는 사용자 표현을 따른 임시 명칭이며 제품명 변경은 미확정이다.

기존 plugin/CLI/headless, batch, session, repro를 기반으로 목적·자료 감사,
과학 recipe, 입출력 계약을 검증하는 DAG planner와 단계별 실행을 설계한다.
없는 기능을 자동 실행 가능으로 표시하지 않고, 자료 부족·상한·검증 실패도
결과에 기록한다. ICL–위성 공동 연구는 이 도구의 장기 과학 recipe다.
현재는 설계 문서만 등록했으며 planner·신규 UI·공동 분석 엔진은 미구현이다.

### 8.1 같은 날 추가 구현: Desktop / Web 순서도와 저장

위 문서 작성 이후 두 저장소에 11개 분야의 순서도, 단계 도구 연결, 검토 메모,
중간 저장·복원을 구현했다. 사용법은 [science_pipelines.md](docs/science_pipelines.md).
Desktop은 `plugins/workflow_designer`, `ogfkit/science_workflow.py`이며
Measure → Workflows에서 연다. Web은 별도 구현으로 동일한 선언적 recipe 데이터만
공유하고 Desktop 코드를 import하지 않는다. 모든 레이아웃에서 Ctrl-K → Science
pipelines로 열 수 있다.

Desktop은 완료 산출물과 카탈로그를 새 스냅샷에 복사하고 SHA-256을 검증한다.
원본 영상은 경로·크기·수정시간 참조이므로 별도 보존한다. 자동 기록되지 않은
결과 파일은 사용자가 연결한다. 복원은 작업 사본을 만들어 저장본을 보존한다.
Web은 계정 세션의 저장 시점까지 완료 산출물을 ZIP으로 내려받고, 서버의 세션
저장 지점으로 복원한다. 누락·변조 파일은 저장 성공으로 처리하지 않는다.

Desktop 저장/복원 6건, Web API 3건과 브라우저 선택·메모·ZIP·특정 저장 지점
복원 1건이 통과했다. 타입 검사·Web 빌드·Desktop manifest 검사도 통과했다.
검사는 CPU로 lageunha에서 수행했다. Desktop Tk 화면의 실제 클릭 검사는
아직 하지 않았다. 자연어/입력 감사 기반 DAG 자동 생성, 미연결 분석 adapter,
모든 분야의 최종 과학 검증은 후속 개발이며 현재 UI의 `reviewed`는 수동 검토다.
이 작업은 커밋·푸시하지 않았다.

### 2026-10-13 사용자 UI 수정: 본창 파이프라인

사용자는 팝업 대신 본창의 분석 메뉴를 클릭하면 왼쪽 위 PIPELINE 단계 목록이
해당 분야에 맞게 바뀌도록 지시했다. Web은 11개 분야를 분석 메뉴에 직접 넣고,
선택 단계의 설명·메모·저장을 중앙 분석 영역으로 옮겼다. 일반 도구는 기본
단계 목록으로 연결한다. Desktop도 별도 toplevel을 제거하고 기존 분석 패널에
순서도와 선택기를 넣었다. Web의 이전 캡처 4장은 같은 파일 경로에 본창 캡처로
갱신했다. 저장/복원 동작은 유지하며 과학 검증 완료와 수동 검토를 구별한다.

### 2026-10-13 추가 구현: 단계 실행과 작업 연결

Web 중앙 단계 영역에 직접 실행과 작업 ID·산출물 목록을 추가했다. 단계별
`runs` 참조는 계정·기록 세션·분석기 일치를 서버가 검사하며 실행 상태는 실제
작업 기록에서 확인한다. 이전 저장 지점은 당시 완료 산출물만 포함한다. 완료
기록이 아직 없는 저장 지점은 완전한 결과 저장으로 표시하지 않는다. 서버에
남은 작업은 재개 시 기존 SSE로 관찰하며 새 작업을 중복 제출하지 않는다.
Desktop은 순서도에서 직접 실행한 도구의 recorder 항목을 단계에 연결하고
write trace/idle callback으로 상태를 갱신한다. 실행 중 스냅샷을 새 프로그램에서
열면 재검토로 바뀌며 계산 프로세스가 복원됐다고 표시하지 않는다. 앞 단계
재실행·입력/설정 변경은 후속 결과를 재검토 상태로 바꾼다. 수동 검토와 과학
검증은 계속 구별한다. 사용법은 `docs/science_pipelines.md`에 추가했다.

이번 추가 구현 검사는 Desktop 7건, Web API 6건, 브라우저 3건 통과다. 실제
합성 영상 검출·산출물 표시·새로고침 재개·설정 변경 재검토를 포함한다. 타입
검사·빌드도 통과했고 CPU 검사는 lageunha에서 수행했다. 본창 캡처 4장은
갱신했다. Desktop 실행 상태는 Tcl로 검사했으며 실제 Tk 클릭 검사는 아직
하지 않았다. 커밋·푸시와 관측 과학 분석 실행은 하지 않았다.

### 2026-10-13 미연결 분석단계 추가 연결

영상 입력 검사와 결과 저장을 실제 연결했다. 원본 공유 코어에 동일 은하단 가중 kpc 프로파일 stack을 추가했다. Web에는 계정 소유 TSV import 및 통합 PSF/LSBG·구경 측광·PSF stamp·검출 주입 회수 adapter를 연결했다. 통합 실행과 개별 검토를 구별한다. Desktop 명시적 산출물도 recorder에 등록한다. 계약과 남은 전용 분기는 docs/science_pipelines.md 마지막 절 참조. 합성자료 검사: Desktop 9, 서버 15, 브라우저 4개 통과. 실제 HSC 계산/다운로드는 수행하지 않았으며 science_validated=false다.

### 2026-10-13 다음 구현: 1D 분광 준비

공유 코어 spectrum_prepare의 단일 1D TSV 입력 검사·오차 가중 Chebyshev 연속광을 양쪽에 연결했다. 웹 계정 소유 import와 Desktop catalog 로드, 산출물·입력 표·공분산 저장/재개를 포함한다. 파장계와 flux 단위는 명시 필수이며 AGN/host/Fe II 분해를 대신하지 않는다. 웹의 미구현 후속 1D 과학 adapter는 버튼을 비활성화한다. Desktop workflow의 close 이름 충돌도 ::close로 수정했다. 입력 계약과 남은 단계는 docs/science_pipelines.md 마지막 절 참조. 실제 관측 데이터 계산은 수행하지 않았다.

1D 분광 최종 검증: Desktop 11, 서버 25, 브라우저 5개(총 41) 통과. Web typecheck/build 통과; Desktop 39 plugins/289 steps/873 params manifest 검사 문제 0. 검증은 lageunha CPU에서 수행했다. 연속광 출력에도 원 MASK를 유지하며 각 단계의 입력 snapshot 파일은 서로 다른 이름으로 보존한다.


### 2026-10-13 후속: 동정된 고립선/적색편이 후보

두 제품의 항성/은하 선 측정 및 은하/AGN 적색편이 노드를 original shared-core `spectrum_lines.py`에 연결했습니다. 사용자 rest 파장·초기 z·반창 폭·방출/흡수 입력, 국소 연속광+Gaussian 공동 적합, covariance/모델/원본/측정표 보존. 여러 accepted 선이 일치할 때만 조건부 z 후보를 보고하며 자동 선 동정·전역 검색 및 AGN 분해는 미구현입니다. 상세 사용법·판정 범위는 docs/science_pipelines.md 마지막 절에 있습니다. 실제 HSC 데이터 분석·다운로드, 커밋·푸시는 수행하지 않았습니다.

추가로 Web `sessionController.ts`의 영상 없는 workflow 복원에서 활성 분석 표를 생략하던 문제를 수정했습니다. 저장→새로고침→후속 선 측정→다시 ZIP 저장 흐름을 브라우저 회귀 검증에 포함합니다.

검증(lageunha, CPU/BLAS 1 thread): 서버 기존 25개 및 새 선 측정 15개, Desktop 12개, 실제 브라우저 5개 통과(총 57개). Web typecheck/build 통과. Desktop manifest 39 plugins/291 steps/878 params, 0 problems. 변경 파일 diff whitespace 검사 통과. 합성 자료 검증이며 실관측 과학 검증은 수행하지 않았습니다.


### 2026-10-13 후속: 초기 z 없는 피크 탐색/적색편이 후보

항성·은하·AGN 분광에 search 노드를 추가(11 recipes/110 nodes). Web ds10_spectrum.search / Desktop ds10core.spectrum-search, original NumPy/SciPy shared core. 사용자 기준선 목록, 구간별 median detrending/피크 대조, 다중 대응과 ambiguity 보존; 후보표·피크표·JSON·입력 원본을 기록/체크포인트/ZIP에 포함. 과학 승인·설정 자동 변경을 하지 않으며 후보를 검토해 다음 정밀 측정 설정에 반영합니다. AGN 분해는 아직 미구현입니다. docs/science_pipelines.md 마지막 절에 계약·한계가 있습니다. 실제 관측 데이터 계산·커밋·푸시는 하지 않았습니다.

검증(lageunha, CPU/BLAS 1 thread): 서버 52개·Desktop 13개·브라우저 5개, 총 70개 통과. Web typecheck/build 통과. Desktop manifest 39 plugins/292 steps/882 params, 0 problems. 초기 z 없는 합성 방출/흡수선 회수, 모호한 대응 보존, 잡음/마스크/스파이크, 실제 owned job 및 저장/재개 검증.


### 2026-10-13 후속: AGN 방출선 공동 적합

AGN lines를 Web ds10_spectrum.joint / Desktop ds10core.spectrum-joint에 실제 연결했습니다. 공통 z에 묶인 여러 방출 Gaussian+국소 선형 연속광, 사용자 observed sigma bounds/fit-range, 결정적 3-start, 차원 정규화 퇴화 검사, 전체 parameter/성분 flux covariance, 판정·잔차·MASK·원본 저장. 이름 broad/narrow의 물리적 해석이나 host/Fe II 템플릿 분해를 대신하지 않습니다. docs/science_pipelines.md 마지막 절에 계약과 한계가 있습니다. 실제 관측 계산·커밋·푸시 없음.

검증(lageunha, CPU/BLAS 1 thread): 서버 66개·Desktop 14개·브라우저 6개, 총 86개 통과. Web typecheck/build 및 Desktop manifest(39 plugins/293 steps/885 params, 0 problems) 통과. 합성 겹침선/동일 중심 좁은+넓은 성분/작은 물리 flux 단위/마스크/잘못된 모형과 계약을 확인했습니다. 실제 owned job, Desktop 체크포인트, Web ZIP 저장→새로고침→동일 단계 재실행도 통과했습니다.


### 2026-10-13 후속: 사용자 host/Fe II 템플릿 분해

AGN decompose를 두 제품에 실행 연결했습니다. 같은 관측 TSV의 TEMPLATE_HOST/FEII, 출처·frame/grid/redshift/resolution template-note 필수. 고정 F_lambda index/pivot + nonnegative AGN/host/FeII 세기, 퇴화·경계 검사와 조건부 공분산. joint continuum-model=templates에서 template 세기와 선을 원본 flux에 함께 재적합하여 전체/선 flux 공분산에 반영합니다. 분해표·원본·JSON은 저장/재개합니다. 템플릿 모양/폭/AGN 지수 자체 탐색·블랙홀 질량은 남아 있습니다. 실제 관측 계산·다운로드·커밋·푸시 없음. 상세는 docs/science_pipelines.md 마지막 절.

검증(lageunha, CPU/BLAS 1 thread): 서버 75개·Desktop 15개·브라우저 7개, 총 97개 통과. Web typecheck/build, Desktop manifest(39 plugins/294 steps/889 params, 0 problems) 통과. 합성 계수/선 회수·작은 물리 flux 단위·퇴화/경계·입력 계약, owned job·체크포인트·ZIP 저장→새로고침→template joint 재실행 검증.


### 2026-10-13 후속: 공동 적합 선비 진단

AGN diagnostic 노드를 두 버전에 연결했습니다. 선행 lines의 joint.json을 자동 입력하며 활성 스펙트럼 유지. 성분 합계 ratio/log10와 전체 flux covariance를 전파한 오차·비 사이 covariance, S/N/품질/행렬 검증, source report snapshot/SHA 기록. Desktop 기록·복원 artifact 경로 연결, Web bound owned job/artifact 연결. 낮은 S/N/flagged 성분에 숫자나 상한을 만들지 않습니다. 자동 분류·소광·광도·블랙홀 질량은 남아 있습니다. 실제 관측 계산·커밋·푸시 없음. 상세는 docs/science_pipelines.md 마지막 절.

검증(lageunha, CPU/BLAS 1 thread): 서버 87개·Desktop 17개·브라우저 7개, 총 111개 통과. Web typecheck/build 및 Desktop manifest(39 plugins/295 steps/892 params, 0 problems) 통과. 선비·log 및 교차 covariance의 해석식/작은 물리 단위·잘못된 행렬/flags/SNR/소유권 검사, Desktop 복원 파일 연결·최신 실패 차단, Web 활성 스펙트럼 유지·ZIP 저장→새로고침→진단 재실행 검증.


### 2026-10-13 후속: 분류·소광·선 광도·블랙홀 질량 통합 해석

남아 있던 네 항목을 AGN inference 노드로 묶어 두 제품에 연결했습니다(11개 목적/111개 단계). 선행 lines/joint.json 자동 연결, Web은 popup 없이 단계 내 JSON 편집기, Desktop은 spectrum-inference-config. 명시한 좁은선 BPT와 모호 경계, 외부 또는 Balmer E(B-V)/사용자 k, 명시한 광도거리·cgs 변환의 선 광도, 확인된 넓은 Hα와 기기 Gaussian 폭의 Greene-Ho 2005 식6 조건부 질량을 계산합니다. 전체 공동 적합 covariance 및 독립 보조입력 오차를 전파하고 측정오차와 사용자 채택 scatter는 분리합니다. 생략은 not_configured, 부적절/부족한 입력은 blocked+이유이며 job 성공은 보고서 생성 성공입니다. science_validated=false/requires_review=true. 네 산출물과 원본 적합·설정은 체크포인트/ZIP 저장·복원에 포함됩니다. 상세 계약·합성 예시·1차 논문 링크는 docs/science_pipelines.md 마지막 절.

검증(lageunha, CPU/BLAS 1 thread): 서버 100개·Desktop 18개·실제 브라우저 8개, 총 126개 통과. Web typecheck/build, Desktop manifest(39 plugins/296 steps/894 params, 0 problems) 통과. 네 결과 값/교차 covariance, 질량오차의 독립 유한차분 검증, 경계 모호성/음의 소광/오류 가정/해상도 차단, owned API·Desktop subprocess와 체크포인트·Web inline 입력→네 결과→ZIP 저장→새로고침→재실행 검증. 실제 관측 계산·HSC 데이터 변경·커밋·푸시 없음. 원본 코드/허용 NumPy·SciPy만 사용, 추가 dependency/외부 템플릿 없음. 확정 AGN 판정/자동 선 식별, cosmology 자동선택, 소광법칙·host/Fe II 모양/virial factor 등 체계오차 posterior는 여전히 범위 밖입니다.


### 2026-10-13 통합 해석 추가 검증과 결함 수정

정상 회귀 126개 재통과 후 독립 합성 주입(z guess=0/.35, truth=guess+.0008, flux scale=1/1e-17)의 참값 회수 및 외부/Balmer 소광 질량오차 유한차분 검증을 추가했습니다. 운동학 보고서 변조 세 경우에서 질량을 차단하지 않는 결함을 발견하여 finite/양의 기준파장·유효 z·fitted shift와 common z의 일치 및 finite uncertainty 검사로 수정했습니다. 최종 서버 108개, Desktop 18개, 실제 브라우저 8개, 총 134개 통과(lageunha CPU/BLAS 1 thread). Web typecheck/build 및 Desktop manifest 0 problems. 범위·허용오차·결함 재현/수정·검증 한계는 docs/validation_spectrum_2026-10-13.md에 기록했습니다. 실제 관측 승인 검증·HSC 데이터 변경·커밋·푸시 없음.

커밋 준비 추가 검증: 커밋 대상 파일만 추출한 lageunha 체크아웃에서도 서버 126개·Desktop 18개·실제 브라우저 8개, 총 152개 및 Web typecheck/build 통과. 파이프라인 직접 의존 morphology/PSF stamp/counts 포함, 기존 런처·요금제·persona·회절 PSF 변경 제외. 상세와 기존 zero-noise 테스트 RuntimeWarning 3개는 검증 보고서 마지막 절.

### 2026-10-13 사용자 AI 도구 Import 사전검사

Tools 메뉴에 AI submission kit 내보내기와 local Import review를 추가했다.
`tools/user_tool_import.py --kit <새 ZIP>` / `--review <제출폴더>` CLI도 제공한다.
Python 3.10 이상 필요, Web 서버나 sibling checkout 없이 동작한다. 검증기는
Web first-party `tools/validate_user_tool.py`의 사본이며 제출 코드를 실행하지
않는다. `tools/user_tool_kit/`에 규약·예제·런타임 정보를 포함했다. 내보내는
kit의 target은 standalone이다. 실행 승인/자동 설치/분석 단계 등록은 별도다.

검증: `python3 regression/test_user_tool_import.py` 4개(CLI, kit/덮어쓰기 거부,
검증기 동일성, Tcl 메뉴 명령 stub 실행). Desktop 바이너리를 재빌드하지 않았으므로
기존 실행본의 메뉴 사용은 Tcl 변경 반영 후 가능하다. 상세 `docs/user_tool_import.md`.
Web 측에서는 해당 기능을 열면 로컬 모드로 자동 전환한다. 커밋·푸시 없음.


### 2026-10-13 사용자 도구의 실행 위치 분리 (최신 지시)

이전 절의 전체 로컬 모드 자동 전환 설명을 대체한다. Web 사용자 도구 Import는
독립된 브라우저 Worker에서만 검사하며 전역 서버/로컬 모드, localStorage 모드,
일반 분석의 실행 위치를 변경하지 않는다. 모드 변경 이벤트를 발생시키지 않고
사용자의 전체 모드 변경에도 검사 Worker를 취소하지 않는다. 후속 등록·실행도
도구별 browser 실행 위치로 구현해야 한다. 현재는 사전검사까지이며 사용자 도구
자동 등록·실행 연결은 여전히 미구현이다. 규약과 standalone 동봉 사본도 갱신했다.

### 2026-10-13 개인 도구 native 실행·분석 단계 편입

앞선 사용자 도구 자동 등록/실행 미구현 설명을 대체한다. Tools의 native Python
등록으로 고정 manifest 해시 버전과 plugin.json을 생성하고, Tools의 selected
pipeline 메뉴로 현재 공식 분석 단계 뒤에 사용자 단계를 넣는다. 입력·매개변수
버튼에서 inputs-json / params-json을 설정한다. 각 실행은 별도 작업 폴더에서
해시를 재검사하고 실행하며 로그·보고서·선언 산출물을 기록한다.
시간 제한과 SIGTERM 시 자식 프로세스 그룹을 정리한다. 네이티브 Python은
사용자 권한으로 실행되며 OS 보안 샌드박스는 아니다. 자동 의존성 설치는 없다.

연구 세션에는 사용자 단계 정의, 원본 코드/manifest, 결과와 실제 실행 입력
사본을 보관한다. 복원은 코드 검증/재등록과 입력 경로 재연결만 수행하며 실행하지
않는다. 등록 버전이 바뀌면 기존 노드는 실행을 막는다. 기존 공식 recipe 파일은
수정하지 않고 사용자 세션 안에서만 순서도를 확장한다.
사용법: docs/user_tool_execution.md. 검사: regression/test_personal_pipeline.py,
regression/test_user_tool_import.py, plugins/workflow_designer/test_workflow.py의
27건이 lageunha CPU에서 통과했다. Web 관련 실제 CSP E2E 10건도 통과했다.
Desktop 바이너리는 재빌드하지 않았으며 Tcl/CLI 소스 수준에서 검증했다.
실제 과학 자료 계산, 커밋·푸시는 하지 않았다.


### 2026-10-13 HSC 50장 검출·마스크·PSF 입력 준비 완료

lageunha CPU에서 검증 표본 50/50장 처리, 오류 0. 검출 156,984개, 임시 PSF 후보
569개를 저장했다. 면적 집합 40장과 마스크 스트레스 10장은 분리 유지한다.
BCG 위치 1 arcsec 내 고유 검출 11장, 검출 없음 39장이므로 후자는 BCG 보존 미확정으로
차단했다. 원본 품질 마스크의 BCG 중심 가림은 18장이다. PSF 후보 4개 미만은
3806(1), 6940(3), 3700(2)이며 대체 Gaussian을 생성하지 않았다.
Kron 2.5/3.5/5.0 마스크, 검출/멤버 표, 원본 PSF 우표와 방사 진단, 전체/BCG 확대
이미지를 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/detection_psf_2026-10-13_01`
에 저장했다. 합성 회귀 3건과 독립 FITS/표/WCS/마스크 산출물 검사 50건 통과.
`review_queue.tsv`에 후속 검토를 기록했다. 과학 승인 마스크·PSF가 아니며
`science_measurement=false`, `mask_ready_for_science=false`다. 원본 영상 변경,
소광 적용, 대조장 차감, 관측 스태킹은 하지 않았다.
상세: `analysis/icl_wgl_membership/hsc_debug/detection_psf_README.md`.
동시 진행 중인 사용자 CAS 수집은 이번 확인 시 running/random_008이었다.
앞선 인증 조회 미실행 기록은 당시 상태이며 현재 상태로 해석하지 않는다.
커밋·푸시는 하지 않았다.


### 2026-10-13 후속: BCG 정점 대조와 멤버 마스크 연결 보완

lageunha CPU에서 50장 검출을 재현하고 deblend contrast 0.005/0.001의 원본·필터
정점과 중심 분할 번호를 대조했다. 두 설정에서 안정된 위치 후보는 35장으로,
기존 centroid 후보 11장을 유지하고 24장을 추가했다. 남은 15장은 중심 품질 마스크
8장과 정점 조건 불충족 7장이다. 고유 BCG 멤버가 미검출 위성으로 재마스킹되는
연결을 보완했으며 실제 이웃 마스크와 원본 품질 마스크는 보존했다.

산출물은 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/bcg_review_2026-10-13_01`
및 `bcg_masks_2026-10-13_01`이다. 후자가 멤버 처리까지 보완한 검토용 마스크다.
Kron 2.5/3.5/5.0에서도 중심이 가려진 영상은 31/32/35장이므로 후보 35장을
사용 가능 과학 표본으로 해석하지 않는다. 모든 과학 승인 상태는 false로 유지했다.
PSF 후보 569개의 외부 대조 좌표를 준비했으나 별 분류 확인은 아직 수행하지 않았다.
회귀 8건과 독립 저장 산출물 검사(마스크 50장, PSF 좌표 569개) 통과.
기존 결과·원본·활성 manifest를 보존했고 실제 스태킹·새 인증 요청·커밋·푸시는 없다.
상세는 `analysis/icl_wgl_membership/hsc_debug/bcg_review_README.md`.
후속은 중심부 이웃 분리/마스크 유효면적 검토, 외부 별 분류 및 넓은 PSF 날개,
주입 회수와 공식 random에 기반한 대조장 검증이다.


### 2026-10-13 후속: 마스크 면적/분리 민감도와 Gaia PSF 후보 대조

lageunha CPU에서 50장의 0–500 kpc, 25 kpc 고리·8방향 면적을 계산했다.
품질 마스크와 분리 2종×Kron 3종의 7,000행, 기존 마스크 200개의 독립 재계수
검사 통과. 신규 합성 테스트 5건 통과. 기준 contrast .005/Kron 2.5에서 BCG 후보와
모든 고리 유효면적≥50%를 동시에 만족한 것은 5724·619·5814 세 장, contrast .001에서는
619·5814 두 장이다. Kron 3.5/5.0 및 스트레스 집합은 통과 0이다. 과학 승인은 아니다.
5724는 BCG 정점/모멘트가 같아도 이웃 검출 차이로 최대 면적 차이 약 72.7%p 발생.
마스크를 임의로 축소하거나 중심을 강제로 비우지 않았으며, 이웃 분리/광량 모형과
주입 회수 검증을 다음 단계로 유지한다.

공개 CDS Gaia DR3 위치 대조 후 ESA source_id 조회로 외부 정보를 확보했다.
569개 후보 중 대응 73개, 보수적 별 근거 충족 32개(23개 영상), 근거 불충분 41개,
미대응 496개다. 입력 좌표/ID·응답 해시·천구 거리 독립 검사를 통과했다.
PSF 승인과 HSC epoch 전파는 미완료다. 지연된 ESA 공간 결합은 ABORTED로 정리했다.
HSC 별도 CAS 조회 161개도 모두 완료: 소광 800개 대상 1,212행, random 후보
800개 대상 160,000행, patch QA 56,994행. 완료 파일 감사 기록을 저장했다.

결과는 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/` 아래
`mask_coverage_2026-10-13_02`, `mask_case_5724_2026-10-13_01`,
`psf_cds_2026-10-13_01/evidence`다. 상세와 그림·검증 범위는
`analysis/icl_wgl_membership/hsc_debug/mask_psf_review_README.md`.
원본·활성 manifest 보존. 과학 스태킹·대조 영상 다운로드·재인증·커밋·푸시 없음.


### 2026-10-13 후속: 이웃 외곽·확산광의 주입/회수 예비 실험

lageunha CPU에서 실제 HSC 품질 마스크/분산/멤버 위치 6개 조건을 사용한 합성
432장면×5방법=2160개 회수 기록을 완료했다. science 픽셀은 합성이고 Gaussian
PSF·고정 기하/형태·독립 잡음을 가정했다. 실제 관측 및 upstream 하늘 전달 검증이 아니다.
결과는 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/neighbor_injection_2026-10-13_01`.
단위 검사 5건, 입력/격자 감사, 독립 열 정규화 SVD 재계산 48건 통과.

면적 집합 5개 조건의 이웃 중심+외곽 공동 모형에서 확산광 상대 오차 중앙값은
같은 모형 대조 −3.82%, 외곽 크기 오차 +11.38%, 합성 PSF 날개 누락 +3.39%,
하늘 곡률 누락 −51.65%다. 이 값은 관측 편향 보정값이 아니다. 스트레스 6940은
따로 보고했고 외곽광 형식 오차가 최대 주입량의 약 33배여서 과학 회수로 승인하지 않는다.
마스크 방식은 조건별 108장면 중 41–45개만 성공했다. 알려진 위치를 사용하는 모형
방식과 성공 사례만 비교해 우월성을 주장하지 않는다. 완전한 진폭 격자만 전달계수 계산.
독립 주입 성분의 오차에도 반상관(면적 집합 약 −0.281)이 생겨 공동 공분산 처리가 필요하다.

외곽 포함 모형은 후속 후보지만 하늘 잔차/외곽 형태/경험 PSF/상관잡음 검증은 남아 있다.
science_measurement=false, observational_validation_passed=false 유지. 원본·활성
manifest 보존, 실제 스태킹·추가 인증·커밋·푸시 없음. 상세:
`analysis/icl_wgl_membership/hsc_debug/neighbor_injection_README.md`.


### 2026-10-13 후속: 외곽 크기·하늘 제약 공동 검증

lageunha CPU/1 thread에서 이전 합성 실험과 같은 입력·잡음으로 864장면,
5방법=4320개 회수 기록, 크기별 적합 25,056건을 완료했다(약 19.2초).
4–16 kpc 공통 외곽 크기 격자, 평면/이차 하늘, 독립 합성 하늘 제약 및
편향된 제약을 비교했다. 실제 대조장 측정이나 경험적 PSF를 사용한 실험은 아니다.

면적 집합 5개 조건에서 외곽 12 kpc+하늘 곡률 동반 시 확산광 상대 오차 중앙값은
기존 고정 크기/평면 −34.42%, 자유 크기/이차 하늘 −4.63%, 합성 외부 제약 −4.02%다.
각각의 절대 상대 오차 중앙값은 34.42%, 17.00%, 9.21%이므로 장면별 정확도 4%를
주장하지 않는다. 편향된 제약은 +7.78%, 모형에 없는 국소 하늘 잔차는 합성 제약에도
+27.96%가 남았다. PSF 날개 오차도 미해결이다. 6940의 720개 방법 기록은
수치적으로 성공해도 가려진 성분의 정보 부족으로 별도 차단 표시했다.

이전 기준 432건 재현(최대 차이 1.93e-10), 독립 SVD 검산 60건(최대 3.23e-9),
신규 6건 포함 테스트 11건과 해시/격자/실패/profile 구간 감사 통과.
결과: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/halo_sky_validation_2026-10-13_01`.
상세: `analysis/icl_wgl_membership/hsc_debug/halo_sky_validation_README.md`.
크기 모호성을 반영한 구간도 보정된 신뢰구간은 아니다. 실제 대조장 하늘 공간 변동,
경험적 PSF/상관잡음과 추가 독립 주입 실험을 다음 조건으로 유지한다.
science_measurement=false, observational_validation_passed=false. 원본/활성 manifest
보존, 관측 스태킹·추가 인증·다운로드·커밋·푸시는 하지 않았다.


### 2026-10-14 후속: 잠정 대조장 좌표와 실제 외곽 잔차

lageunha CPU/1 thread에서 공식 random 160,000행(고유 ID 109,402개), patch QA,
기존 SFD 지도와 전체 보유 CAMIRA 7,319개를 사용했다. 은하단 배제 반경
0.5/1/2 Mpc에서 각각 대상 787개의 비중첩 잠정 대조 좌표를 배정했다. 기본 1 Mpc
대비 좌표 교체는 316/505개다. 13개 대상은 5각분 안의 유효 QA 대용점이 없어 남겼다.
대상 patch는 근방 random의 QA 대용값이며 정확한 patch 식별은 미확인이다.
전체 대조 영상 footprint/마스크/필터/PSF 검증과 대조 영상 다운로드는 아직 없다.
`request_plan.json`은 잠정 좌표 계획이고 과학용 또는 실행용 manifest가 아니다.

기존 검증 영상 50장(면적 집합 40, 스트레스 10)의 600–900 kpc 외곽을
Kron 2.5/3.5/5.0으로 검사했다. 150개 진단, 타일 26,906개, 상관 기록 1,200개다.
면적 집합의 타일 평면 잔차 MAD/명목 평균 오차 중앙값은 2.79/2.58/2.41,
스트레스 집합은 3.15/2.85/2.57이다. 면적 집합 외곽 유효면적은 44.7/37.3/26.0%다.
1픽셀 상관은 약 .054–.055. 은하단 빛/PSF/미검출 소스/마스크 선택도 포함되므로
하늘 보정값, 순수 잡음 커널 또는 최종 오차 배율로 편입하지 않는다.

독립 검산에서 기존 거리 격자 보간과 Astropy의 약 .002–.005% 차이로 고리 경계
픽셀 불일치를 확인했다. 이번 외곽 진단을 FlatLambdaCDM으로 통일한 `_02`가 최종이다.
첫 결과와 공용 코어/과거 주입 결과는 보존했다. 신규 테스트 8건, 입력 해시 190건,
반경별 787개 대조 좌표의 전체 CAMIRA/상호 겹침 검사, 원본 타일 독립 재계산
54건이 통과했고 타일 평균 최대 차이는 0이다.

자료 상위 경로 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13` 아래:
`control_candidates_2026-10-14_01`, `outer_sky_diagnostics_2026-10-14_02`,
`control_sky_audit_2026-10-14_02`. 상세는
`analysis/icl_wgl_membership/hsc_debug/control_sky_README.md`.
원본/활성 manifest 보존, science_measurement=false. 다음은 정확한 patch 연결,
대조장 전체 영상 검증, 경험적 PSF/배경 공분산 및 같은 타일의 마스크 비교다.
추가 인증·새 영상 다운로드·관측 스태킹·커밋·푸시는 하지 않았다.


### 2026-10-14 후속: 공식 patch 연결과 경험적 PSF 진단

lageunha CPU/1 thread에서 공개 tract·patch 다각형 91,449개를 사용해 근방 QA
대용값을 교체했다. 대상 800개 중 단일 유효 patch 607개, 복수 patch 192개,
미연결 1개다. 복수 patch 모두의 화질 조건을 만족하도록 해 비중첩 대조 좌표
798개를 선정했다. 1497은 공개 경계 미연결, 6470은 겹친 네 patch와 동시에
맞는 후보가 없어 남겼다. 실제 primary tract 및 전체 cutout 품질은 미확인이다.

`control_download_manifest.json`은 원래 성공한 `download_pilot.py`를 그대로
사용하는 실행 입력이다. matched_control 역할·대응 대상 ID·자체 좌표·BCG 없음과
미검증 소광을 보존한다. 예상 42.10 GiB이고 실제 새 대조 영상은 0개다.
현재 프로세스에 HSC 인증이 남아 있지 않아 전송하지 않았으며 인증을 재요청하지 않았다.

Gaia 별 근거 32개/23영상에 401×401 native 우표와 마스크 2종×배경 2종을
검사했다. 우표 20개 확보, 12개 경계 실패. 128개 방법 기록 중 진단 가능 30개,
배경 면적 부족 49개, 영상 경계 48개, 중심부 마스크 1개다. 일부 조건 통과 별
11개, 네 조건 모두 통과 별 4개다. 바깥 10.75–16.13초각 고리 면적 통과 기록
25개 중 23개가 배경 차감 후 0 이하다. 중심 광량 안정성을 PSF 날개 승인으로
대신하지 않았고 음수·실패 사례를 모두 보존했다.

신규 테스트 8건, 입력 해시 80건, 접평면 다각형 재검산 2,840건, 798개 전체
CAMIRA/대조장 겹침 검사, 원본 우표 20개·방사 구간 330개·마스크 표본점
3,240건의 독립 재계산을 통과했다. 원본 픽셀 평균 최대 차이는 0이다.

자료는 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13` 아래
`control_patches_2026-10-14_02`, `empirical_psf_2026-10-14_01`,
`patch_psf_audit_2026-10-14_02`. 상세는
`analysis/icl_wgl_membership/hsc_debug/patch_psf_README.md`.
빈 QA 중심을 처리하지 못한 최초 patch 실행과 마스크 HDU 참조 키를 수정하기 전
감사 `_01`은 보존했다. 원본·활성 manifest 보존, science_measurement=false,
psf_approved=false. 다음은 실제 대조 영상의 배경/공분산과 PSF 입력·epoch 검증이다.


### 2026-10-14 후속: 공개 확장 PSF 실제 확보·전달 시험

lageunha CPU에서 Garate-Nuñez et al. (2024)의 공개 i-band PSF 대칭/비대칭 FITS를
실제 내려받았다. Zenodo doi:10.5281/zenodo.14915578, CC-BY-4.0. 두 파일 합계
19,618,560 bytes(18.71 MiB), 배포 MD5와 SHA-256을 확인했다. 4001×4001 배열의
중심·유한성·단위 합을 검사했고 원본 음수도 보존했다. 외부 R 코드 실행 없음.

동일한 마스크와 배경 연산으로 기존 관측 성공 30조건×참조 2종, 방사 구간
660개를 비교했다. 참조 날개의 10.75–16.13초각 평균에서 배경 추정이 제거하는
비율은 10.32–28.39%(중앙값 26.54%)다. 남는 참조 신호는 명목 고리 평균 오차의
.210–.709배다. 관측 23/25개 방법의 0 이하 값은 그대로 두었다. 이 자료로
개별 PSF 날개를 승인하거나 공개 평균을 800장 전체의 참 PSF로 고정하지 않는다.

원본에 참조 신호를 더한 1,980개 차분 기록이 선형 연산과 일치했다.
분산·마스크 고정, photon noise 추가 없음이므로 HSC upstream sky나 재검출
완전성 검증은 아니다. 신규 6개·기존 4개 테스트, 입력 해시 31건, 독립 방사
비교 660개/신호 증분 1,980개 검산 통과. 최대 독립 증분 오차 5.68e-14다.

기존 800장 위치/크기 선별 후 관련 24개 헤더에서 별 좌표 33조합을 WCS로
확인했다. 경계 때문에 잘린 별 12개를 대체할 보유 영상은 없다. 날짜/MJD 키는
0개이고 PSF archive ID는 있어도 모델 본체는 없다. 최초 inventory `_01`의
상대 영상 파일명 해석을 고쳐 `_02`에 실행했으며 첫 폴더를 보존했다.

자료 상위 경로 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13` 아래
`extended_psf_reference_2026-10-14_01`, `extended_psf_comparison_2026-10-14_01`,
`psf_support_2026-10-14_02`, `extended_psf_audit_2026-10-14_01`.
설명은 `analysis/icl_wgl_membership/hsc_debug/extended_psf_README.md`.
대조장 798개 입력은 준비 상태이고 실제 대조 영상은 인증 세션이 없어 0개다.
science_measurement=false, psf_approved=false, reference_approved=false 유지.


### 2026-10-14 사용자 다운로드 승인과 실제 전송 재개

사용자는 새 자료 다운로드 시 대상·수량·예상 용량을 제시하고 명시적 승인을
받도록 지시했다. 이후 **HSC 대조장 남은 794장, 예상 약 41.9 GiB의 다운로드를
명시적으로 승인**했다. 승인한 동일 묶음의 전송/재개에 이 승인을 유지하며,
다른 자료나 수량 확대는 별도 승인을 받는다. 암호 자체는 이 문서에 기록하지 않는다.

이전에 대화로 제공된 인증정보를 실제로 사용한 결과 인증에 성공했고 대조장
4장의 FITS·헤더 검증을 완료했다(240,894,720 bytes). 새 승인 조건을 확인하는
동안 전송을 중단했다가, 명시적 승인 후 기존 `download_pilot.py --resume`으로
lageunha CPU에서 남은 794장을 재개했다. 중단된 다섯 번째 응답은 보존하고 재수신한다.

입력: `sample_800_2026-10-13/control_patches_2026-10-14_02/control_download_manifest.json`.
실행: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/control_download_01`.
실시간 상태는 `download_status.json`, 전체 완료 후에는 `download_manifest.json`이다.
앞 절의 "인증 세션이 없어 0개"는 재개 전 기록이며 현재 상태는 이 실행 파일로 확인한다.
인증정보는 실행 중인 프로세스의 숨김 입력으로만 전달했고 파일에 저장하지 않았다.

### 2026-10-14 대조장 798장 다운로드 완료 확인

기존 downloader가 798/798장 완료 후 정상 종료했다. 2026-10-14 17:23 KST에
승인 입력 SHA-256, 입력/상태/출력 manifest의 전체 고유 ID, 실제 파일 존재·크기,
control 역할·대상 ID·좌표·각척도 보존을 대조했다. 누락·크기 불일치·역할 불일치 0건.
실제 합계는 45,254,626,560 bytes(42.1467 GiB)다. 다운로드 당시 FITS 구조·WCS
검증 기록은 798건이며 헤더 필터는 i 312개/i2 486개다.

완료 자료는 `sample_800_2026-10-13/control_download_01`의 `images/`,
`download_manifest.json`, `download_status.json`이다. 확인 보고서는
`confirmation_20261014T082353Z.json`. 이 확인에서 45 GB 전체 해시를 다시 읽지는
않았으며 다운로드 중 기록한 SHA-256을 보존한다. 전체 픽셀·품질 마스크·유효면적·
대조장 배경/공분산과 실제 passband에 따른 소광 검증은 다음 단계다.
science_measurement=false, credentials_saved=false 유지.

### 2026-10-14 대조장 전체 품질검사와 대응 외곽 배경 비교

승인된 다운로드의 후속 분석을 lageunha CPU에서 완료했다. 추가 다운로드·GPU
사용은 없다. 기존 `hsc_quality` 픽셀 검사와 SEP/외곽 타일 함수를 재사용하고,
대조장 중심에는 BCG 보존 예외를 두지 않는 검증 스크립트를 추가했다.
상세 방법·수치·실행 명령은 `analysis/icl_wgl_membership/hsc_debug/control_images_README.md`다.

`sample_800_2026-10-13/control_quality_2026-10-14_01`:
798/798장의 45,254,626,560 bytes 전체 SHA-256을 다시 읽어 다운로드 기록과
대조했고 구조·WCS 오류 0, 비유한 science 픽셀을 가진 영상 0이다.
26장의 일부 분산 픽셀 1,205개를 제외했다. random 중심 BRIGHT_OBJECT는 0개다.
전체 품질 유효면적 중앙값 47.44%, 밝은별 마스크 면적 중앙값 42.13%다.
품질 마스크만 적용한 모든 0–500 kpc 고리의 50% 유효면적 문턱을 통과한
동일 대조장 표본은 234개이며, 고리별 수는 716개에서 338개까지 달라진다.
이는 천체 확장 마스크 적용 전 기술 진단이며 최종 과학 표본 선정이 아니다.

`control_outer_2026-10-14_01`:
기존 면적 집합 40개·마스크 스트레스 10개의 은하단과 대응 대조장 50쌍을
600–900 kpc에서 비교했다. Kron 2.5/3.5/5.0의 150개 방법 행,
대조장 타일 25,229개, 간격/방향별 상관 1,200행을 완료했다. 50쌍 모두
헤더 i/i2가 일치하고 각 방법에서 양쪽 유효 타일 10개 이상이다.
Kron 3.5의 타일 잔차 MAD/명목 평균 오차 중앙값은 면적 집합 은하단/대조장
2.577/2.687, 스트레스 집합 2.854/2.826이다. 유효면적 중앙값은 각각
37.31%/34.18%, 35.68%/28.03%다. 이는 순수 잡음 배율이나 ICL 신호가 아니다.
마스크 크기별 살아남는 타일도 달라지므로 같은 유효 영역에서 후속 비교가 필요하다.

`control_image_audit_2026-10-14_01/audit.json`:
입력 해시 108건·코드 해시 5건, 반경별 면적 1,000건, 마스크 150개,
WCS 표본점 600개, 원본 타일 450건, 쌍별 광량 환산 150건 검산이 통과했다.
타일 평균 최대 차이 0. 새 테스트 4개와 기존 회귀 7개 총 11개도 통과했다.
비교 그림은 `control_outer_2026-10-14_01/paired_outer_diagnostics.png`다.

다음은 이 50쌍에서 같은 유효 영역의 마스크 민감도, 실제 대조장 신호 주입·회수,
배경 공분산을 검증하는 단계다. PSF/소광/유효 passband 검증도 남아 있다.
원본·기존 결과·활성 science manifest는 보존하며 science_measurement=false,
controls_ready_for_science=false, sky_correction_approved=false를 유지한다.

### 2026-10-15 대조장 비교 확대 재검증

사용자의 재검증 요청에 따라 lageunha CPU에서 기존 테스트 11개와 원본 감사
전체를 다시 실행해 통과했다. 감사 결과는
`sample_800_2026-10-13/control_image_audit_2026-10-15_01/audit.json`이다.
같은 폴더의 `summary_recheck.json`에서 품질 집계와 그룹별 통계 60개도 확인했다.
234개 유효면적 조건 통과, 26장의 분산 문제 픽셀 1,205개 등 보고 수치는 그대로다.

기존 감사의 일부 타일 표본 검사에 더해, 별도 코드
`analysis/icl_wgl_membership/hsc_debug/recheck_control_statistics.py`로
대조장 798장 전체 45,254,626,560 bytes 해시와 양쪽 영상 50쌍의 모든 통계를
재계산했다. 기존 타일 통계 함수를 호출하지 않고 가중 정규방정식과 영상 전체의
평균 차감 픽셀쌍 공분산을 사용했다. 타일 경계 밖의 쌍은 명시적으로 제외한다.

`control_full_recheck_2026-10-15_01/verification.json`:
타일 52,135개·평면 적합 300건·상관 2,400건이 통과했다.
타일 평균·중앙값·명목 평균 오차 차이는 0이고, 잔차/명목 오차 비율의 최대
절대 차이는 1.42e-14, 상관계수는 1.11e-16으로 반올림 수준이다.
별도 계산식 테스트 4개를 포함해 총 15개 테스트가 통과했다.
원본·기존 결과·과학용 활성 상태는 유지했고 추가 다운로드는 없다.
이번 단계는 수치 재현성 검증이며 신호 주입·회수와 과학용 공분산 검증은 다음 단계다.

### 2026-10-15 실제 영상 주입·회수와 스태킹 전 공분산 검증

사용자 지시에 따라 lageunha CPU에서 대조장 798개 전체를 재검출하고
Kron 2.5/3.5/5.0 마스크, 50–500 kpc 여섯 구간 프로필을 추출했다.
기존 50개 검증 대조장에 원형 지수 h=100 kpc와 타원 h=250 kpc/q=0.65,
100 kpc의 관측 SB 1/3/10 nJy/arcsec²를 주입해 총 300장을 재검출했다.
이 단계는 코애드 주입이며 추가 photon noise·HSC 이전 처리·PSF 날개 전달을
검증한 것은 아니다. 기존 은하단 50개의 마스크 기반 프로필도 별도 추출했다.
추출/주입 시간 약 1,254초, 오류 0, 회수 행 10,800개, 집계 432개다.

기준/방법/결과는 `analysis/icl_wgl_membership/hsc_debug/prestack_validation_README.md`다.
Kron 3.5에서 면적 조건 집합 40개의 50–300 kpc는 모형별 유효 13–25개로,
원본 고리 추출 24개 집계의 회수 편향 중앙값이 −3.31%~+3.39%였다.
300–500 kpc의 유효 영상은 6–8개라 미승인이다. 마스크 스트레스 집합은
구간별 0–4개여서 전체가 표본 부족이다. 전체 고리 참값과 열린 면적 참값을
구분해 보조 회수비도 저장했으며 원래 기하학적 고리 기준 판정은 유지했다.

외곽 평면을 빼면 넓은 모형의 200–300 kpc에서 A=1/3 주입에 각각
−16.44%/−10.02% 편향이 생겼다. 평면 차감은 민감도 시험으로만 남기고
기존 원본 고리/대조장 기반 분석에 새 과학용 보정으로 편입하지 않았다.

여섯 구간 모두 면적 50%를 충족한 대조장은 Kron 2.5/3.5/5.0에서 140/40/0개다.
Kron 3.5 대조장 공분산은 40개·28개 5° 하늘 셀, rank=6, 조건수 약 84.1이다.
독립 bootstrap 대각은 S/N 대비 최대 4.72% 차이, 공간 hold-out의 1.96σ 포함률
80%로 정의한 탐색적 기술 기준을 통과했다. 공간 bootstrap 표준오차는 독립
영상 가정보다 약 7–17% 크다. 원점/격자 크기·필터·각척도 구간 민감도도 저장했다.
이는 대조장/잔차 성분의 조건부 공분산이며 순수 잡음 배율이나 전체 과학 오차가 아니다.

기존 은하단–대조장 50쌍 중 양쪽 모두 여섯 구간을 지원하는 것은 Kron 3.5에서
0쌍이다. 구간별로도 8/6/8/7/4/2쌍이다. Kron 2.5의 공통 3쌍은 rank=2여서
6차원 전체 스태킹 공분산에 부족하다. 대조장 공분산을 전체 오차로 대체하지 않는다.
**원본 추출의 지원 집합에서 50–300 kpc 회수는 확인했지만 전체 과학 스태킹은 보류**다.
다음은 양쪽에 공통 면적을 갖는 짝 표본과 마스크/반경 정의 검증이다.

실행 폴더는 `sample_800_2026-10-13` 아래
`prestack_extraction_2026-10-15_01`, `prestack_clusters_2026-10-15_01`,
`prestack_validation_2026-10-15_01`, `prestack_audit_2026-10-15_01`이다.
검토 그림은 `prestack_validation_2026-10-15_01/prestack_validation_review.png`다.

새 테스트 6개(독립 잡음 500회/공간 공유 오프셋 포함)가 통과했다.
독립 감사도 입력·산출물 해시 1,896건, 기존 마스크 150개, 픽셀 프로필 값 1,512개,
공통 면적 증분 1,296개·고정 면적 모형 1,296개, 공분산 84개 설정을 통과했다.
원본 재계산 최대 차이는 2.51e-13 nJy/arcsec² 미만이다.
감사 폴더에 해시 이름의 코드 사본을 보존했다. 추가 다운로드·GPU 사용은 없고,
원본·기존 결과·science_measurement=false·sky_correction_approved=false를 유지했다.

### 2026-10-16 matched pair 공통 면적 민감도

기존 은하단·대조장 50쌍의 저장된 프로필에서 Kron 2.5/3.5/5.0, 면적 바닥
25/35/40/50%, 50–150/200/300/500 kpc와 개별 고리를 조합해 120개 공분산 지원
설정을 계산했다. 신호 주입 회수 민감도 72개 그룹도 재집계했다. 코드는
`analysis/icl_wgl_membership/hsc_debug/evaluate_prestack_pair_support.py`이고,
상세 설명은 `analysis/icl_wgl_membership/hsc_debug/prestack_validation_README.md`다.

결과는 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_pair_support_2026-10-16_02/`다.
원래 Kron 3.5·50% 기준에서는 양쪽이 50–300 kpc 네 구간을 모두 만족하는 쌍이
1개, 50–500 kpc는 0개였다. 탐색 조건 Kron 2.5·40%는 50–300 kpc 24쌍
(모두 area-qualified), 21개 공동 5° 블록, rank 4를 남겼다. bootstrap 대각
차이 6.8% 이내, 공간 hold-out 포함률 100%로 기술 공분산 진단은 통과했다.
공동 블록 표준오차는 [11.30, 2.80, 1.34, 0.892] nJy/arcsec²이며,
전체 800개 표본의 오차나 순수 픽셀 잡음은 아니다.

같은 Kron 2.5·40% 주입 민감도에서 각 6개 모형이 28개 필드의 네 구간을 모두
지원했고 중앙 회수 편향 절대값은 최대 5.90%다. 이것은 coadd smooth-signal
주입의 회수이며 실제 미검출/잔류 천체 외곽빛은 검증하지 않았다. Kron 2.5/40%
조건은 표본 수를 확보하는 후보일 뿐 과학용 측정 정의로 채택하지 않았다.
기존 50개는 전체 800개 대표 표본이 아니고 PSF·상류 sky·i/i2 소광/응답·위성–ICL
광량 누출이 미승인이다. `science_measurement=false`,
`sky_correction_approved=false`를 유지한다. 전체 matched sample의 잔류 천체광과
면적 민감도 평가 후 본 분석 반경/마스크를 고정한다. 다운로드·GPU 사용 없음.

### 2026-10-16 전체 표본 자동 마스크·paired 면적 지원 계산 완료

전체 `download_01/download_manifest.json` 800개와 798개 matched control 연결을
확인했다. 기존 BCG 검토 마스크는 50개에 한정되어 있어 새 스크립트
`analysis/icl_wgl_membership/hsc_debug/prepare_full_cluster_mask_support.py`로
800개 source-mask 후보와 50–500 kpc 구간별 면적을 계산한다. 이 작업은 밝기
프로필을 만들지 않고 자동 중심 검출 후보를 BCG 확정으로 취급하지 않는다.

Slurm job 424800은 800/800 성공·오류 0·46:03 경과로 완료했다. 자동 중심 후보
258개, 중심 품질 마스크 391개, Kron 2.5/3.5/5.0에서 가려진 중심 671/696/725개다.
후보 중 품질/천체 마스크가 중심에 걸리지 않은 수는 129/104/75개다. 이는 자동
위치 후보이며 BCG 검증이나 측광용 표본 승인이 아니다. cluster-only 50–500 kpc
지원 수는 Kron 2.5에서 면적 25/35/40/50% 기준 331/186/131/40개,
Kron 3.5에서 171/50/27/4개다. paired 면적 집계 job 424865도 완료했다.
전체 798개 matched control 중 2개 대상은 대조장이 없어 분석 연결에서 빠졌다.
탐색 조건 Kron 2.5·40%, 50–300 kpc는 공통 면적 138쌍(동일 필터 쌍 137; 대상
필터 i 44/i2 94), 공동 5° 공간 블록 44개로 표본 수 스크린을 통과했다.
같은 마스크/바닥에서 50–500 kpc는 69쌍(동일 필터 68), 블록 37개로 수 스크린을 통과했지만,
그 반경의 주입 회수와 광도 공분산 검증은 완료되지 않았다. Kron 2.5·50%의
50–300 kpc는 41쌍/30블록으로 스크린을 통과하나 50–500 kpc는 11쌍/10블록으로
실패했다. 기준 Kron 3.5·50%는 50–300 kpc 4쌍, 50–500 kpc 1쌍뿐이다.
이 표는 면적 지원만 계산하며 flux covariance, 광도 프로필, 영상 스택을 만들지
않았다. 자동 중심 후보도 BCG 판정이 아니다. 모든 조건에서
`flux_covariance_available=false`, `science_measurement=false`다.
결과는 `prestack_full_pair_area_support_2026-10-16_01`에 있다. 작업 스크립트는
`run_full_cluster_mask_support.sh`, `run_full_pair_area_support.sh`이며 로그는
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_mask_support_logs_2026-10-16_01/`다.
출력은 `prestack_cluster_masks_2026-10-16_01`과
`prestack_full_pair_area_support_2026-10-16_01`이다. mask summary에 입력/출력 해시
800개를 보존했다. 두 작업 모두 GPU 코드를 사용하지 않고 과학용 마스크/측정을
승인하지 않는다.


### 2026-10-20 후속: 전체 Gaia PSF 후보 진단과 weak-lensing 좌표 반출 상태

기존 Gaia DR3 별 근거 후보 480개에 대해 `vari_summary`와
`vari_classifier_result`를 조회했다. 두 표에서 일치 source 행은 0개였으며,
이는 비변광 판정이 아니라 Gaia DR3 변광 처리 기록이 없다는 뜻이다. 전체
HSC 진단은 288개 프레임, 333개 401×401 우표, 1,920개 방법 조합을 기록했다.
137개 후보가 적어도 한 조합의 유효 프로필을 냈고 42개는 네 조합 모두에서
중심/배경 분석 조건을 통과했다. 64–96픽셀 외곽 고리에서는 기하 지원을 통과한
방법 프로필 274개 중 118개 양수·156개 0 이하였고, 후보별 묶음 113개의 지원
설정 부호는 전체 양수 38, 전체 0 이하 58, 혼합 17이었다. 외곽 PSF 날개는
검출/승인되지 않았다. Gaia 변광 분석 자료 부재, 정확한 epoch별 Gaia↔HSC 후보
식별, 영상별 포화·고립도 및 프레임 전달 검증이 남아 `science_measurement=false`,
`psf_approved=false`다.

Gaia DR3 `ref_epoch=2016.0` proper motion을 구면 접평면에서 ±10년까지 전파하는
민감도 분석도 수행했다. HSC 후보 좌표를 고정한 시나리오에서 480개 모두 기존
2″ 연결반경 안에 남았고 최대 분리는 1.572″였다. HSC coadd의 실제 epoch를
대입하지 않은 시나리오이므로 실제 전파나 천체/PSF 검증이 아니다. 결과와
재현기는 `empirical_psf_gaia_epoch_sensitivity_2026-10-20_01/` 및
`summarize_psf_gaia_epoch_sensitivity.py`다.

재현 스크립트와 수치는 `analysis/icl_wgl_membership/hsc_debug/`의
`diagnose_empirical_psf.py`, `query_psf_gaia_variability.py`,
`summarize_empirical_psf_star_level.py`, `extended_psf_README.md`, 통합 연구평가에
있다. 원자료 기반 진단과 Gaia 조회 산출물은 Git 밖의 `/gpfs/kjhan/Astrafek/data/`
아래에 보존한다.

Weak-lensing 대조장용 좌표 export 프로세스 PID 1049265는 2026-10-20 확인 시
`lageunha`에서 살아 있었다. checkpoint에서 SQL batch 000과 001은 완료,
002는 CAS job 1062312에 제출된 상태다. 프로세스를 재시작하지 않는다. Slurm
제어기 DNS SRV 조회가 실패해 `squeue`와 새 `sbatch` 작업을 사용할 수 없다.

### 2026-10-20 후속: 승인 shape 행과 렌즈 단계 재확인

현재 GPFS의 authoritative shape-row run은 실패한 `_02` 시도가 아니라
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_shape_rows_export_2026-10-18_03/`이다.
체크포인트와 요약을 다시 읽어 32/32 batch, 3,069,544 lens–source pair,
909,833,598 bytes, 22개 승인 열, `credentials_saved=false`, 전체 photo-z PDF
미포함, `science_measurement=false`를 확인했다. 각 배치의 SHA-256과 렌즈별
count 대조는 체크포인트 및 `row_summary.json`에 기록되어 있다. 재조회나 추가
반출은 하지 않았다.

이 승인 행과 로컬 HSC Y3 SACC `n(z)`, 좌표로 검증한 raw B-mode mask 극성을 사용한
cluster-only 약한렌즈 스택은 2,153,004 lens–source pair와 609개 렌즈를 포함한다.
52-cell 연결 sky jackknife 및 source-bin 일관성 진단이 Astrafek 스크립트로
저장되어 있다. 0.1–0.7 Mpc source-bin 간 차이가 남고 random 중심 차감 및 대응
random-field 공분산은 아직 없으므로 해당 스택 역시 `science_measurement=false`다.
수치와 해시는
`analysis/icl_wgl_membership/hsc_debug/wl_tomographic_consistency_assessment_2026-10-20.md`
및 `wl_source_reuse_covariance_validation_2026-10-20_02` 검증 산출물을 참조한다.

random-control 좌표 내보내기 PID 1049265는 이번 확인에서도 `lageunha`에서 살아
있었다(상태 `hrtimer_nanosleep`). 좌표 체크포인트는 batch 000·001 완료
(33,514/780,191행), batch 002는 CAS job 1062312 제출 상태로 그대로다. 이 CAS
job을 재제출하거나 살아 있는 프로세스를 재시작하지 않는다. 현재 output은 좌표
내보내기뿐이며 random source shape 행은 포함하지 않는다. 다음 실행 단계는
기존 좌표 체크포인트가 32/32가 된 뒤 행 수와 해시를 검증하고 전체 CAMIRA 조리개
중첩·중복좌표·렌즈별 지원을 다시 screening하는 것이다. 이후 random shape count와
shape-row 반출은 각각 실제 필터링된 중심 수·행 수·전송량을 확정해 별도 승인을
받은 뒤 진행한다.

2026-10-20 21:44 KST의 read-only SSH 재확인에서 PID 1049265는 34시간 37분째
살아 있었고 `hrtimer_nanosleep` 상태였다. 표준입출력은 원격 `/dev/pts/0`에 연결돼
있으며, 실행 중 메모리에 올라간 구버전 submitter는 CAS 상태를 확인한 뒤 60초
주기로 대기한다. 저장된 `coordinate_status.json`은 여전히 2/32 배치이며 최상위
`HTTPError` marker가 남아 있다. 성공 poll 시각을 저장하지 않는 구버전이므로
서버 측 작업 상태는 체크포인트로 판별할 수 없다. 프로세스와 제출 배치를 유지하고
중복 제출하지 않는다.

### 2026-10-20 후속: 전체 Gaia PSF 스택 및 승인 DESI 적색이동 보정

확장 Gaia ledger의 480 후보 중 hash-verified stamp 333개(227 프레임)를
기존 800개 HSC-I 컷아웃으로부터 스택했다. 이 표본에서도 외곽 64–96픽셀은
배경 annulus/source-mask 선택에 따라 부호가 달라져 확장 PSF 날개를 측정하지
못했다. 표본 수 증가로 두 설정은 0과 양립하고 한 설정은 양의 평균을 보였으나,
나머지 두 설정은 음수다. 국소 sky와 PSF 날개를 분리할 수 없어 `psf_approved=false`,
`science_measurement=false`다. 권위 산출물은
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/empirical_psf_frame_weighted_full_gaia_2026-10-20_02/`
이며 summary SHA-256은
`0be06f1a827823bcd04ce23db4113dca1bd368e1a2b7205d0bab0ed76911a7ba`다.

사용자가 DESI DR1/DR2 HSC PDR3 공개 보정 자료 4개만 다운로드하도록 승인했다.
고정 commit `abc70af29910e49fce17a63ef3083160d515fb6c`의 두 scale-cut NPZ와
JSON sidecar만 수신했고 총 549,818 bytes, 파일별 URL/크기/SHA-256은
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/desi_y3_hsc_nz_abc70af2_2026-10-20/manifest.json`에
있다. 4개 bin의 1,000점 `n(z)` mean/median와 점별 lower/upper, magnification
corrected 곡선을 확인했다. joint posterior 샘플은 포함하지 않는다. 계산 재현기는
`audit_wl_desi_nz_distribution_sensitivity.py`; SACC 기준과 두 scale-cut mean
곡선의 동일 shape rows 및 5도 sky jackknife, 점별 envelope stack을 포함한다.
random-lens shape rows가 없어 아직 cluster-only 진단이다. 결합 PSF/DESI Slurm job
430248은 69분 계산 후 pointwise envelope 행에 없는 jackknife-error 열을 요구하는
출력 조합 버그로 종료됐다. 승인된 기존 SACC 기준 summary가 3,069,544행 계획,
각 배치 해시, SACC 해시, 5도 52셀, 반경 구간, 우주론과 정확히 일치함을 확인해
재사용하도록 고쳤다. 반복 렌즈 적색편이의 `Sigma_crit^{-1}` 적분도 lens-z/bin별
한 번만 계산하도록 최적화했다. 회귀검사 356개와 subtest 6개가 통과했다. DESI만
계산하는 Slurm job 430481은 37분 41초 만에 완료했고 독립 validator도 통과했다.
기존 SACC stack과 비교해 DESI central curve는 7개 반경 median `DeltaSigma`를
source bin 1에서 +12.5–14.8%, bins 2–4에서 약 −1–3.5% 바꿨다. 네 central
model의 52셀/56×56 covariance는 rank 51이고 상대 비대칭은 4.1e-17 이하이다.
0.1–0.7 Mpc의 내부 토모그래피 불일치는 낮아지지만 계속 남았다(세 radial raw
4-bin p 범위 약 4.4e-8–1.2e-3). random lens subtraction이 없으므로 모두
cluster-only diagnostic이며 `science_measurement=false`다. 정량표와 해시는
[`wl_desi_nz_distribution_sensitivity_assessment_2026-10-20.md`](analysis/icl_wgl_membership/hsc_debug/wl_desi_nz_distribution_sensitivity_assessment_2026-10-20.md)
및 GPFS의 `spatial_jackknife_2026-10-20_01/summary.json`, `validation.json`에 있다.

### 2026-10-20 후속: 초기 20개 별 우표 PSF 파일럿 (확장 표본으로 대체)

초기 파일럿에서 승인된 기존 HSC 별 우표 20개(15개 프레임)를 Astrafek 재현 스크립트
`analysis/icl_wgl_membership/hsc_debug/stack_empirical_psf_star_profiles.py`로
프로필 스택했다. Kron 2.5/3.5와 120–150/160–190픽셀 배경 설정 모두 48–96픽셀
외곽 net profile이 음수였다. 따라서 현 국소 sky estimator와 별 표본으로는
양의 확장 PSF 날개를 측정하거나 보정할 수 없다. 이는 날개가 없다는 검출이
아니며 배경·ICL 혼입과 epoch 문제도 남는다. 결과는
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/empirical_psf_frame_weighted_stack_2026-10-20_02/`
및 `analysis/icl_wgl_membership/hsc_debug/empirical_psf_frame_weighted_stack_assessment_2026-10-20.md`에
해시와 함께 기록했다. 이 파일럿은 333개 stamp/227개 프레임의 확장 표본 평가에
의해 해석상 대체됐다. 확장 결과와 해시는 위 절에 기록했다. 두 결과 모두
`science_measurement=false`, `psf_approved=false` 상태를 유지한다.

### 2026-10-20 후속: Gaia PSF 별 우표의 matched-control WCS 지원

확장 Gaia ledger의 333개 별 우표(227개 표적 프레임) 각각에 대해 표적 은하단 중심
기준 east/north 각거리 오프셋을 기존 매칭 대조장 중심에 옮겨 WCS 대응을 확인했다.
319개 후보(217개 표적 프레임)는 같은 하위 필터 `i` 또는 `i2`였고 401×401 전체
control 우표가 영상 안에 들어왔다. 14개 후보(10개 프레임)는 target/control 하위
필터가 달라 제외했다. 검사기는 영상 헤더와 기존 star stamp만 읽고 coadd science
픽셀·새 HSC 자료는 읽지 않았다. Astrafek 재현기
`audit_empirical_psf_matched_control_support.py`와 테스트
`test_audit_empirical_psf_matched_control_support.py` 3개가 있다. 산출물은
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/empirical_psf_matched_control_support_2026-10-20_01/`에
있으며 summary SHA-256은 `a6d80adbfe8d01e466930dc68d44355287124b1b5ad8d8221a0c025f4cd7346f`,
행 표 SHA-256은 `fc255bf0704ddfd43e817e6b21e4a66f28c593e6748fda50976be12558f5da15`다.
이는 다음 대조장 잔차 profile의 입력 지원만 확인한다. PSF wing 보정 및
`science_measurement`는 승인되지 않았다. 상세는
`analysis/icl_wgl_membership/hsc_debug/empirical_psf_matched_control_support_assessment_2026-10-20.md`를
참조한다.

후속 paired profile은 기존 표적·대조 영상만 사용해 같은 필터 319개 별 위치(217 프레임)를 세 mask/coverage 설정으로 평가했다. outer-ring 중심 mask 조건을 조정한 뒤에도 최대 지원은 i2 64–96픽셀 구간의 12 프레임이며, sky annulus에 민감하다. control 위치는 별 없는 placebo이고 target에는 ICL 곡률이 남아 PSF wing을 분리할 수 없다. `_03` summary SHA-256은 `7478c905a4e46c0c82f3fe569e4ee1b4cfb53280c9acb635267bc4188a98cd5b`다. 방법·상태별 제외·모든 실행 해시는 `analysis/icl_wgl_membership/hsc_debug/empirical_psf_matched_control_profiles_assessment_2026-10-20.md`에 있다. `psf_approved=false`, `science_measurement=false`.

2026-10-20 22:04 KST read-only 원격 재확인: lageunha PID 1049265는 `Ss+`로 살아 있으며 누적 경과는 34시간 55분이다. 승인된 `wl_random_point_coordinate_plan_2026-10-19_01` 체크포인트는 `running`, 32개 중 2개 완료(33,514/780,191 예상 좌표), batch `002`는 CAS job `1062312` 제출 상태다. 기존 `HTTPError` 표시는 여전히 요청/시각이 특정되지 않아 terminal 실패로 판단하지 않았다. 재제출·프로세스 재시작은 하지 않았다. 이어질 전체 기하 스크리너는 승인된 좌표 32/32와 계획 해시 검증 뒤에만 실행한다.

전체 좌표 기하 스크리너 회귀검사를 보강해 target과 CAMIRA 이웃의 각각의 적색편이에서 계산한 proper-3-Mpc 각반경 합으로 교차/비교차 지점을 검증하고, chunk 크기 1과 10에서 결과 동일성을 확인했다. 좌표 체크포인트 감사기와 함께 focused test 7개가 통과했으며, 실제 전체 좌표 자료가 도착하기 전에는 기하 스크리닝을 실행하지 않는다.

Gaia↔HSC 후보 연결의 epoch 범위를 점검하려고 기존 333개 별 우표의 227개 target coadd FITS header를 조사했다. `DATE-OBS`/`DATE-*`/`MJD-*` 키는 0/227이다. 단, PSF 우표 중심은 Gaia가 아닌 HSC 검출의 `X_IMAGE/Y_IMAGE`이므로 이 누락이 영상 중심정렬을 막지 않는다. ±2/5/10년 Gaia 이동 결과는 실제 epoch 확인이 아닌 cross-match 민감도이며, 480개 모두 2초각 내에 남았다. 명확화한 Astrafek header-only 검사기는 `analysis/icl_wgl_membership/hsc_debug/audit_empirical_psf_hsc_coadd_epochs.py`, 테스트 3개, 상세·체크섬은 `analysis/icl_wgl_membership/hsc_debug/empirical_psf_hsc_coadd_epoch_audit_2026-10-20.md`에 있다. `_02` 요약 SHA-256은 `fb5ff0c7b87d1d17bdeb42ef824980b3a50eb7c7f0a487d7449b1fd9fe88d93c`이며 픽셀과 HSC 서비스는 읽지 않았다. `psf_approved=false`, `science_measurement=false`를 유지한다.

같은 333개 별 우표에 대해 HSC 포화/`BRIGHT_OBJECT` 비트, 검출 이웃, Kron source-mask 지지율을 추가 점검했다. 중심 12픽셀 내 포화는 0/333, 밝은천체 비트는 3/333이고, 12픽셀 이웃 검출은 없다. 중심 source-mask 개방률 99% 기준을 통과한 수는 Kron 2.5에서 290개/201프레임, Kron 3.5에서 242개/179프레임이다. 하지만 48픽셀 내 이웃은 317개, 96픽셀 내 이웃은 333개다. 중심 적격 subset에서도 48–64 및 64–96픽셀 개방률 중앙값은 Kron 2.5에서 0.718/0.694, Kron 3.5에서 0.605/0.579로 외곽 날개 측정에 강한 마스크·혼잡도 통제가 필요하다. 검사기 `analysis/icl_wgl_membership/hsc_debug/audit_empirical_psf_star_stamp_quality.py`, 관련 테스트 3개 통과. 결과는 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/empirical_psf_star_stamp_quality_audit_2026-10-20_02/`에 있고 summary SHA-256은 `bdb176c6ede927369e0c9b37bbc9ec88b0e1ba607bce89dfbab18795fe925bc0`, 행 표 SHA-256은 `ca2b908c717358ca0c3bb581da87ad717af6206b40ef0c6333d462a80ed48187`이다. 중심부 후보 진단이며 `psf_approved=false`, `science_measurement=false`를 유지한다. 해석과 한계는 `analysis/icl_wgl_membership/hsc_debug/empirical_psf_star_stamp_quality_assessment_2026-10-20.md`에 있다.

### 2026-10-20 후속: weak-lensing random 후보 count 전량 감사

이미 승인·수신된 count-only 표의 32개 CSV 해시와 800개 렌즈별 row를 좌표 SQL 계획에 대조했다. 30 arcmin 외곽 및 각 렌즈의 두 배 proper-3-Mpc 제외반경 안에서 HSC quality-filtered random 표본은 총 780,191개였다. 렌즈별 count는 364–1,443(중앙값 1,000)이고 전 렌즈가 계획한 200개를 넘는다. 유효 고리 면적으로 나눈 점밀도 중앙값은 0.566 arcmin⁻²다. raw count와 렌즈 적색편이의 고정표본 Spearman `rho=0.595`는 면적 정규화 뒤 `rho=0.080`으로 줄어, count의 z 의존 대부분이 각도 제외면적과 연결됨을 시사한다. 이는 중첩 제거 전 후보 가용성으로 shape 지원이나 최종 random pool을 보장하지 않는다. 전체 결과·입력 해시는 `analysis/icl_wgl_membership/hsc_debug/wl_random_coordinate_candidate_count_assessment_2026-10-20.md` 및 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_point_candidate_count_audit_2026-10-20_02/`에 있다. 검사기 `audit_wl_random_coordinate_candidate_counts.py`, 테스트 2개 통과. CAS 추가 요청·좌표 반출은 없었다.

현재 read-only lageunha 확인에서 PID 1049265는 여전히 살아 있음(`Ss+`, 경과 1일 11시간 32분)이며 좌표 체크포인트는 `running`, 완료 2/32, 33,514/780,191행, batch `random_coordinates_002` CAS job `1062312` 제출 상태다. 최상위 `HTTPError` 표시는 남아 있고 checkpoint mtime은 바뀌지 않았다. 서버 job의 terminal 상태를 확인하지 못했으므로 살아 있는 프로세스를 재시작하거나 batch를 재제출하지 않는다. 남은 전체 기하 스크리닝은 기존 export가 끝나고 32개 batch의 row count·SHA 검증을 통과한 뒤 수행한다.

### 2026-10-20 후속: 경험 PSF radial coverage와 ICL annulus

기존 333개 Gaia 별 우표/227 프레임의 표적 적색편이와 HSC 픽셀 scale을 사용해 경험 PSF 스택의 radial extent를 proper kpc로 환산했다. 실제 stack은 96픽셀까지이며 물리반경 중앙값 88.4 kpc (q10–q90 76.8–101.4 kpc)다. ICL 분석의 200–300 kpc annulus에는 227개 frame 중 하나도 닿지 않는다. 전체 우표의 circular-inscribed radius도 중앙값 184.1 kpc여서 200 kpc까지 닿는 프레임은 52개, 300 kpc는 0개다. 반면 96픽셀 endpoint는 75 kpc 이상을 덮는 stamp 316/333개로 satellite 외곽 분석과 일부 기하 중첩은 있다. 이것은 측정 PSF가 아니라 coverage audit다. 재현기 `analysis/icl_wgl_membership/hsc_debug/audit_empirical_psf_physical_radius_coverage.py`, 테스트 2개 통과. output `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/empirical_psf_physical_radius_audit_2026-10-20_01/`; summary SHA-256 `6d88022648930be6fd6f5f6785d66271cb399b1f6c4107dd7305453fead5d872`, table SHA-256 `b1b4cfee369a20939dbf3b6b040ff4a08b7efb4239e155829f65fc1bd99a1f32`. 상세 범위와 입력 provenance는 `analysis/icl_wgl_membership/hsc_debug/empirical_psf_physical_radius_coverage_assessment_2026-10-20.md`에 있다. PSF 보정과 `science_measurement`는 승인되지 않았다.

### 2026-10-20 후속: full-frame outer-field PSF-wing 파일럿

기존 333개 hash-verified Gaia 별 우표 중 18개가 target/control 양쪽에서 300 kpc 기하 지원, 동일 하위필터, 800 kpc 외곽장, 중심 Kron 2.5 개방률 조건을 통과했다. 이미 승인·다운로드된 800 target 및 798 control HSC 영상만 읽었다. Slurm job 430734는 18쌍의 25 kpc profile을 처리해 중복 없는 864행과 24개 frame-equal paired stack bin을 만들었다.

독립 검증에서 profile arithmetic·SHA-256이 일치했다. 하지만 radial bin별 paired support는 최대 4쌍이며 250–300 kpc에는 없고 200–225 kpc는 1쌍이다. `_01`에는 control profile 중복 버그가 있었으므로 사용하지 않는다. 수정된 `_02`도 PSF 보정 또는 ICL 과학 측정으로 쓸 수 없다: `science_measurement=false`, `psf_approved=false`. 핵심 산출물은 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/empirical_psf_outer_field_wings_2026-10-20_02/`이며 summary SHA-256 `c8085c8a7459a34ab03e551de6d2d37cccf637a80156bff4f70ba61b0dc46369`, validation SHA-256 `c70a144aae21922044d14a230404c8633a7d9be755580b976d7752a0fad860f6`이다. 재현기·실패 원인과 후속 기준은 `analysis/icl_wgl_membership/hsc_debug/empirical_psf_outer_field_wings_assessment_2026-10-20.md`에 있다.

### 2026-10-20 후속: random coordinate radial-yield

이미 승인된 첫 2/32 coordinate batch(50개 계획 순서 렌즈, 33,514행)만 사용해 inner 2×proper-3-Mpc 경계에서 거리에 따른 geometry-clean 후보수를 나눴다. 30′ outer limit은 각 렌즈에 7.14–8.76′(중앙값 8.10′)의 annulus 폭만 남긴다. 0–5′에서 1,199/19,920, 남은 shell에서 930/13,594가 geometry-clean이며 총 2,129(6.35%)다. 기존 50개는 비대표 표본이다. 30–60′ 추가 고리의 이상적 spherical area는 현재 검색 고리의 6.02–7.15배(중앙값 6.42배)이지만 HSC coverage·depth를 반영한 것은 아니다. 보충 count-only CAS 계획은 아직 제출하지 않았다. output `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_coordinate_radial_yield_audit_2026-10-20_03/`, summary SHA-256 `cd5a30122a2beb2ad875dec3f59bbca35286ee0862cbc1473d021cbb240402ad`. Astrafek 재현기와 해석은 `analysis/icl_wgl_membership/hsc_debug/wl_random_coordinate_radial_yield_assessment_2026-10-20.md`에 있다. `science_ready=false`, `science_measurement=false`.

### 2026-10-20 후속: random-coordinate export live audit

이번 read-only lageunha 확인에서 PID 1049265는 살아 있었다(`Ss+`, elapsed 약 1일 11시간 59분). 기존 체크포인트는 `running`; `random_coordinates_000`·`001` complete, `002`는 CAS job 1062312 submitted이며 top-level `HTTPError` marker가 남아 있다. `audit_wl_random_coordinate_checkpoint.py`로 저장된 CSV 2개를 count plan별 33,514/33,514행 및 SHA-256과 대조해 검증했다. 새 CAS 요청·중복 제출은 하지 않았고, 해당 batch의 실제 CAS terminal state는 여전히 모른다. checkpoint audit `_11`은 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_coordinate_checkpoint_audit_2026-10-20_11/`에 기록했다.

### 2026-10-20 후속: ICL×위성 반경 교차 지원

기존 800개 HSC full-sample 프로필의 여섯 ICL 고리(50–500 kpc)와 일곱 위성 고리(5–75 kpc)의 exact cluster/filter 교집합 84개를 작성·해시 검증했다. 가장 큰 지원은 ICL 100–150 kpc × 위성 5–10 kpc로 i 78개, i2 150개지만 위성 내부 구조를 뜻한다. 관심 외곽 조합인 ICL 200–300 kpc × 위성 50–75 kpc는 i 26개, i2 68개이며, 400–500 kpc × 50–75 kpc는 21/57개다. 최대 지원 조합으로 과학 질문을 바꾸지 않았다. 이번 단계는 지원 감사뿐이며 상관·유의성·선택 보정을 하지 않았다. 이 표본도 미검증 PSF 날개·sky 전달 때문에 `science_measurement=false`다. 재현기 `analysis/icl_wgl_membership/hsc_debug/audit_icl_satellite_radial_support_matrix.py`, 테스트 3개 통과; 보고서 `analysis/icl_wgl_membership/hsc_debug/icl_satellite_radial_support_matrix_assessment_2026-10-20.md`. 결과 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_radial_support_matrix_2026-10-20_01/`, summary SHA-256 `7d00bc9fad2e9066cc152c2485c3f5c301637594ac6acd37212cea1264204126`, 84행 table SHA-256 `1ee4f8e4a9c682296c64781d82cb74a31267acc1952b4814b643aa3b779be32d`.

### 2026-10-20 후속: 전체 ICL×위성 반경 가족 검정

상기 84개 반경/필터 조합의 exact 프로필 교집합과 추출된 대조장 좌표를 사용했다. 420개 고유 표적/대조 cluster-filter 쌍을 공통 5도 공간블록에 배정했고, 각 셀의 `satellite ΔSB ~ ICL ΔSB + z + richness` 부분상관을 CR1 공간블록 t 통계와 함께 계산했다. 84개 셀 전부 10개 행·8개 블록 최소기준을 통과했다. 모든 셀에 공통 Webb 6점 wild-cluster bootstrap 5,000회를 적용해 Holm 및 max-T 완전귀무 가족조정을 수행했다.

미보정 `p<0.05` 8개가 있었으나 Holm 조정 최소 p는 0.588, max-T 조정 최소 p는 0.156으로 84개 중 어느 것도 보정 후 유의하지 않았다. 후보 ICL 200–300 kpc × 위성 50–75 kpc의 i2는 `r=0.293`, 68쌍/25블록, `p_raw=0.0288`이지만 `p_Holm=1.000`, `p_maxT=0.987`; i는 `r=0.079`, 26쌍/9블록, `p_raw=0.621`이다. 표본별 지원 차이와 필터 불일치가 남고 PSF 날개·sky 전달·완전도 미해결이어서 탐색 결과이며 `science_measurement=false`다. 상세·한계는 `analysis/icl_wgl_membership/hsc_debug/icl_satellite_radial_matrix_familywise_assessment_2026-10-20.md` 및 재현기 `analyze_icl_satellite_radial_matrix_familywise.py`에 기록했다. 회귀검사 6개와 독립 계수·p조정·해시 검증을 통과했다. 결과 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_radial_matrix_familywise_2026-10-20_01/`; summary SHA-256 `9e4bf2e8369c2c7b3de025a51eafc74538b560b9c42a9cf2fc68de49e1a8ccfb`, association TSV SHA-256 `99ae3afd7a8845de9fb668cb415dcc2b7f8edded60bba550e2d69b5765014b0d`, max-T null TSV SHA-256 `2fcfe2a763d33a1ee12255c35dc77e5623dfdcf5121f4283597629d3368f4c19`.

### 2026-10-20 후속: CAMIRA 위성 가중을 포함한 168개 가족 민감도

기존 hash-verified member-weighted 프로필의 CAMIRA 열 8 `weight`를 사용해 84개 반경/필터 equal-satellite 분석과 CAMIRA-weighted 분석을 같은 168개 검정 가족으로 결합했다. 두 estimator는 모든 셀에서 동일한 cluster/filter 행과 공간블록을 유지했다. 각 클러스터 안에서만 CAMIRA weight를 적용하고 cluster는 동일 가중했다. 5,000회 공통 Webb wild-cluster bootstrap의 max-T와 Holm 보정 후 유의한 셀은 0개였다(17개 원시 `p<0.05`, 최소 Holm `1.000`, 최소 max-T `0.255`). ICL 200–300 × 위성 50–75 kpc의 i2는 equal `r=0.293`, weighted `r=0.241` (둘 다 68쌍/25블록); 원시 p는 각각 `0.030`과 `0.055`지만 168개 보정 후 둘 다 `p_Holm=1.000`, `p_maxT≥0.999`다. CAMIRA weight는 독립 불편 확률이 아니며 효과 차이에 별도 추론을 하지 않았다. PSF·sky·완전도 문제가 남아 `science_measurement=false`다. 재현기, 독립 validator, 결과·제한은 `analysis/icl_wgl_membership/hsc_debug/analyze_icl_satellite_weighted_radial_matrix_familywise.py`, `validate_icl_satellite_weighted_radial_matrix_familywise.py`, `icl_satellite_weighted_radial_matrix_familywise_assessment_2026-10-20.md`에 있다. 테스트 8개 및 168행/5,000 null·Holm/max-T·입력/출력 해시 독립 검증 통과. Slurm 430809는 4분 16초, max RSS 57,896 KB, stderr 0으로 종료했다. 결과 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_weighted_radial_matrix_familywise_2026-10-20_01/`, summary SHA-256 `0a337df80fda6fb92f178d251e09347f27012f5cc53f088fd967311b31ef8068`, 결과표 SHA-256 `6e6c37a44823ab0d6ac5f1795af87f4d47cc3dda1f65b0872ea9e622828142e7`, max-T null SHA-256 `04386baebe37577b40b96d2c859f9202e5cdc73b5eaf01b7ed6c1cc6b53e56a5`.

같은 200–300 × 50–75 kpc endpoint에서 equal 대비 CAMIRA-weighted 상관 변화량을 같은 cluster/filter 지원과 연결 5도 공간블록으로 10,000회 paired bootstrap했다. i: Δr=−0.003, 95% 구간 [−0.225, 0.092] (26개/9블록); i2: Δr=−0.052, [−0.139, 0.061] (68개/25블록). 두 구간 모두 0을 포함하며, 선택된 탐색 endpoint에 대한 민감도이지 독립 유의성 검정이 아니다. 분석 및 입력과 bootstrap 재계산 검증은 `analysis/icl_wgl_membership/hsc_debug/analyze_icl_satellite_weighting_contrast.py`, `icl_satellite_weighted_radial_matrix_familywise_assessment_2026-10-20.md`에 기록했다. 결과 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_weighting_contrast_2026-10-20_01/`; summary SHA-256 `b1a07b7cd46c8a2dc8fd562663a3297ba1d4f46c3a0b25df9ab35f9b56ffc699`, 결과표 SHA-256 `992bf7e7883b3b21babf88d6ba73c0d54e52d1e2ed8b1a7ff2ed07b294d0a012`. `science_measurement=false` 유지.

2026-10-20 cluster-only weak-lensing cross-shear null을 반경/source-bin별로 분해했다. 전체 28성분 χ²=74.75/28 (p=3.93×10⁻⁶)와 일치하지만 어떤 단일 반경 가족도 Holm 후 유의하지 않다(최소 p=0.102); source-bin family도 모두 비유의(최소 p=0.310)이며, 28개 개별 성분 family의 최소 Holm p=0.966이다. 현 요약통계에서 null 오염의 위치를 특정하지 못했다. 재현기·테스트·결과 SHA는 `analysis/icl_wgl_membership/hsc_debug/wl_tomographic_consistency_assessment_2026-10-20.md`에 있고, JSON은 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_cross_shear_null_decomposition_2026-10-20_01/`에 있다. 추가 CAS 조회 없음. `science_measurement=false` 유지.

2026-10-21 PSF support sweep은 outer 600 kpc 기준 91개 star/control 쌍, 78프레임을 찾았다. 첫 `_01` stack에서 `i`/`i2` 필터 혼합 오류를 발견해 이를 제외하고 필터별 `_02` 48셀 집계만 사용한다. 독립 검증기가 프로필 4,368행과 10,000회 frame-bootstrap을 재현했다. `i` 29쌍/25프레임, `i2` 62쌍/53프레임이다. 엄격 지지 기준에서 `i`의 50–75 kpc 표본은 0, `i2`는 4–6프레임이고, 200–300 kpc `i2` 지원은 0이다. mask 지지 완화에도 bin당 독립 프레임 최대 1–6개뿐이라 PSF 보정은 불가하다. 보고서 `analysis/icl_wgl_membership/hsc_debug/empirical_psf_outer_frame_bootstrap_assessment_2026-10-21.md`; 권위 산출물 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/empirical_psf_outer_field_wings_frame_bootstrap_2026-10-20_02/`, summary SHA `556fe0e828a5f2b0ff1600c4e4f2449147e9b383a115d3757eb4b5ea18faa388`, validation SHA `844337eb4fa89f72648ec1d660f57cc89208df2fa3ff9b44f15b406a3116da6f`. `psf_approved=false`, `science_measurement=false`.

### 2026-10-21 재연결 후 read-only 확인

2026-10-21 01:59 KST `lageunha` read-only 확인에서 승인된 random-center coordinate export PID 1049265가 여전히 살아 있었다(`Ss+`, elapsed 1일 14:48:55). 체크포인트는 `running`, mtime `2026-10-19 11:07:53 +0900`, 크기 936 bytes다. `random_coordinates_000`/`001`은 완료(총 33,514행), `002`는 CAS job `1062312` 제출 상태이며 최상위 `HTTPError` marker가 남아 있다. 이후 batch는 제출하지 않았고 PID를 재시작/종료하지 않았다. 완료된 32개 batch의 row/hash 검증 뒤에만 전체 기하 스크리닝을 진행한다.

같은 확인에서 target/control manifest(800/798개)를 다시 읽어 1,598개 고유 중심 중 승인 PSF 코어가 있는 위치 56개, 미포함 위치 1,542개라는 기존 오프라인 계획을 확인했다. 계획 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_psf_picker_target_control_centers_plan_2026-10-21_01/psf_picker_plan.json`은 `prepared_not_submitted`; `pdr3_wide`, `coadd`, `HSC-I` 중심 모델 PSF 1,542개만 요청하며 image/shape/photo-z는 요청하지 않는다. 좌표 목록 SHA-256은 `78bb8b426ad3aa99df68d757b8981e0893517c9e246667aa15d0c891d071d99c`, mapping SHA-256은 `4c2f060f350b65fc93dbb0a8b7c12421517b9be54ddba09c866b08a401dc0212`다. 이전 PSF 파일 크기 분포로 계산한 예상량은 평균 27.6 MiB, 중앙값/관측 최대 기준 29.6 MiB였다. 사용자가 2026-10-21에 이 1,542개 범위만 승인했다. HSC PDR3 문서상 i/i2 visit은 하나의 HSC-I coadd로 합쳐지므로 요청 band는 HSC-I로 유지하되, 분석 진단의 i/i2 표본은 계속 분리한다.

승인 후 기존과 동일한 공식 HSC PDR3 PSF Picker `downloadPsf.py`를 재사용했다. 내려받은 원본은 SHA-256 `17cee33a86de87a77f3f800049b3fdb056f2e0c50b96c2b8cf7e80c1f62f2298`로 앞선 성공 실행과 일치한다. 요청은 좌표목록 방식으로 `rerun=pdr3_wide`, `type=coadd`, `filter=HSC-I`, centered stamp만 받았으며 정상 종료(exit 0)했다. 데이터는 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_psf_picker_target_control_centers_download_2026-10-21_01/`에 보존했다. 승인된 1,542개 FITS가 좌표·매핑 해시와 일대일 대응하고 누락·중복·추가 파일은 0개다. 총 29,269,440 bytes이며 크기는 17,280–20,160 bytes, 형상은 43×43 911개, 43×41 242개, 41×43 266개, 41×41 123개다. 모든 픽셀은 유한하고 전체 합은 0.9999999999999986–1.0000000000000018, 양의 광량 중심의 최대 편차는 0.2761 pixel이다. 음수 픽셀 699,893개를 포함한 원본 부호를 보존해야 한다. 파일별 SHA-256 및 좌표·FITS 검증 manifest는 같은 폴더의 `psf_picker_audit.json`; SHA-256 `b6079bd47f5dd0a716a5841fc7c1cca01f96a857a4cbd7a4d8d9132b4856e6b4`다. 감사기는 `audit_hsc_psf_picker_target_control_batch.py`로 저장하고 3개 테스트를 통과시켰다. 이 스크립트로 실제 1,542개 전체를 다시 감사해 기존 manifest와 바이트 단위로 같은 해시를 얻었다. 재현 실행은 기존 파일을 덮어쓰지 않도록 `psf_picker_audit_reproducible.json`에 기록했으며 SHA-256도 동일하다. 기존 승인 PSF 56개와 합쳐 800 target + 798 matched-control 중심 1,598곳 모두에 모델 PSF 코어가 확보됐다. 그러나 이는 중심 코어 입력의 완전성만 확인한 것으로 PSF 날개·영상별 하늘 전달·주입회수 성능은 검증하지 않았다. `psf_approved=false`, `science_measurement=false`를 유지한다. 영상·shape·photo-z는 받지 않았고 자격정보도 저장하지 않았다. 이 다운로드는 실행 중인 random-coordinate CAS 작업을 건드리지 않았다.

공개 확장 HSC-i PSF 대칭·비대칭 FITS 2개(각 9,809,280 bytes, 합계 18.71 MiB)는 이미 사용자 승인 후 수신되어 Zenodo 체크섬과 독립 재검증을 마쳤다. 직전 대화에서 이 자료가 아직 미수신인 것처럼 말한 것은 부정확했다. 기존 위치별 PSF-picker 744 코어도 이미 수신됐다. 다만 두 종류 모두 PSF/하늘 전달을 과학적으로 승인하지 않았으며 `psf_approved=false`, `science_measurement=false`를 유지한다. 세부 결과는 `analysis/icl_wgl_membership/hsc_debug/extended_psf_public_reference_assessment_2026-10-21.md`에 있다.

### 2026-10-21 전체 target/control 위치별 PSF 코어 및 ICL 의존성

신규 재현 감사기 `audit_hsc_psf_picker_target_control_batch.py`를 저장하고 테스트했다. 기존 감사 manifest와 공식 도구 해시 `17cee33a…f2298`를 확인한 1,542개 수신 PSF에서 이 스크립트를 다시 실행해 바이트 단위로 동일한 audit SHA를 재현했다. 감사 테스트와 새 target/control 연결 테스트를 합쳐 9개가 통과했다. 신규 1,542개와 기존 승인 위치 PSF 중 target/control 중심과 겹치는 56개를 합쳐 1,598개 중심을 `analyze_hsc_target_control_psf_cores.py`로 분석했다. 중심 table은 1,598행, paired table은 798행이며 800 target 중 control이 없는 것은 2개다. target-control 장축 moment FWHM 차이 중앙값은 −0.0070 arcsec, q16–q84는 −0.0572–+0.0392 arcsec, 중앙 절대차는 0.0319 arcsec다. 해당 요약 SHA는 `cb9244c39e5608932434eec8eb0e7aca0a50ebac585cd769a19f996937843bb0`이다.

이 코어 차이가 기존 ICL 200–300 kpc 초과광과 함께 변하는지 공통 면적 40% 통과 137쌍(44 i, 93 i2)에서 추가 진단했다. 적색편이·풍부도·공통 피복률을 조정한 partial r은 i `0.154`, i2 `0.038`; 중앙 마스크 flag까지 넣으면 i `0.016`, i2 `0.036`이다. 5,000회 연결 5° 공간블록 bootstrap에서 여섯 모델의 모든 2.5–97.5 percentile 범위가 0을 포함한다. 독립 validator가 6개 점추정, 총 30,000개 재표집 상관, leave-one-block-out 범위를 재현했다. 보고서는 `analysis/icl_wgl_membership/hsc_debug/icl_psf_core_mismatch_full_sample_assessment_2026-10-21.md`; 분석·독립검증 코드는 `analyze_full_sample_icl_psf_core_dependence.py`, `validate_full_sample_icl_psf_core_dependence.py`다. 결과 요약 SHA-256은 `f67019b8cb6f3c096b45e8640d489022db19f2b1d29f7b4db1ab89441936c34e`, joined table SHA-256은 `d6e0f0b669e3eb2c0f7a57c37669260da4ad17216daaf66fda2c60af7223fcf5`다. Slurm controller DNS 실패로 이 137행 단일 CPU 작업은 로그인 노드에서 약 6초에 끝냈다. 코어가 target-control 중앙값에서 비슷해 보인다는 결과는 확장 PSF 날개나 하늘 전달 검증을 뜻하지 않으므로 `psf_approved=false`, `science_measurement=false`다.

weak-lensing random-control-coordinate export는 재확인 당시 PID 1049265가 살아 있었으나 체크포인트 mtime은 계속 `2026-10-19 11:07:53 KST`였다. CAS job `1062310`, `1062311` 두 개는 33,514행으로 완료되어 있고 `1062312`는 체크포인트에서 `submitted`로 남아 있다. Astrafek에 추가한 read-only `check_wl_cas_job_status.py`로 CAS에 딱 한 번 상태를 조회한 결과 `1062312`는 실제 `running`이었다. 이 상태 조회는 `query_submitted=false`, `rows_downloaded=false`이며 자격정보를 저장하지 않았다. `audit_wl_random_coordinate_checkpoint.py`도 다시 실행해 2/32 배치, 33,514/780,191행과 파일 해시를 확인했다. 새 audit `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_coordinate_checkpoint_audit_2026-10-21_01/checkpoint_audit.json` SHA-256은 `d5072ee7cecd4ec4798892dc6d3be4acf1d916738e0ff65564666234a3ec8212`, checkpoint SHA-256은 `84cccc6443ec8d2fc8bee1611ab6b9b57a3aabdc90978fcdfc558151d2bb92b2`다. 활성 submitter를 중단·재시작하지 않았고 같은 CAS 쿼리를 재제출하지 않았다. 저장된 완료분과 결합하기 전에는 기존 작업의 terminal status 및 체크포인트 갱신을 확인해야 한다. 상태조회기 테스트 2개 통과.

### 2026-10-21 확장 피복 계획 및 CAS 좌표 작업 재확인

lageunha에서 승인된 `submit_icl_expansion_coverage_metadata.py --check`를 실행했다. 기존 저장 plan의 체크섬·SQL 구조·열 제한이 모두 통과했고, 후보 2,596개·probe 23,364개·배치 52개를 확인했다. 이 계획은 PDR3 `mosaic`/`skymap`/`patch_qa` 메타데이터만 반환하며 image, object/source rows, shape, 개별 photo-z 요청 플래그가 모두 false다. 이 검사는 네트워크 요청을 하지 않았고 아직 `prepared_not_submitted` 상태다.

같은 read-only 확인에서 weak-lensing 좌표 export PID 1049265가 `lageunha`에서 계속 살아 있었다(`Ss+`, `hrtime`, elapsed 1일 15:56). `cas_coordinates/coordinate_status.json`은 2026-10-19 11:07 KST 이후 갱신되지 않았으며 32개 중 2개 완료(33,514/780,191행), CAS job 1062312 제출 상태, `error_type=HTTPError`를 기록한다. 현재 `squeue`는 Slurm controller DNS SRV 조회 오류로 상태를 가져오지 못했다. 이 확인에서는 HSC CAS job 상태를 새로 조회하지 않았고 좌표 worker를 종료·재시작하거나 새 CAS 쿼리를 제출하지 않았다. 다음 단계는 기존 worker/job의 권위 있는 terminal 상태를 확인한 뒤, 이미 승인된 피복 메타데이터 52개 batch를 기존 HSC 세션에서 순차 실행하는 것이다.

### 2026-10-21 전체 표적·대조 PSF 코어/공개 확장 PSF 비교

새 Astrafek 분석기 `analyze_public_psf_target_control_overlap.py`로 800개 표적과 798개 대조장의 모델 PSF 코어 1,598개를 공개 HSC-I 확장 PSF와 1.5–3.3 arcsec에서 비교했다. target filter는 기존 검증 스택, control filter는 해당 대조장의 실제 `coadd_header_filter`에서 가져왔다. 동일 filter 쌍은 i 302개, i2 481개(총 783개); 다른 filter 15쌍은 paired 통계에서 제외됐다. 2개 target에는 control이 없다.

전 창 진폭 scale 중앙값(i target/control, i2 target/control)은 0.713/0.740, 0.787/0.797이다. 중앙 scale은 inner 1.5–2.4에서 outer 2.4–3.3 arcsec로 각각 0.731→0.520, 0.754→0.560, 0.803→0.555, 0.817→0.572로 감소한다. full-window 모양 잔차 RMS 중앙값은 11.7–12.8%다. same-filter target-minus-control scale 차이 중앙값은 i −0.0239 (302쌍, q16–q84 −0.1294…0.0644), i2 −0.0033 (481쌍, −0.1181…0.0988)이다. 이는 위치별 분포/민감도이지 보정계수나 추론 검정이 아니다. 한 진폭으로 public wing을 접합할 수 있다는 증거가 없으며 200–300 kpc wing 및 sky 전달은 계속 미검증이다.

초기 `_01` 실행은 control에 target filter를 상속하여 paired 결과에 사용하지 않는다. 실제 control filter를 조인한 `_02`만 권위 출력이다. 독립 FITS 재계산기가 1,598개 전체의 네 fit window, 동일-filter 783쌍, 15 mismatch를 재현했다. 7개 focused test와 syntax/whitespace 검사 통과. 분석기 `analysis/icl_wgl_membership/hsc_debug/analyze_public_psf_target_control_overlap.py`, 독립 검증기 `validate_public_psf_target_control_overlap.py`, 상세 보고서 `extended_psf_target_control_overlap_assessment_2026-10-21.md`. GPFS 결과 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/extended_psf_target_control_overlap_2026-10-21_02/`; summary SHA-256 `8ba338cf139ab469d1c21a10f4401ea3d6a86e3dd402ea07bc49ca512be2e6d8`, position TSV `ffe0a557d47694f73156b8f084d22f5977d658bb1b6932fe1cb5400dba22c62a`, paired TSV `b6276e52678abc91499ce56a65acbaec0ceaa055599aa56fbecb2a6be1686556`, validation JSON `add0ce56323fe32564147f6ab9ec91d05d8ee3d2810f62c69494ac1b38b607ec`. `psf_approved=false`, `science_measurement=false` 유지. 후속은 기존 image mask에서 splice radius·진폭별 점광원/diffuse injection recovery다.

### 2026-10-21 `i` 필터 공개 PSF 주입회수

Astrafek `extended_psf_injection.py`로 동일 필터 PSF 지원과 K2.5/40% ICL 지원을 만족하는 13개 `i` 프레임에서 HSC-I 대칭·비대칭 public PSF를 각각 주입 민감도에 사용했다. 두 run 모두 13/13, error 0, 각각 약 68초. 독립 `validate_extended_psf_injection_batch.py`가 저장 응답을 재계산했고 최대 산술 오차는 0이다. 기존 public PSF는 새 승인 위치의 동일 byte file(SHA `8897af…2d85`)임을 확인했다.

대칭 모델 회수 중앙값/프레임 범위: smooth diffuse→ICL 200–300 kpc `0.215` (`−2.734…0.271`), 위성 외곽광→ICL `0.353` (`0.024…1.269`), diffuse→위성 50–75 kpc `−0.039` (`−0.915…4.324`), 위성 외곽광→위성 `0.294` (`−3.024…1.895`). 13개 표본은 선택된 지원 집합이며 quantile은 모집단 CI가 아니다. 비대칭−대칭 paired 모델 차이는 네 응답에서 최대 절댓값 `4.77×10⁻⁷`; 동일한 공개 별 표본을 쓰므로 독립 검증은 아니다. 음수와 넓은 분산 때문에 scalar recovery correction을 적용하지 않는다. 실제 science coadd는 수정하지 않았고 source photon noise/coadd covariance는 추가하지 않았다. `science_measurement=false`, `psf_approved=false`.

결과 보고서 `analysis/icl_wgl_membership/hsc_debug/extended_psf_i_injection_recovery_assessment_2026-10-21.md`. symmetric run `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_extended_psf_i_sensitivity_2026-10-21_01/` summary SHA `e9853c6626cace7b630d6db17e3d9a8587a26007d726b884801608f2719be228`, validation SHA `aa1c1b3fd80f044d76f66a25ca2da19517b90c2e511ef48d3d2c3033b43df2d1`. asymmetric run `icl_extended_psf_i_asymmetric_sensitivity_2026-10-21_01/` summary SHA `dd3f40c56c55840a6a268eb9192cfaf7acede0183a138951df978bf63f882527`, validation SHA `b62255417a19e069bbb53879b002bde3932a40d10a98cf2f3f8cba80db382aa8`. Filter-aware paired comparison `variant_comparison_filter_aware.json` SHA `9281dca86f14aaa4ce7e63537b021e2ce861e1c3663f5b4a17b7ff86f5fba575`; earlier `variant_comparison.json` is retained but superseded. Comparator `compare_extended_psf_injection_variants.py`를 i/i2 필터 선택형으로 일반화; focused test 8개 통과. 다음은 여러 splice 반경·진폭의 real-mask 주입과 더 넓은 local wing 검증이다.


### 2026-10-21 고정 마스크 주입 민감도 검증 및 Slurm 위치

`icl_extended_psf_i_fixed_masks_2026-10-21_01/`의 13개 프레임을 독립 validator로 다시 계산했다. 13/13 완료, 오류 0, 저장 산술 최대 오차 0이다. 고정 baseline mask에서 diffuse→ICL 회수 중앙값은 0.22495 (q16–q84 0.22175–0.22717), diffuse→위성 50–75 kpc는 0.00062 (−0.01161–0.01260)다. 위성-envelope→ICL은 0.494 (0.165–0.713), 위성-envelope→위성은 0.278 (−0.174–0.509)으로 여전히 분산이 크다. 재탐지−고정마스크 diffuse→ICL 차이는 중앙값 −0.0115지만 최솟값 −2.955이며 diffuse→위성은 범위 −0.900–4.326이다. 따라서 일부 프레임의 극단값은 재탐지/마스크 선택에 민감하다. 이는 합성 주입 민감도일 뿐 PSF 보정이나 과학 측정이 아니다. 보고서 `extended_psf_i_injection_recovery_assessment_2026-10-21.md`에 추가했다. 결과 SHA-256: summary `a732fbc5fdaf41932f5685a329143f8d794f543885b5a5c0917ea321c3c264c5`, validation `2bcd49fa584c18e2415be56ed4fa75fa4894f2762b839ec9cdf689b2147c0b45`. 관련 focused pytest 11개 통과.

운영 위치를 확인했다. Slurm 작업의 조회·제출은 syntax에서 한다. 단, 1062312는 Slurm job ID가 아니라 HSC CAS job ID이므로 `squeue`/`sacct` 조회 결과는 이 CAS 요청의 상태 증거가 아니다. lageunha 좌표 exporter PID 1049265는 계속 살아 있으며 checkpoint는 2/32 배치, 33,514/780,191행, 마지막 갱신은 2026-10-19 11:07 KST다. CAS 상태를 별도 인증 조회하지 않았고, 프로세스를 종료·재시작하거나 요청을 재제출하지 않았다.


### 2026-10-21 공개 PSF 주입의 재탐지 영향도

Astrafek `audit_extended_psf_injection_reselection_influence.py`와 독립 focused tests 2개를 추가했다. 13개 `i` 프레임의 입력 파일 해시, summary-to-frame 응답값, redetection−fixed-mask 산술을 대조하고 프레임별 지원 변화와 leave-one-frame-out 중앙값을 기록한다. 가장 큰 재탐지 outlier는 cluster 4982다. baseline에서 matched satellite 4개/control 16개였으나 diffuse-only·joint injection에서 2개/8개로 감소했고 ICL 200–300 kpc coverage는 0.5395에서 0.4982로 떨어졌다. 이 프레임의 diffuse→ICL 재탐지−고정마스크 응답차는 −2.955, diffuse→위성은 +4.326이다. cluster 4982 제외 시 diffuse→ICL 재탐지 median은 0.2150→0.2181, 고정마스크 median은 0.22495→0.22522로 유지되지만 satellite-envelope→satellite median은 0.294→0.218로 변한다. 선택된 13프레임의 기술적 leave-one-out일 뿐 모집단 추론이 아니다. 보고서 `extended_psf_i_injection_recovery_assessment_2026-10-21.md`; 산출물 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_extended_psf_i_reselection_influence_2026-10-21_01/reselection_influence_audit.json` SHA-256 `09dc009995b8dc8bbbc29cff0b392812f47b5351e399c57c943a26d8a56cee29`, 도구 SHA-256 `6b9f2a1e9f9f97dcbc1c27a21e1383524a42aed8a023c5a8b5fe086c8129ccad`. `science_measurement=false`, `psf_approved=false` 유지.

기존 random-control 좌표 exporter PID 1049265를 60초 간격 재확인했다. lageunha에서 여전히 `hrtimer_nanosleep`이며 checkpoint는 2026-10-19 11:07:53 KST 그대로(2/32 배치, 33,514/780,191행; batch 002 CAS 1062312 `submitted`)다. CAS 재조회·재제출·프로세스 재시작은 하지 않았다. 1062312는 Slurm ID가 아니라 HSC CAS ID이며 Slurm 작업은 syntax에서 조회·제출한다.

### 2026-10-21 위치별 PSF 코어/공개 날개 하이브리드 구성기

`position_psf.py`의 `build_core_wing_hybrid_psf`는 signed local core와 public extended wing을 smooth radial taper로 결합해 단위 정규화한다. 중앙 코어·음수 픽셀 보존 및 범위 검사를 단위시험으로 확인했다. 이 커널은 synthetic-light transfer 민감도 시험에 적용됐으나 observed coadd를 수정하지 않았고 새 HSC 요청도 없다. 결과는 `psf_approved=false`, `science_measurement=false`다. 재현 및 제한은 `analysis/icl_wgl_membership/hsc_debug/core_wing_hybrid_psf_design_2026-10-21.md`에 있다.

### 2026-10-21 i2 확장 PSF 프록시/위치별 코어 주입 민감도 비교

이미 완료된 HSC-I public symmetric/asymmetric proxy 주입을 위치별 PSF core 주입의 동일 44개 `i2` 프레임과 paired 비교했다. 기존 배치 각각 44/44 프레임, 오류 0이고 저장 응답을 독립 재계산한 최대 산술 차이도 0이다. 프레임·지원·주입 규격을 맞춘 재현 감사기 `audit_extended_psf_core_vs_wide.py`에 `i2` 선택 경로를 추가했고 focused tests 3개가 통과했다. 산출물 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_extended_psf_i2_proxy_variant_comparison_2026-10-21_01/core_vs_wide_audit.json`, SHA-256 `ba1287493d5d3db4d8050f601ac982a659e1ca8588c024c90d8a7814c30e0f97`; 감사기 SHA-256 `8d08096c9d4babab32aa61fc6bef34d710267e65873c7cc1296ac3c36009f625`.

Diffuse→ICL 200–300 kpc 회수율 중앙값은 local core `0.2245`, public proxy `0.2235`, paired 차이 `−0.00036` (q16–q84 `−0.00092–0.00093`)였다. Satellite-envelope→ICL 중앙값은 `0.4904`와 `0.4849`, paired 차이 `0.00846` (`−0.00773–0.02847`)였다. 반면 satellite 50–75 kpc 회수 분산은 두 모델에서 매우 컸다. Satellite-envelope→satellite paired 차이의 q16–q84는 `−0.191–0.081`, 전체 범위는 `−3.551–0.377`이다. 대칭·비대칭 public PSF 모델 차이는 최대 `9.5×10⁻⁷`로 수치 정밀도 수준이며, 두 모델은 같은 별 표본을 써 독립 검증이 아니다. HSC-I public PSF는 i2 전용 보정이 아니므로 프록시 민감도에 한정하고 회수율 보정이나 과학 측정으로 승격하지 않는다. 보고서 `analysis/icl_wgl_membership/hsc_debug/extended_psf_i2_proxy_core_comparison_assessment_2026-10-21.md`; `science_measurement=false`, `psf_approved=false` 유지. 다음 ICL gate는 기존에 설계한 multi-radius/amplitude local-core/public-wing hybrid 점광원·diffuse real-mask 주입 파일럿이다.

### 2026-10-21 위치별 하이브리드 위성 PSF 경로 및 i-band 파일럿 실행

하이브리드 주입에서 위성 envelope에 은하단 중심 PSF 코어를 재사용하던 것을 수정했다. 확산 ICL에는 중심 위치 코어+public wing hybrid를, 각 위성 envelope에는 해당 위성 위치 코어+같은 wing을 적용한다. 결과 파일과 독립 validator는 각 위성 코어의 filename/SHA, 정규화, wing parameters, 승인 플래그를 확인한다. `test_position_psf.py`, `test_extended_psf_injection.py`, `test_validate_extended_psf_injection_batch.py`, 새 grid 테스트를 포함한 focused test 22개 통과했다.

미리 고정한 12개 조합은 splice 1.8/2.4/3.0 arcsec, blend 0.6 arcsec, wing amplitude 0.50/0.65/0.80/1.00이다. cluster 1056의 10개 위성 위치를 이용한 grid job `432792`는 syntax `syn04`에서 13m58s 만에 정상 완료됐다(1 CPU, 32 GB, GPU 미요청; peak RSS 2,866,960 KiB). 12/12 설정과 독립 response 산술 검증이 통과했고 최대 오차는 0이다. Diffuse→ICL 응답은 grid에서 0.22111–0.22128로 안정적이었지만 위성 외곽광→위성은 mask 재탐지 시 0.09297–0.16836, 고정 mask에서 0.31063–0.32470이었다. 단일 cluster라 모집단 범위는 아니다. 재현 보고서 `core_wing_hybrid_i_pilot_assessment_2026-10-21.md`; `grid_audit.json` SHA-256 `18606dcd49e8b918186eda6158e2209f2755a76b767806756e3571d7b4990ee1`. 다음은 13개 supported i 프레임의 frame-to-frame 민감도다. 이 grid에는 point source model profile이 포함되지 않았고 별도 noiseless point-source mask diagnostic만 cluster 1056에서 수행됐다(아래 절). 모든 과학·PSF 승인 플래그는 false다.


### 2026-10-21 추가: position-specific hybrid 단일 프레임/점광원 진단 완료

`lageunha`에서 사전 정의한 단일 프레임 cluster 1056의 새 고유 출력 경로로 splice 2.4 arcsec, blend 0.6 arcsec, wing amplitude 0.65를 실행했다. 앞서 기록한 syntax Slurm array(job 432792)와는 별도 실행이다. 이번에는 Slurm 명령을 실행하지 않았다. 1/1 프레임, 10개 satellite 위치 코어, elapsed 12.12초, 오류 0으로 완료했다. 독립 validator는 모든 저장 응답을 재계산해 최대 오차 0, `validated_exploratory_transfer_only`를 반환했다.

고정 mask 응답은 diffuse→ICL 0.230718, 위성 envelope→ICL 0.506852, diffuse→satellite −0.012496, satellite envelope→satellite 0.314188이다. 동일 cluster의 기존 public PSF fixed-mask 값과 차이는 각각 약 +0.000105, −0.000554, −0.000093, −0.012785다. 별도 noiseless unit point-source 진단에서 50–100 kpc annulus의 mask coverage/마스크 내 flux fraction은 46.3%/0.001996, 200–300 kpc는 55.4%/0.000834다. 이는 단일 프레임 민감도 진단이며 Poisson noise, coadd covariance, sky transfer, 관측원 제거를 반영하지 않는다.

결과 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_hybrid_psf_i_pilot_2026-10-21_r2p4_a0p65_perpos_point_01/`; `summary.json` SHA-256 `8a2fa26dd197aae9f29fafaaaac7e91afeab0bcc9c1d2299b2a73704a367d16e`, `cluster_1056.json` `ddb650883c7e664fb9c026d28a83768c195d7b5ccaa867b6b5c58d8048666711`, `validation.json` `d547fba7bfbae4869f48c6946a316190c860afe4dd711b301c2205e55d03c568`. focused tests 23개 통과; `py_compile`, array script `bash -n`, `git diff --check` 통과. 설계와 다음 grid 단계는 `analysis/icl_wgl_membership/hsc_debug/core_wing_hybrid_psf_design_2026-10-21.md`. 해당 결과도 `psf_approved=false`, `science_measurement=false`다.

2026-10-22 `validate_extended_psf_injection_batch.py`로 위 단일 hybrid+point-profile batch를 독립 재검증했다. 출력 `validation_repro_2026-10-22.json`은 기존 검증과 동일 SHA-256 `d547fba7bfbae4869f48c6946a316190c860afe4dd711b301c2205e55d03c568`, 1/1 프레임, 최대 응답 산술 오차 0이다.

### 2026-10-22 supported i 프레임 전체 hybrid 민감도 grid 실행 중

2026-10-21 cluster 1056 단일 프레임에서 12개 PSF 모델 grid가 끝난 뒤, 13개 supported `i` 프레임으로 동일한 frozen grid를 확대했다. 프레임별 결과에는 redetection/fixed-mask diffuse·satellite 응답과 별도의 noiseless unit point-source profile/mask transmission이 포함된다. 재현 가능한 `--all-supported` driver, updated point-profile grid summarizer, CPU-only Slurm wrapper, 범위 문서 `core_wing_hybrid_i_supported_plan_2026-10-22.md`를 추가했다. 관련 grid/summarizer/validator focused test 14개 통과.

Slurm job `432803`은 `syntax`의 `a10`/`syn04`에서 2:07:02 만에 정상 완료했다(1 CPU, 32 GB, GPU 미요청, peak RSS 4,363,184 KiB). 12개 설정×13 프레임=156개 frame record, 오류 0, 각 저장 응답의 독립 재계산 최대 오차 0이다. 점광원 profile은 별도 noiseless unit source 모델을 saved masks에 통과시킨 것으로 source detection/noise 회수가 아니며 독립 FITS 재구성은 하지 않았다.

중앙 `r2.4_a0.65` 모델에서 fixed-mask diffuse→ICL median 0.2251 (q16–q84 0.2219–0.2272), redetection median 0.2141이며 12개 모델 medians의 envelope는 각각 0.2250–0.2252 및 0.2141–0.2150이다. Satellite envelope→ICL redetection−fixed paired median은 −0.1022이고 PSF grid에 따른 median 범위는 −0.1045…−0.0825다. Satellite envelope→satellite paired 차이는 중앙 모델 median −0.0536이나 q16–q84 −0.3158…+0.2065, 전체 범위 −3.0894…+2.1677 (cluster 430 및 4982 극단)으로 단일 보정값을 쓸 수 없다. Point source 유효 flux fraction은 중앙 모델 median 50–100 kpc 0.00259, 200–300 kpc 0.000957다. 결과는 과학 측정·PSF 승인 또는 ICL–위성 상관의 증거가 아니다.

독립 Astrafek 요약은 `analysis/icl_wgl_membership/hsc_debug/core_wing_hybrid_i_supported_assessment_2026-10-22.md`. v2 audit은 12 configuration 전체가 동일한 13개 고유 cluster ID를 썼는지도 확인한다. `grid_audit_v2.json` SHA-256 `c0545e9f36c1b0d47589744abfcf21264461e8bf467f6f63fd379436ee35ba2c`, checkpoint SHA-256 `91f5aed2dbff64fd5ae07290c30b8f77b69fde15535c881d1822352ba25747b6`. 최초 `grid_audit.json`도 보존했다. 출력 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_core_wing_hybrid_i_supported13_2026-10-22_01/`. 다음 단계는 matched-control sky/PSF-wing leakage와 pixel/noise response를 관련 물리반경까지 검증하는 것이며, 현재 응답 중앙값을 보정값으로 사용하지 않는다. `science_measurement=false`, `psf_approved=false` 유지.
## 2026-10-21 ICL 주입 응답의 마스크 재선택 분리

### 고정 마스크와 재검출 비교

이전 2026-10-19 위치별 HSC PSF 코어 민감도 시험을 이어, PSF를 합성 은하단 영상에
합성곱한 뒤 확산 성분과 위성 성분의 주입 응답을 다시 계산했다. 기존 방식은 각
주입 시나리오마다 source detection을 다시 수행하므로, 주입된 빛의 전달과 검출
마스크·위성 표본의 변경 효과가 한 응답에 섞여 있었다.

이를 분리하기 위해 로컬 작업 트리에서
`analysis/icl_wgl_membership/hsc_debug/extended_psf_injection.py`에 기준 검출 마스크를
고정한 분기와 `calculate_transfer_response`를 추가했다. 프레임 JSON schema v2에는
`measurements_by_scenario_fixed_masks`, `fixed_mask_transfer_response`,
`redetection_minus_fixed_mask_response`를 저장한다. 독립 검사기
`validate_extended_psf_injection_batch.py`는 v1/v2를 읽고 두 분기 및 차이를
재계산한다. 이 코드 변경은 이 인수인계 커밋에 포함하지 않는다. 재개 시 현재 로컬
작업 트리 상태와 테스트를 확인하고 코드 변경은 별도로 검토한다.

### 13개 HSC-I 프레임의 탐색 결과

lageunha CPU 실행은 `icl_extended_psf_i_fixed_masks_2026-10-21_01`에서 13/13 프레임,
오류 0, 72.627초로 완료했다. 독립 검사 결과는
`validated_exploratory_transfer_only`, 13개 검증, 저장 응답과 재계산 응답의
최대 차이 0이다. 요약 SHA-256은
`a732fbc5fdaf41932f5685a329143f8d794f543885b5a5c0917ea321c3c264c5`,
검증 SHA-256은
`2bcd49fa584c18e2415be56ed4fa75fa4894f2762b839ec9cdf689b2147c0b45`이다.
산출물은 실행 데이터의
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_extended_psf_i_fixed_masks_2026-10-21_01`
아래에 있으며 Git에 넣지 않는다.

재검출 마스크 결과에서 확산광→ICL 응답 중앙값은 0.215였으나 범위는
−2.734–0.271이었다. 기준 마스크를 고정하면 같은 응답 중앙값 0.22495,
16–84백분위 0.22175–0.22717, 전체 범위 0.22130–0.23061로 좁아졌다.
고정 마스크 확산광→위성 측정 응답 중앙값은 0.000624
(16–84백분위 −0.0116–0.0126)였다. 이는 재검출 결과의 극단값 상당 부분이
영상 내 광량 전달 자체보다 재검출·측정 대상 변경에서 왔음을 보여준다.

주요 사례인 cluster 4982에서 확산광→ICL은 재검출 −2.734, 고정 마스크 0.2213,
차이 −2.955였고, 확산광→위성 지표는 각각 4.324, −0.00287, 4.326이었다.
재검출 중 위성 표본은 4개에서 2개, ICL 면적 지원율은 0.539에서 0.498로 변했다.
cluster 411은 위성 측정에 usable satellite가 1개뿐이므로 극단 응답을 일반화하지
않는다. 위성→위성 자체 응답도 고정 마스크 중앙값 0.278이지만 표본 간 산포가
커서 보정값으로 사용하지 않는다.

이 결과는 선택된 13개 프레임의 모델 PSF 합성곱 주입 민감도일 뿐이다. 국소 PSF
보정, 확장 날개, coadd 상관 잡음과 주입 광자잡음, upstream sky 전달, 전체 표본
대표성은 검증하지 않았다. 따라서 PSF 보정값·과학 측정·ICL–위성 상관의 증거로
사용하지 말고 `science_measurement=false` 및 PSF 승인 보류를 유지한다.
새 고정 마스크 테스트를 추가하려던 마지막 편집 이후 전체 검사는 다시 실행하지
못했으므로, 앞서 통과한 9개 테스트 결과를 최종 코드 상태의 검사로 간주하지 않는다.

### 작업환경·HSC 반출 상태와 다음 순서

2026-10-21 세션에서는 `exec_command`가 Astrafek 경로, `/tmp`, `/` 및 명시적
`/bin/sh` 실행에서도 셸 프로세스를 시작하지 못하고
`CreateProcess: No such file or directory`를 반환했다. 따라서 이 세션에서
공유 지침, 로컬 `git status`, Slurm 상태 및 실행 데이터 경로를 다시 읽지 못했다.
이는 `/tmp`의 파일 권한 거부를 뜻하지 않는다. 도구가 셸을 띄우기 전에 실패한
것이며 실제 `/tmp` 권한 모드는 확인되지 않았다.

사용자는 추가 HSC PDR3 Wide HSC-I PSF 코어 1,542개 반출을 명시적으로 승인했다.
이 작업 세션에서는 셸·HSC 전송을 시작하지 못했으므로 추가 요청이나 다운로드를
완료했다고 기록하지 않는다. 재개 시 기존 manifest, 완료 상태, 파일 크기·해시를
먼저 확인하고 이미 수신된 파일을 중복 요청하지 않는다. 암호를 다시 요청하거나
문서·로그에 저장하지 않는다.

Weak-lensing 대조장용 좌표 export 프로세스 PID 1049265는 2026-10-20 확인 시
`lageunha`에서 살아 있었다. checkpoint에서 SQL batch 000과 001은 완료,
002는 CAS job 1062312에 제출된 상태다. 이 상태는 마지막 확인 시점의 것이므로
재개할 때 같은 PID와 Slurm 작업 ID의 현재 상태를 먼저 확인한다. 살아 있는 작업을
재시작하거나 중복 제출하지 않는다.

다음 순서:
1. 작업환경 복구 후 `SHARED_CONTEXT.md`, 저장소 로컬 지침, branch/status와 변경 파일을 확인한다.
2. 위 고정 마스크 테스트 추가분을 확인하고 관련 테스트·validator를 재실행한다.
3. 기존 프로세스·checkpoint 및 PSF manifest를 확인한다. 승인 범위를 넘는 반출은 하지 않는다.
4. ICL 및 위성 공동 분석의 전체 표본 마스크·광량 누출·공분산·주입 회수를 계속 검증한다.
5. 관측 weak-lensing 스택과 최종 과학 해석은 별도 검증을 마친 뒤 진행한다.
   이 인수인계 커밋에는 `HANDOVER.md`만 포함한다.

### 2026-10-22 후속: weak-lensing random-control 좌표 export 중간 감사

`syntax`에서 `lageunha`를 read-only 확인했을 때 이전 export PID 1049265는 더 이상
존재하지 않았다. 기존 HSC CAS checkpoint는 여전히 `running`으로 표시되지만, 승인된
32개 배치 중 `random_coordinates_000`과 `_001`만 완료되어 있었다. 이 두 CSV를
Astrafek 재현 감사기 `audit_wl_random_coordinate_partial_export.py`로 렌즈별
count-only 결과, 행 수, checkpoint SHA-256에 대조했다. 15,879행 및 17,635행, 총
33,514행이 모두 통과했고 기존 파일 해시도 일치한다. 승인 계획의 전체 예상은
780,191개 coordinate row다.

배치 `_002`는 CAS job `1062312`로 `submitted` 기록이며 결과 CSV는 없다. 기존 계획과
같은 제출기를 재개해 CAS 상태를 다섯 차례 확인했으며 매번 `running`이었다. 이때 다른
exporter 프로세스는 없었다. 세션 종료 후 원격 polling process를 남기지 않도록 로컬
poller를 Ctrl-C로 중단했으며, checkpoint는 `interrupted`, job `_002`는 여전히
`submitted`와 마지막 CAS 상태 `running`으로 보존됐다. CAS query를 재제출하지 않았고
새 CSV도 다운로드하지 않았다. 나머지 29개 배치는 아직 제출되지 않았다. 오프라인
부분 감사 결과는
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_point_coordinate_plan_2026-10-19_01/cas_coordinates/coordinate_partial_audit_2026-10-22.json`
(SHA-256 `29c24052e09eca4102dcd2837ebdcba701087208767b82ebbe55384218b4d327`)이다.
checkpoint가 `interrupted`인 시점의 감사는 같은 경로의
`coordinate_partial_audit_2026-10-22_02.json` (SHA-256
`06aa349317e3fc7e7da017be36fadbadf01893cfb8e288b5e024a4c55960b049`)이다. 감사기는
HSC에 연결하지 않으며 제출·다운로드·과학 측정을 하지 않는다는 것을 3개 단위시험과
실제 부분 결과로 검증했다. 다음 실행은 job `1062312`의 CAS 상태가 완료 또는 오류로
바뀌었는지 확인한 뒤 같은 checkpoint에서 이어간다. Random-lens subtraction 및
spatial-block covariance가 없으므로 weak-lensing 결과는 계속 `science_measurement=false`다.

### 2026-10-22 후속: 독립 footprint 표본 확장 게이트

800개 parent와 33개 i2 spatial block은 현재 ICL–위성 검정력 목표에 부족하다.
공식 HSC 자료 설명상 PDR3 CAMIRA는 S20A 기반이며, HSCLA2020은 훨씬 넓은 PI 관측을
포함하지만 다중 필터 영역은 성기고 별도 데이터 릴리스다. S19A shape catalog는 별도
PDR3 제품이다. 그러므로 넓은 영상 footprint를 곧바로 weak-lensing 표본으로 세지 않는다.
새 요청 전에 기존 CAMIRA 후보, coadd/control/mask 지원, S19A shape footprint, 독립 sky
block을 로컬 매니페스트로 교차검사하도록
`analysis/icl_wgl_membership/hsc_debug/hsc_independent_footprint_expansion_plan_2026-10-22.md`
에 게이트를 작성했다. 이 문서를 위해 외부 서비스 조회나 신규 데이터 요청은 하지 않았다.

### 2026-10-22 후속: S23B 대조장 후보 오염 민감도

고정된 1,333개 추가 CAMIRA 후보표를 이미 받은 HSC 영상 800개 target과 798개
matched-control의 FITS WCS에만 대조했다. 새로운 영상·카탈로그를 받지 않았고,
픽셀값을 읽지 않았다. target 영상에서 후보 중심은 0개, 대조장 영상 7개에서 서로
다른 후보 중심 7개가 나왔다(필터별 i 3개, i2 4개). 이 결과는 대조장 기준선에
후보 은하단의 확산광이 들어갈 가능성을 표시할 뿐, 검출이나 실제 ICL 광량 오염을
증명하지 않는다. 오프라인 감사 요약은
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_s23b_existing_cutout_support_2026-10-22_01/summary.json`
(SHA-256 `071377813fcfbf15f298a5d48f3e7a3f0b74f60899d88c46eafb7d5bb58ffebb`)이다.

`stack_hsc_icl_target_control_2d.py`에 `--exclude-candidate-control-audit` 입력을
추가했다. 스태커는 후보표 고정 해시, 감사 범위·상태, control manifest/status
해시, target-control 대응, 후보 ID를 검증한다. 후보 7개 중 1개 control은 기존
필터 불일치로 baseline에서 이미 빠져 있어 중복 제외하지 않고, eligible한 나머지
6쌍만 별도 제외한다. 따라서 baseline 783쌍에서 777쌍을 남기는 비교다. 독립
verifier는 감사 JSON과 원 manifest로 두 제외 집합을 재구성해 pair table을 검사한다.
첫 Slurm job 434536은 위 중복 제외를 거부하며 출력 생성 전에 종료됐다. 중복 제외를
분리한 뒤 전체 493개 테스트와 하위 검사 6개가 통과했고, 작업 434537은 19분 10초에
완료됐다. 독립 verifier는 798개 입력 control 중 15개 필터 불일치와 6개 추가 후보
겹침 제외, 777개 pair 행, 40개 radial 행, 두 201×201 map을 검증했다. map은 40,401
pixel을 모두 덮고 `DELTA - (TARGET - CONTROL)` 최대 잔차는 i 0.00391, i2
0.00098 nJy/arcsec²다. 검증 JSON SHA-256은
`f0fb9a4ba6e9048366e481507c73074cebb88cccd6bb87cfa842963bb69535ca`다.

검증 baseline 대비 50–200 kpc 면적 가중 평균은 i +0.357%, i2 −0.020%, 200–300 kpc는
i +1.321%, i2 −1.867% 이동했다. 개별 annulus의 최대 상대 이동은 i 250–275 kpc에서
2.633%, i2 200–225 kpc에서 4.611%다. 따라서 여섯 대조장의 제외는 넓은 구간 평균을
크게 바꾸지 않지만 일부 outer annulus에는 수 퍼센트 변화를 준다. 필터별 이동 부호도
다르며, 이 subset 민감도는 차이의 공분산이나 유의도를 제공하지 않는다. 상세 수치와
해시는 `analysis/icl_wgl_membership/hsc_debug/icl_s23b_control_exclusion_assessment_2026-10-22.md`
에 기록했다. 이 결과는 계속 exploratory이며 `science_measurement=false`; baseline
783쌍 스택은 보존됐다. 영상 반출이나 CAS 조회는 없다.

추가로 새 연관표 중첩 감사기를 포함한 전체 494개 테스트와 하위 검사 6개가
통과했다. 7개 WCS 후보 겹침을 저장된 ICL 200–300 kpc 지원표, 위성 50–75 kpc
프로파일표, 94행 ICL–위성 결합표에 해시 검증 후 대조했다. 6개 신규 제외 control은
모두 94행에 없다. 그중 target 4803/i2는 ICL 지원만 있고 위성 프로파일이 없고,
target 4749/i2는 위성 프로파일 1개만 있고 ICL 지원이 없다. 나머지 4개에는 어느
쪽 지원도 없다. 따라서 이 6개 후보 control flag는 기존 결합 상관표 행을 제거하지
않는다. 재현 감사기
`analysis/icl_wgl_membership/hsc_debug/audit_s23b_candidate_overlap_association_support.py`
의 결과는 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_s23b_association_overlap_2026-10-22_03/summary.json`
(SHA-256 `51c979d9bfe698fbd94036579b0b82997c57cb0d535e20f54fa69495a459b532`)이다.
분석 스크립트 SHA-256은 `f43d3057253132a8308e259ab1f6c6d200849f1ce62d1e0f8ced95d8714781ce`다.
이는 저장된 표만 읽으며 새 영상·측광·연관 적합을 수행하지 않는다.

### 2026-10-22 후속: S23B 후보 HSC-I metadata coverage 조회안

앞선 S23B 공간 장부는 i2 결합표 밖의 후보 target 5도 셀이 최대 67개라고 보였지만,
새 후보들의 실제 HSC 영상 지원은 확인되지 않았다. 이를 확인할 PDR3 Wide metadata-only
계획을 `analysis/icl_wgl_membership/hsc_debug/prepare_s23b_icl_coverage_plan.py`로
준비했다. S23B 후보 1,333개 각각에 중심과 proper 1,000 kpc 원주 8점을 놓아 총
11,997 위치를 27개 SQL 배치로 검사한다. 열은 candidate/probe key, skymap/tract/patch,
HSC-I mosaic 크기·seeing·zeropoint, i-band patch QA depth/seeing 및 image-product
flag 등 14개로 고정했다. 각 probe 최대 4 patch라는 상한으로 예상 metadata 출력은
47,988행이며 영상, mask pixels, object/source rows, shape, photo-z, PSF FITS는 요청하지
않는다.

기존 2,596-candidate S20A 계획은 별도 승인 scope이며 저장 manifest는 여전히
`prepared_not_submitted`다. 이 새 S23B 계획도 `prepared_pending_user_approval`이며
CAS에 연결하거나 질의하지 않았다. SQL 27개 및 총 29개 artifact hash와 metadata-only
allowlist를 검증했고 focused tests 2개가 통과했다. 계획 JSON SHA-256은
`e615ee7c2414109af6d3c156104d71766355d3b4e2a5dba3636a0142995557d7`이다. 상세 범위,
제한과 재현 명령은
`analysis/icl_wgl_membership/hsc_debug/icl_s23b_coverage_metadata_plan_2026-10-22.md`에
있다. 새 S23B CAS 조회는 사용자의 해당 정확한 범위 승인 후에만 진행하며, 이 계획은
이미지 다운로드를 승인하지 않는다. 공통 submitter는 두 계획 모두 `--check`로 검사하며,
S23B 제출에는 정확한 계획 SHA-256을 `--approved-plan-sha256`으로 다시 전달하도록
가드했다. coverage 후처리기도 `S23B:` release-prefixed ID를 지원한다. 변경 후 focused
submitter/plan/aggregation 검사 10개가 통과했다.

### 2026-10-22 후속: S23B ICL–위성 후보의 좌표·공간 셀 장부

기존 S23B 후보 1,333개를 현재 800개 표본의 좌표 및 실제 ICL–위성 분석의
94행(i 26, i2 68)과 비교하는 오프라인 장부기
`analysis/icl_wgl_membership/hsc_debug/build_s23b_spatial_candidate_ledger.py`를
추가하고 실제 자료로 실행했다. 서로 다른 CAMIRA 릴리스의 숫자 ID를 결합 키로
사용하지 않고 좌표를 비교했다. 후보는 100개의 5도 RA/Dec 셀에 있으며 현재 parent
800개 표본과 1초각 이내 좌표 일치는 0개다. 연관표에는 i 13개, i2 33개 target 셀이
있고 후보는 이 중첩을 제외한 셀을 각각 87개, 67개 차지한다. 따라서 i2는 카탈로그
기하 상한으로 현재 33개 target 셀에서 최대 100개까지 확장 가능하다. 이 수치는 독립
공분산 블록이나 usable 표본이 아니며, nominal 92셀 설계선은 기하 상한으로 넘을 수
있어도 보수 180셀 설계선에는 80개가 모자란다. 신규 영상 coadd, mask, depth, PSF,
위성 완전도, shape 및 matched-control 지원은 미검증이다. 기존 target 영상 WCS에는
후보 중심이 없고 matched-control 영상에는 7개가 있어 오염 확인 후보로 남긴다.
재현 결과 TSV SHA-256 `ab8df9fce1435d1d9eca15f742f43356cbd430b2159852a14d9b67963d83b775`,
요약 JSON SHA-256 `f6ac867881d8184d02a580e9a18540b7513a90f530e51e6cdf04205db535e8d1`이다.
세부 방법과 입력 해시는
`analysis/icl_wgl_membership/hsc_debug/icl_s23b_spatial_candidate_ledger_2026-10-22.md`에
기록했다. 단위시험 2개를 통과했다. HSC 질의·다운로드·픽셀 판독은 없고
`science_measurement=false`를 유지한다.

### 2026-11-03 후속: 전체 ICL 영상 sky 제품 출처 감사

현재 HSC 자료가 ICL에 적절한 global-sky 제품인지 전체 매니페스트와 완료기록에서
다시 확인했다. Astrafek 재현 감사기 `audit_icl_sky_product_provenance.py`가 표적
800개와 matched-control 798개 전부 `sky=global-sky`, `product=coadd/bg`,
`product_reviewed=true`임을 검사했다. 파일 1,598개가 모두 존재하고 완료기록의 바이트
크기와 일치한다. 이번 확인에서 파일 전체 SHA나 픽셀은 다시 읽지 않았다. 결과는
GPFS의 `icl_sky_product_provenance_2026-11-03_01/summary.json` (SHA-256
`235c4f5cedd9506eb874c9cfe0ac09507a71b1a5f729b1816779b63da4f2c839`)이다.

HSC PDR3 공식 처리 설명에서 `deepCoadd/` patch는 global-sky coadd, `deepCoadd-results/`
의 `calexp`는 source detection/측광용 local-sky 제품으로 구분된다. 저장된 cutout manifest도
`coadd/bg`를 aggressive background subtraction 이전 `deepCoadd` 제품으로 기록한다. 따라서
현재 ICL 스택은 확장광에 맞는 global-sky 쪽을 사용하고 있으며 local-sky calexp를 쓰지는
않는다. 그러나 global-sky도 노출별 sky 배경 및 잔류 local pattern 처리 후의 coadd이므로
sky-free는 아니다. 이번 감사는 제품 선택만 확인했고 diffuse ICL의 upstream transfer는
측정하지 않았다. 다음 sky 관문은 현재 global-sky 자료에서 target-control 공간 잔류와
고정 ICL/위성 가중치의 공분산을 검증하는 것이다. 관련 테스트 3개 통과; 과학 상태는
`science_measurement=false` 유지. 상세는
`analysis/icl_wgl_membership/hsc_debug/icl_hsc_sky_product_provenance_assessment_2026-11-03.md`.

### 2026-11-03 후속: 전체 표적–대조장 외곽 residual 진단 완료

783개 same-filter matched target/control pair의 기존 global-sky coadd/bg 영상에서
600–900 kpc 타일 통계를 같은 64-pixel plane-fit 방식으로 비교하는 Astrafek 도구
`analysis/icl_wgl_membership/hsc_debug/audit_icl_target_control_outer_sky.py`를 추가했다.
표적은 저장된 K=2.5 mask, 대조장은 기존 SEP/Kron-2.5 절차를 쓰며 입력 파일은 변경하지
않는다. 이 고리는 표적 외곽 ICL/은하광·PSF 날개·cirrus를 포함할 수 있으므로 intercept를
순수 하늘로 보거나 측광에서 빼지 않는다. 결과는 기존 74개 ICL–위성 표본에서 그 paired
outer-residual 진단값을 공변량으로 넣은 부분상관 민감도도 계산한다.

기술 사전실행에서 i2 target 2269 / matched control `control_2269_43150080429003818`은
각각 156/291 유효 타일을 만들었고 단위·마스크 경로가 정상 동작했다. 이는 한 쌍의 입력
검증이지 과학 결과가 아니다. 전체 계산은 syntax Slurm job 448403, CPU 1개·16 GiB로
26:08에 정상 종료했다. 783쌍( i 302, i2 481), 총 1,566 frame row, connected 5도 block
38개다. target-control 600–900 kpc plane intercept median 차이는 i −0.055
[q16,q84 −3.942,+3.816], i2 −0.370 [−4.863,+3.596] nJy/arcsec²다. 타일 잔차 MAD/명목
평균 sigma 중앙값은 양쪽 이미지 모두 약 3.0–3.2다. 이것은 실제 외곽광·PSF·cirrus·잔류
천체와 마스크를 포함할 수 있는 진단량이지 순수 sky나 보정 배율은 아니다. 파일은
`icl_target_control_outer_sky_2026-11-03_01/`에 있다. 신규 HSC 접속·질의·다운로드는 없다.

처음에는 74개 ICL–위성 연관 표본에 전체 783쌍 그래프의 연결 block label을 재사용해
연관 표본의 공간 영역이 일부 합쳐지는 문제가 발견됐다(초기 요약의 i2 20, i 8 block).
pair/frame 영상 진단은 유효하지만 초기 bootstrap 구간은 폐기했다. 표본 내부에서 filter별
block을 다시 계산하자 i2 53쌍/28 block, i 21쌍/13 block이다. 수정 partial r은 i2 기준
0.239→0.252, 차이 +0.013 [−0.034,+0.111]; i 0.112→0.107, 차이 −0.005
[−0.137,+0.106]이다(각각 5,000 block resample). 차이 구간은 0을 포함하고 outer
intercept가 물리적 cluster light를 포함할 수 있으므로 보정·검출로 해석하지 않는다.
정정 결과와 독립 검증은
`analysis/icl_wgl_membership/hsc_debug/icl_target_control_outer_sky_assessment_2026-11-03.md`에
기록했다. 원·수정 산출물 모두 `science_measurement=false`.
완료 산출물용 독립 결과 검증기
`analysis/icl_wgl_membership/hsc_debug/verify_icl_target_control_outer_sky.py`도 추가했다.
저장 manifest에서 same-filter pair 선택·좌표를 대조하고 plane 차분, 연결블록,
deterministic bootstrap을 재계산한다. 정적 컴파일 및 diff 검사는 통과했으며 실제 산출물
대조는 Slurm 계산 완료 후 수행한다. 검증기 commit은 `45c7408a9`다.

### 2026-11-03 후속: 약중력 공간 잭나이프 격자 크기 민감도

기존 800 렌즈 HSC Y3 형상 행과 SACC `n(z)`를 재사용해 equal-area 공간 잭나이프를
2.5도·5도·7.5도 셀 크기에서 비교했다. CPU-only Slurm 작업 448324/448323은 `syntax`에서
각각 28:54/7:09에 정상 종료했다. 2.5도·5도·7.5도 결과는 점유 셀 156/52/36개, 사용
lens-source 쌍 2,153,004개로 동일했다. 중심 `DeltaSigma` 프로파일도 모두 같았고, 대각선
잭나이프 오차는 5도 대비 2.5도에서 0.829–1.167배, 7.5도에서 0.878–1.177배였다. 교차성분
null 진단은 각각 chi2/dof=7.10/7, 9.30/7, 6.21/7 (p=0.418, 0.232, 0.515)다.

따라서 중심 프로파일은 이 격자 변화에 안정적이지만 대각 오차는 약 18%까지 변한다. 격자
셀은 HSC tract/독립 필드가 아니므로 이 비교는 근사 공분산 진단이고 paired
cluster-minus-random 공분산을 대체하지 않는다. raw B-mode 제외값 0, 고정 fiducial `n(z)`,
미적용 reduced-shear 및 랜덤 차감 등 기존 제한을 유지한다. random-center shape-count 신규
집계 쿼리는 별도 명시 승인 대기 상태이고 제출하지 않았다. 과학상태 `science_measurement=false`.
상세 결과와 checksum은 `analysis/icl_wgl_membership/hsc_debug/wl_spatial_jackknife_cell_scale_assessment_2026-11-03.md`.

### 2026-10-22 후속: S20A·S23B 후보 footprint 합집합

기존 승인 S20A wide-redshift candidate plan과 S23B selection-compatible ledger를
좌표 격자에서 합쳤다. 기존 i2 discovery 33 target 셀 기준 S20A pool은 51개 새 셀,
S23B는 67개 새 셀, 둘의 합집합은 69개 새 셀(총 candidate 셀 102개)이다. 따라서
S20A 계획은 S23B에 새 셀을 2개만 더하며, 추가 표본은 BCG spectroscopic selection도
달라진다. r=0.293 단일 고정 endpoint 80% 표본력 92 blocks에 union 상한은 23개 부족하고,
30비교 조정 180 blocks에는 discovery 33셀을 더해도 78개 부족하다. 이는 실제 usable/
independent block이 아닌 기하 상한이다. 계산기
`analysis/icl_wgl_membership/hsc_debug/assess_icl_cross_release_replication_capacity.py`,
평가 `analysis/icl_wgl_membership/hsc_debug/icl_cross_release_replication_capacity_2026-10-22.md`,
focused tests 2개를 추가했다. S20A metadata query도 이 부족을 해결하지 않는다는 점을
확인했으며, 다른 wide-area footprint/candidate source가 필요하다. 새 HSC 요청은 없고
`science_measurement=false`다.

### 2026-10-22 후속: S23B 독립 재현 표본력 관문

S23B candidate pool을 기존 exploratory discovery 표본과 단순 합산해 검정력을 과장하지
않도록 해시 고정 장부에서 재현 표본력 상한을 계산했다. 기존 i2 연관표의 33 target
셀 밖에서 S23B는 최대 67개 새 target 셀을 제공한다. 단일 사전고정 i2 endpoint에서
양측 α=0.05, partial r=0.293을 가정한 Fisher-z 근사 80% 목표는 92 independent blocks로,
S23B 단독 67개 target-cell 상한의 power는 0.662, 최소 25개 block이 부족하다. r=0.20이면
196개가 필요해 최소 129개 부족하다. 기존 discovery 33셀을 더한 100은 독립 재현 표본으로
셈할 수 없다. 영상·마스크·위성·control 필터링과 block 병합은 이 상한을 더 줄일 수 있다.
계산기 `analysis/icl_wgl_membership/hsc_debug/assess_s23b_replication_power.py` 및
평가 `analysis/icl_wgl_membership/hsc_debug/icl_s23b_replication_power_gate_2026-10-22.md`,
focused tests 5개를 추가했다. 새 HSC 요청이나 science measurement는 없고
`science_measurement=false`를 유지한다.

### 2026-11-03 후속: weak-lensing 대조 중심 기하 검사 및 좌표 재사용 민감도

기존에 승인된 random-coordinate export가 32/32 배치, 780,191행까지 완료된 것을
`coordinate_status.json`에서 확인했다. export는 800개 렌즈와 `lens_id, ra_deg,
dec_deg` 세 열만 포함한다. syntax Slurm job 448283이 기존 CAMIRA S20A 카탈로그와
렌즈 적색편이를 사용한 오프라인 3 proper-Mpc 겹침 검사를 31초 만에 정상 완료했다.
기하 표 780,191행 및 렌즈별 800행 집계가 summary와 재계산상 일치한다. summary는
`full_aperture_footprint_verified=false`, `hsc_i_depth_psf_sky_matched=false`,
`science_measurement=false`다.

원래 screen은 random coordinate가 전체 렌즈 연결표에서 단 한 번만 나타나야 후보로
인정했다. 이는 겹침 검사를 각 lens-coordinate 쌍에 적용해도 되는 상황에서 다른 렌즈의
중복을 함께 버리는 과도하게 보수적인 선택 규칙일 수 있다. 재현 감사기
`audit_wl_random_coordinate_reuse.py`는 원래 전역 고유 기준과 렌즈별 겹침 판정 기준을
비교한다. 원 기준은 68,352 후보/544 렌즈(5개 이상은 483렌즈)다. 렌즈 간 같은 위치를
재사용하되 각 렌즈 안 중복은 배제하면 95,847개 연결/667 렌즈(5개 이상 620, 10개 이상
577)가 남는다. 이는 82,366개 고유 좌표이며 11,514개 좌표가 여러 렌즈에 공유되고 최대
5개 렌즈가 같은 지점을 사용한다. 재사용은 별도 독립 표본이 아니므로 최종 추정에서는
공유 위치를 공간 블록 공분산에 넣어야 한다.

예비 count-only 계획을 위해 렌즈별 기하 후보가 5개 이상인 620개 렌즈에서 SHA-256 순위로
렌즈당 5개, 총 3,100개 중심을 선택했다(3,020개 고유 위치; 77개 위치가 여러 선택 렌즈에
공유). HSC S19A의 shape catalog 세 테이블에서 각 3 proper-Mpc 원뿔 안 `COUNT(DISTINCT
object_id)` 정수만 반환하도록 25개씩 124 SQL batch를 생성했다. 열은
`random_center_id, matched_lens_id, n_shape_sources`로 제한된다. 계획은 GPFS의
`wl_random_shape_count_plan_lens_reuse_2026-11-03_01`에 `prepared_not_submitted` 상태로
있으며 HSC CAS 연결/제출은 하지 않았다. source rows, shape, photo-z, 이미지, 새 좌표는
포함하지 않는다. 이는 target 렌즈에 이미 받은 count-only 승인과 별개의 random-center
범위이므로 사용자 명시 승인 전 제출하지 않는다. 계획 SHA-256은
`18aeb1e831d847818b2a671dfe5f75c3122c38637c91b9abbd8fabff5faae625`이다.
재현·한계는 `analysis/icl_wgl_membership/hsc_debug/wl_random_coordinate_reuse_assessment_2026-11-03.md`에,
수정기 및 단위시험은 같은 폴더의 `prepare_wl_random_centers.py`,
`test_prepare_wl_random_centers.py`에 기록했다. 관련 audit, selector, shape-count plan 및
submitter focused tests 30개 통과. 무작위 대조장은
아직 full-aperture coverage, mask/depth/PSF/sky matching, shape selection, spatial covariance를
통과하지 않았으므로 weak-lensing 결과는 계속 `science_measurement=false`다.

### 2026-11-03 후속: supported i2 hybrid-PSF grid 독립 검증

Slurm job 434556의 i2 44프레임×12 splice configuration 완료 상태를 재확인했다.
Astrafek 독립 summarizer가 총 528 frame-model 기록의 응답 산술을 다시 계산했으며,
모든 configuration이 동일한 44 cluster ID를 사용하고 최대 산술 차이는 0이었다.
검증 결과는 GPFS의 `icl_core_wing_hybrid_i2_supported44_2026-10-22_01/grid_audit_repro_2026-11-03.json`
(SHA-256 `eaf08ddc0cbaec3243f3e13935bdd5c960b1acdf450f682ca94f44fcbce13712`)에 있다.

고정 마스크 smooth-diffuse→ICL 중앙 회수율은 PSF grid 전역에서 0.225702–0.225874로
매우 안정적이다. 반면 redetection satellite-envelope→satellite 응답은 대표 설정의
q16–q84가 −0.613–1.187, 전체 범위 −5.92–10.78로 극단적이다. 위성광→ICL redetection
중앙값도 고정 마스크보다 0.088–0.110 낮다. 따라서 마스크 재선택이 위성 관련 전달을
크게 바꾸며 scalar 보정은 불가하다. 200–300 kpc annulus mask coverage 중앙값은 0.620이고,
별도 noiseless 1 nJy 점광원 모델의 unmasked-annulus 광량 비율은 0.000724–0.001411이다.
이는 위치/epoch가 맞지 않는 공개 HSC-I wing proxy와 선택된 44프레임에 대한 합성 주입
진단이지 실제 i2 wing 측정이 아니다.

평가와 재현 조건은 `analysis/icl_wgl_membership/hsc_debug/core_wing_hybrid_i2_supported_assessment_2026-11-03.md`에
기록했다. source Poisson noise, coadd covariance, real sky transfer, matched-control residual,
full-sample completeness는 아직 포함하지 않았다. 다음 ICL 관문은 target/control 동일 마스크·sky
estimator를 통과시키는 noise-aware injection과 공간블록 불확실도이며,
`science_measurement=false`, `psf_approved=false`,
`observational_validation_passed=false`, `flux_covariance_available=false`를 유지한다.

### 2026-11-03 후속: 약중력 CAMIRA 풍부도 분할 진단

기존 800-lens shape-row export를 재검증하고, HSC profile을 보기 전에 부모 sample manifest의
CAMIRA 풍부도와 lens ID로 표본을 400개씩 균형 분할했다(경계 20.490/20.498). 동일한 HSC Y3
`n(z)`, raw B-mode 제외값 0, 7개 반경 bin 및 5도 equal-area 잭나이프를 사용해 공통 52개
공간 셀을 삭제할 때마다 양쪽 stack과 보정을 다시 계산했다. lageunha에서 CPU로 실행했고
신규 HSC 접속·쿼리·다운로드는 없었다.

사용된 lens는 저풍부도 299개/1,043,918쌍, 고풍부도 310개/1,109,086쌍이다. 상위−하위
`DeltaSigma`는 0.1–0.2 Mpc에서 `15.13 ± 21.71`, 0.2–0.4에서 `41.58 ± 12.26`,
0.4–0.7에서 `35.38 ± 8.71`, 0.7–1.0에서 `27.92 ± 6.76`, 1.0–1.5에서
`15.32 ± 4.88`, 1.5–2.0에서 `17.23 ± 3.56`, 2.0–3.0 Mpc에서 `9.25 ± 2.12`
이며 단위는 `10^12 M_sun Mpc^-2`다. 조건부 Gaussian 7-bin 차이 진단은 `chi2=82.66`,
rank 7, nominal `p=3.95e-15`이나 공분산과 approximate sky blocks에 민감하므로 검출
유의도로 해석하지 않는다. 그룹별 부모 적색편이 중앙값은 0.442와 0.427이다. 랜덤 중심
차감, redshift/selection 보정, shear·`n(z)` 불확실도, reduced shear 및 검증된 survey-region
공분산이 빠져 있어 `science_measurement=false`다. random-center count-only 쿼리도 이 실행에서
제출하지 않았다.

상세 평가 `analysis/icl_wgl_membership/hsc_debug/wl_richness_split_assessment_2026-11-03.md`,
계산기 `analysis/icl_wgl_membership/hsc_debug/audit_wl_richness_split.py`; 요약 JSON SHA-256
`20f28c5cdec8a8a6ca9d3041bfcdad352999d1a6e4b3e809eb6bc5e6e1827a33`이다.

### 2026-11-03 후속: ICL–위성 공분산 지원범위 정정과 공유 고리 민감도

기존 63쌍 공분산 실행은 `prestack_flux_profiles` 137행 공통-stack table을 ICL 지원 입력으로
사용해 고정 74쌍 연관표본 중 22쌍을 빼고 그 표본 밖 11쌍을 넣었다. 이는 annulus별 지원표
범위에 맞지 않아 이전 63쌍 평가는 역사적 subset 결과로 표시하고 새 실행으로 교체했다.
지원 선택기의 기본 입력을 `icl_annulus_support_2026-10-18_01`로 수정했다. 실제 200–300 kpc
지원표 224 pair를 40% ICL coverage로 선별하고, 위성 50–75 kpc 양쪽 60% coverage와 대조점
3개 기준을 적용하자 94 pairs/235 profiles가 선정됐다(26 i, 68 i2). 고정된 74 association
pairs(21 i, 53 i2)는 모두 포함되고 누락은 0이며, 바깥 추가 지원 pair는 20개다.

실제 마스크 유효 픽셀을 비교한 206개 association 위성 profile 중 140개(68%)가 ICL ring과
공유 pixel을 가진다(29/41 i, 111/165 i2). 지원 위성 target 50–75 kpc pixel 합집합을
200–300 kpc target ICL ring에서 제외해도 paired 측정값을 직접 재구성했다. 5,000회 연결
5도 블록 재표집에서 i2 공통 coverage 표본은 52개가 남고 부분상관이 0.217에서 0.176으로
낮아진다. 짝지은 변화 95% 구간은 [−0.065, −0.018]이다. 전체 filter subset 53개에서 0.234→
0.193, 변화 구간 [−0.064, −0.021]이다. 이는 동일 픽셀 재사용이 관측 관계 일부를 만들 수
있음을 보여주는 aperture-exclusion 민감도이지 correction/detection은 아니다. ICL coverage
40% 아래로 내려가는 쌍이 filter별 하나씩 있어 공통지원 재분석을 함께 기록했다.

표적+matched-control 분산을 포함한 94행 공분산은 독립 formula/provenance 검증을 통과했고,
새 exact-support verifier는 output pair set이 independently recomputed 94개와 같으며 74개
association pair를 전부 덮음을 확인했다. 분석은 기존 HSC 컷아웃만 썼고 신규 query/download는
없다. 상세 보고서는 `analysis/icl_wgl_membership/hsc_debug/icl_satellite_full_support_overlap_assessment_2026-11-03.md`,
지원범위 검증기 `verify_satellite_icl_covariance_support_scope.py`, noise-error 진단기
`assess_icl_satellite_pixel_noise_bias.py`, shared-shell 민감도기
`measure_icl_satellite_overlap_sensitivity.py`다. 최신 sensitivity summary SHA-256
`733433f2495ef4cf1592d7885a9845de57c726b6b64a2fdf93ca91514a0260e9`; 모든 산출물은
`science_measurement=false`를 유지한다. 다음 ICL 단계는 위성-shell 제외 profile을 포함한
joint injection/recovery와 큰 규모 sky·PSF wing 공간 블록 공분산이다.

### 2026-11-03 후속: 약중력 풍부도 차이의 적색편이 매칭 민감도

기존 CAMIRA 풍부도 상·하위 400개 cluster-only HSC Y3 스택의 적색편이 중앙값은
0.442/0.427이었다. 부모 표본에서 shear profile을 보지 않고 `|Δz|≤0.02` 1:1 최대매칭을
구성해 392쌍을 얻었다. 승인된 shape-row export에서 양쪽 렌즈 모두 source 행이 있는 쌍만
남기자 251쌍으로 줄었다(141쌍은 한쪽 또는 양쪽 export 지원 부족). 매칭된 실제 표본의 평균
`|Δz|=0.00720`, 중앙값 0.00253, 최대 0.01999이며 그룹별 평균 z는 0.4259/0.4198이다.

공통 49개 5도 등면적 하늘 셀 paired jackknife에서 상−하 풍부도 `DeltaSigma` 차이는
0.1–0.2 Mpc에서 `35.28±24.90`, 이후 0.2–0.4, 0.4–0.7, 0.7–1.0, 1.0–1.5,
1.5–2.0, 2.0–3.0 Mpc에서 각각 `40.73±14.38`, `32.79±8.74`, `24.67±6.63`,
`15.31±5.56`, `17.61±3.57`, `9.59±2.28`이다(단위 `10^12 M_sun Mpc^-2`).
0.2–3 Mpc의 차이는 원래 400/400 split과 유사하다. 이는 lens z 한 변수만 맞춘 cluster-only
robustness 진단이며 random-lens subtraction, survey 연결영역 공분산, selection/shear/`n(z)`
불확실도 및 reduced shear가 빠져 있다. nominal profile Wald 값 `χ²=62.66`, rank 7,
`p=4.44e-11`은 검출 유의도가 아니다. 새 CAS 요청·자료 다운로드 없이 lageunha CPU에서
실행했으며 `science_measurement=false`를 유지한다.

재현기 `analysis/icl_wgl_membership/hsc_debug/audit_wl_redshift_matched_richness_split.py`,
평가 `analysis/icl_wgl_membership/hsc_debug/wl_redshift_matched_richness_assessment_2026-11-03.md`.
summary는 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_redshift_matched_richness_2026-11-03_01/summary.json`,
SHA-256 `f050c2c90236a4215d951cb939bb65c145b98d6526dea90a88007beff2f152fb`이다. random-center
count-only plan은 여전히 `prepared_not_submitted`; 본 분석은 제출하지 않았다.

### 2026-11-03 후속: ICL–위성 공유 고리 제외 주입 전달

고정 74개 ICL–위성 연관표본과 joint-annulus 지원 94쌍, 승인된 compact PSF 위치를 교차해
정확한 association/PSF 교집합 48쌍만 실행하도록 batch runner를 확장했다. 실제 검증된
74행 association table SHA-256 `483f81540543fc46ec1cf4333bbd3e238799b493e48fe03235da3038f5714917`
및 independent audit를 selector가 확인한다. 48프레임(12 i, 36 i2)의 관측 mask/segmentation에서
지원된 위성 50–75 kpc target pixels를 합쳐 ICL 200–300 kpc에서 제외하고, diffuse, satellite,
joint synthetic component를 local HSC-I model PSF core로 주입했다. 실제 coadd는 변경하지 않았다.

Diffuse→ICL 회수 median은 i `0.22564→0.22491`, i2 `0.22594→0.22623`으로 안정적이다.
위성 envelope→ICL marginal median은 i `0.46709→0.46624`, i2 `0.68125→0.65981`이다.
하지만 각 프레임의 짝지은 response 차 q16–q84는 i `[-0.01361, +0.00287]`, i2
`[-0.06558, +0.01211]`로 0을 포함하고 프레임 분산이 넓다. 고리 제외 coverage 손실 중앙값은
i 2.46%, i2 3.70%다. 따라서 주입 전달 scalar correction을 쓰지 않는다. 57-frame broader
joint-support batch에는 association 밖 9개가 있어 이 association 결과의 기준으로 삼지 않았다.

재현 옵션은 `run_joint_icl_satellite_batch.py --association-pairs .../association_rows.tsv
--exclude-supported-satellite-shells-from-icl`이다. 단일 프레임 계산은
`joint_icl_satellite_injection.py`에서 실제 supported profile 픽셀 수와 기존 segmentation을
대조한다. batch summary `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_shell_exclusion_association_psf48_2026-11-03_01/summary.json`
SHA-256은 `2f6c3d2ae393dd14532d640a1f9710ddee3064f8daf4feaca8f63c57aaf1821d`다. 검토와 한계는
`analysis/icl_wgl_membership/hsc_debug/icl_satellite_shell_exclusion_transfer_assessment_2026-11-03.md`에
기록했다. PSF compact core만 사용했고 source Poisson noise, 확장 PSF wing, coadd covariance,
상류 sky, 마스크 재선택과 matched-control/공간공분산을 포함하지 않았다. 결과는 계속
`science_measurement=false`, `psf_approved=false`, `observational_validation_passed=false`,
`flux_covariance_available=false`다.

### 2026-11-03 후속: 확장 PSF wing + ICL 위성 고리 제외 민감도

앞의 정확한 association/PSF 교집합 48프레임(12 i, 36 i2)에 공개 HSC-I 확장 wing을 compact local
PSF core와 혼합해 같은 위성 50–75 kpc shell-exclusion 민감도를 실행했다. 매개변수는 splice 2.4″,
blend 0.6″, wing amplitude 0.65다. 두 배치 모두 오류 0이며 독립 집계기가 프레임별 SHA-256,
74행 고정 association membership, 48개 유일 target/filter key를 확인했다.

기존 고정 마스크 대비 diffuse→ICL paired 변화 중앙값은 i −0.00018, i2 0.00000이고,
위성 envelope→ICL은 i −0.00022, i2 0.00000이다. 네 조합 모두 프레임 차이 q16–q84가 0을
포함해 상수 응답 보정을 뒷받침하지 않는다. 단위 점원의 ICL ring 내 모델 wing flux fraction은
중앙값 i 0.000994→0.000902, i2 0.000970→0.000896으로 낮아졌다. 다만 유효 ring coverage도
중앙값 2.46 및 3.70 percentage points 감소한다. i2는 필터별 wing이 없는 상태에서 HSC-I proxy를
쓴 탐색 민감도다. 매칭된 target-specific PSF, source Poisson noise, full coadd/sky/control covariance,
관측 마스크 재선택은 아직 없다.

결과 해석·모든 수치는 `analysis/icl_wgl_membership/hsc_debug/icl_satellite_shell_exclusion_hybrid_psf_assessment_2026-11-03.md`,
독립 해시 및 표본 감사기는 `summarize_icl_hybrid_shell_exclusion.py`에 있다. i/i2 summary SHA-256은
각각 `ee3c1ae31b6ae2cef099702f071899e6b54a652935d8456f22893b37f476ec9e`와
`56e628cfcf39fbdf80d499376d158b0056bbafeab7de8a648c9e0b327b7ce7b3`이다. 새 CAS 요청·다운로드는
없었다. 계속 `science_measurement=false`, `psf_approved=false`,
`observational_validation_passed=false`, `flux_covariance_available=false`다.

### 2026-11-03 후속: 약중력렌즈 무작위 중심 형상원 수 사전 집계 완료

사용자가 무작위 중심 3,100개 count-only 조회를 명시 승인했다. 검증된 계획은 620개 matched lens,
lens당 5개 중심, proper 3 Mpc, 25개씩 124 SQL batch다. 반환 열은 `random_center_id`,
`matched_lens_id`, `n_shape_sources`이며, 고유 S19A 형상 catalog object 수만 집계한다. shape row,
개별 좌표, photo-z, 이미지 다운로드는 포함되지 않는다. 실행 전 `--check`는 network 사용 없이
계획/SQL/checksum을 확인했다. plan SHA-256은
`18aeb1e831d847818b2a671dfe5f75c3122c38637c91b9abbd8fabff5faae625`, random-center source manifest SHA-256은
`728c8e7ae25b961a6d68f4dfff1ac37913535ef20ee074b48b999b6a9bb2f133`이다.

lageunha에서 기존 `submit_wl_random_shape_counts.py --submit --approved-random-shape-count-query`로
제출했다. 실행기는 CAS 상태를 확인하며 기존 SQL과 제출 job ID를 바꾸지 않았다. 124개 batch 전부
완료됐고, 각 결과 파일은 지정된 25개 중심 수와 SHA-256으로 검증됐다. 최종 checkpoint 경로는
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_shape_count_results_lens_reuse_2026-11-03_01/count_status.json`이다.
checkpoint 최종 status와 124개 job 상태는 모두 `complete`다. 과거 단말 중단 때 기록된
`error_type=KeyboardInterrupt` 필드는 checkpoint에 남아 있지만 완료된 배치의 결과 상태를 바꾸지 않는다.
승인한 3,100개 count row가 반환됐으며 HSC 비밀번호는 checkpoint나 결과 파일에 저장되지 않았다.
source shape row, 개별 좌표, photo-z, 이미지는 반환되지 않았다. 별도 승인 전 원자료 행 반출은 하지 않는다.

기존 승인 렌즈 위치 count-only preflight(800개, 3,069,544 후보)와 전체 random-center 결과를
해시 검증해 paired support를 비교했다. 전체 800개 사전검사 중 lens 위치 count가 0인 렌즈는 168개지만,
random 5개가 모두 완료된 matched lens 620개 안에는 150개가 포함된다. 이 150개 중 131개는 다섯
random 중심도 모두 0이었다. lens count-positive 470개 중 385개는 다섯 random 중심 모두 양수,
35개는 전부 0이었다. 470개에서 random 5개 평균/lens 위치 후보 수 비의 중앙값은 0.878,
q16–q84는 0.419–1.110이다. Random 중심 3,100개 중 1,011개(32.6%)는 0 count이고 합계는
9,010,324이다. 이는 후보 수 비교일 뿐 shape 품질·tomography·B-mode·방사상 절단 이후의 지원,
매칭 보정 또는 렌즈 신호 측정이 아니다. 따라서 raw random-zero 비율을 단독 random-control 실패로
간주하지 않고, 이후 동일한 source selection을 두 endpoint에 적용해 다시 평가한다.

최종 감사 결과 경로는
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_shape_count_support_audit_2026-11-03_01/lens_paired_support_audit_ratio_v25.json`,
SHA-256 `76c7b87f5c099ad2770b56b39aeac0e324ea293e3e644b7495d244fa4bf0c1c4`다. 재현 감사기와 해석은
`analysis/icl_wgl_membership/hsc_debug/audit_wl_random_shape_count_support.py` 및
`analysis/icl_wgl_membership/hsc_debug/wl_random_shape_count_support_assessment_2026-11-03.md`에
기록했다. `science_measurement=false`는 유지한다.

무작위 중심 shape-row export plan은 `submit_wl_random_shape_rows.py --check`로 오프라인 검증한 뒤
사용자가 정확한 계획 SHA-256 `fc883444d2b23663d5f6034d773b76db948add29f1cda0c23384d17cfed2996a`와 범위를
승인해 제출했다. 범위는 620개 matched lens, lens당 5개 무작위 중심, proper 3 Mpc, 124 SQL batch,
예상 9,010,324 `(random_center_id, object_id)` 행이다. 제안 23개 열은 source 좌표, ellipticity,
shape weight/calibration, aperture magnitude, resolution, PSF shape, raw B-mode mask, HSC tomographic bin을
포함하고 개별 photo-z/PDF 및 영상은 제외한다. 비압축 크기 추정은 2.15–4.03 GiB다.

2026-11-03 현재 lageunha CAS 실행은 진행 중이다. batch 000–034 (jobs 1068173, 1068175, 1068176, 1068177, 1068181, 1068182, 1068183, 1068184, 1068185, 1068186, 1068187, 1068188, 1068189, 1068190, 1068191, 1068192, 1068195, 1068202, 1068207, 1068214, 1068215, 1068216, 1068217, 1068218, 1068219, 1068225, 1068276, 1068277, 1068278, 1068279, 1068280, 1068281, 1068282, 1068283, 1068284)은 각각
58,162, 68,595, 21,091, 46,205, 41,681, 80,090, 107,855, 99,798, 109,409, 64,912, 75,494, 91,881, 93,850, 97,237, 104,882, 101,084, 145,275, 51,720, 108,181, 88,754, 58,370, 109,504, 127,802, 67,391, 53,630, 90,493, 83,412, 91,082, 21,775, 37,938, 19,427, 25,185, 16,150, 33,505, 51,059행을 반환해 planned per-center count 대조를 통과했다.
누적 2,542,879행이다. 순서대로 파일 SHA-256은
`97f8dabcce114e2b182a7fa351fd3549ad92cbf4eaf5b6b058de28779a2a37e8`,
`82c45ca2f1dbedbafd5a4aba68c015274a291e1701dfe552c180e24714d98808`,
`f49467d3ab8e140d46e3e5ae137dc3564ac4af620630e06ba1e8082a064750b4`,
`31da89520d8e8b22027844beb348fd55b894e140476b2287a0b1b945dac46ca6`,
`a5b5c28f12fa158f1863f018415d4e6138076f21b9046024861c96f83d97c7cf`,
`b7d00359aad94325151bf6a58084ace301698634c4fbcd24bc88e409d274dcd8`,
`747c7cb1ca744abf6d1a10d17784ca2dfecddffb41c0828d3f72920363489ba3`,
`d51e19b1e6f111d59a70169a00ae17e55a685e0bc3897a00e545630bf25ad01c`,
`936c0dff19a7a0d30a1540be36b0fcea290e012420f86bcfa5737aa6834eea62`,
`28e57d51964a04691b28a37f46b53f62221d710dca1264c1fd0b99988eb0c72c`,
`3179df66db0737754aec22c50bfd79873fe752f063d90d446db0bccfd7cb2815`,
`123922011b04487310af35d89d496f68f27aecd3476a655a34b63d955fbc08e2`,
`dd34988086f5a0a835cefad186a464b65a40375fd0816c4381d287de5d91d1b1`,
`3dba0c0bbc0aba29f4da0f80364f120c98751569e12838600fda078e4ddede0c`,
`795dd7c8484a13dc7fdfb42e554df8f86edec2dfc50fce75c879e3712b73d757`,
`cabc885374d06a717376eca11033b4d748797ed9d7de943a360eec9d0bc0c864`,
`e6dbc23165d3793a3cf7923ddbb960af7c08d81dbb515d0b678d0cafe0d0ea69`,
`3d4ae3921a943b031a4a6f66da506e1298a76391eb2a753d93c439a7f11baf63`,
`c2089218f7160d0177cc08f3f272a6bf38cb1aa964b6d26ad7ed95311f28db6b`,
`0be812934fbc2c95f9ebafae9c6b25cb252b29f47096daaddea393fc0270c219`,
`d72762323140542e626753df6dc0db43f468ff9271a41447842425934ae90389`,
`ea923cf61d97540cd5c172a5d276c6a2560f7e2ed5aac24efed5abdf2f46897f`,
`67c55b67b1ce62518ff2e1ad2fcbe2d6b3789cdc98831cc612df3b77216d700d`,
`7cc8b8066309c29586c8d162f86cc8a648f07fc1e000530b7182a451b4ec6468`,
`e75aa1ff106929117a58b5c0c70e1c5bd3415e434623d03588fca252c9cb49d7`,
`0ef934308ceed24b4dcbe5b71981734f2ad33ded1f016d7f1c69b72866b43243`,
`4380c9d93d5406a8994872c8b4b6481915e8b7d0fcec1e06aa99a2f7679f0c40`,
`944303fd361e380a36aa7175657ee0c2fc5f1ef6d29a1d8279acb6cf6c2550c3`,
`51c7b91936690d38e50de1d55de119f6b9cfad08ca8fa6441272ea8faa31d752`,
`1f9bdce19812bfd380d77de8a9e78df7eccaf0636bbea6434544133b23216bfb`,
`b1ec2011efbd94357ac9edb8dc5beb6a43037aaddef07bc147df8fee380747ae`,
`be5f5c31e6688fd6aa5b2ccaee93a24e633873cdbc4f609ab9f8f9e61ef43cb7`,
`95d696df02d42c4f475e1ef73dcde76a73d402fd585a54f094a24a39c7b0593e`,
`83bc75dc2a9187fbb5414c6376d47552f28d917ee3dd0017f1516a38c366db45`,
`23091616982940c92cf44565596440cd8d0ed01fd156cbe7535995e68830e362`이다.

2026-11-04 재개 후 batch 035–040 (jobs 1068285–1068290)은 각각 38,478, 0, 0, 51,488,
36,726, 0행을 반환했고 전부 preflight center count 대조를 통과했다. 누적 완료는 41/124 batch,
2,669,571행이다. batch CSV SHA-256 순서는
`1c02eb1b2222d0f9c497684779eabead60922fd009f407afbd46c921634a845d`,
`5cab5731ae6624e7c0b5362a9c08cb08e3ab9358935aede924cc96736b594522` (0행),
`5cab5731ae6624e7c0b5362a9c08cb08e3ab9358935aede924cc96736b594522` (0행),
`bb6a89c97df9cbc37e19725cecf04e2d995fb7a0d57a8a953e3ef30bde56d71b`,
`0c988313dce3fba8b508cd751940f3464addb6335871d9f3026f0295222075eb`,
`5cab5731ae6624e7c0b5362a9c08cb08e3ab9358935aede924cc96736b594522` (0행)이다. Batch 041–043
(jobs 1068291–1068293)은 각각 38,732, 12,477, 42,984행을 반환하고 preflight 대조를 통과했다.
현재 누계 44/124 batch, 2,763,764행이다. 세 CSV SHA-256은
`c4b70c6dbe924093808d296f79aa16e2cdc626991f4cba4955296248d82a404b`,
`f410aeb763fdd7499bf57a1fe1de46b656c7a483b79ad1ecfb035e8c6d850c80`,
`8f644d55e42fe254002273003de87454da33dc1255e64f30b43248f632381f02`이다. 이후 batch 044–049
(jobs 1068294–1068303)은 각각 83,154, 49,733, 146,988, 147,935, 75,023, 73,306행을 반환해
preflight 대조를 통과했다. batch 050 (job 1068304)도 114,673행으로 검증됐다. 누계는 51/124 batch,
3,454,576행이다. batch별 SHA-256은 checkpoint의 jobs 항목에 기록되어 있으며 CSV와 최종 대조한다.
Batch 051 (job 1068305)이 실행 중이며, checkpoint는
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_shape_rows_lens_reuse_2026-11-03_01/results/row_status.json`이다.
shape 원자료 행 반출은 기존 count-only 승인과 별도이며 이번 제출은 위에 적은 승인 범위와 일치한다.
실행 완료 뒤 모든 batch의 행 수, center별 preflight count, 파일 해시 및 실제 총량을 확인해야 한다.

Plan 디렉터리는
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_shape_rows_lens_reuse_2026-11-03_01/`,
manifest SHA-256은 `fc883444d2b23663d5f6034d773b76db948add29f1cda0c23384d17cfed2996a`다.
shape row 제출 스크립트는 `submit_wl_random_shape_rows.py`이며 명시적 승인 플래그가 없으면 제출을
거부한다.

완료된 첫 68 batch만으로 한 과거 중간 audit에서는 340개 matched lens에 1,700개 중심이
확인됐고, 합계 source-count는 4,749,572였다. 중심 600개(35.3%)는 count 0이며, 107/340 lens는
5개 중심 모두 0, 211/340은 모두 양수였다. 이는 manifest 순서의 전반부만 본 비최종 결과라
공간 순서/표본 편향 가능성이 있고 전체 지원율로 일반화하지 않는다. 해시 확인된 CSV만 사용했다.
재현기 `analysis/icl_wgl_membership/hsc_debug/audit_wl_random_shape_count_support.py`, 평가와
snapshot 경로는 `analysis/icl_wgl_membership/hsc_debug/wl_random_shape_count_support_assessment_2026-11-03.md`에
기록했다. 위 전체 124-batch 감사가 최종 count-only 후보 지원 비교다.

### 2026-11-03 후속: ICL–위성 sky·coadd 상관잡음 공동 민감도

기존 고정 ICL–위성 연관 74쌍(i2 53, i 21)에 association 표본 내부에서 재계산한 5도 연결 블록을
사용해, 적색편이·풍부도만 조정한 기준모형과 외곽 600–900 kpc target−control plane residual,
local 9×9 lag-window coadd sigma-inflation target−control 차이를 각각 및 함께 추가한 네 모형의
partial Pearson 상관을 비교했다. 입력은 저장된 outer-sky 진단과 독립 산술검증을 통과한 전체 표본
coadd residual-covariance frame 표다. 필터별 5,000회 paired block bootstrap을 수행했고, 74행 join 및
모든 모델 점추정치를 별도 직접 최소제곱으로 재현했다.

i2의 기준 r=0.239 [0.022, 0.423], outer residual 조정 r=0.252 [0.036, 0.456], coadd proxy 조정
r=0.211 [−0.006, 0.403], 두 변수 모두 추가한 r=0.223 [0.009, 0.440]이다. 두 proxy를 함께 넣은
기준 대비 paired 변화는 −0.016 [−0.091, +0.083]이다. i는 기준 r=0.112 [−0.229, 0.514], joint 조정
r=0.098 [−0.374, 0.727], 변화 −0.014 [−0.300, +0.272]로 불확실성이 크다. 대괄호는 spatial-block
bootstrap percentile 95% 구간이며 검출 유의도가 아니다.

따라서 선택된 i2 표본의 point estimate 방향은 두 진단 공변량 조정에 대체로 안정적이지만, 이는
ICL–위성의 물리적 연관을 확증하지 않는다. 국소 lag-window proxy는 방사 ICL estimator covariance가
아니며 outer annulus도 순수 sky가 아니다. 상류 HSC sky 전달, target-specific extended PSF,
위성 완전도, full image-pixel covariance 등은 여전히 미보정이다. 분석기
`analyze_icl_satellite_joint_sky_covariance.py`, 평가
`analysis/icl_wgl_membership/hsc_debug/icl_satellite_joint_sky_covariance_assessment_2026-11-03.md`.
summary `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_joint_sky_covariance_2026-11-03_01/summary.json`
SHA-256 `947d9ecd8bb728e5a3bb1ec193f88687a893752403c9c7a4d409a3e6752888e2`; 계속
`science_measurement=false`다.


### 2026-11-03 후속: ICL 200–300 kpc 고유 aperture lag-covariance

실제 저장된 target/control aperture 유효 픽셀 가중치에 per-frame 대칭 9×9 coadd lag kernel을
수축해 finite-lag aperture 분산을 계산했다. 지원된 224쌍 전부(i 71, i2 153), 448프레임을 처리했고
재구성된 aperture coverage와 평균 표면밝기를 저장 profile과 대조했다. `i`의 독립 프레임 pixel-noise
stack sigma는 대각만 가정할 때 0.06274에서 finite-lag 0.07869 nJy arcsec⁻²로 1.254배,
`i2`는 0.06006에서 0.07356으로 1.225배다. 해당 값은 기존 spatial-block-bootstrap 구간 halfwidth의
각각 7.56%, 6.94%다. 이는 추가 pixel-noise floor 진단이며 전체 estimator covariance도 아니고
관측 ICL profile 보정값도 아니다. Target/control median inflation은 각각 i 1.2588/1.2546,
i2 1.2576/1.2441로 유사하다.

9×9 kernel의 aperture 내 국소 정상성 가정, 4픽셀 밖 상관, 상류 sky 전달, source Poisson noise,
확장 PSF wing, mask uncertainty, 프레임 간 covariance는 반영되지 않았다. 과학 측정과 보정은 적용하지
않으며 `science_measurement=false`다. 재현기는 `analyze_icl_annulus_aperture_covariance.py`, 평가는
`analysis/icl_wgl_membership/hsc_debug/icl_aperture_covariance_200_300_assessment_2026-11-03.md`다.
summary 경로 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_aperture_covariance_200_300_2026-11-03_01/summary.json`,
SHA-256 `82b0f2ee2d364b204552933f60dfcc0853e5c9869d014957a07801dfe090fb03`; frame table SHA-256
`d7632bf7b43a978d31045a1bb7919bddb7237580b7664e51243f598b7188605a`다. 독립 해시 manifest와 재계산이
summary 및 두 TSV 파일 해시와 일치했다.

### 2026-11-04 후속 진행: 무작위 중심 weak-lensing source row export

기존 50/124 checkpoint 이후 batch 050–054 (jobs 1068304–1068308)가 각각 114,673, 117,799,
93,429, 127,958, 97,884행으로 완료되어 center별 count-only preflight 대조를 통과했다. 총 55/124
batch, 3,891,646행이다. SHA-256은 차례대로
`473c8f57f8d1eb7c87b9adce1d3edb1165aa896550ba442e9bde29351eef3c2f`,
`8b929f9912c0db548633088bbddba829f933e55652fdaf406d5774bb9fb0aa06`,
`fdc6b20d7eb311cc892842424b6a9daaf16d2fb34f5f74429939b4d65ba897bc`,
`3bc6b00dce4b62a8e6e87c5236d69a7677888b44654eae6d3943781fd186c3ad`,
`128e9ebfc772f9186d049938c3ede9fbcabba2871150f9ad7175c23a2cb0c738`이다. Batch 055 (job
1068309)이 실행 중이며 최신 상태는 `row_status.json`이다. 기존 user 승인 범위와 checkpoint flags는
유지된다: credentials 저장 없음, 개별 photo-z PDF와 이미지는 미포함, `science_measurement=false`.

2026-11-04 추가 진행: batch 055–059 (jobs 1068309–1068313)은 85,314, 124,773, 131,901,
84,962, 43,949행을 반환해 per-center preflight 대조를 통과했다. 누적 60/124 batch,
4,362,545행이다. SHA-256 순서는
`6709975594ceb6017a490ef59f24b57ba486567e46b05aad183bf01ccf63a6ae`,
`930e24befebc9f240e8d4d537a662d0a93b3f04feb5d3efe8e462c9cf4146849`,
`de614d1da55e6daa089cd61c49db4c333324546abced4ab9e77df620ea38d314`,
`9c5b2706796d45a84d1e7955d146fd10ebb4655ee8796256afbd0afa153e5324`,
`0b99e058b76a3b6d423110cebcd96e8391d5a262ba026bdcd5b85dfeba8cb001`이다. Batch 060 (job
1068314)이 실행 중이다.

2026-11-04 진행 갱신: 승인된 random-center shape-row export가 63/124 batch 완료,
누적 4,574,956행·1,504,916,757 bytes다. 최신 `row_status.json`의 plan SHA-256은
`fc883444d2b23663d5f6034d773b76db948add29f1cda0c23384d17cfed2996a`이며,
batch 063 (job 1068317)이 제출 상태다. batch 061 (job 1068315)은 65,200행,
SHA-256 `30d6e51844c6956c867029826e8b570bf5add0d2d8d54b3487ca01a9f8bc03d8`;
batch 062 (job 1068316)은 76,846행,
SHA-256 `2fb4d954dee7ed88b0ed7d904df80058ad5ee1ff210e9d1f6196b436bfa98f44`로
완료되어 각각 center별 count-only preflight와 대조됐다. data는 계속
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_shape_rows_lens_reuse_2026-11-03_01/results/`에만
저장된다. 모든 124 batch가 완료되기 전에는 random stack을 실행하지 않는다.

batch 063 (job 1068317)은 90,879행으로 완료되어 per-center count-only preflight와
일치했고 CSV SHA-256은 `385c432ab147a50ddc6abac929c0cec10329a5c51acbf2ce38acd2c1224e40bb`다.
이에 따라 현재 누계는 64/124 batch, 4,665,835행, 1,534,853,155 bytes다.
batch 064 (job 1068318)은 실행 중이며 제출기 PID 1619431이 살아 있음을 확인했다.

2026-11-04 즉시 갱신: batch 064 (job 1068318)은 83,737행, CSV SHA-256
`2640632164525a6628add6ea3996790650fe560c50bde1911aa7052cba6253ce`로 완료되어
전 center count와 일치했다. 누계는 65/124 batch, 4,749,572행, 1,562,533,263 bytes다.
Batch 065 (job 1068319)이 다음으로 제출됐다. 실행기는 계속 살아 있으며 checkpoint는
동일 경로의 `row_status.json`에서 권위 확인한다.

batch 065 (job 1068319)는 0행·229 bytes의 빈 결과로 완료됐고, 해당 center별
count-only 사전 집계와 일치했다. CSV SHA-256은
`5cab5731ae6624e7c0b5362a9c08cb08e3ab9358935aede924cc96736b594522`이다. 누계는
66/124 batch, 4,749,572행, 1,562,533,492 bytes다. batch 066 (job 1068320)이 제출됐다.

batch 066 (job 1068320)도 동일한 0행·229 bytes 결과로 완료됐고 center별 사전 집계와
일치했다(동일 빈 CSV 해시). 현재 누계는 67/124 batch, 4,749,572행, 1,562,533,721 bytes다.
batch 067 (job 1068321)이 제출 상태이며 동일 실행기 PID 1619431이 계속 살아 있다.

batch 067 (job 1068321)도 0행·229 bytes로 완료됐고 center별 count-only 결과와 일치했다.
동일 빈 결과 해시를 기록했다. 현재 68/124 batch, 4,749,572행, 1,562,533,950 bytes이며
batch 068 (job 1068322)이 제출됐다. credentials 저장·이미지·photo-z PDF 반출 없음,
`science_measurement=false` 플래그는 그대로다.

batch 068 (job 1068322)은 47,790행, 15,773,332 bytes로 완료되어 center별 사전 집계와
일치했고 CSV SHA-256은 `55fbae44d5f13d450b05c6fa5af53ca5b0825084fc722a633f05287ce63e4c41`다.
누계는 69/124 batch, 4,797,362행, 1,578,307,282 bytes다. batch 069 (job 1068323)이
제출됐고 submission process PID 1619431은 살아 있다.

후속 분석을 자동 연결하기 위해 Astrafek에
`run_wl_random_shape_stack_after_export.py`를 추가하고 Lageunha에서 watcher PID 2161773을
시작했다. watcher는 승인 제출기 PID 1619431이 남긴 checkpoint만 감시한다. export가
완료되면 기존 해시 검증 스택을 CPU 낮은 우선순위로 실행하며 CAS에 접속하거나 재제출하지
않는다. 예정 출력은
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_shape_stack_lens_reuse_2026-11-04_01/summary.json`이며,
아직 export checkpoint는 `running`이고 stack 출력은 없다.

갱신: batch 069 (job 1068323)은 107,318행, 35,380,876 bytes로 완료되어 count-only
preflight와 일치했다. CSV SHA-256 `faa7957abcef5012b26a71e3c69ed783dcfeec0369a402a8ab27abae71cf7cac`.
누계는 70/124 batch, 4,904,680행, 1,613,688,158 bytes이며 batch 070 (job 1068324)이
제출됐다. export는 계속 실행 중이고 watcher는 checkpoint 완료를 기다린다.

batch 070 (job 1068324)은 70,975행, 23,360,449 bytes로 완료되어 count-only 사전
집계와 일치했다. CSV SHA-256 `8e3c31a014140dd2619ce677a5d00504fdda08ead22f2e5484d5593057d496b2`.
누계 71/124 batch, 4,975,655행, 1,637,048,607 bytes다. batch 071 (job 1068325)이
제출됐으며 live submitter와 completion watcher가 계속 실행 중이다.

batch 071 (job 1068325)은 82,165행, 27,078,305 bytes로 완료되어 center별 count-only
preflight와 일치했다. CSV SHA-256 `c06a4e944ada34b6d22a8c9a929971d209a4fdb6165bebfe39b9fb2cda0484e5`.
누계는 72/124 batch, 5,057,820행, 1,664,126,912 bytes다. 다음 batch 072는 CAS job
1068327로 제출됐다.

2026-11-04 상태 갱신: batch 072 (job 1068327)은 44,268행, 14,539,325 bytes로
완료됐고 count-only center별 사전 집계와 검증됐다. CSV SHA-256은
`c9a4f9311ec450abae26d2459bfea770aff65f154c0d91c5321c5256e6d45a91`이다. 누계는
73/124 batch, 5,102,088행, 1,678,666,237 bytes이며 batch 073 (job 1068336)이
실행 중이다. lageunha에서 submitter PID 1619431과 후속 스택 watcher PID 2161773이
살아 있는 것을 확인했다. `row_status.json`에는 `status: running`과 함께 이전 실행의
`error_type: KeyboardInterrupt` 필드가 남아 있다. 프로세스와 상태를 다시 대조하기 전까지
이를 종료 상태로 해석하지 않는다. random stack은 아직 시작되지 않았다.

추가 갱신: batch 073 (job 1068336)은 계획의 124,115행을 반환했고 per-center count-only
preflight와 일치했다. CSV는 40,926,920 bytes, SHA-256
`f3cf370505cc53ed958f557c1f065196bd58b32f72c3968131e0513ba34daa8a`이다. 누계 74/124 batch,
5,226,203행, 1,719,593,157 bytes다. 다음 batch 074 (job 1068337)이 제출됐으며 제출기 PID
1619431과 completion watcher PID 2161773은 계속 실행 중이다. 전체 export 전 random stack은
아직 시작하지 않는다.

batch 074 (job 1068337)은 54,239행, 17,907,763 bytes로 완료되어 per-center count-only
preflight와 일치했다. CSV SHA-256은
`9cb4e88a814faabd6b753642282ad756fa8279dc48a6a91cf9bcbeb4d105a3ba`다. 누계는 75/124 batch,
5,280,442행, 1,737,500,920 bytes다. batch 075 (job 1068338)가 제출됐고 CAS 상태는
`running`이다. 기존 제출기와 스택 watcher가 계속 실행 중이다.

batch 075 (job 1068338)은 67,646행, 22,279,477 bytes로 완료됐고 count-only center별
preflight와 일치했다. CSV SHA-256은
`cf8519fc15935730652fa532eb1a6afc3950e134a92878069f7bda6ae03f398e`다. 누계 76/124 batch,
5,348,088행이다. batch 076 (job 1068339)이 제출돼 CAS `running` 상태다.

### 2026-11-04 후속: ICL coadd 주입 포아송 잡음 입력 감사

Astrafek `audit_icl_poisson_noise_inputs.py`가 HSC 다운로드 상태 manifest에 지정된 표적 800개와
matched control 798개 전체의 FITS 헤더를 읽었다. 1,598개 모두 manifest의 파일 크기 및
IMAGE/MASK/VARIANCE HDU index와 일치하고, 1,598개 모두 variance plane을 갖는다. 그러나
primary와 image-plane headers에서 `GAIN`, `GAIN1/2`, `EFFGAIN`, `COADDGAIN`, `EXPTIME`,
`EXPOSURE`, `EFFEXPT`, `SUMEXPT`, `NEXP`, `NINPUT`를 하나도 찾지 못했다. 따라서 현재 coadd
variance plane과 `FLUXMAG0`만으로 합성 광원의 효과적 coadd gain/shot-noise 분산을 복원할 수
없다. 노이즈가 작다고 결론내린 것이 아니며, 기존 ICL injection에 포아송 variance를 추가하지
않았다.

요약 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_poisson_noise_input_audit_2026-11-04_01/summary.json`
SHA-256 `f75e704de3ae8dbdfc6c2f52700bebb9407e13532d0f1dd10eb6b71990fd533c`, 재현기 SHA-256
`bd42fc675aaa7c4fd2c33c22d9bbbcb35309e6518fd38351ac301e7770c04e2b`다. 세부는
`analysis/icl_wgl_membership/hsc_debug/icl_coadd_poisson_noise_assessment_2026-11-04.md`에 있다.
픽셀 데이터 접근·network·새 HSC 자료 수신은 없었다. HSC PDR3 문서는 coadd가 visit 이미지들의
결합이며 개별 visit `warp`와 pixel별 visit count `nImage`가 제품으로 존재함을 설명한다. 이들은
후속 effective-noise 재구성 검토 후보지만 visit count만으로 coadd weight/gain을 정하지 못한다.
기존 ICL injection은 `science_measurement=false`로 둔다.

HSC PDR3 database 문서의 `mosaic`, `frame`, `warped`는 image metadata tables이며 Wide `warped`
table은 약 1.58M warped-image metadata rows를 보유한다. 다음 분석 전에는 target/control coadd patch를
contributing visit metadata에 연결하고, visit-level variance 및 coadd weights로 synthetic-source shot
noise를 계산할 수 있는지 확인해야 한다. 이 metadata query는 아직 제출하지 않았다. 근거는
`analysis/icl_wgl_membership/hsc_debug/icl_coadd_poisson_noise_assessment_2026-11-04.md`의 HSC 공식 문서
링크다.

추가 갱신: weak-lensing random-center batch 076 (CAS job 1068339)은 113,357행,
37,437,519 bytes로 완료되어 per-center count-only 사전 집계와 일치했다. CSV SHA-256은
`cf3800bec37fbf8860c9ca9b591bd253f84963f35f8c13e2eb9b24d24d11f0c4`다. 누계는 77/124 batch,
5,461,445행, 1,797,217,916 bytes이며 batch 077 (job 1068340)이 실행 중이다. 제출기 PID 1619431과
completion watcher PID 2161773이 살아 있고, 전체 row export가 완료되기 전 random stack은 시작되지 않는다.

batch 077 (job 1068340)은 29,684행, 9,515,117 bytes로 완료되어 per-center count-only 사전 집계와
일치했다. CSV SHA-256은 `ae82fcdb1cff03104559c168f1703a9d1242e163dfbc2fb43360572e75b9ae3b`다.
누계 78/124 batch, 5,491,129행, 1,806,972,033 bytes다. batch 078 (job 1068341)이 제출돼 CAS
`running`이며 제출기 PID 1619431과 watcher PID 2161773이 계속 실행 중이다.

batch 078 (job 1068341)은 count-only 예상과 같은 71,116행, 23,298,516 bytes로 완료됐다.
CSV SHA-256은 `260f33dbfd88fff2faf3d34246195ecee6b017bfb83078a905c0a428c52200fb`다. 누계 79/124
batch, 5,562,245행, 1,830,270,549 bytes다. batch 079 (job 1068342)가 실행 중이고 기존
submitter와 completion watcher가 살아 있다.

batch 079 (job 1068342)은 count-only 기대치와 동일한 34,758행, 11,405,309 bytes로 완료됐다.
CSV SHA-256은 `11713e8881c817b977b52c08b83714e9826d5f284edde3613bd8aa3b830f9feb`이다. 누계는
80/124 batch, 5,597,003행, 1,841,675,858 bytes다. batch 080 (job 1068343)이 제출돼 실행 중이다.

추가 갱신: 2026-11-04 현재 batch 083 (job 1068350)까지 84/124 batch가 완료됐다. 검증된
누계는 5,784,264 source-pair 행, 1,903,446,536 CSV bytes이며 batch 084 (job 1068351)이
제출됐다. submitter PID 1619431과 completion watcher PID 2161773을 재확인했다. status는
`running`이다. export 도중 watcher가 실행할 `run_wl_random_shape_stack.py`를 수정해, random
center가 있는 cluster 620개와 3,100개 무작위 중심만으로 matched estimator와 paired spatial
jackknife를 계산하도록 했다. 무작위 중심은 각 620개 cluster ID에 정확히 5개씩 연결된다.
나머지 180개 cluster에는 적격 중심 5개가 없다. 기존 저장 스택은 800개 전체이므로, 새 runner는
먼저 검증된 원본 source rows로 800-cluster 저장 프로필 재현을 확인하고, 이후 공통 620개 ID의
cluster subset으로 random subtraction 및 jackknife를 수행한다. 두 표본을 맞추지 않은 800 대
620 차감은 실행하지 않도록 했다. runner 구문 AST 확인은 통과했으며, 완료된 HSC row export 및
후속 stack은 아직 없다. 출력은 계속 `science_measurement=false`를 유지한다. 상세한 분석 규칙은
`analysis/icl_wgl_membership/hsc_debug/README.md`의 random-center 절에 기록했다.

같은 날 batch 086 (job 1068353)가 추가로 검증 완료되어 현재 87/124 batch, 5,901,175행,
1,941,776,117 CSV bytes다. batch 087 (job 1068354)가 제출됐고 submitter와 watcher는 계속
살아 있다. ICL 결과 추적을 위해 `build_icl_association_evidence_report.py`와 hash-provenance
보고서를 만들었다. 현재 74쌍 고정 association 표본의 i2 부분상관은 raw `r=0.234`, finite-lag
pixel-noise 진단 후 `r=0.249`이고, 지지된 satellite-shell 픽셀을 ICL aperture에서 제외하면
`r=0.193`이다. 배제 전후 paired 변화의 95% 공간 bootstrap 구간은 `[-0.064,-0.021]`이지만
28개의 근사 공간 블록에 조건부인 민감도 결과다. i-band는 21쌍으로 불확실성이 넓다. 외곽 잔차
공변량을 넣었을 때 i2 변화는 `+0.013`, 구간 `[-0.034,+0.111]`로 0을 포함한다. 따라서 i2 신호는
추적할 근거는 있으나 위성/ICL 픽셀 중첩에 민감하며 확정 검출은 아니다. ICL 합성광원 포아송잡음,
PSF 날개, sky transfer, completeness 및 전체 공분산 문제는 남아 있다. 보고서 파일은
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_association_evidence_report_2026-11-04_01/report.md`,
재생성 스크립트는 `analysis/icl_wgl_membership/hsc_debug/build_icl_association_evidence_report.py`다.
후속 점검에서 보고서의 5개 입력 JSON SHA-256과 builder SHA-256이 현재 파일과 모두 일치했고,
5개 원본 진단의 `science_measurement`도 false임을 확인했다.

HSC frame metadata 경로를 추가 조사했다. 공식 PDR3 [mosaicframe schema](https://hsc-release.mtk.nao.ac.jp/schema/schemata/pdr3.pdr3_wide.mosaicframe.json)는 coadd `(skymap_id, band)`에 포함된 frame/visit/CCD를 제공하며, [frame schema](https://hsc-release.mtk.nao.ac.jp/schema/schemata/pdr3.pdr3_wide.frame.json)는 `exptime`과 `gain1`–`gain4`를 정의한다. 기존 800개 target 및 798개 control 위치의 count-only SQL plan v1은 16개 batch, 1,598개 요약 행, SHA-256 `7d8b04df45489518b24095910bf201a93687f9e317aa2bf350ef0747c1c97570`였다. 승인된 v1과 v2 SQL 모두 CAS HTTP 406을 반환했고 job ID는 발급되지 않았다. v2 첫 batch의 응답 본문은 `column "patch_area" does not exist`였다. 각 결과 checkpoint의 `query_submitted=false`, 자격증명 저장·이미지·shape·photo-z 반출도 false다. 공식 PDR3 [mosaic schema](https://hsc-release.mtk.nao.ac.jp/schema/schemata/pdr3.pdr3_wide.mosaic.json)는 `areacube`와 WCS 중심 좌표를 제공한다. 따라서 v3 계획은 `public.skymap.patch_area` 대신 `mosaic.areacube && radec_to_coord(ra, dec)` 공간 조건을 사용한다. v3는 1,598개 중심·16 batch·13개 count-only 열·173,398 SQL bytes, plan SHA-256 `3dac4db631d0dcbee5c5f4509e639f2bd90ca2af99fa7cf5f66c2a9f60bd9014`, SQL aggregate SHA-256 `924c478ec02b5e703027f0446b834c381cf479f7e4ac2ab7dcf27b3af9d999d6`다. 이 SQL hash는 별도이므로 제출 전 v3 명시 승인이 필요하다. 계획은 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_frame_metadata_count_plan_2026-11-05_03/plan.json`에 있다. 카운트 결과 뒤 실제 gain/exposure metadata row export는 별도 범위와 승인을 받아야 한다. Per-pixel coadd weights가 없어 frame gains만으로 포아송 보정을 확정할 수는 없다.

weak-lensing 반출도 진전됐다. batch 088 (job 1068355)까지 89/124 batch가 완료되어 6,008,508 source-pair 행, 1,976,901,478 CSV bytes로 확인됐다. batch 089 (job 1068356)이 제출됐다. submitter PID 1619431과 post-export watcher PID 2161773은 살아 있고 최종 stack 결과는 아직 없다.

ICL–위성 표본 선택 원인을 추가 분해했다. `audit_icl_satellite_outer_coverage_loss.py`가 224개 ICL 지원 cluster/filter 프레임의 후보 2,691개를 재검사해, unrounded 50–75 kpc coverage 40% 기준을 1,171개가 통과하고 1,520개가 실패함을 확인했다. 실패 단계는 base quality mask 465개, SEP segmentation 81개, 이웃 검출원 마스크(Kron 타원 또는 유효하지 않은 모양의 원형 마스크) 973개, 미검출 CAMIRA 멤버 원형 마스크 1개다. image footprint에서 처음 기준 아래로 내려간 후보는 없다. 따라서 손실은 이미지 경계보다 source masking에 주로 연결되지만, mask 완화의 정당성을 뜻하지 않으며 flux correction이나 ICL–위성 검출도 아니다. 3-decimal saved coverage가 경계 후보 4개를 `0.400`으로 보이게 해 최종 원인 수는 미반올림 픽셀 재구성으로 분류했다. 출력은 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_outer_coverage_loss_2026-11-05_03`; 평가 문서는 `analysis/icl_wgl_membership/hsc_debug/icl_satellite_outer_coverage_loss_assessment_2026-11-05.md`에 있다. `science_measurement=false` 유지.

표본선택 진단도 진행했다. `analyze_icl_satellite_mask_selection.py`가 같은 2,691개 후보를 cluster ID 단위 5,000회 bootstrap으로 비교했다. 50–75 kpc profile 유지율은 cluster-centric radius 최상위 사분위에서 최하위보다 i `+22.7%p` (95% block interval `[+10.6,+35.1]%p`), i2 `+20.6%p` (`[+12.7,+28.2]%p`) 높았다. i에서 richness 최상위 사분위는 최하위보다 `−19.9%p` (`[−29.5,−9.2]%p`) 낮았고 i2 차이는 구간이 0을 포함했다. 이는 retained profile이 반경에 따라 달라지는 선택표본임을 보여주며 보정값은 아니다. 다음에는 cluster-level multivariable selection을 정하고 injection/recovery로 보정 가능성을 평가해야 한다. 평가 문서 `analysis/icl_wgl_membership/hsc_debug/icl_satellite_mask_selection_assessment_2026-11-05.md`, 결과 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_mask_selection_2026-11-05_01`.

다변수 selection 진단도 완료했다. standardized radius·log(S/N)·redshift·log(richness)를 포함한 cluster-robust logistic model 4개가 수렴했다. 다른 변수를 조정한 뒤에도 profile 유지의 radius OR/SD는 i `1.55` (`95% CI 1.25–1.92`), i2 `1.45` (`1.28–1.64`)였고, richness OR/SD는 각각 `0.64` (`0.51–0.81`), `0.81` (`0.69–0.94`)였다. 이는 기술적 선택 모델이지 completeness 보정이나 과학적 가중치가 아니며, injection/recovery 검증 전 association에 적용하지 않는다. Astrafek 실행기 `analysis/icl_wgl_membership/hsc_debug/analyze_icl_satellite_mask_selection_adjusted.py`; 요약 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_mask_selection_adjusted_2026-11-05_02/summary.json` (SHA-256 `9fe70bf1da239c05a4bf72d6675918e1864b7186a458dcf0a7ec6f0ea74800af`).

2026-11-05 support threshold 감도를 기존 60% PSF 주입/재검출 결과와 비교했다. 새 40%와 이전 60% 결과의 동일한 57개 프레임(i 13, i2 44)에서 이미지·위치별 PSF 해시, 재검출 개수와 마스크가 일치함을 확인했다. 프레임당 median satellite 수는 i `2→5`, i2 `2→6`; median control profile 수는 `7→26`, `9→34.5`였다. satellite unit-injection response 차이의 cluster bootstrap 95% 구간은 i `−0.033 [−0.768,0.487]`, i2 `+0.098 [−0.264,0.526]`로 모두 0을 포함했다. 이는 지원/추정량 민감도이며 모집단 completeness 보정이나 과학적 검출이 아니다 (`science_measurement=false`). 보고서 `analysis/icl_wgl_membership/hsc_debug/icl_satellite_support_threshold_comparison_2026-11-05.md`; 재생성기 `compare_icl_satellite_support_thresholds.py`; 출력 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_support_threshold_comparison_2026-11-05_02`. Runner SHA-256 `55b8a89b40fdc9191412f30602d48218015d63c73157431aada0c114d84f66d6`, paired table SHA-256 `c86368c4dfb8aa4811f1d32bf9f133b341c2d14a3189a31fdc0baab084ca0add`, summary SHA-256 `c2722311dc4c5dbfaf6592751a6ae58358152baa0376e27a3cb88f8e2529ffe2`.

사용자가 승인한 HSC coadd frame count-only SQL plan v2 (SHA-256 `b92a3c4931de25cf560abf5f6b0cd39af86b4e6d8cd3001a19cbeaa681f4316a`)를 제출 시도했으나 CAS가 `patch_area` 열 오류로 HTTP 406을 반환했다. job이 생성되지 않았고 checkpoint `query_submitted=false`; 결과 행이나 자격증명은 저장되지 않았다. 기존 SQL 재제출은 하지 말고 공식 스키마의 올바른 영역/WCS 열을 확인해야 한다. 수정 SQL은 새 plan hash에 대한 별도 승인이 필요하다. 현 승인 범위는 count-only로 영상·shape·photo-z를 요청하지 않는다.

224개 ICL-supported pair의 위성 PSF coverage를 늘리는 오프라인 요청계획을 준비했다. 기존 744개 PSF 중 688개는 57프레임 위성 주입 위치, 기존 추가 1,542개는 target/control 중심을 덮는다. 기존 800 HSC 이미지, ICL support table, satellite prep에서 224 pair의 위성 주입 위치 2,693개를 재구성했으며 기존 파일·매핑 해시와 대조해 2,005개 위치가 빠진 것을 확인했다. `prepare_icl_psf_picker_supported_pairs_plan.py` (SHA-256 `5d97f940360220164908b96b3e63bf553782d373fb3484eb9ae1e042519258ee`)가 coordinate-only 계획을 만들었다. 계획 해시 `1d7c0189721600fc31e8524a4258dc70a5c0f0d436df45b57b9a0c3dd0058b65`; 좌표목록 SHA-256 `d71f38e67028d0e258e0b17321229f0f0102c528bdaf23c2e141c4d7d23b5ed6`; mapping SHA-256 `1dc1625cf88b86e95736c54e1a5977c6a787e4d9f8ffbe5d7dc96fe693186deb`. 요청 범위는 HSC-I PDR3 Wide centered coadd model PSF FITS 2,005개, 예상 40,420,800 bytes (38.55 MiB)다. 이미지·shape row·photo-z PDF는 포함하지 않으며, 요청은 제출되지 않았다. PSF 수신은 사용자의 별도 승인 후 진행한다. 수신되면 파일 및 좌표 감사를 통과한 뒤 injection/recovery 검증 표본을 57에서 최대 224개 supported pair로 늘릴 예정이다. 이는 과학적 검출이 아닌 선택/시스템 오차 검증이다.

추가 다운로드 없이 PSF-supported 57프레임에서 위성 주입 진폭 0.5·1·2 nJy arcsec⁻²를 짝지어 비교했다. i2 normalized satellite-profile recovery median은 각각 `0.593`, `0.572`, `0.560`; 두 amplitude 대비의 paired 95% cluster-bootstrap 구간 모두 0을 포함한다. 다만 위성 밝기에 따라 재생성 detection mask가 변했다(i2에서 satellite-only mask pixel 중앙 증가는 약 `230`, `631`, `1,014`). 따라서 median transfer는 이 주입 범위에서 대체로 일정해 보이나 selection은 진폭 독립이 아니다. i 표본은 13개뿐이고 산포가 커 검출 또는 completeness 보정으로 해석하지 않는다. 평가 `analysis/icl_wgl_membership/hsc_debug/icl_satellite_injection_amplitude_grid_assessment_2026-11-05.md`; 분석기 `analyze_icl_satellite_injection_amplitude_grid.py` (SHA-256 `fa64fd9d554e9b44b8306859bc1ce3a1a2a7679db84898fa8eda079498561852`); 출력 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_injection_amplitude_grid_2026-11-05_02`; paired table SHA-256 `dbd5489eb3c799bb6a7fb53f9251a800923dd2fd0cf30f679589a5b421cab5e5`, summary SHA-256 `733bcb90a7471be061b0280259e6cc9bc87fb05f1761f500d25cd91727d22d63`.

약중력 stack의 zero-signal 검증을 5°/2.5° paired endpoint-jackknife full covariance로 추가했다. 435 common lens의 tangential ΔΣ profile은 0.1–3 Mpc에서 각각 χ²/7 dof `479.81` (51 blocks), `642.83` (154 blocks)이며, inner 0.4 Mpc를 제외해도 `460.37/5`, `537.23/5`다. chi-square approximation의 nominal p는 매우 작지만 block partition, random-center depth/mask/PSF/sky, coadd covariance, n(z) calibration, reduced shear가 미완료라 이를 최종 과학 검출로 승격하지 않는다. Cross null은 p=`0.157`/`0.399`; covariance 모두 full-rank positive definite이며 point profile는 block width에 따라 동일하다. Astrafek 분석기 `analysis/icl_wgl_membership/hsc_debug/analyze_wl_random_shape_stack_detection.py` (SHA-256 `e0934796e204acaec678a2a7dde5ae971920d098dae24fc1539659c1dc93bed1`); 평가 `wl_random_shape_stack_detection_assessment_2026-11-05.md`; 결과 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_shape_stack_detection_audit_2026-11-05_01/summary.json` (SHA-256 `ed6040a61157df74532fccf3e67e23c6f8231a33fe953045fd19ae35acc1cffb`). `science_measurement=false` 유지.

2026-11-05 후속 무결성 감사에서 승인된 HSC random-center shape 행 반출의 124/124 CSV를 체크포인트 및 5°·2.5° stack provenance와 대조했다. 총 9,010,324행, 2,962,522,918 bytes이며 파일 hash/bytes, 계획 hash, stack의 batch hash가 모두 일치했다. row checkpoint에는 `status=complete`와 `error_type=KeyboardInterrupt`가 공존한다. 모든 job과 row summary가 complete이므로 stale error 필드로 기록하되 기존 체크포인트를 수정하지 않았다. 스택은 620개 계획 렌즈에서 cluster+random 공통지원 435개, random 중심 2,175개를 사용한다. 이 감사는 파일 무결성만 확인하며 random aperture depth/mask/PSF/sky 일치와 covariance 공간분할을 검증하지 않는다. 분석기 `analysis/icl_wgl_membership/hsc_debug/audit_wl_random_shape_export_integrity.py` (SHA-256 `1fb7ccf2122a118d0be373ec0fbfb70f486413cd125eb13bbb74c0897d51d60a`); 평가 `analysis/icl_wgl_membership/hsc_debug/wl_random_shape_export_integrity_assessment_2026-11-05.md`; 결과 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_shape_export_integrity_2026-11-05_01/summary.json` (SHA-256 `b7a371a6ff5e4c8d546cbc88ec81a2289f57c09816fd5874eef135e2cc7da27d`). `science_measurement=false` 유지. 새 CAS query, image download, stack 실행은 하지 않았다.

추가 PSF 다운로드를 쓰는 경로도 사전 보강했다. 기존 2,005-position 요청 계획을 바탕으로 `SupportedPairsPositionPSFBank`와 `run_joint_icl_satellite_batch.py`, `verify_supported_pairs_position_psf_bank.py`가 이제 정확한 명시 승인 계획 SHA-256을 인수로 요구하며, 실행 결과에는 요청 승인 증거와 과학적 PSF 승인 상태를 따로 기록한다. 기본 승인 해시는 `1d7c0189721600fc31e8524a4258dc70a5c0f0d436df45b57b9a0c3dd0058b65`이나 해당 PSF 수신은 아직 승인·실행되지 않았다. 수정 파일은 `analysis/icl_wgl_membership/hsc_debug/position_psf.py`, `run_joint_icl_satellite_batch.py`, `verify_supported_pairs_position_psf_bank.py` 및 hsc_debug README다. 분석·다운로드를 실행하지 않았다.

2026-11-05에는 승인된 자료만으로 S23B 후보 표본의 기여 상한을 다시 계산했다. 현 800 parent에 catalog-only 후보 1,333개를 더하면 2,133개다. observed i2 satellite support yield를 그대로 가정하면 50–75 kpc pair는 60% 기준 181개, 40% 기준 373개 예상이지만 40%의 `r=0.249` 95% block 구간은 0을 포함하며 estimator 변경도 있다. S23B 후보 ledger는 100개 5° target 셀, 그중 discovery와 겹치지 않는 셀 22개를 보인다. 따라서 S23B 전체도 단일 endpoint 계획 92 independent blocks 및 30 비교 familywise 180 blocks의 공간상한에 미달한다. 이는 이미지·mask·PSF·depth·control 지원을 알기 전 상한일 뿐 실제 표본이 아니다. 재현기 `analysis/icl_wgl_membership/hsc_debug/assess_s23b_support_yield_and_blocks.py` (SHA-256 `5fb4c1a4354546ff972995a861e83a8f87aa92704611669d2e3138442f27b026`); 평가 `analysis/icl_wgl_membership/hsc_debug/icl_s23b_support_yield_and_blocks_assessment_2026-11-05.md`; 산출물 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_s23b_support_yield_and_blocks_2026-11-05_01/summary.json` (SHA-256 `b4f95bc248351a4e05446ce4fdb258485b863004dce55c811727c3af1f5b6951`). HSC query/download 없음, `science_measurement=false`.

224개 ICL-supported cluster/filter pair 중 50–75 kpc 위성 프로파일이 있는 94개( i 26, i2 68 )의 표본선택 민감도를 계산했다. 모든 224개에서 이용 가능한 redshift·log richness·ICL annulus coverage·중심 mask flag로 pooled logistic support 확률을 추정하고 stabilized IPW를 적용했다. selection subset과 나머지의 ICL coverage 표준화 평균 차이는 i `+0.44`, i2 `+0.70`; richness는 `−0.29`, `−0.34`다. 선형 richness와 redshift를 통제한 상관은 i `0.079→−0.011`, i2 `0.293→0.227`이다. BFGS에서 “precision loss”로 보고된 216개 bootstrap 적합은 모두 gradient norm 기준을 충족했다(최대 `1.40×10⁻⁶`); 수정 실행은 허용기준 기반으로 0개 모형 적합 실패를 기록했다. 중심 mask 제외·coverage 제곱항·필터별 slope 대안에서도 IPW 부분상관 점추정은 i `[-0.014,0.021]`, i2 `[0.226,0.231]` 범위였다. 1,994개 i 및 2,000개 i2 연결 5° sky-block bootstrap 유효 표본의 IPW 95% 구간은 i `[-0.566,0.581]`, i2 `[-0.034,0.416]`이다. i의 6개 표본은 최소 6개 지원쌍 미달로 제외했다. 유효 표본수는 22/26, 60/68이다. IPW는 가정된 조건부 무작위 결측의 민감도일 뿐, 그 가정은 검증되지 않았고 검출이나 completeness 보정으로 해석하지 않는다 (`science_measurement=false`). Astrafek runner `analysis/icl_wgl_membership/hsc_debug/analyze_icl_satellite_joint_support_ipw_sensitivity.py` (SHA-256 `00131be3d6f91246d7e1733ee2e06845525359dd183d8464b87db26d822deb85`); 평가 `analysis/icl_wgl_membership/hsc_debug/icl_satellite_joint_support_ipw_sensitivity_2026-11-05.md`; 결과 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_joint_support_ipw_2026-11-05_06/summary.json` (SHA-256 `8dbfae2169361c9fb6ce5ccfcce684b486162dbf4820e531112cdfc4f39fdbc7`). IPW 분석은 HSC 네트워크를 사용하지 않았다.

2026-11-05 S23B CAMIRA 후보 선택을 6개 redshift/BCG-spec 시나리오로 오프라인 비교했다. baseline은 기존 1,333 후보 ID를 정확히 재현했다. BCG-spec을 유지한 z 확장(`0.25–0.65`, `0.20–0.70`)은 후보 1,983/2,339개, 신규 5° 셀 25/26개. BCG-spec 없는 최광 선택은 후보 4,137개, 신규 셀 31개이나 렌즈 redshift 보정이 필요해 주표본으로 채택하지 않는다. 카탈로그 후보/셀은 영상 지원이나 독립 공분산 블록을 보증하지 않는다. 결과는 `analysis/icl_wgl_membership/hsc_debug/icl_s23b_camira_selection_scenarios_2026-11-05.md`, runner는 `audit_s23b_camira_selection_scenarios.py`; 기계 요약 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_s23b_camira_selection_scenarios_2026-11-05_02/summary.json` (SHA-256 `e0fddda7ecc973327842796278dd01b8e2ae3de33f2d07244590cfb085507ff6`). 신규 HSC query/download는 없고 `science_measurement=false`.

2026-11-05 HSC Y3 네 source-z bin을 분해한 paired weak-lensing 진단을 추가했다. 공통 435 lens/2,175 random center에서 전체 bin profile과 pair count 합계는 기존 스택과 재현됐다. zbin1 exact estimator-weight boost는 0.1–0.2 Mpc `4.605±0.274`, 2–3 Mpc `1.260±0.049`; zbin2–4는 약 0.97–1.16이다. zbin1 대 각 high bin의 21 반경 대비가 5°/2.5° paired delete-one jackknife에서 Holm 보정 뒤 모두 남는다. zbin1 제거 profile도 양의 tangential/0-consistent cross 신호를 유지하지만, 안쪽 amplitude가 유의하게 변한다(동일 jackknife covariance 기준). 이는 저적색편이 source selection/photo-z overlap과 일치할 수 있으나 member contamination을 판별하지 않으며 boost correction을 정당화하지 않는다. 전체 random aperture depth/mask/PSF/sky, `n(z)` 불확실성 및 reduced-shear 관문은 미완료다. 재현기 `analysis/icl_wgl_membership/hsc_debug/analyze_wl_tomographic_weight_boost.py` (SHA-256 `497b1ac799f261beae91dd24bb8a2a61da65817f2943ba5411d95173d3e5d1673`); 평가 `analysis/icl_wgl_membership/hsc_debug/wl_tomographic_weight_boost_assessment_2026-11-05.md`; 산출물 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_tomographic_weight_boost_2026-11-05_04/summary.json` (SHA-256 `47419ac229ffec9aab2a001f629252627126a456812b01cf11908ad44122fe14`). HSC 신규 query/download 없음, `science_measurement=false`.

ICL–위성 94 지원쌍을 central-mask eligibility로 나눈 5,000회 block bootstrap에서 적격 표본은 i 8쌍/r=`−0.421`/95% 구간 `[-1,1]`, i2 15쌍/r=`+0.307`/구간 `[-0.069,0.792]`였다. 부적격은 i 18쌍/r=`+0.234`/`[0.014,0.802]`, i2 53쌍/r=`+0.312`/`[-0.070,0.551]`. i2 점추정은 유사하지만 구간이 넓고 0을 포함하며, 이 분리 bootstrap은 eligible/ineligible 차이를 직접 검정하지 않는다. 작은 적격표본은 mask effect 결론이나 mask 완화 근거가 되지 않는다. 보고서 `analysis/icl_wgl_membership/hsc_debug/icl_satellite_center_mask_sensitivity_assessment_2026-11-05.md`, runner `analyze_icl_satellite_center_mask_sensitivity.py` (SHA-256 `839b14e3ce4fabd504fe5398a7d6d20fc9fb0141810cb0d65fa1030904a66fd7`), 산출물 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_center_mask_sensitivity_2026-11-05_01/summary.json` (SHA-256 `f779e417d5a5ebaf9c8e688f3d089a19e9363f15ba65375491717ee938fe91cc`). 새 HSC query/image read 없음; `science_measurement=false`.


2026-11-05 HSC PDR3 frame-count v3 조회(승인 plan SHA-256 `3dac4db631d0dcbee5c5f4509e639f2bd90ca2af99fa7cf5f66c2a9f60bd9014`)가 이미 16/16 batch 완료되어 1,598개 중심 모두 patch에 연결됐다. 표적·대조 프레임 수 중앙값은 각각 42개이고 쌍별 차이 중앙값은 0이다. 이 결과는 촬영 횟수 차이가 큰 편향은 아님을 보일 뿐 픽셀별 분산/가중치를 주지 않는다. 후속 gain/exposure 행 반출을 오프라인 준비했으며, runner `analysis/icl_wgl_membership/hsc_debug/prepare_icl_frame_metadata_rows.py`가 기존 count 결과를 검증해 67,071 center-frame 행·21열·16 batch 계획을 만들었다. 정확한 plan SHA-256 `49aa4b7601dd42243209d7b3b394a4dd3b4cbf00d87d28cbb9467553c29f4ce3`, SQL aggregate SHA-256 `4270f3f0d290ac427b24694f7561f46cdee71cc8c3621a57c3fcf3f034583c76`. 계획 경로 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_frame_metadata_row_plan_2026-11-05_01/plan.json`; 상태 `prepared_not_submitted`, 새 CAS 제출 없음. 별도 사용자 승인이 있어야 제출한다. 이 frame-level 정보로도 per-pixel coadd weights/amplifier map이 복원되지 않으므로 ICL 주입 포아송 잡음이 해결되는 것은 아니다. `science_measurement=false` 유지.


승인된 HSC frame-count 표를 기존 200–300 kpc ICL 대비와 연결해 frame-count 민감도를 추가 평가했다. 224 지원쌍(i 71/18 connected 5° blocks, i2 153/40 blocks)에서 ICL target-control 표면밝기 차이를 종속값으로 두고 target-control unique frame/visit/CCD 수 차이를 각각 시험했으며 redshift, log richness, common aperture coverage를 조정했다. 5,000회 block bootstrap partial-r은 i `−0.085, −0.061, +0.078`, i2 `+0.017, −0.041, −0.008`; 모든 95% 구간이 0을 포함한다. 정수 프레임 수 불균형이 이 선택 표본의 ICL 대비를 설명한다는 증거는 없지만 픽셀 깊이·분산·coadd weight의 동등성을 보이지는 않는다. Astrafek runner `analysis/icl_wgl_membership/hsc_debug/analyze_icl_frame_count_icl_sensitivity.py` (SHA-256 `de0ae9afbe5c6261b24009a9b9b6cdc0f9b24c6d4832e47f26d1eea9968b6bd4`); 평가 `analysis/icl_wgl_membership/hsc_debug/icl_frame_count_icl_sensitivity_assessment_2026-11-05.md`; 결과 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_frame_count_sensitivity_2026-11-05_01/summary.json` (SHA-256 `dd00410649e230359df6a6ea6582be2df6881e247bf81d29347bdca33d569cca`). 기존 저장 표만 사용, image pixels/HSC 추가 질의 없음, `science_measurement=false`.

2026-11-05: Astrafek의 `analysis/icl_wgl_membership/hsc_debug/fit_wl_random_stack_nfw_mass.py`
(SHA-256 `b8999e9e534def7eeeff9f29921e9ed479d57eb31e935d21c1deaba2c271c3e2`)가 저장된 435-lens
random-subtracted weak-lensing 스택에 단일 NFW 질량 선별을 수행했다. source bin 2–4,
proper 0.2–1.5 Mpc에서 고정 `c=3`의 유효 `M200c`는 약 `1.78e14`/`1.80e14 Msun`
(5°/2.5° covariance)이고, `c=5`는 `1.54e14`/`1.44e14 Msun`이며 적합도가 나쁘다.
농도 가정에 따른 약 20% 이동 때문에 질량 보정 결과로 쓰지 않는다. reduced shear,
magnification, 중심 오차, two-halo, photo-z, 무작위 aperture 매칭 및 covariance 오차가
미해결이다. 자세한 평가 `analysis/icl_wgl_membership/hsc_debug/wl_random_stack_nfw_mass_screen_assessment_2026-11-05.md`,
산출물 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_stack_nfw_mass_2026-11-05_01/summary.json`
(SHA-256 `d3ca78a24fc030f64694247fd45d33f607d89c18fb102b6b5a6879d53c52cd4b`). 새 HSC 조회는
없고 `science_measurement=false` 유지.

2026-11-05: reduced-shear 보정 민감도를 Astrafek runner
`analysis/icl_wgl_membership/hsc_debug/fit_wl_random_stack_reduced_shear_sensitivity.py`
(SHA-256 `92394e9ba11bc49582c7aab23a79aeb01d8d7798d1656a6bf6a42d80a8a3c556`)로 계산했다.
동일한 435개 random-subtracted 렌즈와 5°/2.5° 공분산에서 HSC Y3 `n(z)` 모멘트 및
8차 reduced-shear 급수를 적용하자 NFW 유효 질량은 weak-shear보다 1.7–5.5% 낮아졌다.
급수 꼬리 상한은 0.012–0.022%지만 고정 농도에 따른 차이가 약 20%이고 `c=5` 적합은
나빠 질량 보정이나 ICL–위성 검출로 사용할 수 없다. 평가
`analysis/icl_wgl_membership/hsc_debug/wl_random_stack_reduced_shear_sensitivity_assessment_2026-11-05.md`,
summary `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_stack_reduced_shear_sensitivity_2026-11-05_03/summary.json`
(SHA-256 `310a03e29069b7799239f21a9ae3c0d0ef0abda5e2d0da93c9e8bdc2b50bc286`). 새 HSC 요청 없음;
`science_measurement=false` 유지.

2026-11-05: 204쌍 random-subtracted richness contrast에 저/고 그룹 NFW 질량을 동시에 적합했다.
5° shared covariance와 0.2–1.5 proper-Mpc 네 bin에서 고/저 유효 `M200c` 비는 고정 `c=3,4,5`에
대해 source bin1 포함 `1.89, 1.40, 1.21`, bin1 제외 `2.85, 1.73, 1.39`다. profile 비 구간은
매우 넓어 upper limit이 최대 28이고 절대 질량은 농도·tomographic 선택에 크게 변한다. 따라서
고 richness가 더 큰 유효 질량을 갖는 방향은 보이나 CAMIRA richness–mass 보정으로 사용하지 않는다.
runner `analysis/icl_wgl_membership/hsc_debug/fit_wl_richness_split_nfw_mass_contrast.py`
(SHA-256 `25e9c57d59cd3c61568bbbf5478f9726bc31263f49ad442931c7909a5dd22f1d`), 평가
`analysis/icl_wgl_membership/hsc_debug/wl_richness_split_nfw_mass_contrast_assessment_2026-11-05.md`.
All-bin summary SHA-256 `c0e4c92e5ef00fc77ffef7e26e918d3da312049e6d5ee2877edbe4385c32217b`, source-bin1-excluded
summary SHA-256 `6d293f41c12a2d44fc4dd7303be238c9bfbcca27570db180f0bf93e3224eaa94`. 공분산·매칭영역·PSF·
sky 및 질량분포 관문이 열려 `science_measurement=false`; 새 HSC 요청 없음.

2026-11-05: CAMIRA richness 차이를 기존 random-subtracted 스택에서 다시 계산했다.
처음 214쌍 계산에서 각 그룹 10개 렌즈가 target tomography에 기여하지 않는 것을 확인해
그 산출은 폐기 대상이며, estimator-usable target/random 공통지원으로 재계산했다. 최종은
415개 유효 공통지원 렌즈(211 low/204 high) 중 204개 redshift-matched pair, 2,040개 예정
random aperture(1,911개에 shape rows), 50개 공유 5° block이다. high−low random-subtracted
`DeltaSigma` 차이는 proper 0.2–3 Mpc에서 모두 양수이며 0.2–0.4 Mpc에서
`+60.16 +/- 18.09 × 10^12 Msun Mpc^-2`이다. 명목 Wald `p=4.04e-12`는 approximate
equal-area covariance 때문에 검출 유의도가 아니다. 매칭 깊이/마스크/PSF, selected support,
`n(z)`·shear 보정 및 mass distribution 미해결. runner
`analysis/icl_wgl_membership/hsc_debug/analyze_wl_random_subtracted_richness_split.py`
(SHA-256 `6ecd77d5c1ce81b0595a4e04b8aa1b09b118929921d1e4926500adc25004cadb`), 평가
`analysis/icl_wgl_membership/hsc_debug/wl_random_subtracted_richness_split_assessment_2026-11-05.md`,
summary `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_subtracted_richness_split_2026-11-05_02/summary.json`
(SHA-256 `b38a78f6005dc17d8a12ddd6117debd63a005d0b21ff8dbda05961bd31eb65b7`). 새 HSC 요청 없음;
`science_measurement=false` 유지.

## 2026-11-05 weak-lensing reduced-shear richness contrast

Added and ran Astrafek `analysis/icl_wgl_membership/hsc_debug/fit_wl_richness_split_reduced_shear.py` on the saved 204-pair redshift-matched richness split. The same random-subtracted profile, shared 50-block covariance, exact source rows/weights, and HSC Y3 source `n(z)` were retained. With order-16 reduced-shear modeling, high/low effective NFW mass ratios for fixed `c=3,4,5` are 1.168/1.200/1.214 using source bins 1–4 and 1.150/1.172/1.198 excluding source bin 1. Weak-shear fits had 1.894/1.396/1.209 and 2.847/1.727/1.391 respectively. The series order 8-to-16 ratio change is ≤0.3%; maximum conservative per-pair truncation bound is 6.94% for c=3 with source bin 1 excluded. Conditional Δχ²=1 intervals use the approximate saved covariance and are not detections. Mass-distribution, photo-z/calibration, miscentering, magnification, matched-random depth/mask/PSF/sky, and covariance limitations remain; `science_measurement=false`. Full limits and hashes: `analysis/icl_wgl_membership/hsc_debug/wl_richness_split_reduced_shear_assessment_2026-11-05.md`. Summary JSONs: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_richness_split_reduced_shear_2026-11-05_09/summary.json` (SHA-256 `1b1fbbafc7aaf90e57b1b651dc29979b955bb593cf2e8e5e99dbd89dd4c91be4`) and `_10/summary.json` (SHA-256 `c3885d73ec8dced7bb9a87fd96510a00b075de2dfdcd6d54156c11c601af4044`). No HSC query or image read.

## 2026-11-06: approved 2,005-position PSF batch and ICL-supported injection

The user approved the exact HSC PDR3 Wide centered HSC-I coadd model-PSF plan
SHA-256 `1d7c0189721600fc31e8524a4258dc70a5c0f0d436df45b57b9a0c3dd0058b65`.
All 2,005 requested PSFs were received (38,088,000 FITS bytes) and passed
coordinate/FITS audit; no images, shape rows, or photo-z were requested. The
combined position-dependent bank verified 224 cluster/filter frames, 2,693
satellite roles, and 2,917 unique PSF files (summary SHA-256
`9c0002344fc7642a5a224ef108b209bb65ca468948ed0309420e8e12e4fac5ac`).

The CPU redetection injection selected 137 of these 224 pairs by the measured
200–300 kpc ICL common-coverage floor of 40% (44 `i`, 93 `i2`); 87 pairs did
not pass that profile-support threshold. It completed all 137 with zero errors
in 671.5 s. Injections used position-dependent compact PSF cores, SEP mask
regeneration, and the `mask_all` policy for unresolved BCGs. No original FITS
or variance plane was changed. Independent audit verified all 137 frame hashes,
PSF role hashes, code provenance, and support-gap accounting; 79/137 frames
have at least one satellite-profile gap.

The diffuse ICL response median is 0.226 in both filters (16th–84th percentiles
0.178–0.249 in `i`, 0.194–0.257 in `i2`). Satellite envelopes contribute a
median 0.342/0.449 to the measured ICL annulus in `i`/`i2`. Satellite-profile
response is unstable and support-limited: only 14 `i` and 44 `i2` frames have a
valid satellite-injection response, with broad ranges spanning negative values
to around or above unity. These are sensitivity diagnostics, not a correction
matrix or ICL–satellite detection. Injected Poisson noise, upstream HSC sky
transfer, extended PSF wings, and calibrated flux covariance remain absent;
`science_measurement=false` and `psf_approved=false` remain.

Assessment:
`analysis/icl_wgl_membership/hsc_debug/icl_supported_pairs_psf_redetection_assessment_2026-11-06.md`.
Batch summary SHA-256 `be3c2e8cb22185c06768d9022ae3e6192072f180fa54d55e21b85920e5572750`;
independent audit SHA-256
`b91783ee4e1582c24676d9a9bd85aacec2cebafa46d76980e909fcbfecbc0bd1`.
The focused auditor is
`analysis/icl_wgl_membership/hsc_debug/audit_icl_supported_pairs_psf_redetection_batch.py`.

## 2026-11-06: 0.5/1/2 satellite-injection amplitude sensitivity on 137 pairs

Expanded the existing paired amplitude grid from 57 to all 137 ICL-coverage-
selected pairs. The user-approved position-dependent PSF bank, cluster images,
and baseline detections are identical in all runs. The 0.5, 1, and 2 nJy
arcsec⁻² runs each completed 137/137 with zero frame errors and passed the
dedicated PSF-role/hash audit. At least one satellite support gap remains in
78, 79, and 82 frames, respectively.

Median satellite-profile unit response is `i` 0.285/0.349/0.564 (valid n
15/14/12) and `i2` 0.612/0.409/0.411 (n 44/44/43). Paired cluster-bootstrap
intervals for amplitude shifts cross zero in both filters; response scatter
is wide. Mask growth is amplitude dependent: the median satellite-only added
mask pixels rise from 151 to 871 in `i2` and 349 to 1,423 in `i`. This is a
selected-sample pipeline sensitivity, not an astrophysical result or a
completeness correction. Poisson variance for injected light, upstream sky
transfer, extended PSF wings, and validated covariance remain unmodeled;
`science_measurement=false`.

The Astrafek analyzer now accepts `--expected-frames` (default 57), and its
paired bootstrap omits unavailable responses instead of allowing missing
support to propagate as NaN. Full tables, confidence intervals, hashes, and
limitations are in
`analysis/icl_wgl_membership/hsc_debug/icl_supported_pairs_satellite_amplitude_assessment_2026-11-06.md`.
Paired summary SHA-256
`1191414120941dc2022001ed260f277601a4af038367482383f7226cf55707f3`;
analyzer SHA-256
`2291744adf04d752c3f9ffa957541d2ed78c9ca6196ca28925a3928e11529cbb`.

## 2026-11-05: ICL-strength split in random-subtracted weak lensing

Added an Astrafek analysis joining saved 200–300 kpc ICL target-minus-control
surface-brightness contrasts to the verified random-subtracted HSC shape rows.
The `i2` ICL/WL common sample has 70 lenses; using source tomographic bins 2–4,
the median ICL split produced 35 low and 35 high candidates. Maximum-cardinality
matching with `|Δz|<=0.03` and `|Δln(richness)|<=0.25` retained 29 pairs (58
lenses), with 26 approximate equal-area delete-one blocks. Median richness was
18.97/19.47 and median ICL contrast −0.886/+6.360 nJy arcsec⁻² for low/high.
The high-minus-low `DeltaSigma` profile varies in sign with radius. Its
conditional seven-bin tangential Wald result is `χ²=27.77`, nominal `p=2.42e-4`;
cross-component `χ²=4.54`, nominal `p=0.716`. Neither nominal p-value is a
calibrated detection significance. The sample is selected and small; ICL
Poisson-noise/PSF-wing/sky/completeness, random-aperture matching, and validated
contiguous-region covariance remain unresolved. The `i` overlap has 38 lenses
and was not analyzed here. `science_measurement=false` remains.

Follow-up block-width sensitivity used the identical 29 matched pairs, source
rows, and tangential point estimates at 2.5° blocks. It yielded 54 blocks and
`χ²=12.05`, nominal tangential `p=0.099`; cross-component `χ²=5.22`, nominal
`p=0.633`. Per-bin error ratios (2.5°/5°) are 0.89–1.26. The apparent 5° joint
contrast does not pass this sensitivity, so no ICL-strength/WL relation is
claimed. The sensitivity summary is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_subtracted_icl_split_2026-11-05_03/summary.json`
(SHA-256 `809dc2362adc8a3b345fe469abdb50339c9aeb7baf1b8b45e3859a8f43464722`);
plot `icl_split_profile_2p5deg.png` SHA-256
`64aed42cef186cd7499aac16749be92d4332c75c97023f874362e12a89a19d3e`. The 5°
figure is `icl_split_profile_5deg.png` SHA-256
`1566f601e065eb70fa3ee93d63a09f67c5a003686831c41980722ac3cb7cf39d`. The
plotter hash is now `e008ee7150a7da21a67f433b3eaf27ca236a6f56834dd65165a8c55f908f90bd`.
The Astrafek comparison runner `analysis/icl_wgl_membership/hsc_debug/compare_wl_icl_split_blockwidth.py`
(SHA-256 `e9b669dd980d83362f4a08875467e2a24a9d899c397a67f8783735af1119ab4b`)
asserted identical matched pairs and bitwise-identical point profiles. Its
comparison summary is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_subtracted_icl_split_blockwidth_2026-11-05_01/summary.json`
(SHA-256 `6d9be0508c1e89aecb4a2b7f87e69764af21a9534d4449eb9cc028ab8bf0d442`).

## 2026-11-06: ICL–satellite finite-lag pixel-noise correction, 94 pairs

Extended the existing pixel-noise residual-covariance correction to the full
94-pair joint-support sample (26 `i`, 68 `i2`). An independent scope audit
verified covariance coverage for all pairs; the nested 74-pair raw values
match the prior table. Corrected partial correlations controlling for
redshift and richness are `0.081` in `i` (5° block interval
`[−0.203,0.494]`) and `0.306` in `i2` (`[0.045,0.516]`), versus raw values
`0.079` and `0.293`. The finite-lag pixel-noise correction is small. The
intervals come from approximate connected equal-area cells and are not
calibrated detections; support selection, block-width sensitivity, sky
transfer, source-Poisson noise, extended PSF, calibration, and astrophysical
covariance remain unresolved. `science_measurement=false`; no HSC query,
download, or image read occurred.

Runner `analysis/icl_wgl_membership/hsc_debug/assess_icl_satellite_pixel_noise_bias_supported_94.py`
(SHA-256 `9d38b6b692fa15f69201e9e8eae0263bef6a59e235d92e317fe89d4c8eefc8b1`),
assessment `analysis/icl_wgl_membership/hsc_debug/icl_satellite_pixel_noise_bias_supported_94_assessment_2026-11-06.md`.
Summary `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_pixel_noise_bias_supported_94_2026-11-06_02/summary.json`
(SHA-256 `e467e95192cabd411803b8bda350d5492301f508ecda4ee988d85df099119f54`).

Astrafek files: `analysis/icl_wgl_membership/hsc_debug/analyze_wl_random_subtracted_icl_split.py`
(SHA-256 `705d6389c70f1a7d6dce6ff2c07bc1517ede631b0a507757b81ef39417bb4a3f`),
`plot_wl_random_subtracted_icl_split.py` (SHA-256
`e008ee7150a7da21a67f433b3eaf27ca236a6f56834dd65165a8c55f908f90bd`), and
`wl_icl_split_assessment_2026-11-05.md`. Result summary:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_subtracted_icl_split_2026-11-05_02/summary.json`
(SHA-256 `1da415abdc45720294056e7379ceff049ecb2baff66f67921a69ec8ddc1ed495`);
figure `icl_split_profile_5deg.png` in the same directory (SHA-256
`1566f601e065eb70fa3ee93d63a09f67c5a003686831c41980722ac3cb7cf39d`). No HSC
query/download or image-pixel read occurred. README and this handover record
the reproducible path and caveats. The preliminary result directory `_01` is
superseded by `_02`, which includes the cross-component null and final runner
hash.

## 2026-11-05: ICL–satellite association block-width sensitivity

Added `analysis/icl_wgl_membership/hsc_debug/analyze_icl_satellite_supported_94_blockwidth.py`
to recompute the nested 74/94-pair partial correlations under connected sky-cell
widths 2.5°, 5°, and 10° (5,000 block-bootstrap draws per case). The full 94-pair
`i2` estimate remains `r=0.293`; 95% intervals are `[0.002,0.504]`,
`[0.048,0.498]`, `[0.056,0.496]`. The narrower 74-pair `i2` interval crosses
zero at 2.5°. The paired 94-minus-74 scope-change intervals include zero for
both filters at all widths. This is a selected-sample sensitivity, not a
population detection; equal-area cells are not validated contiguous HSC
regions. `science_measurement=false`; no HSC query, image read, or download.

Assessment: `analysis/icl_wgl_membership/hsc_debug/icl_satellite_supported_94_blockwidth_assessment_2026-11-05.md`.
Runner SHA-256 `dfd845145c3e5c752aa9bf8ab0e0f77d7e99ed0bcd8364b197de51812388945b`;
plotter SHA-256 `5115445514c12ff47d32edbcc4cb54f863c632bf83b436ff6c2b44ceb86e83cf`.
Summary `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_supported_94_blockwidth_2026-11-05_01/summary.json`
SHA-256 `7aad81051d46b88804db2923cf0de38e85e62384fa8d82455839b827e96c7227`;
figure `icl_satellite_94_blockwidth_clean.png` SHA-256
`b01765d7a6aab05bf3e60e5ea9f518f65317bad7c128b8cb476890fb7ec45a93`.

2026-11-06 follow-up: extended the 94-pair finite-lag pixel-noise correction
to coordinate-connected 2.5°, 5°, and 10° spatial blocks. The corrected
partial-correlation point estimates remain `i=0.081` and `i2=0.306`; the
95% intervals are `i: [-0.359,0.610], [-0.195,0.511], [-0.501,0.524]` and
`i2: [-0.004,0.523], [0.050,0.519], [0.062,0.512]` at those widths. The
smallest-width `i2` interval marginally includes zero, so this is
block-sensitive near the boundary and is not a detection. The bootstrap
uses approximate equal-area connected cells, not validated contiguous HSC
fields; `science_measurement=false`. No new HSC query, download, or image
read occurred. Runner SHA-256
`4bd020428ba8ff2dc8ccb0a5bfa961e044f71985bdb79e31ff011da7138af7da`; summary
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_pixel_noise_bias_supported_94_2026-11-06_04/summary.json`
SHA-256 `965d0e81e78ad8d27ebfc7726a9ae170828b824fab41d29292398fbdf2940725`.

## 2026-11-06: frame-level calibration metadata export gate

The coadd header audit still finds no explicit gain/exposure keys in 1,598
target/control images, so a physically motivated source-Poisson term cannot be
derived from the coadd headers alone. An existing count-only-approved parent
query maps all 1,598 centers to patches and yields an exact estimate of 67,071
contributing CCD-frame membership rows. The separately prepared metadata plan
returns 21 fields (frame/visit/CCD, original filter, exposure, four amplifier
gains, FLUXMAG0, sky RMS, MJD, and center/patch identity), in 16 SQL batches,
estimated 13–27 MB. It excludes images, shape rows, photo-z, and object rows.

Added
`analysis/icl_wgl_membership/hsc_debug/submit_icl_frame_metadata_rows.py`.
Read-only `--check` passed against exact plan SHA-256
`49aa4b7601dd42243209d7b3b394a4dd3b4cbf00d87d28cbb9467553c29f4ce3`, 16 SQL
hashes, and the parent count audit. Submission has not occurred; this separate
metadata export requires explicit approval of that exact plan. Even if
returned, frame-level gains and sky statistics do not provide per-pixel coadd
weights or amplifier-region mapping and cannot alone calibrate the source
Poisson variance. `science_measurement=false` remains. Official HSC references:
PDR3 [processing](https://hsc-release.mtk.nao.ac.jp/doc/index.php/processing__pdr3/),
[database](https://hsc-release.mtk.nao.ac.jp/doc/index.php/database__pdr3/),
and [available data products](https://hsc-release.mtk.nao.ac.jp/doc/index.php/available-data__pdr3/).

## 2026-11-06 continuation: diffuse-injection linearity and approved frame metadata

The user approved submitting the exact per-center HSC frame-metadata plan
(SHA-256 `49aa4b7601dd42243209d7b3b394a4dd3b4cbf00d87d28cbb9467553c29f4ce3`).
On `lageunha`, the first batch returned original frame filter values `i` and
`i2`, while the submitter incorrectly expected `HSC-I`/`HSC-I2`. The first
job was checkpointed; no further SQL was submitted during the stop. The
validator was fixed to accept the original filter codes and the exact plan was
resumed from its saved job ID. All 16 approved batches completed: 67,071 rows,
1,598 centers (800 target, 798 control), with no nulls in exposure, gain1–4,
FLUXMAG0, sky RMS, or MJD. No images, shape rows, photo-z, object rows, or
credentials were returned/saved. The independent offline auditor rechecked
all 16 batch hashes, center/patch identities, unique frame membership and
parent row counts. In the 798 matched pairs, median target-minus-control
differences are zero for frame, visit, CCD counts, summed exposure, median
sky RMS, median FLUXMAG0, and each amplifier gain. This is a patch-level
summary only; no per-pixel coadd weights or amplifier maps exist in this
export. It cannot yield source-Poisson variance or promote `science_measurement`.

Metadata submitter:
`analysis/icl_wgl_membership/hsc_debug/submit_icl_frame_metadata_rows.py`
(SHA-256 `bacf170f85b3f2fc47310ec10457201981f7c8706677455036f36da44ecc43bb`).
Offline auditor:
`analysis/icl_wgl_membership/hsc_debug/audit_icl_frame_metadata_row_results.py`
(SHA-256 `d0340160ed0f7c2ed7c8e7570b2687ed1009917380d8c8dda1625607d5922336`).
Audit summary:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_frame_metadata_row_audit_2026-11-06_01/summary.json`
(SHA-256 `0f361fcb496fbd7be7c1059f638414f6900535f73a9ac48696030480f400a237`).
CAS result summary:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_frame_metadata_row_results_2026-11-06_01/summary.json`
(SHA-256 `3925f74be34b19cf16ea0183055ba4d74d91568f66150c569e2d48b6bd2f7b53`).

Separately, diffuse-ICL injection recovery was measured on the same 137
PSF-convolved frames at input amplitudes 0.5, 1.0, and 2.0 nJy arcsec⁻² at
100 kpc, with satellite amplitude fixed at 1.0 and detection masks regenerated.
All runs completed 137/137 with zero frame errors. Median unit response is
`i: 0.240, 0.226, 0.224` and `i2: 0.225, 0.226, 0.226`; paired 5° connected
block-bootstrap changes from the 1.0 run include zero in both filters. Absolute
median recovery remains about 23%, and masks change in many frames. This
supports approximate amplitude stability within the tested range but not an
absolute flux correction. Injected Poisson noise, coadd weights, upstream sky
transfer and extended PSF wings remain absent. Keep
`science_measurement=false`.

Analysis tool:
`analysis/icl_wgl_membership/hsc_debug/analyze_icl_diffuse_injection_amplitude_grid.py`
(SHA-256 `8eef6823016e35538a059b0c125d12e88b98e9f8307c52cc21f88ceb04af7448`).
Assessment:
`analysis/icl_wgl_membership/hsc_debug/icl_diffuse_injection_amplitude_assessment_2026-11-06.md`.
Paired output summary:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_diffuse_injection_amplitude_grid_supported137_2026-11-06_01/summary.json`
(SHA-256 `538656e4b3c5141282ddd24c16ad8ce0d20454692228f38b437ad996b29fcd96`).

The approved frame metadata was also deduplicated for a potential pixel-level
variance follow-up. It contains 12,839 distinct patch/filter/visit tuples and
56,731 patch/filter/visit/CCD tuples; those are not cutout request counts. The
authenticated PDR3 DAS Cutout manual documents `type=warp` tar downloads with
image, mask, and variance layers. A bounded 2′×2′, one-center, nine-visit pilot
is estimated at about 55 MB, but no query or download was made; it needs its
own approval. Details and the scope are in
`analysis/icl_wgl_membership/hsc_debug/icl_warp_product_scope_assessment_2026-11-06.md`.
Even warp variance plus the approved frame metadata does not yet establish a
calibrated injected-source Poisson model. Keep `science_measurement=false`.

## 2026-11-06: i-band ICL-split weak-lensing follow-up

The existing `analyze_wl_random_subtracted_icl_split.py` was also run on the
saved `i` common-support subset at block widths 2.5°, 5°, and 10°. It has 37
usable common lenses; the same ICL median split and redshift/richness matching
retains 10 pairs. The exact pair list and seven-bin point profile are identical
across all three outputs. Tangential conditional Wald χ² is 15.79/161.92/101.28
and cross-null χ² is 24.70/90.97/105.50 at 2.5°/5°/10°. Only 9–24 approximate
blocks support the 10-pair comparison; covariance and the nominal diagnostics
are too unstable for a detection claim. The profile is radially mixed and
cross-null is large at all widths. This `i` result does not strengthen the
ICL–lensing association; `science_measurement=false` remains.

Assessment:
`analysis/icl_wgl_membership/hsc_debug/wl_icl_i_band_blockwidth_assessment_2026-11-06.md`.
All results are under
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_subtracted_icl_split_i_2026-11-06_{01,02,03}/`.
Summary SHA-256 values for 5°/2.5°/10° are
`61bbf1d605ed15b7666aee3490ed91051fa22f7452956e839b3d4b41535733ff`,
`5c62835fcae76a8e9251f05218c03c73629579bf15ef8460257aabfc1a1d5baf`, and
`ad38c8b5d7747d264199eaf5a420e48206be8d5c1054b82971fafc98d822007f`.
No new HSC request or image-pixel read occurred.

An additional offline region-balance audit mapped the i and i2 common/matched
lens IDs to the approximate tract components from the existing 800-lens parent
audit. The i common sample has 38 lenses split 18/20 across two components;
the 10 matched high/low pairs are balanced 3/7 and 4/6, with five same-region
and five cross-region pairs. The i2 sample has 70 lenses split 48/18/4 across
three components; 29 matched pairs have 17 same-region and 12 cross-region
pairs. The i cross-null therefore is not explained by a gross high/low lens
count imbalance, but the audit does not assign source/random positions and
cannot identify its cause. Region labels are approximate, not official HSC
fields. Astrafek runner
`analysis/icl_wgl_membership/hsc_debug/audit_wl_icl_split_tract_region_balance.py`
(SHA-256 `768e18989a2907911697d6c25d01a53022bc3572a9a4dc2d06d7ec3ad1b5cbc1`);
summary hashes i `8cba95a740c5d332b8c313c5687e58289652042497a5cc977061032203657130`,
i2 `adb313d390c9def350255dd8719d46ea9c536b1864f142b73283455421c0f6d1`.
Both outputs reside under `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/`.

To locate the i cross-null, a second Astrafek script reapplied the exact saved
HSC Y3 source cuts and stacked the 10 matched low/high ICL groups' target and
five-random-center rows separately by the approximate tract component of the
lens/center. The selection reproduces 53,399 used target pairs and 241,258 used
random pairs. All used sources are assigned to the same component as their
lens/center. The high-minus-low random-subtracted cross profile is oscillatory
in both occupied components and has same-sign, large first two bins; neither
cross-component source leakage nor one component alone explains the pattern.
No regional covariance is available, and each component has only 3–7 lenses
per group; do not infer the source of the cross signal.

Runner:
`analysis/icl_wgl_membership/hsc_debug/audit_wl_icl_split_source_region_cross_null.py`
(SHA-256 `961874ffdeb80efa9abe09c519adac680ea35adf7bdc83adceff5043fbfd4fdb`).
Summary:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_icl_split_source_region_cross_null_2026-11-06_01/summary.json`
(SHA-256 `fb90c982da784dff60b270efdfe986d0fbe5b8f3a656e15dd5bc14629810862c`).
No new HSC request or image-pixel read occurred; `science_measurement=false`.

The saved endpoint profiles were then decomposed into target and random-center
high-minus-low terms. In the inner two bins, target cross differences are
large in both components (component 0: −142.1/+120.2; component 1:
−114.4/+159.5, in `10^12 M_sun Mpc^-2`); random terms modify these and reinforce
component 1. A post-hoc `R>=1 Mpc` radial cut gives conditional cross-null
nominal p values 0.307/0.303/0.107 and tangential-null values
0.0093/6.63e-5/5.27e-5 for 2.5°/5°/10° blocks. The outer tangential contrast
is a candidate for further confirmation only: the cut was post-hoc, there are
10 matched pairs and 9–24 approximate blocks, and the covariance is not
validated. In the larger i2 29-pair split, the same `R>=1 Mpc` cut gives
tangential conditional nominal p=0.172/0.114 and cross p=0.639/0.550 for
2.5°/5° blocks; the outer i tangential candidate is not reproduced.
`science_measurement=false` remains.

Astrafek endpoint decomposition:
`analysis/icl_wgl_membership/hsc_debug/audit_wl_icl_split_target_random_cross_decomposition.py`
(SHA-256 `051bc67f3be1f23f0134d7509287fa201b2da7b2689d86302fbbf1b3e2b4f42c`),
summary SHA-256 `ffc453ac4b78e7f25a354d19469b31411856cc6d4e0ab41274f6940f331fe5af`.
Radial-cut audit:
`analysis/icl_wgl_membership/hsc_debug/audit_wl_icl_split_radial_cut_nulls.py`
(SHA-256 `00a2f24b6a4c7e63b13f84c393d4778224db03b07a4ce5a6316f3b3f70717cf2`),
summary SHA-256 `583426fef8b60c011ea392dd80c3a7783d1572fa06006d72b729bebea1eade45`.

The next ICL/WL step tested continuous ICL against per-lens random-subtracted
weak-lensing profiles for the exact 70-lens i2 common-support sample. The only
nominal tangential-bin result is a negative 2–3 proper Mpc slope
`−12.3×10^12 M_sun Mpc^-2` per 1-SD ICL contrast (`HC3 p=0.021`); its
seven-bin Bonferroni/BH-adjusted minimum is `0.149`. Cross slopes are
consistent with zero. This outer-bin candidate is not reproduced by the prior
i2 median split. Radial support ranges from 58 to 70 lenses and covariance
uses 27 approximate equal-area cells, so retain `science_measurement=false`.
No HSC query was made; only already-approved local source/random rows were
read.

Astrafek runner:
`analysis/icl_wgl_membership/hsc_debug/analyze_wl_continuous_icl_relation.py`
(SHA-256 `9edaec5bf87da2e8739788d62275892ec904f255046bf2081b49b4d10523e6c8`).
Assessment:
`analysis/icl_wgl_membership/hsc_debug/wl_continuous_icl_relation_assessment_2026-11-06.md`.
Results:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_continuous_icl_relation_i2_2026-11-06_02/`
(summary SHA-256 `12c9456799e2671e3a0452c09a52465bdc6434c5b7546659c4bf21ac340f5fea`).
The radial tangential/cross regression figure is in the same folder as
`wl_continuous_icl_relation_i2_v2.png` (SHA-256
`f8225b95111ab29a040edf97826eb3a6469fa0f87ed66c81b27078273c82ec08`) and
`wl_continuous_icl_relation_i2_v2.pdf`; reproduce it with
`analysis/icl_wgl_membership/hsc_debug/plot_wl_continuous_icl_relation.py`.

An offline HSC warp pilot plan is now prepared from the already approved
frame-metadata export. Scope is one `pdr3_wide/HSC-I` `type=warp` request at
target `6966`, tract/patch `9703/306`, 2′×2′, image+mask+variance; metadata
maps it to 43 frame rows, nine visits and 40 CCDs. Estimated full-square
uncompressed planes are 52.65 MiB; the download runner has a 70 MiB hard cap.
Plan SHA-256 is
`a436df245a33bda09cd5599e564b2b9302e5d0d4c006b59fcf6fab1c2c33a65e` at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_warp_pilot_plan_2026-11-06_02/plan.json`.
No warp request has been submitted. `download_icl_warp_pilot.py` requires a
separate explicit approval of that exact plan hash and performs no request by
default. Larger warp transfers would need their own approved plan.

The warp downloader remains plan-only unless the user separately approves the
exact plan SHA. Its offline FITS audit now checks image/mask/variance shapes,
per-HDU WCS alignment, requested-center coverage, 0.168″/pixel scale, finite
image fraction, mask occupancy, and finite/positive variance fractions. It
labels a visit-complete response with failed checks
`complete_with_qa_warnings`; it has not been run against a downloaded warp.
Current downloader SHA-256:
`766dbedb2b76951307b6c10a9d4be309f96ba1f509d6e09073c0716a6c71e291`.

The 27-lens i2 triple-overlap diagnostic is now recorded in
`analysis/icl_wgl_membership/hsc_debug/icl_satellite_wl_triple_overlap_assessment_2026-11-06.md`.
Its Astrafek runner is
`analysis/icl_wgl_membership/hsc_debug/analyze_icl_satellite_wl_triple_overlap.py`;
metrics and joined rows are under
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_wl_triple_overlap_2026-11-06_01/`.
ICL–satellite partial `r=0.067` (block-bootstrap 95% interval
`[-0.144, 0.293]`); all ICL–WL intervals cover zero. Two unadjusted
satellite–WL intervals exclude zero among 29 reported relationships, so the
result remains exploratory and `science_measurement=false`. No new HSC query
was submitted.

The existing supported-94 target/control coadds were also remeasured in ten
radial annuli for delivered VARIANCE-plane balance. Runner:
`analysis/icl_wgl_membership/hsc_debug/analyze_icl_variance_radial_balance.py`;
assessment:
`analysis/icl_wgl_membership/hsc_debug/icl_variance_radial_balance_supported94_assessment_2026-11-06.md`;
outputs:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_variance_radial_balance_supported94_2026-11-06_01/`.
Three of 20 unadjusted target/control radial comparisons have pointwise
block-bootstrap intervals excluding equality; all ICL–variance correlations
include zero. This does not estimate source-Poisson noise or the radial ICL
estimator covariance. It read existing coadd variance/mask planes only, made
no HSC request, and leaves `science_measurement=false`.

The variance comparison was expanded to all same-filter acquired pairs
(783 total, 302 `i`, 481 `i2`) and received a joint max-t correction over 20
filter/radius contrasts. Three `i2` simultaneous intervals are below equal
variance at 75–100, 100–150, and 200–300 kpc (ratios 0.961/0.955/0.961;
intervals 0.924–0.999/0.915–0.996/0.934–0.988); none remain in `i`. The
resampling partition has 38 approximate blocks and one large 107-pair block,
so this is still only a delivered-variance diagnostic. Astrafek record:
`analysis/icl_wgl_membership/hsc_debug/icl_variance_radial_balance_full783_assessment_2026-11-06.md`.
Results are in
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_variance_radial_balance_full783_2026-11-06_01/`
and corrected simultaneous intervals in
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_variance_radial_balance_familywise_2026-11-06_02/`.
The `_01` familywise directory is explicitly invalidated for pooled-filter
metric construction. No HSC query/download; `science_measurement=false`.

### 2026-11-06 variance-adjusted ICL–satellite association

Added offline sensitivity runner `analysis/icl_wgl_membership/hsc_debug/analyze_icl_variance_adjusted_satellite_assoc.py` and assessment `analysis/icl_wgl_membership/hsc_debug/icl_variance_adjusted_satellite_assoc_assessment_2026-11-06.md`. On the saved 94 supported pairs, the `i2` redshift/richness partial correlation is 0.293 (paired 5° connected-block 95% interval 0.043–0.498); adding log target/control variance ratio changes it by −0.008 (paired change interval −0.078–0.034). Adding a coarse mean-depth proxy changes it by −0.032 (−0.148–0.027). The `i` sample has N=26 and broad intervals. These selected-sample model changes include zero and do not establish that the delivered-variance difference explains the ICL–satellite association. The coadd variance plane is not total noise and this analysis does not estimate source-Poisson noise or ICL-estimator covariance. Results: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_variance_adjusted_satellite_assoc_2026-11-06_01/`; summary SHA-256 `8731105a3d5a9c0430ae2f688442beea640277d06af7d70ff03b1f48e9eb0644`. No HSC query/download; `science_measurement=false`.

### 2026-11-06 satellite-support IPW with delivered-variance covariates

Added `analysis/icl_wgl_membership/hsc_debug/analyze_icl_satellite_ipw_variance_sensitivity.py` to test whether target/control delivered VARIANCE balance predicts which of the 224 ICL-supported pairs has a 50–75 kpc satellite profile. The 94 selected pairs were compared under three support-propensity models using the same 2,000 connected 5° block resamples. For `i2`, the IPW partial r is 0.227 baseline, 0.266 with log variance ratio, and 0.193 with the ratio plus mean log variance; paired changes are +0.039 (95% interval −0.069–0.162) and −0.034 (−0.163–0.092), both spanning zero. The richer model reduces `i` effective N to 13.5/26, so that band is unstable. IPW missing-at-random is unverified and delivered variance is not total noise; no correction or science promotion follows. Assessment: `analysis/icl_wgl_membership/hsc_debug/icl_satellite_ipw_variance_sensitivity_assessment_2026-11-06.md`. Results: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_ipw_variance_sensitivity_2026-11-06_01/`; summary SHA-256 `a1e2972d728bb88590dff3f3ff5197e1f10ea3bbeeec7a1aa0d93f6718fc6b2c`. No new HSC request/download; `science_measurement=false`.

### 2026-11-06 joint IPW and finite-lag pixel-noise association

Added `analysis/icl_wgl_membership/hsc_debug/analyze_icl_satellite_ipw_pixelnoise_joint.py` to combine the support-selection weights with the independently scope-verified finite-lag pixel covariance for the same 94 ICL/satellite pairs. The baseline `i2` corrected IPW partial correlation is 0.238 (2,000 connected 5° block interval −0.023–0.446); adding variance ratio gives 0.278 (−0.018–0.495), while adding variance ratio plus mean depth gives 0.201 (−0.054–0.410). Both paired model-change intervals include zero. The `i` interval is very broad, and the richest model leaves effective N=13.5/26. The covariance table covers all 94 pairs, but its estimation uncertainty is not propagated, source Poisson/effective gain and sky/PSF covariance are missing, and IPW missing-at-random is unverified. No science promotion: `science_measurement=false`. Assessment: `analysis/icl_wgl_membership/hsc_debug/icl_satellite_ipw_pixelnoise_joint_assessment_2026-11-06.md`. Output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_ipw_pixelnoise_joint_2026-11-06_01/`; summary SHA-256 `761e3448fd5d01b2937c62cc3668bc8a2f7e3af0bcdc74a3a2c821c193d3d01d`. No new HSC request/download.

### 2026-11-06 ICL-strength split satellite profiles

Added a cluster-level seven-bin high/low ICL satellite-profile comparison for
the existing 94 supported ICL/satellite pairs (26 `i`, 68 `i2`). The OLS
high-minus-low contrast is adjusted for redshift and log richness; 10,000
restricted Webb wild-cluster bootstrap-t replicates use one joint 5-degree
connected-pair partition across both filters, with Holm and max-T correction
across 14 contrasts. All CR1 95% intervals include zero, Holm p=1.0 throughout,
and minimum max-T p=0.520. Positive i2 profile differences remain unresolved;
the i profile is noisy and changes sign with radius. The analysis is
post-hoc and selected, and it does not establish an ICL–satellite relation.
Runner `analysis/icl_wgl_membership/hsc_debug/analyze_icl_strength_split_satellite_profiles.py`
(SHA-256 `2f6208e7f34504f3dc7a5c76d2d202b4c34ae2bdf2ed05ec9e655d8e6fac512b`),
plotter `analysis/icl_wgl_membership/hsc_debug/plot_icl_strength_split_satellite_profiles.py`
(SHA-256 `dfed877d3140fabffa5c3c64dad537293248f98a71c723bd7f3effe0f338f036`),
assessment `analysis/icl_wgl_membership/hsc_debug/icl_strength_split_satellite_profiles_assessment_2026-11-06.md`.
Machine results and figure are under
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_icl_strength_stacked_profiles_2026-11-06_01/`.
No HSC request/download; `science_measurement=false`.

### 2026-11-06 joint radial ICL–satellite relations

Added `analysis/icl_wgl_membership/hsc_debug/analyze_icl_satellite_radial_relations.py`
for continuous one-SD ICL slopes and median-split contrasts in each of seven
satellite-profile radii and both filters. It uses the existing 94-pair selected
sample, adjusts for redshift and log richness, and applies shared-block Webb
bootstrap with Holm/max-T correction across all 28 tests. The smallest raw
bootstrap p is 0.0369 (continuous i2, 50–75 kpc) but max-T p=0.770; minimum
max-T p=0.190 and Holm p=1.0 throughout. Three negative pointwise CR1 intervals
in inner i bins do not pass the joint bootstrap. No ICL–satellite relation is
established; `science_measurement=false`. Assessment:
`analysis/icl_wgl_membership/hsc_debug/icl_satellite_radial_relations_assessment_2026-11-06.md`.
Plotter: `analysis/icl_wgl_membership/hsc_debug/plot_icl_satellite_radial_relations.py`.
Results: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_radial_relations_joint_2026-11-06_01/`.
No new HSC request/download.


## 2026-11-06 continuation: approved frame metadata audit and weak-lensing boost sensitivity

The exact approved HSC frame-metadata plan (`49aa4b7601dd42243209d7b3b394a4dd3b4cbf00d87d28cbb9467553c29f4ce3`) had already completed all 16 CAS batches before this continuation. The user again approved this exact plan; inspection confirmed the prior completed submission, so no duplicate CAS jobs were sent. An independent offline audit completed at `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_frame_metadata_row_audit_2026-11-06_03/summary.json`: 67,071 rows, 1,598 centers (800 target, 798 control), 798 matched pairs, and zero nulls in exposure, four gains, FLUXMAG0, sky RMS, and MJD. The auditor made no network request. This export has no pixel coadd weights or amplifier map and therefore cannot calibrate source-Poisson noise; `science_measurement=false` remains.

The existing weak-lensing boost analysis was extended with a paired-jackknife sensitivity profile, `B_W(R) * DeltaSigma_lens-random`, using the exact response-adjusted S19A estimator-weight ratio, not the raw count boost. Across seven bins, the candidate profile is 12.5–16.3% above the saved profile; central values agree under 2.5- and 5-degree spatial blocks. This is not an applied physical correction because member dilution is not separated from magnification, depth, mask, PSF, sky, or source-selection effects. The saved profile is unchanged; `boost_correction_applied=false` and `science_measurement=false`. Assessment: `analysis/icl_wgl_membership/hsc_debug/wl_boost_profile_sensitivity_assessment_2026-11-06.md`; runner: `analysis/icl_wgl_membership/hsc_debug/analyze_wl_source_density_boost_blocks.py`; output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_source_density_boost_spatial_2026-11-05_06/`.


## 2026-11-06: forced CAMIRA-position coverage sensitivity

The forced-profile runner `analysis/icl_wgl_membership/hsc_debug/measure_forced_satellite_profiles.py` now accepts a `--coverage-floor` option; default remains 0.60. A separate run unmasked other unmatched CAMIRA positions but retained the 60% floor. This gave no support increase: the 50–75 kpc candidate count remained nine; across the 20 retained position-bin measurements, maximum target-coverage change was zero and maximum control-coverage change was 0.00184. Lowering the floor to 0.20 (still requiring at least three rotated controls) yielded 65 forced profiles at 50–75 kpc (25 `i`, 40 `i2`) from 3,197 candidates; 123 candidates passed any bin and 263 position-bin rows were retained. All 721 eligible target frames completed. Forced-only `i2` support is 21 pairs/17 blocks, `r=0.147`, 95% block interval [−0.485, 0.425]; `i` has only eight pairs. Combined support is 29 `i` and 74 `i2`, with `i2 r=0.257` and unpaired spatial-block interval [0.092, 0.435]. Paired comparison on common `i2` support gives 20%-minus-60% change +0.003 [−0.070, 0.126]; adding the six pairs gives −0.025 [−0.127, 0.079]. Neither comparison resolves a change. The positive unpaired interval is post-hoc, unadjusted, and not a physical detection; the CAMIRA coordinates are not confirmed galaxies and low coverage may be non-random. Keep `science_measurement=false`.

Assessment: `analysis/icl_wgl_membership/hsc_debug/icl_satellite_forced_coverage_sensitivity_assessment_2026-11-06.md`; paired analysis: `analysis/icl_wgl_membership/hsc_debug/analyze_icl_satellite_forced_coverage_sensitivity.py`; runner: `analysis/icl_wgl_membership/hsc_debug/measure_forced_satellite_profiles.py`. Results: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_forced_coverage_sensitivity_2026-11-06_02/summary.json` (SHA-256 `14367972f95c554c989509f7f1152f2d7fffe09e28c8bcb45dd448cf595ca7a2`) and `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_forced_outer_profile_unmatched_camira_unmasked_cov20_2026-11-06_01/summary.json` (SHA-256 `8de336e6d22da5ef1679d60a73ffd15e299f1c6f5303ced6bd6d5a1545e5852c`). No HSC query or download occurred.


## 2026-11-06 forced-position radial-support selection

Audited all 3,197 CAMIRA positions without SEP detections against their measured target coverage and the 60%/20% forced-profile outputs. At the 20% floor, 50–75 kpc profiles are retained for 25 `i` positions and 40 `i2` positions, but none is from a candidate within 100 kpc of the host cluster center. Forty-one of the 65 profiles are at 300–500 kpc. Parent position-radius medians are 200/195 kpc (`i`/`i2`) versus retained-profile medians 388/339 kpc. Any profile is retained for 18/282 candidate-bearing `i` clusters and 32/439 `i2` clusters. This is a positivity failure for correction to the full parent radial distribution. Do not pool forced profiles into the primary detected-satellite stack; keep the 20%-cut `i2 r=0.257` as an uncorrected sensitivity result. Assessment `analysis/icl_wgl_membership/hsc_debug/icl_satellite_forced_position_selection_assessment_2026-11-06.md`; runner `analysis/icl_wgl_membership/hsc_debug/audit_icl_satellite_forced_position_coverage_selection.py`; output `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_forced_position_coverage_selection_2026-11-06_01/summary.json` (SHA-256 `6d6650a0ff81d203d11ad012ba82867b9150904af38f511212d1e25931ce9e28`). No HSC request/download; `science_measurement=false`.

## 2026-11-06 continuation: contributing-member count sensitivity

Added Astrafek runner `analysis/icl_wgl_membership/hsc_debug/analyze_icl_satellite_membercount_adjustment.py` and assessment `analysis/icl_wgl_membership/hsc_debug/icl_satellite_membercount_adjustment_assessment_2026-11-06.md`. The fixed 94-pair sample was checked against the per-member 50–75 proper-kpc profile table. The baseline partial correlation controls redshift and linear CAMIRA richness, matching the primary 94-pair analysis. Adding log contributing-member count changes the `i2` estimate from `r=0.293` to `0.295`; the paired shared-block change is `+0.002` with 95% interval `[-0.023, 0.044]` from 9,989 valid draws out of 10,000 requested. For `i`, the estimate changes from `0.079` to `0.088`, with change `+0.009` `[-0.081, 0.192]`. Member count does not show a resolved partial association with ICL or satellite excess in either filter. This is a selected-sample precision sensitivity, not a completeness correction. Input tables were local; no HSC request or image-pixel read occurred, and `science_measurement=false` remains. Results: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_membercount_adjustment_2026-11-06_02/`. Runner SHA-256 `13a6c6b6168a65e5dc0f71e40a963f0298c62c04a14266683aae3b466a39f16b`; summary SHA-256 `12420feadba8d83badf3c1cb8446652c4bf3c3041bcbc70eae84e2f94011b485`.

## 2026-11-06 continuation: long-lag ICL aperture covariance and random-point auditor

Current data inspection confirms that the 17×17 residual-covariance kernel has already been contracted through the exact 200–300 proper-kpc masks for the same 224 supported target/control pairs. The independent-frame noise estimates are 0.0793 nJy arcsec⁻² in `i` and 0.0745 in `i2`, only 0.8% and 1.3% above the exact-aperture 9×9 estimates. The finite-lag terms are 7.6% and 7.0% of the saved spatial-block interval half-widths. These results remain a stationary-kernel noise floor, not total uncertainty or a flux correction. Runner `analysis/icl_wgl_membership/hsc_debug/analyze_icl_long_lag_aperture_covariance.py`; assessment `analysis/icl_wgl_membership/hsc_debug/icl_long_lag_aperture_covariance_assessment_2026-11-06.md`; summary `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_aperture_covariance_200_300_2026-11-06_01/summary.json` (SHA-256 `ebb78c5202ff60c5484b5b6aff69d61050bf57d11b005f548e10e8fdc6fb1be4`). No HSC query or image download occurred; `science_measurement=false` remains.

Added `analysis/icl_wgl_membership/hsc_debug/audit_icl_random_pixel_variance_results.py` to validate a future completed random-point result against the pinned exact plan, checkpoint, all batch hashes, annulus membership, returned flags, center identities, and target/control mapping. The tool writes 1,598 center summaries and 798 paired summaries. Syntax parsing passed; no HSC results exist because plan SHA-256 `eeef98ef2210337f714e36d8a4afcd352fe312a536c90e8e236ff9f843f53e52` remains unsubmitted pending the user's explicit approval. Auditor SHA-256 `05352d72efbd95d97ce540e5008970b4ff75fcb7fe121efb467b59e913bb80cc`. This export still cannot provide a coadd weight map or complete source-Poisson/correlated-noise correction.


## 2026-11-06 continuation: host-centric satellite-position mix sensitivity

On the fixed 94-pair ICL–satellite sample, adding the median host-centric radius of contributing members to the redshift- and richness-adjusted model changes `i2` partial `r` from 0.293 to 0.322 (paired change +0.030; shared 5-degree block-bootstrap 95% interval −0.019–0.101). In `i`, `r` changes from 0.079 to 0.029 (change −0.050; interval −0.330–0.044). All position-mix and model-change intervals include zero. This sensitivity does not explain or validate the selected-sample `i2` association and does not correct for population incompleteness. Source Poisson, PSF wings, sky transfer, and full image covariance remain unresolved; `science_measurement=false`.

Runner `analysis/icl_wgl_membership/hsc_debug/analyze_icl_satellite_clustercentric_position_adjustment.py`; assessment `analysis/icl_wgl_membership/hsc_debug/icl_satellite_clustercentric_position_adjustment_assessment_2026-11-06.md`; output `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_clustercentric_position_adjustment_2026-11-06_01/` (summary SHA-256 `78317905c3d0b2b0ada48c8b39a27ae564caaeb0bd744f365c12fedc0a7d8971`). The analysis used existing local tables; no HSC request, download, or image-pixel read occurred.


## 2026-11-06 continuation: S20A coverage metadata query rejected by CAS

The already approved metadata-only coverage request is 2,596 candidates × nine probes (23,364 positions), 52 batches, and 14 mosaic/patch-QA columns. The first submission returned HTTP 406. A missing `ON` condition on the mosaic join was repaired, with candidate centers and probe inventory unchanged byte-for-byte; the corrected plan passed offline validation but its first batch also returned HTTP 406. Both checkpoints show the first batch pending with no CAS job ID, no result files, and no later batches. The server rejection cause remains unknown. Do not retry until CAS job state and the rejection reason are checked. No images, source rows, shapes, or photo-z were requested; `science_measurement=false`. Assessment: `analysis/icl_wgl_membership/hsc_debug/icl_expansion_coverage_submission_assessment_2026-11-06.md`. Corrected plan SHA-256 `974d4e3db450703481ac2b004d37860806586900097e626c5b456702a9e0346f`; generator and submitter hashes are recorded in the assessment.


## 2026-11-06 continuation: consolidated ICL and weak-lensing evidence

Created `analysis/icl_wgl_membership/hsc_debug/build_icl_wl_current_evidence_report.py` and generated `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_03/icl_wl_current_evidence_report.md`. It validates and hashes 13 latest Astrafek component summaries. On 94 ICL/satellite common-support pairs, i2 `r=0.293` raw and `r=0.306` after finite-lag noise; support-weighted plus noise `r=0.238` has a 95% interval containing zero. On the 435-lens random-subtracted weak-lensing stack, deleting the largest of three occupied tract components shifts inner bins by 26.8–28.6%. The ICL–lensing slope family has no adjusted detection. Both analyses remain exploratory and `science_measurement=false`. No HSC access or image-pixel read occurred. Report SHA-256 `fb061edb26f6ea25666eefcc644b464bc4b54eb68239f649a920ca8145c4439f`; summary SHA-256 `1b6fd47fea010fbb8622def7326e9ca62e467aff8e98dc1f713b9f241ea6644e`; runner SHA-256 `bc0a9dd6b85bb87eb59eaa6eee5ded8437fe3b40442cca2e8c9db8abf8c6041d`. The separate S20A metadata expansion query remains CAS HTTP 406 with no returned results.

## 2026-11-06 continuation: independent hybrid-PSF shell-exclusion audit

Re-ran the Astrafek offline summarizer against the already completed hybrid-PSF satellite-shell-exclusion injections. It verified the pinned 74-row association table, all per-frame checksums, 12/12 `i` and 36/36 `i2` frames, and 48 unique target/filter keys. No HSC request, download, or observed-image read occurred. The diffuse-to-ICL paired response change has q16–q84 `[-0.00167,+0.00102]` (`i`) and `[-0.00110,+0.00132]` (`i2`); neither supports a scalar correction. The satellite-envelope response is broad across frames and likewise does not support a common correction. Excluding supported satellite-shell pixels removes 62,858 (`i`) and 257,890 (`i2`) pixels; median 200–300 kpc coverage falls from 0.681 to 0.618 and from 0.652 to 0.611, respectively. This support loss accompanies the small decrease in the synthetic point-source wing fraction. This is an independent audit of model sensitivity, not observational validation; `science_measurement=false`, `psf_approved=false`, `observational_validation_passed=false`, and `flux_covariance_available=false` remain. Assessment: `analysis/icl_wgl_membership/hsc_debug/icl_satellite_shell_exclusion_hybrid_psf_audit_2026-11-06.md`. Runner: `analysis/icl_wgl_membership/hsc_debug/summarize_icl_hybrid_shell_exclusion.py` (SHA-256 `f94191ddee4b9ad163a0ba02587cb67d885e65d359a2f9cb024baad57a2865b8`). Summary: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_hybrid_shell_exclusion_audit_2026-11-06_01/summary.json` (SHA-256 `599006d4f9d2d0f114df413defe3f9040a0c8bfa96fadf08ccf72509534c9aa9`).

## 2026-11-06 continuation: consolidated ICL/WL report includes shell-exclusion audit

Extended `build_icl_wl_current_evidence_report.py` to require the checksum-audited 48-frame hybrid-PSF shell-exclusion summary, validate its filter and observational guard states, and include the measured paired diffuse-response and mask-coverage changes in the ICL discussion. Moved output-directory creation after all source and query-attempt validation so a failed audit does not leave a misleading empty report directory. Rebuilt the consolidated report with 20 science-component summaries plus the coverage-query audit; it remains exploratory and makes no HSC request or image-pixel read. Report: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_14/icl_wl_current_evidence_report.md` (SHA-256 `769a30f6e37285dae04f7958e144b348f3be3c4276298061f2a9c812fa32ad1a`); summary SHA-256 `ebd7b8a09cf0d8e5e4b524af683c0430717079349b7ffee411cbd3ed3926b467`; builder SHA-256 `44b4150bc67e51b7052561dda83e74930cc5abb77a4b9b847b79488e23149215`. `science_measurement=false` remains.


### 2026-11-06: approved HSC warp pilot and coadd-variance consistency check

The user approved one HSC PDR3 Wide HSC-I 2′×2′ warp request at target 6966
(tract 9703 / patch 306): image, mask, and variance for nine specified visits.
The archive has exactly those nine FITS files (55,511,040 bytes; 52.92 MiB),
with 715×715 aligned planes and warp-to-coadd WCS residual below 7e-10 pixel.
Visits 77296 and 77304 have no finite image/variance coverage; 77330 is partial
(59.5%). The other six cover 86.9–100%.

Added `analysis/icl_wgl_membership/hsc_debug/analyze_icl_warp_pilot_coadd_variance.py`
to compare delivered coadd VARIANCE against `1/sum(1/warp_VARIANCE)` using only
the downloaded local data. At 200–300 proper kpc, primary coverage is 99.71%
and the coadd/predicted variance ratio has pixel q16/median/q84
1.00372/1.00748/1.01506. A sensitivity using only mask planes common to the
warps and coadd retains 48.83% of the annulus and has q16/median/q84
1.00312/1.00656/1.04093. These pixel quantiles are descriptive, not confidence
intervals. The ratio near unity is a one-position product-consistency check,
not independent noise calibration; source Poisson/effective gain, resampling
covariance, sky-transfer covariance, PSF wings, and full ICL-estimator
uncertainty remain open. Keep `science_measurement=false`.

Results and pinned input hashes are under
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_warp_pilot_coadd_variance_2026-11-06_02/`.
Assessment: `analysis/icl_wgl_membership/hsc_debug/icl_warp_pilot_coadd_variance_assessment_2026-11-06.md`.
The earlier `_01` draft output is superseded by `_02` because development
corrected the strict-mask comparison and footprint overlap denominator. No GPU
or Slurm job was used.

The consolidated report builder now validates and includes the warp pilot summary as
an additional input (22 inputs total), without rereading image pixels or accessing
HSC. Updated report:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_15/icl_wl_current_evidence_report.md`
(SHA-256 `381d20917d57c12d3977713d9269ec370c37082b67fb7a66f29c35812b5e05d9`);
summary SHA-256 `28fdd07bf7b80c575419e251bd5ba92abc81e3e98067a8de1ab0695481be1a46`.
The report remains exploratory and records `science_measurement=false`.
Runner SHA-256 `dc67da5b7e6c88a6202f3d07de67e52bd7765cafe34339e3b0eb705ad48f3695`;
report-builder SHA-256 `799f6b8f81fb0fce98bde8e208565843c053a0c4ea8957752420a0d9d8d72b5e`.

### 2026-11-06: supported-sample warp noise plan prepared; awaiting exact-scope approval

After the one-position 6966 pilot, used the approved 67,071-row frame metadata
and the hash-pinned 94-pair ICL–satellite support table to inventory every
matched target/control center. All 188 positions resolve to metadata. The new
plan covers all 94 pairs (26 `i`, 68 `i2`), with 188 HSC-I 2′×2′ requests and
1,471 expected visit products. Pilot-scaled expected transfer is 9,072,971,094
bytes (8.45 GiB); the
runner stops at 160 MiB per request or 12 GiB total. Only image/mask/variance
warp products are in scope. No HSC request was submitted while preparing this
plan. `science_measurement=false` remains.

Plan: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_supported94_warp_noise_plan_2026-11-06_01/plan.json`
(SHA-256 `bd936a94bcdbff83a607acbfca66616c396db526ada50cc89c9e3e312d37007a`).
Planner `analysis/icl_wgl_membership/hsc_debug/prepare_icl_supported94_warp_noise_plan.py`;
submission runner `analysis/icl_wgl_membership/hsc_debug/download_icl_supported94_warp_noise.py`.
Both are plan-only unless an exact plan hash is approved and passed to `--submit`.
The pilot assessment is `analysis/icl_wgl_membership/hsc_debug/icl_warp_pilot_coadd_variance_assessment_2026-11-06.md`.
Planner SHA-256 `475cb7a3436415f8f920ab2fe04e9c20ee28f2565373ddd2f598c8e4fc2ca059`;
downloader SHA-256 `1fe8720f87ae31c951eadd40d83988402d89971f4f2d37ea1bc08da8c2da8da53`.


### 2026-11-06 continuation: approved 94-pair HSC warp download and variance audit

The user approved the exact 188-position 2′×2′ HSC-I warp plan (plan SHA-256 `bd936a94bcdbff83a607acbfca66616c396db526ada50cc89c9e3e312d37007a`). The sequential download completed 184 archives (9,939,640,320 bytes; 9.26 GiB, within the 12 GiB cap). Four approved positions returned HTTP 404 twice: i2 target 3010, i2 control `control_7281_43773774104897710`, i control `control_418_36411590274189484`, and i control `control_6229_41567062857876122`. Every returned archive includes all metadata-expected visits; 161 additional visit files came back within the same cutouts. Download summary SHA-256 `0039321462f8158d1c006a1df64bb3e96ae988744ed780c845f74a60c79b95b4`.

Added `analysis/icl_wgl_membership/hsc_debug/analyze_icl_supported94_warp_noise.py`. It maps each visit WCS onto its coadd using a locally validated affine fit (maximum allowed residual 0.05 pixel), nearest-neighbor samples variance/masks, and compares delivered coadd VARIANCE with `1/sum(1/visit_VARIANCE)`. 184 centers passed; 88 matched pairs remain (24 i, 64 i2). In the 200–300 kpc annulus, the median of center-level ratio medians is 1.00669; strict common-mask coverage has median 81.9%. Target-minus-control medians are −0.000905 (i) and +0.000559 (i2); sky-block bootstrap 95% intervals include zero for both. This is a variance-product consistency diagnostic, not exact HSC resampling/covariance, source-Poisson calibration, PSF/sky validation, or an ICL association result. `science_measurement=false` remains. Assessment: `analysis/icl_wgl_membership/hsc_debug/icl_supported94_warp_noise_assessment_2026-11-06.md`. Output directory `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_supported94_warp_noise_analysis_affine_03_2026-11-06/`, summary SHA-256 `991274bb24525442b44d457540eaa7bd0b0e479f540b43e66d58cbd16d3fd031`. No GPU/Slurm, commit, or push.


### 2026-11-06 continuation: 184-center warp-variance audit in consolidated evidence

Extended `build_icl_wl_current_evidence_report.py` to hash-check and include the approved 94-pair warp-variance audit. The consolidated report now carries 23 pinned summaries, notes the four repeated-404 gaps and 9.26 GiB acquisition, and reports the 200–300 kpc variance ratio and block-bootstrap target/control comparisons. The report remains exploratory (`science_measurement=false`), read summaries only during this synthesis, and made no HSC request or image-pixel read. Report: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_17/icl_wl_current_evidence_report.md` (SHA-256 `5a5d074965d37cf360e42c7dbd984c38ce758663f0b6d15b8381b498ae1102ba`); summary SHA-256 `c6fe422210f825cc2171a85d8235507726a8d29c8a0fdcce06fdd21c4e3cbf64`; builder SHA-256 `c1ae82d3677ed77f84b3ae093b5e31b0973470f16bfae1d2a739ca9233fe6a75`. Next noise gate is kernel-aware warp-to-coadd resampling/covariance, followed by source-Poisson and PSF/sky-transfer checks. No commit or push.


### 2026-11-06 continuation: strict warp-mask ICL association sensitivity

Added and ran `analysis/icl_wgl_membership/hsc_debug/analyze_icl_supported94_warp_mask_association.py`. It verifies exact reconstruction of the prior target/control 200–300 kpc aperture flux and coverage, then recomputes the ICL aperture under the strict warp-valid mask while keeping the saved K=2.5 source masks and original 40% coverage floor. Of 90 pairs with both warp cutouts, 81 remain after strict-mask coverage (23 i, 58 i2). On those identical retained rows, partial r changes +0.0020 in i and −0.0048 in i2; spatial-block paired 95% intervals include zero in both filters. The 13-pair loss relative to the original 94 makes this a selected-subset sensitivity, not a parent-sample correction. `science_measurement=false`. Assessment `analysis/icl_wgl_membership/hsc_debug/icl_supported94_warp_mask_association_assessment_2026-11-06.md`; output `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_supported94_warp_mask_association_2026-11-06_02/`, summary SHA-256 `f7d62ae88864157636d9d7be64ad09d4907060a2cefd7d9ce9b57faf09ef69ad`.

### 2026-11-06 continuation: updated combined ICL/WL evidence report

Extended the hash-validated report builder to include the strict warp-mask association sensitivity (24 component summaries total). Latest report `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_18/icl_wl_current_evidence_report.md` SHA-256 `8bdbcf46e018820fb3819372b69e71162b085fda3f1fef0cb8731093a14636fc`; summary SHA-256 `e76c2b88905ffd411c56b6b56e989e142ddf82d52b5ab3c732572f9e75d6aa22`; builder SHA-256 `8febd51821440a0dd5b8d4ab605fd640c63f8b7a059de4ce82209b80e8466a29`. The synthesis reads JSON only and does not access HSC or image pixels. No commit/push.


### 2026-11-06: approved random-point variance export not submitted

The user approved the exact HSC PDR3 random-point plan for the paired
200–300 proper-kpc annuli at 1,598 target/control centers (16 SQL batches,
estimated 46,424 rows/about 10 MB; plan SHA-256
`eeef98ef2210337f714e36d8a4afcd352fe312a536c90e8e236ff9f843f53e52`). The
scope is limited to 27 random-point, position, variance, sky, input-count, PSF,
and mask columns; it excludes images, shape rows, photo-z, and object-catalog
rows. The plan checker passed. Submission did not reach HSC: the hidden-
credential execution session exited before a CAS job ID was recorded. The
checkpoint remains pending with zero jobs, `query_submitted=false`, and
`network_used=false`. The approval remains valid for this exact plan; no
re-approval is needed. Results path:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_random_pixel_variance_results_2026-11-06_01/`.
No query or download occurred; `science_measurement=false` remains. No commit
or push.


### 2026-11-06 continuation: 16-pixel ICL covariance sensitivity

Added an Astrafek runner extending the 224 supported target/control pairs from
±8 to ±16-pixel variance-normalized residual covariance. All 448 frames
completed; the 129,472 lag cells overlapping the previous ±8 result are exactly
identical. Contracting the 33×33 kernels with the exact saved 200–300 kpc
annulus masks gives finite-lag stack sigma 0.06311 nJy arcsec⁻² in `i` and
0.06290 in `i2`, versus 0.07931/0.07449 for 17×17. Paired 5° block-bootstrap
5,000-replicate intervals for the sigma change are −21.3% to −19.1% (`i`) and
−21.5% to −9.4% (`i2`). These conditional intervals do not include uncertainty
in the lag estimates.

The 16-pixel square-kernel sum is below one because negative long-lag shells
dominate after affine-plane removal and separate per-lag patch medians; this
assembled median kernel is not guaranteed positive-semidefinite. Treat the
change as estimator sensitivity, not a physical noise correction. The finite-
lag term is only 6.1%/5.9% of the saved spatial 95% interval half-width in
`i`/`i2`. Profiles remain unchanged and `science_measurement=false`.
Assessment `analysis/icl_wgl_membership/hsc_debug/icl_covariance_16px_aperture_assessment_2026-11-06.md`.
Kernel, aperture, and block-bootstrap summaries are under `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_supported_covariance_to_16px_2026-11-06_01/`,
`icl_aperture_covariance_200_300_16px_2026-11-06_01/`, and
`icl_aperture_covariance_16px_block_bootstrap_2026-11-06_01/`. No HSC access,
additional download, GPU, Slurm, commit, or push.


### 2026-11-06 continuation: affine-plane covariance-null calibration

Added `analysis/icl_wgl_membership/hsc_debug/analyze_icl_plane_subtraction_covariance_null.py`
to quantify lag anticorrelation induced by the local residual-plane fit. For an
unmasked, homoscedastic 64×64 white-noise patch, the exact OLS projection
kernel yields square-lag sums 0.8087 at radius 8 and 0.3629 at radius 16. A
100-realization simulation through the actual clipped estimator gives means
0.848 and 0.358, consistent within its broad finite-simulation spread. The
observed HSC target/control medians are 1.56–1.66 at radius 8 and 0.74–0.83 at
radius 16. The ideal null confirms that plane removal creates negative tails,
but its predicted 8-to-16 decrease is smaller than the observed ≈−0.8 shell;
it does not explain the entire HSC result. Masks, heteroscedasticity, sky
structure, source masking, and the per-lag median remain entangled. No ICL
correction or science promotion follows; `science_measurement=false`.
Assessment: `analysis/icl_wgl_membership/hsc_debug/icl_plane_subtraction_covariance_null_assessment_2026-11-06.md`.
Runner SHA-256 `150fedeb2a8e00ec011dea525f0bf3c71c31fe9d66f31dc1bc0463724f9a9aee`;
summary SHA-256 `bf93d7c0b9189782c9a510cffee6964a45ab87b49bb7d5630ecc6822419b2760`.
No HSC request or download.


### 2026-11-06 continuation: same-visit split-stack covariance pilot

Added `analysis/icl_wgl_membership/hsc_debug/analyze_icl_visit_split_covariance_pilot.py`
and ran a target/control pilot for supported i2 cluster 101 using already
downloaded warp cutouts. The same seven visit IDs are alternated into two
inverse-variance half-stacks on the coadd grid. One expected visit has no valid
pixels per center; usable split counts are 3/2 for the target and 4/2 for the
control. Clean common support in 200–300 proper kpc is 67.5%/58.3%. Normalized
split-difference RMS before plane fitting is 1.093/1.166, while the locally
plane-fitted zero-lag covariance is 1.182 in each center. Cumulative symmetric
lag sums through ±8 pixels are 0.948/1.187. This is one matched pair, and
nearest-neighbor warp mapping does not reproduce HSC resampling covariance;
residual sky/source variation and variance calibration remain unresolved. No
correction or association claim follows; `science_measurement=false`. Assessment
`analysis/icl_wgl_membership/hsc_debug/icl_visit_split_covariance_pilot_assessment_2026-11-06.md`;
runner SHA-256 `1f28f66ed9477abbef03bb53188f0cb14cd9354498a94eaf43989802be68c929`;
summary SHA-256 `09e557de604d4161cd0281065a59f32669744ab480b180029d192f79ba192359`.


### 2026-11-06 continuation: supported-sample visit-split covariance

Scaled the approved warp-split diagnostic to 90 complete supported pairs (24
i, 66 i2; four center downloads were 404). The initial `_01` aggregate had a
zero-support lag-cell issue and is superseded. Corrected `_02` requires valid
data in all 289 lag cells for the ±8 sum; its common complete-case set is 21 i
pairs and 64 i2 pairs. Within-filter four-metric max-T intervals do not resolve
target/control differences in split RMS or ±8 covariance. The i2 zero-lag
difference is +0.0060 with simultaneous 95% interval +0.0010 to +0.0110; clean
annulus coverage differs by +0.074 (+0.006 to +0.141), so a small support/noise
matching difference remains. Split-RMS medians are 1.082/1.082 in i and
1.085/1.079 in i2 (target/control); these include sky/source residuals and are
not a detector-noise calibration. Target/control visit-set Jaccard medians are
0.80/0.845. Nearest-neighbor resampling remains approximate. No correction or
science promotion; `science_measurement=false`. Assessment
`analysis/icl_wgl_membership/hsc_debug/icl_visit_split_covariance_sample_assessment_2026-11-06.md`;
sample summary SHA-256
`82f7ebaa199dca4d278152969fe2a37bf5e8b4a3380bc999a74e6a276ef4f35c`;
joint-bootstrap summary SHA-256
`973d2def87d97e8c5ef5a1527ddd2e66a76fb8bdeb0dc52ac3c6e7fca0a68f26`.


### 2026-11-06 continuation: white-noise propagation through split estimator

Added `analysis/icl_wgl_membership/hsc_debug/simulate_icl_visit_split_covariance_null.py`
and propagated 100 heteroscedastic white-noise realizations per center through
the supported i2 cluster-101 target/control warp masks, exact visit split,
nearest-neighbor map, weighted stacks, plane fit, clipping and lag estimator.
The null gives raw normalized RMS ≈1.000, fitted zero-lag covariance ≈0.997,
and ±8 cumulative kernel ≈0.715 at both centers. Observed split values are
RMS 1.093/1.166, zero lag 1.182/1.181, and ±8 sum 0.948/1.187 (target/control),
all above the corresponding Monte Carlo 16th–84th ranges. Plane removal
accounts for a negative long-lag tail, but independent white noise and the
supplied variance planes do not explain the observed residual covariance.
Correlated warp noise, residual sky/source structure, and variance calibration
remain confounded. No correction or science promotion; `science_measurement=false`.
Assessment `analysis/icl_wgl_membership/hsc_debug/icl_visit_split_covariance_null_assessment_2026-11-06.md`;
summary SHA-256 `d9dc7090b1cee14fdb88a020adc5164a6da610bcaa76d9d190142b8a2ad86890`.


### 2026-11-06 continuation: annular coadd covariance comparison

Compared the 94 supported same-filter target/control coadds in the 200–300
proper-kpc ICL annulus and the 350–500 kpc outer ring. The estimator uses
64×64 fully annular patches, existing target/control source masks, affine-plane
removal, clipping, variance normalization, and ±8-pixel lag kernels. Five-degree
paired sky-block bootstrap intervals for zero lag, cumulative lag sum, coverage,
and radial difference-in-differences all include zero in both `i` and `i2`.
Strict inner-annulus support is only 8 `i` and 16 `i2` pairs; this is
underpowered and cannot establish equivalence or justify transferring the outer
kernel inward. The outer ring may contain cluster-envelope light, not blank sky.
No correction; `science_measurement=false`. Assessment
`analysis/icl_wgl_membership/hsc_debug/icl_annular_coadd_covariance_assessment_2026-11-06.md`;
runner SHA-256 `16ad75b7f9ac3f2e0cd2e1759a84037f4c9519307e982a03e32174526b80f7d3`;
summary SHA-256 `c03b1c60fcc3415aee615c8685700b42633fa769df0ac4a8173013bafe2f93b0`;
results `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_annular_coadd_covariance_2026-11-06_01/`. No HSC request or image download occurred.

The user reiterated download approval, but the only pinned outstanding HSC
export is the already-approved random-point variance plan
`eeef98ef2210337f714e36d8a4afcd352fe312a536c90e8e236ff9f843f53e52`. Its
prior credential-gated session exited before a CAS job ID; checkpoint remains
pending with zero jobs, `query_submitted=false`, and `network_used=false`. No
new or duplicate request was submitted in this continuation.


### 2026-11-06 continuation: familywise max-T audit of ICL–satellite–WL overlap

Added and ran `analysis/icl_wgl_membership/hsc_debug/adjust_icl_satellite_wl_triple_overlap_max_t.py` on the saved 27-lens i2 common-support table. A restricted Webb six-point wild-cluster bootstrap with shared multipliers across 20 connected 5-degree sky blocks tested all 29 conditional associations together (10,000 valid replicates). The original WL-cross/satellite pointwise pattern at 0.4–0.7 Mpc has unadjusted wild-cluster p=0.043 but max-T p=0.672; the tangential 2–3 Mpc pattern has p=0.070 and adjusted p=0.923. No association survives familywise adjustment. All 29 partial-r point estimates match the source output exactly. This is a low-power selected-sample diagnostic, not evidence of absence or presence of a physical relation; `science_measurement=false`. Assessment `analysis/icl_wgl_membership/hsc_debug/icl_satellite_wl_triple_overlap_max_t_assessment_2026-11-06.md`; runner SHA-256 `162c90585875d0811f6921407c77ec72837f6ecfb5fd9e571bf2c5bbde5e156d`; summary SHA-256 `d5e2d3aa05be32b07be5fa62b8de467aea5314411daba49cb733444850126b0e`; output `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_wl_triple_overlap_maxt_2026-11-06_01/`. The analysis used only saved data, with no HSC request or download.


### 2026-11-06 continuation: consolidated evidence report and pending approved query

Extended `analysis/icl_wgl_membership/hsc_debug/build_icl_wl_current_evidence_report.py` to hash-check and include the direct annular coadd covariance comparison and 29-way ICL–satellite–WL max-T analysis. Latest report is `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_19/icl_wl_current_evidence_report.md` (SHA-256 `5a4d2f7517e3fa78625c19ec96e51dbe8920c0b3c1bcf6ef3293edb6a7bf5426`); summary SHA-256 `ff2d720b09d6f0112670c956dad9053d440b02508f6ce2f4484d6a28b12d941f`; builder SHA-256 `af028a910fe9e5308d4b382029544f9968818e9dd77c97ca08ba1cbd3298418a`. It hash-pins 28 inputs and reads no image pixels or HSC endpoints.

Rechecked the exact approved random-point variance export plan `eeef98ef2210337f714e36d8a4afcd352fe312a536c90e8e236ff9f843f53e52` and its checkpoint on LagEunha. Checkpoint remains `pending`, zero jobs, `query_submitted=false`, and `network_used=false`; no submitter process or `HSC_USER`/`HSC_PASSWORD` environment variables were present. No duplicate request or credential workaround was attempted.


### 2026-11-06 continuation: i2 ICL–satellite–WL support selection audit

Added and ran `analysis/icl_wgl_membership/hsc_debug/audit_icl_wl_overlap_selection.py` on the pinned 68-lens i2 ICL–satellite sample. WL eligibility falls through 60 matched-parent lenses, 35 with target shape rows, 31 with positive random counts, and 27 verified common-support lenses. These occupy 33 connected 5-degree blocks, but only nine blocks contain both selection groups; those mixed blocks have 10 selected and 16 excluded systems. The remaining 17 selected and 25 excluded systems lie in single-status blocks. Spatial positivity fails, so do not use IPW to transport the 27-lens relation to all 68. Within-block max-T balance finds no adjusted redshift/richness/ICL/satellite difference, but is low power with only nine mixed blocks. The target-pair count difference (median 4,342 vs 0) is expected from the shape-support eligibility rule. `science_measurement=false`. Assessment `analysis/icl_wgl_membership/hsc_debug/icl_wl_overlap_selection_i2_assessment_2026-11-06.md`; runner SHA-256 `b2d474d5ecc234e3b3fcb9bea21739f72af720d895d730017c1dfba111895a94`; summary SHA-256 `fef9cdbfdd595f770d1a1a933d29618443524896a318298266a6a1134cb92a98`; outputs `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_overlap_selection_i2_2026-11-06_03/`.

Extended `build_icl_wl_current_evidence_report.py` to verify and include this selection audit alongside the annular covariance and familywise association audit. Latest report `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_20/icl_wl_current_evidence_report.md`, SHA-256 `dae6062a10ca70f0b2eeb6f2fdfbcae76e52665a1bf2ffac31eea9cfa9a7826a`; summary SHA-256 `c665590ba37c5ecbbd9ba08d291434d5f98fdc7d5d6dfd1d8c3fef5c4cea6c2a`; builder SHA-256 `507c0d07e6d545c5afe8f34f6f59d3ab5260129e94aadb2be7245017dd5bba80`. It pins 31 inputs and reads no image pixels or HSC endpoints.


### 2026-11-06 continuation: i2 overlap selection and consolidated report

Created and ran `analysis/icl_wgl_membership/hsc_debug/audit_icl_wl_overlap_selection.py` against the pinned ICL–satellite support and WL eligibility tables. The i2 sample falls from 68 ICL–satellite-supported systems to 27 with verified random-subtracted WL support (60 in matched parent, 35 with target shape rows, 31 with positive random counts, 27 common). The 68 systems occupy 33 connected 5-degree blocks, but only nine blocks contain both groups: 10 selected and 16 excluded lenses. The other 17 selected and 25 excluded lenses are in single-status blocks, so spatial positivity fails and IPW cannot transport the 27-lens result to the parent sample. A four-variable within-block max-T found no adjusted balance difference in redshift, richness, ICL, or satellite contrast; only nine mixed blocks make this low-power. Target shape-pair counts (median 4,342 vs 0) are a direct eligibility mechanism. `science_measurement=false`. Assessment `analysis/icl_wgl_membership/hsc_debug/icl_wl_overlap_selection_i2_assessment_2026-11-06.md`; runner SHA-256 `b2d474d5ecc234e3b3fcb9bea21739f72af720d895d730017c1dfba111895a94`; final result `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_overlap_selection_i2_2026-11-06_03/`; summary SHA-256 `fef9cdbfdd595f770d1a1a933d29618443524896a318298266a6a1134cb92a98`.

Extended the consolidated report builder to hash-validate and include this selection audit alongside annular covariance and the 29-way max-T result. Latest report `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_20/icl_wl_current_evidence_report.md`; SHA-256 `dae6062a10ca70f0b2eeb6f2fdfbcae76e52665a1bf2ffac31eea9cfa9a7826a`; summary SHA-256 `c665590ba37c5ecbbd9ba08d291434d5f98fdc7d5d6dfd1d8c3fef5c4cea6c2a`; builder SHA-256 `507c0d07e6d545c5afe8f34f6f59d3ab5260129e94aadb2be7245017dd5bba80`. It pins 31 inputs and reads no image pixels or HSC endpoints.


### 2026-11-06 continuation: approved CCD-frame metadata status and Lanczos-5 sensitivity

The user approved the exact HSC CCD-frame metadata export plan SHA-256 `49aa4b7601dd42243209d7b3b394a4dd3b4cbf00d87d28cbb9467553c29f4ce3`. It was already complete before this approval reply: all 16 CAS batches returned 67,071 rows (12,196,059 CSV bytes) for 1,598 centers; the latest offline audit confirms all eight calibration columns have zero nulls, exact row totals, and no images, shape rows, photo-z, object catalog, credentials, or pixel coadd weights. No duplicate CAS submission or download was sent. This export is metadata only and does not set a pixel-noise correction; `science_measurement=false`.

Ran `analysis/icl_wgl_membership/hsc_debug/analyze_icl_supported94_warp_noise_lanczos5.py` on the existing approved visit-warp downloads. It maps diagonal warp VARIANCE planes to coadd pixels with a secondary separable Lanczos-5 interpolation: 184/188 centers analyzed, four repeat HTTP 404 gaps. The 200–300 proper-kpc median coadd/prediction ratios are 1.00707/1.00760 in i target/control and 1.00639/1.00611 in i2. Paired target-control block intervals include zero in both filters. Relative to nearest-neighbor mapping, median center-level ratio changes are −0.000245 (i) and −0.000161 (i2); strict shared-mask coverage medians range from 69.0% to 71.9%. This does not model upstream CCD-to-warp covariance or reproduce exact coadd weights, clipping, and mask interpolation. No noise correction or science promotion; `science_measurement=false`.

Assessment `analysis/icl_wgl_membership/hsc_debug/icl_supported94_warp_noise_lanczos5_assessment_2026-11-06.md`; result directory `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_supported94_warp_noise_lanczos5_2026-11-06_02/`; summary SHA-256 `16c8afbcc8287444860bacdcfd01038105bd51d5e82d61e44d8b82d8853e9f74`; table SHA-256 `a0322cc88e63768d172ccc8da52ce5148bef1a4407232b53e93e00bc9ae0f06f`; runner SHA-256 `0df14f6c6be1d858a38437c3ae7300adce81fba859047c8d29d8fe4d07232db9`.

The consolidated report builder now hash-validates the Lanczos-5 summary and both per-center tables. Report `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_22/icl_wl_current_evidence_report.md`, SHA-256 `2ffcb4383e9923f4a904f4b6b684a5c9c03d3e3fc52d5bafe23c8e1f2b0a6389`; summary SHA-256 `51ac24ee822eedd5b8af107ef67ab9b3b893686d422a5604b59411c6482b40cf`; builder SHA-256 `cdd70be6227f272adefc1343c1c19c5d0911fb3f26718763fd3c7c8b0c57eaef`. The report reads no image pixels and makes no HSC requests. No tests, commit, or push were performed.


### 2026-11-06 continuation: full-783 annular covariance scale sensitivity

Expanded the existing delivered-coadd lag-covariance diagnostic to all 783
same-filter target/control pairs. The 64×64-pixel run supports 17/302 `i` and
48/481 `i2` complete pairs; none of three within-filter familywise intervals
excludes zero, and support is strongly low-redshift selected. A 32×32 sensitivity
with up to 64 patches, a minimum of 8 per center, 80% valid pixels, and ±8-pixel
lags raises support to 215/302 `i` and 402/481 `i2`, but redshift SMD remains
−0.65/−0.39. Familywise intervals show a positive `i` lag-sum difference and a
small negative `i2` zero-lag difference; coverage differs in both filters.
Because these results depend on patch scale and valid-mask support, and the lag
sum may contain astrophysical residuals, they are delivered-coadd diagnostics,
not pure detector-noise calibration. No correction; `science_measurement=false`.

The 32×32 covariance output is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_annular_coadd_covariance_full783_32px_2026-11-06_01/summary.json`
(SHA-256 `6e9c786f10116de6973bdf70dda7c0d335d4cc71df717d5ab207d8d42ab25bd2`).
The corrected support audit is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_annular_covariance_full783_32px_support_audit_2026-11-06_02/summary.json`
(SHA-256 `8ad6bcba35b9ce9e7e87a014b97608b4f095b85cebab1e4639b76e0e0357343a`).
Assessment: `analysis/icl_wgl_membership/hsc_debug/icl_full783_annular_covariance_32px_assessment_2026-11-06.md`.
The main summary's numeric patch-size field is 32; its descriptive label had a
stale 64×64 literal. The runner and support auditor now render the description
from the numeric configuration, and the corrected support audit records the
actual 32×32 method. No raw inputs or covariance tables were changed.

### 2026-11-06 continuation: approved random-point variance export is running

The user reaffirmed approval for the exact existing random-point plan SHA-256
`eeef98ef2210337f714e36d8a4afcd352fe312a536c90e8e236ff9f843f53e52`. Resumed
its pinned Lageunha submitter with one hidden interactive credential prompt.
The checkpoint now records submitted CAS job IDs and verified batches; credentials
are not saved. At this handover update, 10/16 batches were complete and the
next batch was running. Returned rows remain limited to the 1,598 approved
centers and the 200–300 proper-kpc annulus, with random-point variance, local
sky/input-count/PSF-moment and mask metadata. No coadd image pixels, shape rows,
photo-z, or object-catalog rows are requested. The post-download integrity and
paired-coverage audit must run after all 16 batches finish; no science correction
or promotion follows from this export.


### 2026-11-06 continuation: approved random-point variance export complete

The user approved and the Lageunha submitter completed the exact HSC PDR3
random-point plan SHA-256
`eeef98ef2210337f714e36d8a4afcd352fe312a536c90e8e236ff9f843f53e52`.
All 16 CAS batches returned 40,048 rows (11,138,863 CSV bytes). The query
covered 1,598 target/control centers in the 200–300 proper-kpc annulus and
returned points at 1,587 centers (791 targets, 796 controls). Credentials were
not saved; the export contains no image, shape, photo-z, or object-catalog rows.
The independent auditor verified every CSV hash, row schema, center identity,
and radial boundary. Both ends have a valid random-point `i_pix_variance`
summary in 787/798 matched pairs. Per-center counts are q16/median/q84 11/22/39.
This remains sparse local variance/sky/PSF metadata, not a coadd weight map or
full correlated-noise/source-Poisson model; no correction or science promotion.
`science_measurement=false`.

Assessment `analysis/icl_wgl_membership/hsc_debug/icl_random_pixel_variance_export_assessment_2026-11-06.md`;
run results `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_random_pixel_variance_results_2026-11-06_01/`;
independent audit `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_random_pixel_variance_audit_2026-11-06_01/`;
audit summary SHA-256 `96245d565375635a994c93ec3294fc581d007161c4479232d5091bbd3cf52389`.

The consolidated report was rebuilt after the variance export and full-783
32-pixel sensitivity. Report:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_23/icl_wl_current_evidence_report.md`
(SHA-256 `14ebf8e0965e12f9dfc5aabad4bc359d8e670428de2919c643ca3830f2d1936d`);
summary SHA-256 `5fe134a256ea35e4a7780a796f9b6c38d6cc460432dcddad956c314b318a5cfc`;
builder SHA-256 `adc7c9d3e86823720643644ea23fe332de32e2602062a3dd94aaa4d3ec08b87d`.
It hash-pins 67 inputs, keeps `science_measurement=false`, and the report build
itself made no network request or image-pixel read. Changes are uncommitted and
unpushed.


### 2026-11-06 continuation: random-point conditions in the ICL–satellite relation

Added `analysis/icl_wgl_membership/hsc_debug/analyze_icl_random_point_variance_satellite_relation.py`.
It joins the independently audited PDR3 random-point medians to the fixed
94-pair support table (26 `i`, 68 `i2`) and tests target/control log ratios and
pair-average logs for local pixel variance, sky RMS, input count, and PSF trace.
Each term is added separately to the redshift/richness partial-correlation
model. One shared max-T family covers all 16 comparisons across filters with
5,000 connected 5-degree block resamples. The baseline r=0.079/0.293 (`i`/`i2`);
`i2` adjustments range 0.260–0.295, all simultaneous change intervals include
zero, and minimum adjusted p=0.951. This does not show that noise/selection has
no effect: the 94-pair sample remains selected, `i` has only 26 pairs, and
random-point median uncertainty is not propagated. No stack correction;
`science_measurement=false`.

Assessment `analysis/icl_wgl_membership/hsc_debug/icl_random_point_variance_satellite_relation_assessment_2026-11-06.md`;
runner SHA-256 `738703bc9cadd561595440eaa9e395a783d1b67ffcabfb4c68e0515085a4b826`;
summary `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_random_point_variance_satellite_relation_2026-11-06_01/summary.json`
(SHA-256 `0c75520ad6b04dd500ba66ae10982e30d6fa9758379615a4372e7aa129569f0f`);
joined table SHA-256 `e9d84e1ee6e78d483685893ade99849850019e1fd66713907b0a627d064207cf`.

Rebuilt the consolidated report with hash verification for the completed
random-point export and this 16-test sensitivity. Report
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_24/icl_wl_current_evidence_report.md`
(SHA-256 `bb45e9d3ed817222debbe0d5aff2bad54582eb3b8d71b974234e58e933a0d122`);
summary SHA-256 `87a7bed6f574ea12a1946f16d2b9da7654b258f02c310c1e532612cc30ddb78f`;
builder SHA-256 `f79589e34c903152237ad04a7fb1607aa6001f819f93f3b24a3d14b54be07574`.
The report pins 69 artifacts and remains exploratory with `science_measurement=false`.


### 2026-11-06 continuation: HSC expansion-coverage query field audit

Diagnosed a likely cause of the repeated HTTP 406 responses in the S20A/S23B
coverage-metadata SQL: the query referenced `pdr3_wide.mosaic.band`, while the
official HSC PDR3 mosaic example uses `filter01`. The CAS response body was
suppressed, so this remains a likely cause rather than a confirmed server
diagnosis. Updated both plan generation and submitter validation to select
`m.filter01 AS band` and filter on `m.filter01 = 'HSC-I'`, preserving candidate
IDs, probe coordinates, geometry, output fields, and data scope.

The corrected S20A plan has 2,596 candidates, 23,364 probes, 52 SQL batches,
and SHA-256 `96ad3b9078beff082d3dce4de1f9bb7608791915ffc7a46a656d9e6cf30e1149`.
The corrected S23B plan has 1,333 candidates, 11,997 probes, 27 batches, and
an upper estimate of 47,988 metadata rows; plan SHA-256
`ee6d2b11893dab0e9ceb0c0d2a3c80c2dee27b600683839d2be294441578ef05`. Offline
scope/hash checks pass for both; no HSC request was sent. The S23B plan remains
approval-gated at that exact hash. These metadata-only requests exclude images,
object/source rows, shape rows, and photo-z. Any later image acquisition needs
its own reviewed scope and approval. Assessment:
`analysis/icl_wgl_membership/hsc_debug/icl_expansion_coverage_filter_column_assessment_2026-11-06.md`.
`science_measurement=false` remains.

## 2026-11-06 continuation: forced-profile coverage selection and ICL association

Fixed `analyze_icl_satellite_forced_coverage_both_sides.py` so its selection
ledger includes all 224 ICL-supported cluster/filter pairs, not only those with
a forced-profile association row. The 20/20 coverage selection retains 50/71
`i` and 110/153 `i2` pairs. Selected pairs have higher median ICL 200–300 kpc
residual than excluded pairs (3.060 vs 1.116 nJy arcsec⁻² in `i`; 4.446 vs
1.853 in `i2`), with similar median redshift and richness. This signals that
coverage support is related to the ICL outcome and the selected subset cannot
be treated as representative.

In paired connected-5° block resampling, all intervals for the full-vs-detected
partial-correlation changes include zero. For `i2`, adding the newly supported
pairs shifts the point estimate by −0.164 (95% interval −0.339 to +0.012),
while adding forced measurements within the detected common sample shifts it
by −0.003 (−0.161 to +0.166). For `i`, those shifts are −0.007 (−0.396 to
+0.432) and −0.175 (−0.451 to +0.113), respectively. This remains a post-hoc
coverage sensitivity with unconfirmed CAMIRA positions, not a full-sample
correction or a science detection. Keep `science_measurement=false`.

Assessment: `analysis/icl_wgl_membership/hsc_debug/icl_satellite_forced_coverage_both_sides_assessment_2026-11-06.md`.
Runner: `analysis/icl_wgl_membership/hsc_debug/analyze_icl_satellite_forced_coverage_both_sides.py`.
Machine output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_forced_coverage_both_sides_2026-11-06_03/`
(summary SHA-256 `b8b79a65e1c838e757fc4ee6c1e1fdf61340d72e989bf035ab29bc0a7a6a314a`; runner SHA-256 `7e9ec1672e0ac000e5e28cf9810c9f98b1e538e60465664853b1dc9a0fe5609f`).
The user reaffirmed approval for frame-metadata plan `49aa4b7601dd42243209d7b3b394a4dd3b4cbf00d87d28cbb9467553c29f4ce3`; it was already complete at all 16 batches / 67,071 rows, so no duplicate request was sent.

### Latest consolidated ICL/WL report: forced-coverage selection included

Rebuilt `build_icl_wl_current_evidence_report.py` to validate and include the
20/20 target/control forced-profile sensitivity. The synthesis verifies the
224-row ICL-supported comparison ledger, the 50/71 `i` and 110/153 `i2`
support counts, all six paired intervals, and hashes for the summary, paired
table, and six underlying inputs. The selected 20/20 pairs have higher median
ICL residual than excluded pairs; the relation remains exploratory and
`science_measurement=false`. The report builder reads saved summaries/tables
only and made no HSC request or pixel read.

Report: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_27/icl_wl_current_evidence_report.md`
(SHA-256 `af88c214eadad945c5177e701f4d7496f95222b297ad84e0be2ba91f595f8b0f`);
summary SHA-256 `42780af71c7392f90f80d53057a8efefe5f9062a84e8f08c595dbbc045eb89b5`;
builder SHA-256 `6a253f424e4e3e3688bbc4698f73c2744a41a6c76413f092944b20d66e8f182a`.
It pins 77 artifacts.

### 2026-11-06 continuation: spatial support overlap in the 20/20 ICL selection

Added `analysis/icl_wgl_membership/hsc_debug/audit_icl_forced_coverage_selection_blocks.py`
to join the 224-pair 20/20 sensitivity to target/control coordinates and
measure selected/excluded overlap in connected paired 5-degree sky blocks.
Mixed blocks are 9/18 in `i` and 21/40 in `i2`; the equal-block selected-minus-
excluded ICL residual mean is +1.80 nJy arcsec⁻² in `i` (95% block-bootstrap
interval −3.26 to +6.66) and +3.73 in `i2` (+0.55 to +6.93). The `i2`
difference remains positive within the observed mixed-block support. The `i`
result is unresolved. This is descriptive, unadjusted for redshift/richness,
post-hoc at the 20% threshold, and does not validate transport to blocks with
only selected or only excluded pairs. No correction or science promotion;
`science_measurement=false`.

Assessment: `analysis/icl_wgl_membership/hsc_debug/icl_forced_coverage_spatial_selection_assessment_2026-11-06.md`;
runner SHA-256 `35888cca6115619e77b23a465ce769aa61e6af495bdd4c6bacd45382a73d3ae8`;
summary SHA-256 `3eb63f02507fa83cf1eea377ebc3574bd2b66e86f204bb76face995dcf5f5742`.
The current 80-artifact integrated report is `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_28/icl_wl_current_evidence_report.md`
(report SHA-256 `33e41f660a70d3a81114b11801b47793fbdad161b8087da87c3a2458d1417fc3`;
summary SHA-256 `bc97fe8ba720c1edb90538d09cbdb918cfb9ac976f5b8cf8a4001a75028351f8`;
builder SHA-256 `e3ade262b4f2b5cffd07deac2fc63645f186565af329d536df3e1579a5f05999`).
All 80 input hashes and the report hash were independently checked. No HSC query
or image read occurred.

### 2026-11-06 continuation: covariate-adjusted support selection

Added `analysis/icl_wgl_membership/hsc_debug/estimate_icl_forced_coverage_selection_adjusted.py`
to estimate the selected-minus-excluded ICL residual only within mixed 5-degree
blocks, with block fixed effects, redshift, and log richness. The pair-weighted
contrasts are +1.47 nJy arcsec⁻² in `i` (9 blocks; 95% block-bootstrap interval
−2.10 to +4.79) and +2.13 in `i2` (21 blocks; −0.64 to +5.91). Both intervals
include zero; the unadjusted i2 contrast is attenuated after covariate control.
The small mixed-block counts and single-status blocks preclude population
transport. This is exploratory support diagnosis only; no correction or
science promotion, `science_measurement=false`. Assessment:
`analysis/icl_wgl_membership/hsc_debug/icl_forced_coverage_adjusted_selection_assessment_2026-11-06.md`;
runner SHA-256 `5f81094c8a8f4a0b37567a6c05cc9e2b482fcd786588a308ae734b3bb2530cf0`;
summary SHA-256 `7e3c57cbd95b32986a7f52b707429f682a6d4e358e12724f9866f14f69b6f07c`.

Integrated report version 29: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_29/icl_wl_current_evidence_report.md`
(report SHA-256 `8350d58462b4eabd228920d970de9852d5b66272e58b629a52577a280dff8d10`;
summary SHA-256 `a3bc38d6191c84c99d58d9391285903670e7b5723674ddfc9627fa74287031a1`;
builder SHA-256 `1968cf3bad3f932d2238601bae668b6fb3833386de76410ecf17848a5720c36e`).
All 82 input hashes and report hash verified. Report build made no HSC request or
image-pixel read.

### Forced unmatched-position radial sensitivity (2026-11-06)

Added `analysis/icl_wgl_membership/hsc_debug/analyze_icl_forced_position_radial_relations.py`
and the assessment
`analysis/icl_wgl_membership/hsc_debug/icl_forced_position_radial_relations_assessment_2026-11-06.md`.
It finds no familywise-resolved relation among the nine supported radial tests;
the smallest max-T p is 0.059 (i, 50–75 kpc). Forced positions have no unique
SEP identity, so the measurement is unmatched-position light rather than
confirmed satellite light. The result is selected and omits measurement,
noise, PSF, and sky-transfer uncertainty; `science_measurement=false`.

An attempt to retrieve the separately approved S23B coverage metadata was
rejected by HSC CAS with HTTP 401 before the first job was created. The exact
approved plan hash was `ee6d2b11893dab0e9ceb0c0d2a3c80c2dee27b600683839d2be294441578ef05`.
No metadata rows were returned, and credentials were not persisted.

Integrated report version 30 is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_30/icl_wl_current_evidence_report.md`
(SHA-256 `62a41b863e1ee2d2150c8647829626972fb3a99d5b123cd993289cc130a65515`;
91 validated inputs; `science_measurement=false`). It includes the forced
unmatched-position radial sensitivity and pins the S23B HTTP 401 attempt audit.
Its summary and all input hashes were independently checked.

### S23B HSC-I metadata coverage screen

The exact previously approved scope completed: 1,333 S23B candidates,
11,997 center/eight-boundary probes, 14 metadata columns, 27 batches, 13,575
result rows. The CAS route and SQL geometry were corrected to the established
`hsc-release.mtk.nao.ac.jp` endpoint and `mosaic.areacube` operator. Independent
Astrafek audit: 1,048 centers matched; 1,043 candidates match all nine points
and carry image/background-product metadata flags at each point. This is not
image-level WCS/mask validation or matched-control validation; no images or
source/shape/photo-z rows were downloaded. Assessment:
`analysis/icl_wgl_membership/hsc_debug/icl_s23b_coverage_metadata_assessment_2026-11-06.md`.

Integrated evidence report version 31 is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_31/icl_wl_current_evidence_report.md`
(SHA-256 `3d6aebd01cb2727ca21d2617ef2de1fd17a85c526199315510710f167a1acf15`,
124 verified inputs). `science_measurement=false`.

### 2026-11-06 S23B matched-control candidate gate

The HSC-I footprint metadata audit identified 1,043 S23B CAMIRA candidates with
metadata matches at the center and all eight 1-Mpc probes. To test whether these
candidates can be paired with usable controls, Astrafek prepared an offline
PDR3 random-position query plan: up to 200 points per candidate within 30′,
29 position/ID and scalar I-band QA/depth columns, maximum 208,600 rows, 105
SQL batches, estimated 80 MB. It requests no images, detected-source rows,
shapes, or photo-z. The plan generator is
`analysis/icl_wgl_membership/hsc_debug/prepare_s23b_icl_control_candidate_plan.py`;
plan directory is `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_control_candidate_plan_2026-11-06_01/`;
plan SHA-256 is `143a1bd0a3a3c341576ce420f77ba0089fbfe80fdf58fbfd7afe30ea8e3f6714`.
All 105 batch SQL hashes were checked locally. The corresponding submitter is
`analysis/icl_wgl_membership/hsc_debug/submit_s23b_icl_control_candidates.py`;
its `--check` mode passed without network, while `--submit` requires the exact
approved SHA-256 and validates returned rows. No CAS query was submitted; this
new metadata scope awaits explicit user approval. Any later HSC image retrieval
requires its own specific approval. The data will support candidate matching,
not a science measurement by themselves; keep `science_measurement=false`.

### 2026-11-06 radial profile contrast

Astrafek added
`analysis/icl_wgl_membership/hsc_debug/analyze_icl_radial_outer_baseline_contrast.py`
and ran it against the saved full-sample profile archives. The six-annulus
50–500 proper-kpc complete-case sample contains 20 `i` pairs/10 connected
5-degree blocks and 48 `i2` pairs/30 blocks. The 200–300 kpc target-minus-control
surface brightness exceeds the area-weighted 300–500 kpc profile by 2.780
nJy arcsec⁻² in raw `i` (95% block interval 0.825–5.115) and 2.494 after
outer-plane correction (0.640–4.774). The `i2` contrasts are 2.359
(−0.110–4.219) and 2.159 (−0.296–4.047). Intervals come from 20,000 connected
paired-sky block resamples. The outer profile can include real cluster light.
The result is a relative radial contrast, not a pure ICL measurement. Sample
selection, source-Poisson, exact coadd covariance, upstream sky transfer, and
extended PSF wings remain unresolved. Assessment:
`analysis/icl_wgl_membership/hsc_debug/icl_radial_outer_baseline_contrast_assessment_2026-11-06.md`.

Integrated report version 32 includes 933 hash-verified inputs at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_32/icl_wl_current_evidence_report.md`
(SHA-256 `6a76bec4c8549d7cd887a5d2b115a3559d1ec8e03e10b719391813ebde83d06a`).
All input hashes were independently checked. No HSC request or image-pixel
read occurred. `science_measurement=false`.

### 2026-11-06 ICL–satellite outer-profile relation

`analyze_icl_satellite_outer_profile_relation.py` joins the existing ICL radial
contrast with the detected-member satellite-centered 50–75 proper-kpc
target-minus-control profile, adjusting for redshift and log richness. The
partial correlations are −0.140/−0.126 in `i` (20 pairs/10 connected 5°
blocks) and +0.417/+0.391 in `i2` (47/30). All four 95% block-bootstrap
intervals include zero. One i2 pair without supported satellite profile stays
missing. This remains exploratory: the satellite sample is detection- and
coverage-selected, and measurement/source-Poisson/coadd covariance/sky-transfer/
PSF-wing uncertainty is not propagated. Report:
`analysis/icl_wgl_membership/hsc_debug/icl_satellite_outer_profile_relation_assessment_2026-11-06.md`.

The integrated evidence report was rebuilt as version 33 with this relation and
940 hash-verified inputs:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_33/icl_wl_current_evidence_report.md`
(SHA-256 `1ef9acc5bff699fed8e1d75131192cc03e88d48832e8ea485918b70c94832742`).
An independent pass verified all 940 input hashes. `science_measurement=false`.

The approved S23B random-control-candidate metadata plan
(`143a1bd0a3a3c341576ce420f77ba0089fbfe80fdf58fbfd7afe30ea8e3f6714`)
was submitted from LagEunha after the user's approval on 2026-11-06. It is a
maximum of 208,600 random-position rows across 1,043 candidate centers, with
29 position/ID and I-band QA/depth fields in 105 batches (estimated 80 MB).
Images, source detections, shapes, and photo-z are excluded. Batches 000 and
001 returned and passed validation with 2,000 rows each; batch 002 was running
at the latest checkpoint. Resume
only the existing checkpoint in
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_control_candidate_results_2026-11-06_01/`;
do not resubmit completed batches. The interactive HSC password is not saved.
Image retrieval remains a separate approval scope.

Follow-up check (2026-11-07): live session 3123 remains active. The checkpoint
now records batches 000–061 complete and verified (124,000 rows); batch 062 is
submitted. `credentials_saved`, `images_returned`, `source_detections_returned`,
`shape_rows_returned`, `photoz_rows_returned`, and `science_measurement` remain
false. Continue polling the same session and checkpoint; do not restart it or
request credentials again while it is live. The export is limited to the
approved random-position and scalar I-band QA/depth metadata plan; image
retrieval remains a separate approval scope.

The local-only postprocessor `analysis/icl_wgl_membership/hsc_debug/audit_s23b_icl_random_control_candidates.py`
is now prepared for the completed export. It verifies the approved batch hashes,
screens depth/seeing/SFD extinction and CAMIRA cluster contamination, and emits
only provisional one-per-target control candidates. Full cutout validation and
any subsequent image retrieval remain outstanding.

### S23B target image receipt (2026-11-07)

The user approved image receipt. `prepare_s23b_icl_target_image_manifest.py`
prepared the exact target scope from the 1,043 nine-probe-complete S23B CAMIRA
centers: HSC PDR3 `pdr3_wide`, HSC-I-family `coadd/bg`, IMAGE/MASK/VARIANCE,
1,000 proper-kpc half-width. The target manifest SHA-256 is
`e487582ba40cd3e58da67704242db4d7219e381fb08f4b33fed47ce4ff268b44`; estimated
payload is 58,311,064,944 bytes (54.31 GiB), with 75,804,384,428 bytes of
30%-margin free space required. The existing `download_pilot.py` now accepts
up to 5,000 frames for this expansion while retaining its per-file cap,
streaming validation, exclusive output, and resume checkpoint behavior.

Download session 5633 is active at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_target_image_plan_2026-11-07_01/target_download_01/`.
At the latest checkpoint, 121/1,043 FITS files (6,826,320,000 bytes) passed
FITS/HDU/WCS validation. `credentials_saved=false`, `science_measurement=false`.
The separate S23B random-control metadata session 3123 is also active: 80/105
batches and 160,000 rows verified; batch 080 submitted. No control images are
included yet; they will be selected from the approved random-point metadata.

`analysis/icl_wgl_membership/hsc_debug/audit_s23b_icl_target_images.py` is
prepared for after the transfer. It independently rechecks each FITS checksum,
IMAGE/MASK/VARIANCE planes, HSC-I filter, WCS alignment, and requested footprint
before any stack is attempted.

Live checkpoint follow-up (2026-11-07): approved target download session 5633
remains active and has validated 373/1,043 images (21,290,526,720 bytes).
Random-control metadata session 3123 has now completed all 105 batches; its
final audit and separate control-image plan are recorded below. The completed
query confirms no saved credentials, images, detected-source rows, shapes, or
photo-z, and `science_measurement=false`.

S23B random-control metadata completed in session 3123 on 2026-11-07: all 105
batches and 208,600 rows passed the submitter's hash/schema checks. The local
audit `audit_s23b_icl_random_control_candidates.py` then verified the batches,
catalogue and SFD input hashes, metadata QA/matching, and cluster exclusions.
At the 1 proper-Mpc CAMIRA exclusion, it produced 1,041 provisional
non-overlapping controls for 1,041 of the 1,043 targets; the 0.5/2-Mpc
sensitivity yields 1,041/1,037 controls. There are 164,984 unique random
positions, 148,682 center-QA passes, and 132,886 metadata matches. Audit output
is at `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_s23b_control_candidate_audit_2026-11-07_01/`;
`candidate_audit.csv` SHA-256 is
`ac9cb72ca769b19ace23e2b90de3f6dd5676c9f3813a525e8f8319f6fae6cb14` and
`provisional_controls.csv` SHA-256 is
`6a9d80db375f5a1266ab386a618849bebfe66f8723ddca5e693142526210ee78`.
The 1-Mpc control image manifest was prepared locally at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_control_image_plan_2026-11-07_01/`;
SHA-256 `4e182a3160055f1dfec4171035bceb1f4a8b574d10f0b12660f67a39dd037c35`.
It specifies 1,041 HSC PDR3 `pdr3_wide` HSC-I `coadd/bg` IMAGE/MASK/VARIANCE
cutouts, each 1,000 proper-kpc half-width. Maximum estimated payload is
58,209,418,068 bytes (54.21 GiB); 30%-margin free-space requirement is
75,672,243,489 bytes. The user approved this separate control-image scope.
Session 81219 stopped after 14/1,041 FITS files (876,274,560 bytes) passed
validation. The HSC cutout API returned HTTP 404 at control `CTRL_S23B_456`,
random-table ID `42740933254447653`, RA 8.823223439694, Dec 2.837294350839,
tract 9718/patch 507. No partial file remains. The first non-TTY attempt
stopped at getpass before network access; `--resume` under TTY used the same
checkpoint. Credentials are not saved.

The offline replan script
`analysis/icl_wgl_membership/hsc_debug/replan_s23b_icl_control_manifest_after_404.py`
prepared v2 while preserving the 1,041-control scope and the other 1,040
positions. For `S23B_456`, it selected the next-ranked metadata-matched,
1-Mpc CAMIRA-clean, globally non-overlapping point: random ID
`42740795815504337`, RA 8.842668812326, Dec 2.786488591937, tract 9718/patch
407. This changes only that cluster entry. V2 manifest SHA-256 is
`ae9eb406cfc179d9980031bb67f85f8ee71c4580ba59e97713278823db5c06df`; the
58,209,418,068-byte estimate (54.21 GiB) and 75,672,243,489-byte 30%-margin
space requirement are unchanged. The user approved this exact v2 hash. The v1
14-file failure checkpoint is preserved as
`download_status_v1_404_checkpoint.json`; the resume checkpoint now records the
approved v2 hash. The 14 existing FITS files match unchanged v2 rows and were
revalidated on resume. Control-download session 56526 is active; 86/1,041
files (5,233,340,160 bytes) passed validation at the latest checkpoint,
including the replacement position. Credentials are not saved.

Approved target download session 5633 remains active with 850/1,043 target
images (47,507,889,600 bytes) validated at the latest checkpoint. The partial
sample currently contains 685 coadds with header filter `i2` and 165 with `i`;
this split is provisional until all images arrive.
The latest `/gpfs` capacity check reports 157,141,126,086,656 bytes available,
above the approved target and control-plan 30%-margin requirements.
Prepared `analysis/icl_wgl_membership/hsc_debug/audit_s23b_icl_control_images.py`
for an independent post-download checksum/FITS/HDU/filter/WCS and target-pair
audit. It pins the approved v2 manifest and requires all 1,041 images; it has
not run on the partial transfer. The legacy 2D stacker is specific to the
previous 800-target/798-control sample, so the S23B files will not be passed to
it as though they were that sample. Both image audits and the centering/member
catalog gates remain prerequisites for a scientific ICL–satellite result.
Centering note: the S23B method paper says a bright galaxy is selected as the
BCG and richness/redshift are recomputed around it, while the released seven-
column cluster table labels its fields only as RA and Dec. The available
catalog does not provide a separate BCG coordinate pair. Thus the earlier v1
manifest phrase “not a BCG centroid” is not independently established; read it
as “not yet verified against a separate BCG position.” The approved manifest
and download coordinates remain unchanged. The manifest generator now uses
this more cautious wording for future plans. Sources:
https://arxiv.org/abs/2512.13954 and
https://github.com/oguri/cluster_catalogs/tree/main/hsc_s23b_camira.

Local S23B CAMIRA catalog check: `catalogue_s23b_camira_2026-10-20_01/camira_s23b_wide_sm.tbl`
is 827,339 bytes with 12,162 data rows (SHA-256
`5af057028feab7976cd3ed389a7246cbe73add8bbb96891101f0f692b1de4d71`). Its
published columns are cluster ID/name, RA/Dec, photometric redshift, richness,
and BCG spectroscopic redshift; it contains no member positions or membership
rows. The S23B image expansion therefore cannot measure the satellite-light
relation on its own. Obtain and validate member/source membership data, and
verify the ICL center definition against a BCG center before interpreting a
satellite association. No such data were requested in the approved image scope.
The reproducible local audit is
`analysis/icl_wgl_membership/hsc_debug/audit_s23b_camira_catalog_content.py`
(SHA-256 `8bbf3c4ea8de9df9902528c4cada49ff6dbcfea9d0bd6507a60c2f7b0bd631ad`);
it verified 12,162 well-formed seven-field rows against the pinned file hash.
Summary: `/gpfs/kjhan/Astrafek/data/s23b_camira_catalog_content_audit_2026-11-07_01/summary.json`
(SHA-256 `28f7836a0520f6b59a32241531a8e3ce5fc2e0bb7b71248dfad05f8870743ac8`).
The public `oguri/cluster_catalogs` S23B directory lists the cluster
catalog files and points separately to photometric LRG catalogs; it does not
list a CAMIRA member table. See
https://github.com/oguri/cluster_catalogs/tree/main/hsc_s23b_camira.
The associated paper reports more than 600,000 probable CAMIRA member galaxies,
but describes its photometric LRG catalog as independently constructed rather
than the CAMIRA member list; the public GitHub directory does not identify a
downloadable member table. Track the membership source as unresolved and do not
substitute the LRG catalog as confirmed cluster members. Paper:
https://arxiv.org/abs/2512.13954.

## 2026-11-07 S23B image transfer and pixel-support checkpoint

The approved target transfer completed: all 1,043 S23B CAMIRA HSC PDR3
`pdr3_wide` HSC-I `coadd/bg` IMAGE/MASK/VARIANCE FITS files are present. The
target manifest SHA-256 is
`e487582ba40cd3e58da67704242db4d7219e381fb08f4b33fed47ce4ff268b44`; total
file size is 58,373,470,080 bytes. Independent local checksum/FITS/HDU/filter/
WCS audit verified 1,043/1,043 frames (header filter `i`: 251; `i2`: 792).
Audit summary:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_target_image_audit_2026-11-07_01/summary.json`
(SHA-256 `8580d64ba76e1c1f64331f7ce29b8d6f9dde4b56914b699e9e09e57f339e4b25`).
The audit confirms file integrity, layout, filter evidence, and WCS only; it
does not certify pixel usability. `science_measurement=false` and
`images_ready_for_stacking=false` remain in force.

Automatic source-mask and radial-area preparation is running at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_mask_support_2026-11-07_01/target/`.
Its 90-frame checkpoint has 89 successes and one failure (`S23B_10647`, “Too
few valid pixels”). A separate pixel diagnostic found zero unflagged pixels:
all 4,420,504 pixels carry both `MP_INTRP` and `MP_NO_DATA`; this frame must be
excluded from any stack. The preparation is support-only and does not approve
an ICL measurement. Do not rerun into the existing output directory; inspect
its final summary after the current run finishes.

The separately approved control-image v2 transfer is still running from its
existing resume checkpoint. At this checkpoint it has 308/1,041 files; the
approved manifest SHA-256 is
`ae9eb406cfc179d9980031bb67f85f8ee71c4580ba59e97713278823db5c06df`. Status is
at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_control_image_plan_2026-11-07_01/control_download_01/download_status.json`.
No credentials are saved. Do not start the target-control stack until the
control transfer and its independent audit finish and filter/usable-area
matching is reviewed.

For satellite positions, the one-file member-catalog intake plan remains
separate from the image approval. Plan SHA-256:
`724860324ee069e1b07c1e324ca48c14fde9c6278343b4f645044c002e4eeb56`;
the candidate is `camira_s23b_wide_sm_member.dat.gz` (10,626,590 bytes; Drive
ID `1DEnjGj9P1ilUooEm2U5_A-qXUrxeiJdE`). Download approval is pending. External
sharing visibility and file contents remain unverified; do not download until
the user explicitly approves this exact file. No ICL–satellite science stack
has been produced from the S23B sample.

Follow-up checkpoint on 2026-11-07: target mask/area preparation reached
430/1,043 frames (427 successful, 3 excluded for zero usable pixels):
`S23B_10647`, `S23B_11786`, and `S23B_432`. Pixel checks for all three show
every pixel flagged `MP_INTRP` and `MP_NO_DATA`; each has zero unflagged,
finite, positive-variance pixels. The reproducible audit is
`analysis/icl_wgl_membership/hsc_debug/audit_s23b_icl_zero_support_frames.py`;
its output is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_zero_support_audit_2026-11-07_01/summary.json`
(SHA-256 `ffa87816d6b92da75ded4fb4ef79628e50661d7fa03a60c53b6f43a62da3bbd5`).
The approved control download reached 376/1,041 files with no transfer error;
keep using its existing v2 resume session and checkpoint.

The target mask/area run has now finished with `partial` status: 1,039/1,043
frames processed, and the four failed IDs are exactly the zero-pixel frames
listed above. Post-exclusion header-filter counts are 251 `i` and 788 `i2`.
The descriptive support assessment is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_mask_support_assessment_2026-11-07_01/summary.json`
(SHA-256 `c140d6e02862b6d7654494a9c243eba5e664fb04234d25d5ae3e18e2ecb27c47`).
Filter-separated diagnostic figure:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_mask_support_assessment_2026-11-07_01/radial_unmasked_area_by_filter_v4.png`.
For K=2.5, median unmasked annular area fractions are 0.350 at 50–100 kpc and
0.417–0.423 across 200–500 kpc; 144/1,039 frames have zero support at
50–100 kpc, while 331/1,039 reach at least 50% there. The center pixel is in a
base quality mask in 446 frames and in the K=2.5 automatic source mask in 824.
These are support diagnostics, not adopted cuts; no frame-selection threshold
has been applied. The output remains `science_measurement=false` and
`images_ready_for_stacking=false`. Control-image v2 transfer has reached
533/1,041 files at the latest check, with no transfer error.

Latest control-transfer checkpoint on 2026-11-07: 567/1,041 files, status
`downloading`, no error, and the checkpoint still pins approved v2 manifest
SHA-256 `ae9eb406cfc179d9980031bb67f85f8ee71c4580ba59e97713278823db5c06df`.

Later on 2026-11-07 the approved v2 transfer reached 627/1,041 and stopped at
`CTRL_S23B_8547` with `HSC mask-plane mapping differs from the reviewed probe`.
The complete response remains at
`.../control_download_01/images/cluster_CTRL_S23B_8547.fits.part` (66,349,440
bytes). Its FITS layout, image shape (2335×2367), and HSC-I filter are readable.
The named mask planes are the established map except that `MP_NOT_DEBLENDED`
and `MP_UNMASKEDNAN` are interchanged between bit positions 11 and 12. The
other 627 control headers and all 1,043 target headers use the original map.
HSC's pipeline workshop explicitly says the bit/event correspondence can vary
and to consult the MASK header; see
https://hsc.mtk.nao.ac.jp/pipedoc/DataWorkshop2022/_downloads/2682edf0a494c8b005acf4efa51952ae/hsc_kaiseki_20221128.pdf
(section “マスクとフラッグ”). Accordingly, `download_pilot.py` now validates
all 17 known named planes and a unique 0–16 bit assignment; mask preparation
and control stacking derive bit selections from each image's own header. The approved v2 manifest
and 1,041-control scope are unchanged. Its original resume checkpoint is being
revalidated before the remaining files are fetched; no new catalog rows or
coordinates have been added. The S23B-specific exploratory stack driver is
`analysis/icl_wgl_membership/hsc_debug/stack_s23b_icl_target_control_2d.py`;
it pins both manifests and independent audits, excludes exactly the four
zero-support target frames, reports `i`/`i2` separately, and keeps
`science_measurement=false`. It must wait for the resumed transfer and its
independent audit.

The exact v2 resume then recovered `CTRL_S23B_8547` from the saved response
using its header-declared mask map and continued the approved transfer. Latest
live checkpoint: 734/1,041, status `downloading`, no error, manifest SHA still
`ae9eb406cfc179d9980031bb67f85f8ee71c4580ba59e97713278823db5c06df`; process
is running in the existing LagEunha session. No new HSC files outside the
approved manifest were requested. Two pre-stack artifacts are now prepared:
`analysis/icl_wgl_membership/hsc_debug/prepare_full_cluster_mask_support.py`
accepts the pinned 1,041-control download plus audit, and
`analysis/icl_wgl_membership/hsc_debug/audit_s23b_icl_pair_common_support.py`
will compare warped target/control/common mask area without reading science
image pixels. The flux stack driver now requires that support audit and a
complete control-mask preparation, and it excludes any pair with zero common
mask support inside 1,000 proper kpc. These audits and the stack remain unrun
until the transfer finishes.

## 2026-11-07 S23B approved control download, paired support, and exploratory stack

The exact approved v2 control download completed on LagEunha: 1,041/1,041
HSC PDR3 `pdr3_wide` HSC-I `coadd/bg` IMAGE/MASK/VARIANCE cutouts;
58,271,158,080 bytes; zero failed or partial files. Manifest SHA-256 remains
`ae9eb406cfc179d9980031bb67f85f8ee71c4580ba59e97713278823db5c06df`. The
independent checksum/FITS/HDU/filter/WCS audit verified all 1,041 frames and
found 248 `i` and 793 `i2` headers. Audit summary:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_control_image_audit_2026-11-07_01/summary.json`.

One valid control MASK header interchanges `MP_NOT_DEBLENDED` and
`MP_UNMASKEDNAN` between bits 11 and 12; the remaining 1,040 controls and all
1,043 targets use the other mapping. The complete header mapping is stored per
downloaded frame and all source-mask/quality-mask logic uses each image's own
header. HSC pipeline documentation says to read MASK-header bit/event
correspondence:
https://hsc.mtk.nao.ac.jp/pipedoc/DataWorkshop2022/_downloads/2682edf0a494c8b005acf4efa51952ae/hsc_kaiseki_20221128.pdf.

CPU-only control mask preparation completed for 1,041/1,041 frames with zero
errors at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_mask_support_2026-11-07_01/control/`.
The mask-only paired support audit completed at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_pair_common_support_2026-11-07_02/`.
Its summary SHA-256 is
`376d9c076a54166bf8bc0f30f1cdb26ce65be9c2c43743bf5b60d0367867f91c`.
It reports 1,022 common-filter pairs (`i`: 241; `i2`: 781), 17 filter
mismatches, two controls paired to zero-support targets, two other targets
without matched controls (`S23B_10647`, `S23B_8226`), and no pair with zero
common mask pixels across 0–1,000 proper kpc. This audit reads WCS and prepared
masks only; no science-image pixels or flux were read. Its first output version
omitted per-pair mask paths and could not pass the stacker's provenance join;
version 02 adds those paths, reproduces the same counts, and is the version
used by the stack.

Coverage is a significant limitation. In the innermost 0–25 kpc annulus, only
56/241 `i` and 195/781 `i2` pairs have any common support; median common area
is zero. For 0–50 kpc, median common area is 0.7% (`i`) and 1.3% (`i2`). In
50–100 kpc the medians are 13.1% and 18.5%; across 300–800 kpc they are 11.1%
and 16.5%. The share of pairs with any annular support reaches 90% only at
about 175 kpc in `i` and 150 kpc in `i2`. No coverage threshold was selected
after viewing the data.

The S23B exploratory flux stack completed at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_target_control_stack_2026-11-07_01/`.
It stacks 1,022 pairs on a 401×401 grid spanning ±1,000 proper kpc at 5 kpc
spacing and reports `i` and `i2` separately. The output summary SHA-256 is
`c8bac8dbdd0a886604c216b2410725eb97883b0cdc1a91fd47e7a7cfcaaaa20b`; both
FITS maps pass their stored SHA-256 and FITS-checksum checks. Maps have median
per-pixel pair counts of 30 (`i`) and 134 (`i2`), maxima 59 and 214. The
descriptive equal-pair target-minus-control profile is positive at 50–300 kpc
and trends toward zero outside 300 kpc. For example, 50–100 kpc is
30.7±3.1 nJy arcsec⁻² in `i` (189 pairs) and 37.3±2.0 in `i2` (633 pairs);
100–300 kpc is 4.47±1.14 and 4.77±0.95, respectively. These are between-pair
SEM only, not full uncertainties. Central values have sparse support and are
dominated by compact central light.

Quality assessment:
`analysis/icl_wgl_membership/hsc_debug/s23b_icl_stack_quality_assessment_2026-11-07.md`.
Diagnostic plot:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_stack_quality_assessment_2026-11-07_01/stack_quality_diagnostic.png`.
The stack is useful for pipeline validation and as a candidate 50–300-kpc
diffuse-light signal only. It is not a science result: Galactic extinction,
PSF-wing matching, residual-sky transfer, coadd covariance, spatial jackknife,
and manual BCG/mask review remain. `science_measurement=false` and
`images_ready_for_science=false` remain in the stack and both FITS products.
No commit or push has been made.

## 2026-11-07 S23B exploratory ICL spatial-block uncertainty

Added `analysis/icl_wgl_membership/hsc_debug/analyze_s23b_icl_spatial_blocks.py`
to estimate spatially clustered uncertainty from the saved pair-level profiles.
It verifies the pinned target/control manifests and completed stack, reads no
image pixels, and makes no HSC request. It resamples connected RA-Dec cells at
2.5, 5, and 10 degrees with two half-cell phases; these approximate blocks are
not official HSC field polygons.

The 5-degree baseline gives 50–100 proper-kpc target-minus-control excesses of
30.68 nJy arcsec⁻² (`i`, 189 supported pairs/25 blocks; 95% block interval
25.39–35.29) and 37.26 (`i2`, 633/38; 33.81–40.45). For 100–300 kpc the
corresponding means are 4.47 (`i`, 232/26; 1.31–6.83) and 4.77 (`i2`, 756/39;
3.58–6.17). All six width/phase layouts keep both 50–100 and 100–300 kpc
intervals above zero. At 300–800 kpc intervals include zero; the 800–1,000 kpc
`i2` residual changes zero-crossing with block layout. This is a spatial
resampling sensitivity, not a calibrated survey covariance or ICL detection.

Assessment and plot:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_spatial_blocks_2026-11-07_02/assessment.md`
and `s23b_icl_spatial_block_uncertainty.png`. Summary SHA-256:
`423bb5e3ad52f3b349dedc74b00a2c05ebb4589e20402abb68b4e99d17199afc`.
The stack remains `science_measurement=false`; residual sky, Galactic
extinction, PSF wings, coadd-pixel covariance, and manual center/mask review
remain open. Satellite association still requires the separate CAMIRA member
file; that one-file download approval is pending. No commit or push has been
made.


## 2026-11-07 S23B paired ICL radial contrasts

The exploratory paired inner-minus-outer radial contrast is implemented in
`analysis/icl_wgl_membership/hsc_debug/analyze_s23b_icl_relative_radial_contrasts.py`
(source SHA-256 `1e1f385e32aa67150f7e666ad7e36a80e189be81e5bfdd04b6f1acc47977dbba`).
For each saved matched target/control pair, it subtracts the common-support
300–800 proper-kpc contrast from the 50–300, 100–300, or 200–300 proper-kpc
contrast. Uncertainty is resampled by connected RA-Dec sky blocks at 2.5, 5,
and 10 degrees with two grid phases (5,000 draws per layout). It reads saved
profile/audit tables only; it made no HSC requests and reads no image pixels.

At the primary 5-degree layout, the 50–300 and 100–300 kpc contrasts are
positive for both filters. In HSC-i, they are 6.56 [3.84, 8.67] and 4.53
[1.90, 6.54] nJy arcsec⁻² (95% block intervals; 232 pairs, 26 blocks).
In HSC-i2, they are 7.61 [5.72, 8.88] and 4.69 [3.24, 5.71]
nJy arcsec⁻² (756 pairs, 39 blocks). The 200–300 kpc contrast is weaker:
3.25 [0.498, 5.44] in i and 1.91 [0.401, 3.04] in i2; at least one of the
six grid layouts crosses zero for each filter. The sign of the wider 50–300
and 100–300 contrasts remains positive across all six layouts.

Outputs are in
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_relative_radial_contrasts_2026-11-07_02/`.
Summary SHA-256: `a47c6679ad318abcbee5e5ab46e1a984976bd01c6f1dcbbe3d6de99b89306a8`;
assessment SHA-256: `c5c4e5f3e5e275b9d5a41cb6777adc29d51c88944615c2ae7e61eb3fc752a8fa`.
This is a relative radial excess, not an absolute ICL measurement: the
300–800 kpc reference annulus may contain cluster light. Sky gradients,
Galactic extinction, extended PSF wings, coadd-pixel covariance, and manual
center/mask review remain unresolved. Keep `science_measurement=false`.
Satellite association still awaits the separately approved CAMIRA member-file
download; it has not been downloaded. No commit or push has been made.


## 2026-11-07 S23B approved control-image transfer complete

The approved v2 HSC PDR3 `pdr3_wide` HSC-I `coadd/bg` control transfer is
complete: 1,041/1,041 IMAGE/MASK/VARIANCE FITS, 58,271,158,080 bytes
(54.27 GiB), with 248 `i` and 793 `i2` headers. The download checkpoint and
manifest both pin SHA-256
`ae9eb406cfc179d9980031bb67f85f8ee71c4580ba59e97713278823db5c06df`.
The offline FITS/header/WCS/file-integrity audit is complete at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_control_image_audit_2026-11-07_01/`:
1,041 unique control IDs map one-to-one to 1,041 targets; verified bytes total
58,271,158,080. Summary SHA-256:
`3210513093e363fa7a63c98669ab2cb332f96711295f01936dd121ab64708510`;
verified table SHA-256:
`6e600066dc9d97d1448c8538b5bd3b2e38e1e56797485765296c4fcc44d5574c`.
No credentials were saved and the audit used no network. The target-image
audit is also complete: 1,043/1,043 FITS files. Mask support succeeded for
1,039 targets; four failed the automated mask step for too few valid pixels.
The mask-only paired common-support audit retained 1,022 same-filter pairs and
found zero pairs with zero common support inside 1,000 kpc. The exploratory
2D stack is already complete for those 1,022 pairs (17 filter mismatches and
two zero-support targets excluded; `science_measurement=false`). Summary
SHA-256 is `c8bac8dbdd0a886604c216b2410725eb97883b0cdc1a91fd47e7a7cfcaaaa20b`.
This does not establish sky, extinction, PSF-wing, covariance, or member
completeness corrections; `images_ready_for_science` remains false. CAMIRA
member positions had not yet been imported at the time of this mask-only
stack; the later approved intake below supersedes that status. No commit or
push has been made.


## 2026-11-07 S23B ICL contrast versus CAMIRA richness

Added `analysis/icl_wgl_membership/hsc_debug/analyze_s23b_icl_richness_relation.py`
(source SHA-256 `46eb768c301544802508fce7f15493504c272189a2711e94c71f49e0f2486cf5`).
It models each saved pair's inner-minus-300–800 kpc relative ICL contrast
against log10 CAMIRA richness, adjusted for redshift and target-minus-control
SFD E(B−V) mismatch. It uses six connected RA-Dec block layouts, 5,000
resamples per layout/filter/contrast; no network requests or science-image
pixels. E(B−V) enters only as a nuisance covariate; no extinction correction
was applied.

At the 5° phase-zero primary layout, richness slopes per dex (95% block
intervals; n pairs / blocks) are: `i`, 50–300 kpc 6.14 [−13.66, 22.14]
(232/26), 100–300 4.13 [−15.56, 19.47] (232/26), and 200–300 4.51
[−10.70, 19.96] (231/26) nJy arcsec⁻²; `i2`, 50–300 9.51 [−3.01, 22.32]
(756/39), 100–300 7.82 [−4.77, 18.04] (756/39), and 200–300 2.25
[−10.89, 11.90] (752/39). All primary intervals include zero. The 50–300
and 100–300 `i2` intervals are positive only in the 10° phase-zero layout;
the paired phase-one layouts cross zero. This relation is therefore sensitive
to spatial partitioning and does not establish a richness trend. It also does
not measure satellite positions or satellite radial structure: CAMIRA richness
is only a cluster-scale proxy. Do not interpret a zero-crossing interval as
proof of no relation.

Outputs:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_richness_relation_2026-11-07_01/`.
Summary SHA-256 `144da86644e8fadfff95d18662af4744d1b4eef0c558e5750d7f7f5200e45f16`;
assessment SHA-256 `4190a885f57fa53567edef8e513426caba75505c6ca4b0cdc6d0979a5a52abe2`.
Keep `science_measurement=false`; residual sky, extinction correction, PSF
wings, covariance, centering/mask uncertainty, and member completeness remain
open. At the time of this richness-only analysis, the CAMIRA member file had
not yet been received; later intake and direct-analysis entries below
supersede that status. No commit or push has been made.


## 2026-11-07 S23B richness–coverage sensitivity

Added `analysis/icl_wgl_membership/hsc_debug/analyze_s23b_icl_richness_coverage_sensitivity.py`
(source SHA-256 `4a0f6be1a6eadba647cf135d365414a79e8f68746f4579a545422e795aed81cb`).
It repeats the pinned richness regression with an added log common-mask pixel
count and uses identical connected-block bootstrap draws for the primary and
support-adjusted slopes. This is a sensitivity model; conditioning on coverage
is not promoted as a causal correction.

At the 5° phase-zero layout, the 50–300 kpc slope changes from 6.14 to 5.98
(95% paired block interval for change −5.15 to 6.15) in `i`, and from 9.51 to
9.38 (−1.97 to 1.03) in `i2`, in nJy arcsec⁻² per dex richness. The other
annuli show similarly small changes; all paired slope-change intervals include
zero. The correlation between log richness and log common-pixel support is
small (|r| ≤ 0.017) in the primary layout. This weakens the concern that a
linear richness–support correlation alone explains the slope, but does not
resolve coverage selection or the other science gates. All outputs remain
exploratory, `science_measurement=false`.

Outputs:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_richness_coverage_sensitivity_2026-11-07_01/`.
Summary SHA-256 `7a4df0afee18db5af23c525b9df2184c98b4e08c8b8d49ad45b0f043bffc09ae`;
assessment SHA-256 `be7de62b92fe09296bfc057d283e19423d85ef4f663f5c8a12b5f7b2b62a4b3c`.
Direct satellite analysis still requires the separately approved CAMIRA member
file, whose download remains pending. No commit or push has been made.


## 2026-11-07 S23B ICL extinction sensitivity

Added and ran `analysis/icl_wgl_membership/hsc_debug/analyze_s23b_icl_extinction_sensitivity.py`
(SHA-256 `746228c31f8e80c8b183ef7b50e91b356afec4463001126c12bf547176cc7b2c`).
It applies the nominal hscPipe coefficients `i=1.641`, `i2=1.617` separately
to saved target/control annular means on the existing common support. The
runner pins input hashes, reproduces the shared RA grid origin and connected
block counts, and exactly reproduces the prior raw relative-contrast results
in all 180 checked summary fields. It reads no image pixels and makes no HSC
requests.

For the 50–300 kpc minus 300–800 kpc contrast, the 5° phase-zero result changes
from 6.561 [3.844, 8.667] to 6.930 [4.043, 9.168] nJy arcsec⁻² in `i` (232
pairs, 26 blocks), and from 7.611 [5.722, 8.877] to 8.055 [6.117, 9.385] in
`i2` (756 pairs, 39 blocks). The paired correction is +0.369 [0.190, 0.516]
and +0.444 [0.337, 0.569], respectively; its interval stays positive across
all six spatial layouts. This is a roughly 5.6–5.8% shift, not a science
promotion: rescaling also affects residual sky, and PSF wings, sky transfer,
source-Poisson/coadd covariance, mask/center review, and member completeness
remain unresolved. `science_measurement=false` remains.

Assessment:
`analysis/icl_wgl_membership/hsc_debug/s23b_icl_extinction_sensitivity_assessment_2026-11-07.md`.
Output:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_extinction_sensitivity_2026-11-07_04/`
(summary SHA-256 `5dfb730f26043f2073790d21326156fd7e3d79730f457e7c0bc738745b32e8fe`).
At the time of this extinction-only sensitivity, the separately approved
member-file intake had not yet occurred; later sections record its receipt
and use. No commit or push has been made.
## 2026-11-07 S23B CAMIRA member intake and direct ICL–satellite relation

This is the current status and supersedes earlier entries in this handover that
say the member-file download or direct satellite analysis is still pending.
The user approved and received only camira_s23b_wide_sm_member.dat.gz
(10,626,590 compressed bytes; source SHA-256
66529e02c6ae11ebdccc3fb1d4756b639364f372f9d0f09a721786a19906945a). No
other file was downloaded and no credential was saved.

The importer matched all 12,162 member-file cluster groups uniquely to the
pinned S23B CAMIRA cluster catalog and joined all 1,043 target clusters,
yielding 45,634 member-position rows. WCS checks place 39,166 positions
(85.83%) inside the target cutout pixel-center footprints; 6,468 are outside.
This is a geometric footprint check, not proof of usable or unmasked pixels.
For supported profiles, all catalog candidates in both analyzed annuli are
inside the image footprints. The final import summary SHA-256 is
3d1442420aeae051f5a970c6b717a93c69fdaa57cebb096063c46d8afe5f4d4e.

An offline exploratory association joined those positions to the 1,022
filter-matched target/control profile pairs. It regresses the 50–300 minus
300–800 proper-kpc relative ICL contrast against the unweighted CAMIRA
candidate-member surface-density contrast, adjusting for richness, redshift,
target/control SFD E(B−V) mismatch, and common-pixel support. The 5°
phase-zero connected-block results are:

| Filter | Pairs / blocks | Slope (nJy arcsec⁻² dex⁻¹) | 95% block interval |
| --- | ---: | ---: | ---: |
| HSC-i | 232 / 26 | 5.72 | [−0.06, 13.63] |
| HSC-i2 | 756 / 42 | 4.24 | [0.99, 7.80] |

The point estimates are positive, but the interval excludes zero in only 3/6
spatial layouts for i and 4/6 for i2. At a 10% minimum common-support floor
in both annuli, the i2 interval crosses zero; the i result remains unresolved.
Treat this as a first direct, exploratory association, not a detection. Common
support is sparse, catalog candidates are not spectroscopically confirmed,
and the 300–800 kpc reference can contain cluster light.

Residual sky, extinction correction, extended PSF wings, full coadd
covariance, center/mask sensitivity, and member completeness/purity remain
open. Keep science_measurement=false and images_ready_for_science=false.
The analysis read saved tables and member positions only; it made no network
requests and read no science-image pixels.

Code and assessment:
analysis/icl_wgl_membership/hsc_debug/import_s23b_camira_member_catalog.py
analysis/icl_wgl_membership/hsc_debug/analyze_s23b_icl_satellite_relation.py
analysis/icl_wgl_membership/hsc_debug/s23b_icl_satellite_relation_assessment_2026-11-07.md
The relation summary SHA-256 is
fab29d7aef464707c31e2616fe1db9fc3bad9d281ad6d5dfa6d63f2b484bbb39.
Products are under
/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_satellite_relation_2026-11-07_03/.
No commit or push has been made.

## 2026-11-07 S23B raw member-field-8 weighting sensitivity

A sensitivity runner tested whether the direct ICL–CAMIRA candidate relation
changes when each member position is weighted by the member file's eighth raw
field. Its meaning is unverified: values lie in [0,1] (median 0.373253), and
per-target field sums divided by catalog richness have q10/median/q90
0.830/0.945/1.071, but this does not prove the field is the published CAMIRA
membership weight. Keep it labeled raw field 8.

The adjusted 50–300 minus 300–800 proper-kpc relative-ICL slope per dex of
weighted candidate-density ratio, at 5-degree phase zero, is 5.13
[−0.80,13.05] nJy arcsec⁻² for HSC-i (232 pairs, 26 blocks) and 5.01
[2.44,7.62] for HSC-i2 (756 pairs, 42 blocks). The i2 interval excludes
zero in all six spatial layouts; i does so in two of six. At a 10% minimum
common-pixel floor in both annuli, both filters' intervals include zero.
This is an exploratory support-sensitive relation, not a detection. No HSC
request or science-image pixel read occurred. Keep science_measurement=false
and images_ready_for_science=false.

Assessment:
analysis/icl_wgl_membership/hsc_debug/s23b_icl_satellite_field8_weight_sensitivity_assessment_2026-11-07.md
Runner:
analysis/icl_wgl_membership/hsc_debug/analyze_s23b_icl_satellite_field8_weight_sensitivity.py
Summary SHA-256 `ec2f6a325ac0a6e349b6f0f52cd106c544d3d811594c35568f5712e3273b2500`.
Outputs are under
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_satellite_field8_weight_sensitivity_2026-11-07_01/`.
No commit or push has been made.

## 2026-11-07 S23B CAMIRA member-position mask audit

Audited all 45,634 member positions against the existing target mask products.
The 1,039 successful mask hashes and image-hash joins were verified; saved-mask
WCS positions reproduce the imported coordinates within 1.08e-7 pixel. Four
targets lack successful masks after “Too few valid pixels”; they contain 292
member rows, including 255 positions inside the image bounds. The remaining
1,039 targets provide 38,911 in-footprint mask samples; 6,431 of their member
positions are outside. Overall, 6,468/45,634 positions lie outside footprints.

CAMIRA centers fall on BASE_MASK pixels in 48.9% of HSC-i and 41.4% of HSC-i2
samples. The corresponding median whole-image BASE_MASK areas are 55.2% and
46.5%; this comparison is not spatially matched. Combined Kron source masks
cover 98.4–99.1% of member centers, likely including each candidate's own
source. These masks omit SEGMENTATION labels, so the own-source mask cannot be
safely separated from neighboring sources. Do not apply them unchanged to
satellite outer-light profiles. Keep science_measurement=false and
images_ready_for_science=false.

Assessment:
analysis/icl_wgl_membership/hsc_debug/s23b_icl_satellite_member_mask_audit_assessment_2026-11-07.md
Runner:
analysis/icl_wgl_membership/hsc_debug/audit_s23b_camira_member_mask_overlap.py
(SHA-256 `c494bda32a606d17bee12780c7271f916a201e0267e29f9353f4a24d18657750`).
Summary SHA-256 `98f4d4e5abe043457a19cba229c2b7c0a7f31e3e07ceeb25d7b5e84d1d713710`.
Output directory:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_satellite_member_mask_audit_2026-11-07_02/`.
Next gate: preserve source segmentation, match member positions to detections,
and create member-specific masks that preserve each candidate while excluding
unrelated sources. No commit or push has been made.

## 2026-11-07 S23B member-centered outer-light profiles and ICL relation

The source-aware exploratory stack and spatial-block relation have now been
completed, superseding the mask-audit section's pending next gate. It used only
the already approved CAMIRA member catalog and HSC images; no data was
downloaded this turn. The all-eligible stack contains 979 target/filter frames
(226 HSC-i, 753 HSC-i2), 16,503 CAMIRA candidate positions at 50–500 proper
kpc, and 13,269 retained candidate-by-radial-bin profile measurements. Each
matched SEP source is preserved while unrelated detections and bad/bright
pixels are masked. Azimuthally rotated controls are required to have at least
60% annular coverage, as are the target annuli, with at least three usable
controls.

The median candidate-minus-rotated-control signal declines steeply with
member-centered radius. At 50–75 proper kpc, the cluster median is 0.68
nJy arcsec⁻² in HSC-i (q16–q84 −3.77 to 5.64; 62 clusters) and 1.46 in
HSC-i2 (−4.23 to 8.94; 365 clusters). These spreads are descriptive, not
confidence intervals.

The adjusted 50–75 kpc outer-light versus 50–300 minus 300–800 kpc relative
ICL slopes use richness, redshift, ICL common support, and paired member
coverage as covariates. Five thousand connected sky-block resamples were run
for each of six spatial layouts. At 5° phase zero, HSC-i is 0.086
[−0.098, 0.334] (62 clusters, 19 blocks); HSC-i2 is −0.041 [−0.145, 0.090]
(365 clusters, 43 blocks). No layout interval excludes zero. The support-floor
sensitivity (5%, 10%, 20%) also has no available interval excluding zero;
the point-estimate sign varies with the floor. This is an unresolved
exploratory association, not evidence for a physical positive or negative
relation.

Keep `science_measurement=false` and `images_ready_for_science=false`.
Uncalibrated member purity/completeness, residual sky, PSF wings, extinction,
source Poisson noise, coadd covariance, centering/mask selection, and possible
cluster light in rotated controls or the outer ICL reference still limit
interpretation. The next gates are authoritative verification of CAMIRA member
fields; member-preserving injection/recovery including PSF wings; calibration
of sky and correlated-noise effects; and center/mask/control sensitivity.

Assessment:
`analysis/icl_wgl_membership/hsc_debug/s23b_camira_member_outer_profile_relation_assessment_2026-11-07.md`.
Runners:
`analysis/icl_wgl_membership/hsc_debug/measure_s23b_camira_member_outer_profiles.py`
and
`analysis/icl_wgl_membership/hsc_debug/analyze_s23b_member_outer_icl_relation.py`.
The full stack and block-adjusted outputs are under
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_camira_member_outer_profiles_full_2026-11-07_01/`
and
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_camira_member_outer_icl_relation_2026-11-07_01/`.
No commit or push has been made.

## 2026-11-07 S23B raw member-field-8 outer-profile sensitivity

The saved source-aware profiles were reweighted by the eighth raw field of
the headerless member file. This field remains **unverified** and must not be
called `w_mem` or a membership probability. The runner joined all 13,269
profile-bin rows to imported member positions; no image pixels or network were
used, and no data was downloaded.

At 50–75 proper kpc, the raw-field-8 weighted cluster median is 1.03
nJy arcsec⁻² in HSC-i (62 clusters; q16–q84 −3.87 to 6.36) and 1.99 in
HSC-i2 (365 clusters; −4.43 to 10.89). The adjusted ICL-versus-outer-light
slopes at 5° phase zero are 0.017 [−0.247, 0.334] in i and −0.013
[−0.118, 0.129] in i2. Both intervals include zero; weighted-slope
intervals include zero in all six spatial layouts. The paired
weighted-minus-unweighted i2 slope is marginally positive in 3/6 layouts and
includes zero in the other three; i changes include zero in all six. Treat
this only as a field-8 sensitivity, not as a detection. Keep
`science_measurement=false` and `images_ready_for_science=false`.

Assessment:
`analysis/icl_wgl_membership/hsc_debug/s23b_camira_member_outer_raw_field8_sensitivity_assessment_2026-11-07.md`.
Runner:
`analysis/icl_wgl_membership/hsc_debug/analyze_s23b_member_outer_raw_field8_sensitivity.py`.
Summary SHA-256 `7156397a478932a9d1fb133b88a758a814bfc0543708719cd36f2a21c1b14875`.
Outputs are under
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_camira_member_outer_field8_sensitivity_2026-11-07_01/`.
Next gate: obtain an authoritative file-specific column definition before
labeling the field, then calibrate member-profile recovery and residual
systematics. No commit or push has been made.

## 2026-11-07 WL/S23B CAMIRA positional-overlap audit

The existing 800-lens S19A weak-lensing sample was compared against all 1,043
S23B CAMIRA target centers and the 427 targets in the member outer-profile
relation. The 435-lens common WL support was reconstructed from the saved
count-only results; its ID hash exactly matches the existing WL stack
provenance. No pair is within 60 arcsec. The nearest S23B center is 364.846
arcsec away. There are 66 equal numeric IDs (29 among profile targets), but
none is positionally matched; do not join across these catalog namespaces by
ID alone.

The current WL and S23B member-profile products cannot form a same-cluster
three-way weak-lensing–ICL–satellite measurement. A same-sample continuation
requires either (a) a new S19A count-only support plan for the S23B profile
targets plus separately validated matched random centers, or (b) authoritative
CAMIRA member positions for the existing WL sample. Any new CAS query or
source-row export needs separate explicit user approval.

Runner: `analysis/icl_wgl_membership/hsc_debug/audit_wl_s23b_camira_crossmatch.py`.
Assessment: `analysis/icl_wgl_membership/hsc_debug/wl_s23b_camira_crossmatch_assessment_2026-11-07.md`.
Output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_camira_overlap_audit_2026-11-07_01/`.

A local HSC S19A count-only plan is prepared for 427 S23B member-profile
targets (18 batches, 3 proper-Mpc cones, integer counts only; 122,198 SQL
bytes). It has not been submitted. Plan SHA-256:
`b713eebded489c280fd9ce1b5840bc8d7d3e42d029035f8dbd414b72fb31e984`.
Separate explicit approval is required before submitting it. Runner:
`analysis/icl_wgl_membership/hsc_debug/prepare_s23b_wl_shape_count_plan.py`.
Submitter: `analysis/icl_wgl_membership/hsc_debug/submit_s23b_wl_shape_counts.py`.
Its `--plan` mode verified the local plan without network use. The submit mode
requires an interactive `lageunha` terminal and the exact explicitly approved
plan SHA-256. No query was submitted.
Keep `science_measurement=false`. No commit or push has been made.

## 2026-11-07 S23B raw-field-8 weighted-density sensitivity

The approved CAMIRA member file was already present and validated; this turn
did not fetch it again. An offline, image-pixel-free sensitivity fit compared
raw-field-8 weighted and unweighted member-density ratios against the existing
relative ICL contrast. Because the public sources do not specify this file's
headerless column order, field 8 remains `field8_raw` (likely `w_mem`), not a
verified membership probability.

For the primary 5-degree phase-zero layout, the HSC-i weighted slope is 5.13
[−0.80, 13.05] versus unweighted 5.72 [−0.06, 13.63]. HSC-i2 is 5.01
[2.44, 7.62] versus 4.24 [0.99, 7.80]. Paired weighted-minus-unweighted
shifts include zero in all six spatial layouts. With the 10% common-pixel
floor, the primary i2 interval includes zero; the relation is not robust to
the support cut. Keep `science_measurement=false`.

Runner:
`analysis/icl_wgl_membership/hsc_debug/analyze_s23b_icl_satellite_field8_weighted_density.py`.
Assessment:
`analysis/icl_wgl_membership/hsc_debug/s23b_icl_satellite_field8_weighted_density_assessment_2026-11-07.md`.
Output:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_satellite_field8_weighted_density_2026-11-07_02/`.

The separate 427-target S19A count-only CAS plan remains unsubmitted; approval
for the member catalog did not authorize that new query. No commit or push
has been made.

## 2026-11-07 CAMIRA member-field schema evidence

An offline audit joined all 637,146 rows of the approved S23B member file to
the 12,162-cluster Wide star-mask catalog. Raw field 8 ranges from 0.001000069
to 1.0 and every row is at or above the published `w_mem >= 10^-3` member
threshold. This strongly supports raw field 8 as probable CAMIRA `w_mem`, but
the official public README does not specify the headerless member-file column
order. Keep it stored as `field8_raw` until a file-specific layout is
confirmed. The earlier field-8 weighting remains an exploratory sensitivity.

Runner: `analysis/icl_wgl_membership/hsc_debug/audit_s23b_camira_member_raw_field_schema.py`.
Assessment: `analysis/icl_wgl_membership/hsc_debug/s23b_camira_member_raw_field_schema_assessment_2026-11-07.md`.
Output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_camira_member_schema_audit_2026-11-07_04/`.
Keep `science_measurement=false`. No new data was downloaded; no commit or push has been made.

## 2026-11-07 S23B PSF-candidate field pilot

Used the already received target HSC images and verified CAMIRA member
positions to run the existing CPU-only `prepare_detection_psf.py` on a
stratified 50-frame subset (25 primary-filter `i`, 25 `i2`). All 50 image
sizes and SHA-256 values matched the local receipt. Forty-nine frames
prepared; `S23B_10647` had too few valid pixels. The 49 frames yielded 492
morphology-only PSF candidates: 44/49 frames have at least four candidates
(22/25 `i`, 22/24 `i2`). The median candidate-stamp valid-pixel fraction is
0.863. Candidates remain unconfirmed stars; stamps are raw 101×101 pixels,
not recentered, normalized, or sky-subtracted, and their 8.4-arcsec radius
does not constrain extended PSF wings.

Cross-matching against the previously approved local S19A Gaia subset found
no candidate within 2 arcsec; the nearest old Gaia source is 523 arcsec away.
This only confirms that the saved Gaia subset does not cover these S23B
positions. A separate, exact Gaia DR3 plan is prepared for 492 S23B candidate
positions in two batches (2-arcsec cones; 18 astrometry, quality, and
classification columns). After the user explicitly approved that exact
492-position scope, the query completed with 77 unique matched positions. No
images, CAS rows, HSC shape rows, or photo-z values were requested; no
scientific measurement was made.

Assessment:
`analysis/icl_wgl_membership/hsc_debug/s23b_detection_psf_pilot_assessment_2026-11-07.md`.
The local Gaia plan SHA-256 is
`394083a8cdaf0770376959e71c8bd46bd004c762dcb3e77a6256bac656ee52f0`.
Keep `science_measurement=false`. No commit or push has been made.

## 2026-11-07 S23B Gaia-supported PSF-wing pilot

The exact approved Gaia DR3 query for 492 morphology-only S23B target
candidates completed in two batches and returned 77 match rows (response
SHA-256 `1f1f8b947a3b23e1dad4675447247137d66d5f2586d6aea3a73989df3cb88e6b`).
An offline audit found 48 conservative diagnostic candidates; this is not a
calibrated star classifier. No S23B coadd header among the 50 inspected frames
contains DATE, MJD, or EPOCH, so Gaia positions were not propagated to an
assumed HSC epoch.

Using already received target cutouts, 36 of the 48 candidates supplied
complete 401×401-pixel stamps (16 `i`, 20 `i2`; 22 frames); 12 were blocked by
image edges. Four frame-weighted profile stacks used 5,000 frame-bootstrap
draws. Only 15.5–22.6% of profile rows passed geometry and sky-coverage gates.
The `i2` outer-500-kpc wing is positive in several bins but changes by factors
of about 2–3.9 when the background annulus moves outward; its farthest-bin
interval crosses zero. The `i` profile changes sign with radius/annulus.
Only 2–5 frames support the outer bins, and there are no matched-control star
profiles. The result therefore does not calibrate extended PSF wings or
correct the ICL stack.

Assessment:
`analysis/icl_wgl_membership/hsc_debug/s23b_gaia_empirical_psf_wing_assessment_2026-11-07.md`.
Gaia audit output:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_gaia_audit_2026-11-07_01/`.
Stamp output:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_wing_stamps_2026-11-07_01/`.
Stack outputs:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_wing_stack_{i_all,i_outer500,i2_all,i2_outer500}_2026-11-07_02/`.
The effective coadd epoch, control-field PSF sample, sky-annulus transfer,
coadd covariance, and member-preserving injection/recovery remain unresolved.
Keep `science_measurement=false` and `psf_approved=false`. No commit or push
has been made.

### 2026-11-07 801-pixel extension

An offline WCS geometry screen found room for 801×801 stamps for 8 `i` and 14
`i2` strict Gaia candidates. The local CPU pass produced 22 stamps across 13
frames; 26 candidates fell too close to image edges. Within the >=500-proper-
kpc subset, only 4 `i` stars from 4 frames and 10 `i2` stars from 8 frames
remain. Profiles extend to 200 pixels (33.6 arcsec) and use sky annuli
260–300 and 340–390 pixels. Beyond 96 pixels, `i2` bins use 1–2 frames and
are consistent with zero under the frame bootstrap; the farthest `i` bins
have only one supported frame. The apparent positive inner `i2` wing from the
401-pixel run is not confirmed by this nonpaired, smaller subset. Sky-annulus
transfer is still unresolved, so no ICL PSF correction is supported.

Expanded stamps:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_wing_stamps_2026-11-07_02/`.
Expanded profile stacks:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_wing_stack_{i_all,i_outer500,i2_all,i2_outer500}_2026-11-07_03/`.
Assessment:
`analysis/icl_wgl_membership/hsc_debug/s23b_gaia_empirical_psf_wing_assessment_2026-11-07.md`.
The original 492-position target and 458-position pilot-control Gaia scopes
are complete. The later 4,299-position expansion query was also approved and
completed. Audit and profile results appear in the final expansion section.

#### Matched-control PSF candidate checkpoint

Local candidate preparation is complete for the 49 already received matched
target/control images. The controls yielded 458 morphology-only candidates in
48 fields (180 `i`, 278 `i2`); the target side has 492 (229 `i`, 263 `i2`).
Candidate coordinates have no target/control overlap within 2 arcsec; the
nearest is 58.12 arcsec. No two control candidate coordinates are within
2 arcsec. These are morphology counts, not Gaia-confirmed stars or approved
PSF models.

A separate Gaia approval was obtained for exactly the 458 control positions.
The two coordinate-specific batches (250 and 208 cones, 2 arcsec) completed
with 58 match rows (estimated 72; 14,326 response bytes estimated). The plan
pins the exact ADQL, TAP endpoint, request parameters, returned fields, and
candidate list; its SHA-256 is
`cbf9bb94c73738b06cd5f616ba1c190f4f0e11ec71a5d914333b4b291d0e7daf`.
No images, shape rows, photo-z values, or CAS data were requested. The audit
found 38 strict diagnostic candidates, 20 matched without clean astrometry,
and 400 Gaia non-matches. Together with the target audit, this yields 48
strict target and 38 strict control candidates; these thresholds are not a
calibrated stellar/PSF classifier.

The control sample produced 26 usable 401×401 stamps from 38 strict
candidates; 12 were edge-blocked. Same-grid control stacks were run for `i`,
`i2`, and the center-distance-≥500-kpc subsets. In the paired target/control
profile comparison, every radial bin has at most one matched frame pair with
valid profiles. A paired PSF-wing difference and its uncertainty cannot be
estimated from this pilot. For controls, center distance is from the random
sky position. Keep `science_measurement=false` and `psf_approved=false`.

Outputs:

- Gaia audit: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_control_gaia_audit_2026-11-07_01/`
- Control stamps: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_control_wing_stamps_2026-11-07_01/`
- Four control stacks: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_control_wing_stack_{i_all,i_outer500,i2_all,i2_outer500}_2026-11-07_01/`
- Paired support comparison: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_target_control_profile_comparison_2026-11-07_03/`

Detailed assessment:
`analysis/icl_wgl_membership/hsc_debug/s23b_gaia_empirical_psf_wing_assessment_2026-11-07.md`.

### 2026-11-07 S23B matched PSF expansion and Gaia completion

The local screen selected 100 new matched pairs in each actual HSC-I
subfilter, `i` and `i2`, across five target-seeing strata. The selection
excluded 49 pilot pairs and 17 pairs with target/control subfilter mismatch.
The target screen used 7,505 CAMIRA member positions inside the selected
images. Local CPU preparation completed for 199/200 target frames and all 200
controls. `S23B_11786` had fewer than 100 valid pixels. Three additional
target frames had no compact candidate. The candidate list contains 2,069
target and 2,230 control positions. Per-filter counts are target/control
1,016/1,000 in `i` and 1,053/1,230 in `i2`.

The user approved the exact Gaia DR3 plan for all 4,299 positions. The plan
SHA-256 is
`c00eda5f896f8160022a7e4198b3f8ce69803f0e30adc8f560ce04c8beba0e7d`. The
18 batches returned 706 candidate-source rows. The aggregate response SHA-256
is `151f7c3f227ab7a16b1a6015563afd46fd53ccf4067d44d0264e27be75f665fb`.
Local audit reproduced all associations and verified batch SQL, raw response,
and match-table hashes. The 706 matches map to unique candidates and unique
Gaia sources. No multiple-source cone or reused Gaia source was found. The
maximum separation is 0.2876 arcsec. The conservative diagnostic screen found
448 candidates, with 225 target and 223 control candidates. The screen is not
a calibrated star or PSF classifier. The query returned no images, CAS rows,
HSC shape rows, or photo-z values.

The local CPU stamp pass produced 168 target and 166 control 401×401-pixel
stamps from existing approved HSC images. Target exclusions were 57 edge
crossings. Control exclusions were 56 edge crossings and one known variable.
The code generated eight frame-weighted stacks for target/control, `i`/`i2`,
and all-star/outer-500-kpc subsets. Each stack used 5,000 frame-bootstrap
draws. In the paired comparison, `i` had no radial bin with two valid matched
pairs. The `i2` all-star profiles had at least two pairs in 43/44 bins, with a
maximum of 14. The `i2` outer-500-kpc profiles had at least two pairs in 30/44
bins, with a maximum of eight.

The `i2` all-star target and control frame-weighted means are negative from
32 to 96 pixels under both tested sky annuli. The paired control-minus-target
difference at 24–32 pixels changes from −6.7 to +1.7 × 10⁻⁶ core flux per
pixel when the sky annulus moves from 120–150 to 160–190 pixels. The 16th–84th
percentile bootstrap intervals include zero for both estimates. A few
outer-500-kpc bins have nominal 95% intervals outside zero with only three
matched pairs. Annulus sensitivity and low pair counts prevent a PSF-wing
calibration. Keep `science_measurement=false` and `psf_approved=false`.

The effective coadd epoch remains unknown. The previously approved 67,071-row
CCD metadata export covers 800 S19A lens centers and 798 controls. The prior
coordinate audit found no match within 60 arcsec between the S19A lenses and
S23B CAMIRA centers. The saved MJD rows cannot assign coadd epochs to the 200
S23B matched pairs. Any new metadata query needs separate exact approval.

Outputs

- Gaia audit runner `analysis/icl_wgl_membership/hsc_debug/audit_s23b_psf_expansion_gaia_matches.py`
- Gaia audit results `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_expansion_2026-11-07_01/gaia_audit_2026-11-07_01/`
- Target and control stamps `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_expansion_2026-11-07_01/{target,control}_wing_stamps_2026-11-07_01/`
- Eight profile stacks `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_expansion_2026-11-07_01/{target,control}_wing_stack_{i,i2}_{all,outer500}_2026-11-07_01/`
- Paired comparison `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_expansion_2026-11-07_01/target_control_profile_comparison_2026-11-07_01/`
- Detailed assessment `analysis/icl_wgl_membership/hsc_debug/s23b_gaia_empirical_psf_wing_assessment_2026-11-07.md`

## 2026-11-07 S23B coadd residual covariance

CPU-only diagnostic completed for the exact 1,022 exploratory target/control
pairs (241 HSC `i`, 781 HSC `i2`), using existing coadds, VARIANCE planes,
and saved K=2.5 masks. All 2,044 frames were read; 2,043 had at least eight
eligible outer-field patches. The median zero-lag residual variance is about
1.14 times the delivered variance plane for both targets and controls. The
finite 9×9-lag sum gives approximate variance factors 1.586/1.580 for `i`
target/control and 1.572/1.555 for `i2`; corresponding median standard-
deviation factors are about 1.25. Paired target-minus-control factor medians
are +0.023 and +0.029, with descriptive frame-pair 16th–84th spreads crossing
zero. This is consistent with similar short-lag residual covariance in both
populations, not a target-specific covariance detection.

The finite lag kernel was not propagated through Lanczos reprojection or the
signed annular estimator; outer patches can retain sources and diffuse light.
Do not use these factors to correct profile errors. Source Poisson, residual
sky, PSF wings, selection effects, and covariance-estimator uncertainty remain
unresolved. Keep `science_measurement=false` and `images_ready_for_science=false`.

- Assessment: `analysis/icl_wgl_membership/hsc_debug/s23b_icl_residual_covariance_assessment_2026-11-07.md`
- Runner: `analysis/icl_wgl_membership/hsc_debug/estimate_s23b_icl_coadd_noise_correlation.py`
- Results: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_residual_covariance_2026-11-07_01/`
- Result summary SHA-256: `8034a140deb4e52cb9c47fff7a2475dcaec19f0410305c5c0520995ab9eedbe1`
- No new data/query, commit, or push was performed for this diagnostic.

No commit or push has been made.

## 2026-11-07: S23B relative-contrast covariance propagation

The local CPU run completed the conditional finite-lag covariance propagation
for all 1,022 matched target/control pairs (241 HSC `i`, 781 `i2`). It used the
same bilinear metric-grid warp, saved K=2.5 masks, pair-common pixels, and
50–300, 100–300, or 200–300 proper-kpc means minus the 300–800 proper-kpc
reference annulus. Replayed means and all six propagated pair counts match the
saved profile and 5-degree block-support products. Supported counts are
231–232 in `i` and 752–756 in `i2`. The conditional finite-lag pixel-noise σ
for the equal-pair mean is 1.089–1.098 times the delivered-VARIANCE diagonal
estimate. It remains smaller than the existing block intervals, but is not a
complete uncertainty or a correction. Outer-field kernels were assumed
stationary in cluster annuli; long-range sky covariance, source Poisson, PSF
wings, and other systematics remain. One missing `i2` target kernel was
median-imputed only in three rows that lacked common support and did not enter
the summary. Keep `science_measurement=false` and
`images_ready_for_science=false`.

- Runner: `analysis/icl_wgl_membership/hsc_debug/propagate_s23b_icl_relative_contrast_covariance.py`
- Assessment: `analysis/icl_wgl_membership/hsc_debug/s23b_icl_relative_covariance_propagation_assessment_2026-11-07.md`
- Results: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_relative_covariance_propagation_2026-11-07_02/`
- Summary SHA-256: `9c1ef1740490f7cbc4811a9bd2573d69405f5482d7ce754a09642108507fc94e`

The separately approved Gaia expansion plan for 4,299 PSF-candidate positions
was already completed and audited: 18/18 batches, 706 candidate-source rows,
aggregate response SHA-256
`151f7c3f227ab7a16b1a6015563afd46fd53ccf4067d44d0264e27be75f665fb`.
Its existing audit is complete and local-only; it does not approve PSFs or
promote science measurements. No duplicate Gaia query was submitted.

## 2026-11-07 CAMIRA member outer-profile injection recovery

Completed a CPU-only measurement-transfer check using one seeded-random,
already-profiled member in each of 427 S23B target/filter frames (62 HSC `i`,
365 `i2`). The saved baseline profile replays exactly in all 427 frames.
Across three exponential-halo amplitudes, 1,221/1,281 injected rows retain
adequate support. Median recovered/input ratios are 0.838–0.896 in `i` and
0.807–0.819 in `i2`; all paired amplitude-difference intervals include zero.
Detection counts or segmentation footprints change in many frames. This is
conditional transfer calibration for previously detected members, not a
completeness estimate or an ICL–satellite detection. No PSF convolution or
source Poisson noise was included, and rotated controls were not injected.
Keep `science_measurement=false` and `images_ready_for_science=false`.

- Runner: `analysis/icl_wgl_membership/hsc_debug/run_s23b_member_outer_profile_injection_recovery.py`
- Assessment: `analysis/icl_wgl_membership/hsc_debug/s23b_camira_member_outer_injection_recovery_assessment_2026-11-07.md`
- Results: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_camira_member_outer_injection_recovery_2026-11-07_01/`
- Summary SHA-256: `e9d4cced59f0b0f8ea3f189989886f28af74a339cc5a9b3421f8bca187defc77`

No commit or push has been made for this checkpoint.

## 2026-11-07 S23B member-profile PSF request checkpoint

Prepared a local, unsubmitted request for position-specific HSC PDR3 Wide
coadd model-PSF FITS cores at the 427 CAMIRA member positions used by the
source-aware injection recovery (62 `i`, 365 `i2`). Existing audited PSF
banks have no exact coordinate matches; the nearest position is 69.6869 arcsec
away. The size estimate is 8,608,320 bytes based on the largest audited PSF
file and is a hard review cap; stop if a response exceeds it. This scope asks
for no images, shape rows, or photo-z. The plan was not submitted, and
`science_measurement=false` remains unchanged.

- Planner: `analysis/icl_wgl_membership/hsc_debug/prepare_s23b_member_profile_psf_plan.py`
- Assessment: `analysis/icl_wgl_membership/hsc_debug/s23b_member_profile_psf_request_assessment_2026-11-07.md`
- Prepared plan: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_camira_member_profile_psf_plan_2026-11-07_03/`
- Plan SHA-256: `03f974ad24ead32eabd103ef7a6e0d6b52317aa2647023a5f7402003a2e3b25d`
- Coordinate list SHA-256: `7327d411bacab15954538c64cb44deae9141955377325bc459a1e6a162db7b31`

The existing HSC PSF batch auditor and member-profile injection runner have
been extended for this plan schema. The future run requires the approved plan
SHA-256 and records all PSF file checksums; it preserves signed PSF pixels.
Neither retrieval nor the PSF-convolved run has occurred.

No request, commit, or push was made for this checkpoint.

## 2026-11-07 기존 S19A 약중력 원천행의 S23B 재사용 감사

기존 800개 S19A 렌즈용 3,069,544개 lens-source 행을 32개 배치의 해시,
바이트 수, 행 수, 중복 source ID 필드 일관성, 원래 3 proper-Mpc 원뿔
포함 여부까지 로컬 검증했다. 고유 source ID는 2,929,303개이며 반복
행 사이 원천 필드 충돌은 없다. S23B member-profile 중심 427곳 가운데
기하학적으로 기존 쿼리 원뿔과 겹치는 곳은 43곳, 실제 기존 source 행이
있는 곳은 21곳뿐이다. 단일 기존 원뿔이 완전히 덮는 S23B aperture는
없고, 기존 원뿔 합집합의 면적 커버리지 추정은 중앙값 0, 90백분위
0.00713, 최댓값 0.604다. 따라서 이 반출 자료로 새 S23B 표본의 shape
지원을 대체할 수 없다.

별도 427중심 count-only CAS 계획은 준비되어 있지만 미제출 상태이며,
Gaia 승인과 다른 자료 조회다. 정확한 계획 SHA-256
b713eebded489c280fd9ce1b5840bc8d7d3e42d029035f8dbd414b72fb31e984에 대한
별도 승인이 있어야 제출한다. count 결과 이후에도 matched random-center
지원과 별도 shape-row 계획 검증이 필요하다. science_measurement=false.

- 감사기: analysis/icl_wgl_membership/hsc_debug/audit_wl_s23b_reuse_existing_s19a_shape_rows.py
- 상세 평가: analysis/icl_wgl_membership/hsc_debug/wl_s23b_existing_export_reuse_assessment_2026-11-07.md
- 결과: /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_existing_export_overlap_2026-11-07_01/
- summary SHA-256: 61026a5d82d47303724f576922415c451589f2516f39337adee3abe021c2c54a
- 이 감사에서 원격 조회·자료 전송·커밋·푸시는 하지 않았다.

### 2026-11-07 S19A 렌즈·무작위 중심 원천행 합집합 재사용 감사

기존 승인된 3,100개 random-center shape 행 9,010,324개(124 배치)를
저장된 바이트 수·행 수·SHA-256·계획 schema와 대조하고, 800개 S19A 렌즈
원천행과 `object_id` 단위로 합쳐 427개 S23B member-profile aperture를
로컬 검색했다. 두 반출 사이 반복 ID의 source 필드 충돌은 없었다.
기존 자료의 행은 표적 32곳에서만 발견됐고, 395곳에는 행이 없다. 렌즈
반출만 21곳, random-center 반출만 11곳, 양쪽을 포함해 32곳이다.

두 반출 원뿔 합집합의 4,096점 Sobol 기하 커버리지 추정은 중앙값 0,
p90 0.310, 최댓값 0.994이며, 95% 이상은 2곳, 99% 이상은 1곳뿐이다.
이는 실제 catalog 깊이·완전성의 상한 진단이지 보정된 약중력 측정이
아니다. B-mode polarity, n(z), lensing efficiency는 보정하지 않았고
`science_measurement=false`를 유지한다. 427중심 CAS count-only 계획은
여전히 미제출이며 별도 승인이 필요하다.

- 감사기: `analysis/icl_wgl_membership/hsc_debug/audit_wl_s23b_reuse_existing_s19a_sources.py`
- 평가 문서: `analysis/icl_wgl_membership/hsc_debug/wl_s23b_existing_export_reuse_assessment_2026-11-07.md`
- 결과: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_existing_source_union_2026-11-07_01/`
- summary SHA-256: `876bed400f93a55c391c89128aa9e9214ed489e664ea0e618ea3df6920f91ed7`
- 이 로컬 감사는 네트워크 사용, 새 자료 전송, commit, push 없이 완료했다.

## 2026-11-07 Integrated S23B ICL–satellite–PSF–WL checkpoint

An offline report builder now consolidates the S23B ICL stack, satellite
relations, injection recovery, Gaia PSF-candidate audit, and the mismatch with
the current S19A weak-lensing sample. The completed Gaia expansion matches the
approved 4,299-position plan: 18/18 batches and 706 one-to-one Gaia sources;
the 448 strict astrometric candidates remain diagnostics, not PSF approval.

The stack's 50–300 minus 300–800 kpc contrast is positive in `i` and `i2`, but
common-mask support and unmodeled systematics still limit interpretation. The
unweighted `i2` satellite-density relation weakens to an interval including
zero with a 10% support floor; both adjusted member-centered outer-light
relation intervals include zero. The prior S19A sample has no center within
60 arcsec of an S23B target, and reused lens/random rows cover only 32/427
S23B apertures. The S23B-specific count-only query is still unsubmitted and
requires its own approval. No source-row query is authorized by that plan.

- Assessment: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_s23b_wl_evidence_2026-11-07_01/assessment.md`
- Summary SHA-256: `5684f7e788c72e482ec5f6bbff9e205ffae67a124c0f36fac89ee0fd1055d0fe`
- Builder: `analysis/icl_wgl_membership/hsc_debug/build_s23b_icl_wl_evidence_report.py`
- Pending count-only plan SHA-256: `b713eebded489c280fd9ce1b5840bc8d7d3e42d029035f8dbd414b72fb31e984`
- The report build was local-only; it made no network request or data transfer.
- Keep `science_measurement=false`, `images_ready_for_science=false`, and
  `psf_approved=false`.

## 2026-11-07 follow-up: S23B satellite association inference

A 10,000-draw Webb six-point restricted wild-cluster bootstrap-t reused the
988 pinned S23B pair/filter records. It used a shared 5-degree connected
target/control partition (49 union blocks) across both filters and four
common-support definitions. Three of eight raw p-values were below 0.05; one
`i2` result at >=20% support survived Holm (slope 4.21 nJy arcsec^-2 dex^-1,
raw p=0.0057, Holm p=0.0456), but none passed complete-null max-T (minimum
p=0.2666). Because the support cuts reuse the same images and the relative ICL
measurement retains unresolved sky/PSF/mask/covariance terms, treat this as
exploratory sensitivity evidence and keep `science_measurement=false`.

- Assessment: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_satellite_relation_wild_2026-11-07_01/assessment.md`
- Summary SHA-256: `d95f894c6cf54e92d906063fdc920c24e80bdf9ad4194d2ff63dacd8bdd42494`
- Runner: `analysis/icl_wgl_membership/hsc_debug/analyze_s23b_icl_satellite_relation_wild.py`

## 2026-11-09 S23B weak-lensing count-only query and audit

The user-approved S19A count-only plan (SHA-256
`b713eebded489c280fd9ce1b5840bc8d7d3e42d029035f8dbd414b72fb31e984`) completed
all 18 batches. The resumed run retrieved job `1070884`, completed jobs
`1071485` and `1071508`, and verified 427 integer counts. The local audit
verified every batch hash and center ID against the approved plan and joined
the counts to the saved target manifest. Credentials were not saved; no source
rows, coordinate rows, shape values, photo-z, or images were returned.

There are 356,191 lens-source candidate pairs summed over the 427 cones, but
329 centers (77%) have zero catalog entries within 3 proper Mpc. Only 82 have
at least 1,000 entries, 62 at least 3,000, and 25 at least 5,000; the median is
zero and the maximum is 8,959. Descriptive Spearman coefficients are 0.060 for
count versus redshift and -0.145 for count versus CAMIRA richness; no p-values
were computed. A local cross-tab by the already saved target coordinates shows
strong sky-region variation in the zero counts, so the next analysis must
resolve HSC footprint/support selection before defining a lens sample.

These are availability counts from the released S19A weak-lensing shape subset,
which already applies the release weak-lensing flag and star mask. The project
has not applied source-redshift, tomographic, B-mode, or radial cuts. Counts
are lens-source candidate pairs, not unique sources or a lensing measurement.
Keep `science_measurement=false`.

- Audit: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_shape_count_audit_2026-11-07_01/assessment.md`
- Audit assessment SHA-256: `f2ce0c5f22c2f6f894386f0c682d08aaf21236d5dd56188cf896b65064e6d672`
- Per-center table SHA-256: `d650564b4870d4af93ce04b9ebd82475a99366aa68c2dfe236e7f070826c0f7f`
- Audit runner: `analysis/icl_wgl_membership/hsc_debug/audit_s23b_wl_shape_counts.py`

This local sample and row-plan stage was completed on 2026-11-10. Submitting
the new 400-center count-only query needs separate approval tied to its exact
plan hash. Shape-row export also requires separate approval after the observed
counts and transfer size are known.

## 2026-11-10 S23B weak-lensing common-support sample and export plan

The approved 427-target, 3-proper-Mpc S19A count-only results were audited
locally. Counts are zero for 329 targets; 82 have at least 1,000 entries. These
are raw availability counts from the released weak-lensing shape subset, before
project-specific shape, tomographic, B-mode, and radial cuts. They are not
effective source counts or a lensing measurement.

Of the 82 targets above the provisional 1,000-count screen, 80 have five
qualifying controls in the saved candidate pool. The provisional design has 80
lenses and 400 HSC random-table control positions (five per lens). Control
positions are globally unique, pass the saved center-flag, depth, seeing, and
SFD-EBV criteria, have zero in the earlier CAMIRA 2-Mpc overlap screen, and
lie more than 6 proper Mpc from their paired lens so the two 3-Mpc query cones
do not overlap. The CAMIRA screen is incomplete and does not remove all
environmental or large-scale shear contributions.

The selected lenses have 340,256 verified count-only lens-source pairs; the
per-lens minimum/median/maximum is 1,040/4,323/8,959. Previously selected ICL
controls are unsuitable for weak-lensing control: 114/427 lie within 3 Mpc of
their paired lens and 280/427 within 6 Mpc. Other S23B CAMIRA centers are within
projected 3 Mpc of 110 of those controls at any redshift, and 23 at
|dz|<=0.05. This is a center-proximity diagnostic, not a mass/shear model.
The old S19A lens/random source-row union has rows in only 29 of the selected
80 apertures and none reaches 95% estimated query-cone area coverage.

The prepared next count-only plan covers 400 control centers in 16 batches of
25 and returns only `random_center_id`, `matched_lens_id`, and the integer
`n_shape_sources`. It passes the existing random-count plan validator but is
still `prepared_not_submitted`. It must not be submitted without separate
approval for exact plan SHA-256
`fe19ac426272d4005e3ee74a902dfdb52d1fe83ecb2cc5ec901c10b880d28e8f`.

The later shape-row plan has 80 lens centers in 4 batches and 400 controls in
16 batches. The lens-side expected row total is 340,256 based on the approved
counts; the control row total is unknown until the count query completes. A
proxy from a different completed sample gives 1.05/1.83/2.23 million combined
pairs (q16/median/q84) and 0.25/0.55/1.00 GiB uncompressed CSV (low/point/high).
These are planning estimates, not observed S23B counts. Source-row query and
transfer require another separate approval.

Lens rows have 22 columns and control rows 23, preserving shapes, calibration
terms, PSF moments, the raw B-mode mask, and the HSC Y3 source-redshift bin.
Full object-level photo-z PDFs and images are excluded. Identical shape,
source-redshift, B-mode, and radial cuts must be applied to lenses and controls;
verify B-mode polarity and calibrated n(z)/lensing efficiency before a
measurement. Keep `science_measurement=false`. Neither count nor shape-row
queries were submitted and no network was used in this stage.

- Planner: `analysis/icl_wgl_membership/hsc_debug/prepare_s23b_wl_common_support_plan.py`
- Assessment and plans (outside Git): `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_wl_common_support_plan_2026-11-10_01/`
- Row-plan file SHA-256: `3ec3c8e77cf9cf38fdd02320718c6e38ee3d23e8948584255cfdc27bfe90cfee`
- This stage used local analysis only; no remote query or data transfer occurred.

### 2026-11-10 승인된 대조 count-only 실행 및 원천행 계획 갱신

별도 승인된 count-only 계획 SHA-256
`fe19ac426272d4005e3ee74a902dfdb52d1fe83ecb2cc5ec901c10b880d28e8f`를
lageunha에서 실행했다. 16/16 batch와 400개 결과 행의 ID, 정수 count,
영수증 해시를 확인했다. 결과는 무작위 중심 ID, 짝지은 렌즈 ID, 정수
count뿐이며 shape 값·좌표·photo-z·영상은 반환되지 않았고 자격 증명은
저장되지 않았다.

대조군 400개 원시 count 합은 1,572,453개이며 18개 중심은 0, 382개는
양수다. 80개 렌즈 count 합 340,256개와 합친 source-pair 예상치는
1,912,709개다. 렌즈별 대조 평균 count/렌즈 count 비율의 q16/중앙/q84는
0.690/0.936/1.295다. 렌즈 7곳은 대조 5곳 중 적어도 하나가 0이고, 5곳
모두 0인 렌즈는 없다. 이는 릴리스 원시 count이며 프로젝트 source 선택
후의 유효 공통 지지나 lensing 측정은 아니다.

count 결과로 갱신한 원천행 계획은 렌즈 80개 4 batch, 대조 중심 400개
16 batch다. 예상 총 행 1,912,709쌍, 무압축 CSV는 0.46/0.58/0.86 GiB
(낮음/중앙/높음)이다. byte 수는 과거 export의 행당 크기를 사용한 추정으로
실제 다운로드 용량은 아니다. 이 원천행 계획은 여전히
`prepared_not_submitted`이며 별도 hash 승인 전에는 제출하지 않는다.

- Finalizer: `analysis/icl_wgl_membership/hsc_debug/finalize_s23b_wl_common_support_plan.py`
- Count audit와 갱신 계획: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_wl_common_support_plan_2026-11-10_01/`
- 원천행 계획 파일 SHA-256: `7bfb998149f30a50b5ee7a537bb703ed302b7ea269a3296c35204a836e75aabc`
- 원천행 query/download 미실행; `science_measurement=false` 유지.

### 2026-11-10 원천행 계획 승인 및 실행기 사전검증

사용자가 원천행 계획 파일 SHA-256
`7bfb998149f30a50b5ee7a537bb703ed302b7ea269a3296c35204a836e75aabc`의
제출과 수신을 승인했다. 범위는 렌즈 80개와 대응 대조 중심 400개,
렌즈 source-pair 340,256개와 대조 1,572,453개로 총 1,912,709개다. SQL은
렌즈 4개, 대조 16개 batch이며 출력 열은 각각 22개와 23개다. 영상과 전체
photo-z PDF는 포함하지 않는다. 무압축 CSV 예상량은 낮음/중앙/높음
0.46/0.58/0.86 GiB이며 실제 전송량은 아직 측정 전이다.

새 실행기
`analysis/icl_wgl_membership/hsc_debug/submit_s23b_wl_common_shape_rows.py`의
lageunha `--check`는 네트워크 없이 통과했다. 불변 계획 해시, count-only
영수증, 중심별 예상 수, SQL 내용과 해시를 확인했다. 이 기록 시점에는
아직 소스 행 query와 다운로드를 시작하지 않았다. 다운로드 후에는 배치별
CSV hash/크기/행 수와 count-only 대조를 기록하고 `science_measurement=false`를
유지한다.

### 2026-11-10 승인된 S23B weak-lensing 원천행 수신 완료

정확한 계획 SHA-256
`7bfb998149f30a50b5ee7a537bb703ed302b7ea269a3296c35204a836e75aabc`로
lageunha에서 원천행 query 20개(렌즈 4, 대조 16)를 순차 실행했다. 모든
batch가 count-only 결과와 일치해 완료됐다.

- 렌즈 80곳: 340,256 lens–source 쌍, 102,103,795 bytes
- 대조 400곳: 1,572,453 control–source 쌍, 523,062,067 bytes
- 합계: 1,912,709 center–source 쌍, 625,165,862 bytes (0.582 GiB)
- 파일: 각 렌즈 22열, 대조 23열. 영상과 전체 photo-z PDF 없음.
- 출력: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_wl_common_support_plan_2026-11-10_01/source_rows_post_count_01/`

배치별 count-only 행 수, `(center_id, object_id)` 유일성, 열 스키마,
수치 유한성, 좌표 범위, mask/z-bin 정수 및 shape weight 범위 검사를
통과했다. 별도 재감사에서 20개 CSV 모두 실제 크기와 SHA-256이 checkpoint와
일치했다. 자격 증명은 저장되지 않았다. 행은 center–source 쌍이므로
전역 고유 source 개수로 해석하면 안 된다. 프로젝트별 최종 source 선택 전
릴리스 원천행이며 `science_measurement=false`다.

다음 작업은 다운로드 자료만 사용하는 로컬 source 품질·공통 지지 감사다.
렌즈와 대조에 동일한 문서화된 source cut을 적용해 aperture별 생존 수와
유효 weight를 요약하고, B-mode mask 의미를 확인한다. 보정 source n(z),
lensing efficiency, PSF/공분산 검증 전에는 lensing 측정으로 승격하지 않는다.

## 2026-11-10 S23B weak-lensing 로컬 QA 완료

승인된 로컬 QA를 lageunha CPU에서 완료했다. 이미 받은 20개 CSV와 기존 HSC
Y3 SACC n(z)만 읽었고, CAS 조회·네트워크 요청·추가 다운로드·이미지 읽기·shear
stack은 없었다. 모든 batch의 SHA-256, byte 수, count-only 중심별 행 수,
중심/원천 쌍 중복, 중심 ID, 80개 렌즈와 400개 대조 연결을 확인했다.

렌즈 340,256쌍 중 239,565쌍, 대조 1,572,453쌍 중 1,105,779쌍이 공통 source
cut, 양의 fiducial n(z) 효율, 0.1–3 proper-Mpc, calibration-input 유효성
검사 후 남았다. 렌즈 2곳은 생존 쌍이 0이고 대조 중심 24/400곳은 0이다.
적어도 하나의 대조가 0인 짝은 9/80이며, 다섯 대조가 모두 0인 짝은 없다.
렌즈 지지가 양수인 짝에서 대조 평균/렌즈 유지 수의 q16/중앙/q84는
0.690/0.925/1.311이다. 형상 및
weight 중앙값은 양쪽에서 비슷하지만, 0 지지 짝을 어떻게 다룰지 사전 규칙이
필요하다.

GAMA09H 공개 사각영역 좌표 대조는 정확히 일치했다. raw mask 0인 48,765쌍은
모두 영역 안에 있고, raw mask 1인 1,863,944쌍은 모두 밖에 있었으며 반대
분류는 0건이다. 따라서 이 export에서는 raw 0을 제외값으로 유지한다. 다만
HSC 공식 문서는 mask의 목적은 설명하지만 정수 방향은 지정하지 않는다.

이는 catalog-row 지지/분포 QA 통과이지 lensing 측정은 아니다. SACC n(z)
불확실성, per-object 보정·responsivity·selection-bias, PSF 공간 일치,
대조 오염, 반복 원천/공간 공분산, injection/recovery가 남았다.
`science_measurement=false` 유지.

- 실행기: `analysis/icl_wgl_membership/hsc_debug/audit_s23b_wl_common_shape_rows.py`
- 산출물: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_wl_common_support_plan_2026-11-10_01/quality_audit_2026-11-10_01/`
- summary SHA-256: `9ef472c1ae0888e187fce05f3449b1614d1da7c99e7a400403b92f8aa7e2dc54`
- 다음 단계 제안: 0 지지 짝 처리 기준, 보정 random-subtracted estimator, 공간
  covariance를 명시한 계획을 만들고 승인받은 뒤 shear profile을 계산한다.

## 2026-11-11 반경·tomographic 지지 및 0 지지 대조 민감도

기존 승인된 S23B 원천행 CSV 20개와 기존 HSC Y3 SACC n(z)만 사용해 7개
반경 구간과 4개 source-z bin으로 지지를 세분화했다. CAS 조회, 추가 다운로드,
영상 pixel 읽기, shear profile 계산은 없었다. 13,440개 중심×반경×z-bin 셀은
이전 반경 집계와 정확히 일치하며 공통 선택 후 렌즈 239,565쌍, 대조
1,105,779쌍이다. source가 aperture 간 반복되므로 이 수는 전역 고유 source
수로 해석하면 안 된다.

가장 약한 곳은 0.1–0.2 proper Mpc다. z-bin 1/2/3/4별 지지 렌즈 수는
65/63/59/52 of 80, 대조 수는 285/303/293/262 of 400이다. 바깥 2–3 Mpc에서는
렌즈 78/78/78/78, 대조 376/376/376/376곳에 지지가 있다.

렌즈 지지가 양수인 78개 짝 중 양의 지지를 가진 대조를 최소 3/5개 요구하면
74개 짝(92.5%), 대조 중심 367개, 렌즈 233,064쌍이 남는다. 4/5 기준도 동일하다.
5/5 기준은 71개 짝(88.8%), 대조 355개, 렌즈 224,375쌍이다. 대조/렌즈 유지
count 비율 중앙값은 3/5와 4/5에서 0.940, 5/5에서 0.950이다. 이는 지지
민감도 비교이지 최종 선택 기준은 아니다.

전 구간을 정량적으로 다시 쓴 TSV와 assessment를 보존했다. 결과는 QA이지
lensing 측정이나 질량 추정이 아니다. n(z) 불확실성, 보정 response, PSF 공간
일치, 대조 오염, 반복 source/공간 공분산, injection/recovery가 해결될 때까지
`science_measurement=false`를 유지한다.

- 실행기: `analysis/icl_wgl_membership/hsc_debug/audit_s23b_wl_common_shape_rows.py`
- 산출물: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_wl_common_support_plan_2026-11-10_01/quality_audit_2026-11-11_01/`
- summary SHA-256: `eb73aa0365813d880d3edf60b12e4039f2575f8a97eab8868fc918121a34d445`
- assessment SHA-256: `35369bbc74fa576570860af09d95d29c33c32cf8b1dfb96c30bb5949d5ddcf95`
- per-center radial/tomographic TSV SHA-256: `f8d604e91c43705cb24152c7e98553efd8838933034a8657bcc603c2f8c11bed`
- aggregate radial/tomographic TSV SHA-256: `ae84048ee707ffe531ecff80ede270544a7475846079db7520ee59bd62706911`
- zero-support policy TSV SHA-256: `ba89948681cb0560408030a4df66a232a2989cd5415c34ca1cb99c59d8bc26da`
- 다음 단계 제안: 동일 데이터만 사용해 보정 random-subtracted profile의
  구현 계획과 0 지지 기준, 공간 covariance block 정책을 구체화하고, 해당
  계획 승인 후 탐색 profile 계산을 수행한다.

## 2026-11-11 S23B matched-sample random-subtracted 약중력 프로파일

이전 QA에서 승인된 80개 렌즈/400개 matched-control HSC S19A shape 행과
고정된 HSC Y3 SACC `n(z)`만 사용해 `lageunha` CPU에서 탐색 프로파일을
계산했다. CAS 조회·네트워크 요청·추가 다운로드·영상 pixel 읽기는 없었다.
4/5 positive-control 기준으로 74개 렌즈 세트와 367개 대조 중심이 남았다.
주 프로파일 source bins 2–4는 0.1–0.2 proper Mpc에서
`23.845 ± 3.432`, 2–3 Mpc에서 `0.762 ± 0.400` (`10^13 Msun Mpc^-2`)였다.
5도 equal-area endpoint jackknife의 명목 tangential zero-profile 진단은
`chi2=221.85`, rank 7, `p=2.67e-44`다. 근사 블록 공분산에 조건부이며 검출
확률로 볼 수 없다.

교차 shear null은 4/5 기준에서 `p=0.0741`, 5/5 기준에서 `p=0.0341`,
2.5도 블록 민감도에서 `p=0.3182`로 변한다. 대조 선택과 블록 폭 민감도가
남아 profile을 검출로 주장하지 않는다. 대조는 pooled stack이어서 매칭 세트별
동일 가중 결과와 비교가 필요하다. `n(z)` 오차, contiguous survey-region
공분산, random 위치의 PSF/depth/mask/sky 일치, reduced shear, magnification,
member dilution/boost, selection-response 오차도 미해결이므로
`science_measurement=false`를 유지한다.

- Astrafek 실행기: `analysis/icl_wgl_membership/hsc_debug/analyze_s23b_wl_random_subtracted_profile.py`
- 상세 평가: `analysis/icl_wgl_membership/hsc_debug/wl_s23b_random_subtracted_profile_assessment_2026-11-11.md`
- 산출물 (Git 외): `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_random_subtracted_profile_2026-11-11_01/`
- summary SHA-256: `2c899b0090dd94a34f1fd6b0cb111e3d1ada696da2178f98e9fd5b659e699edf`
- assessment/TSV/PNG SHA-256: `d152c22f21a244e9848f453523c0a346120e30f4897b9a427d65776aa804ba2a` / `3dd1ef2b7f4ac7aa33d6b046ef9b528ae6c8818d1f3992037b6b70aa6e5548b2` / `4c3593236d5a46dd504d0f54b3f2396e0a925460c27f3f879bc6fef7bef92ce9`
- 다음 단계: 렌즈 세트별 동일 가중 random 민감도와 보유 S19A 전용 `n(z)` 대안을 평가한다. 새 자료 반출은 포함하지 않는다.
## 2026-11-11 S23B 가중 및 DESI `n(z)` 민감도

기존 승인된 S23B HSC S19A shape 행, HSC Y3 SACC, DESI/HSC `n(z)` 자료만으로
오프라인 민감도 분석을 완료했다. CAS·네트워크·다운로드·영상 pixel 접근은
없었다. Runner는 `analysis/icl_wgl_membership/hsc_debug/analyze_s23b_wl_weight_and_desi_nz_sensitivity.py`이고
SHA-256은 `af42929b6c392a2ccf98cf11780de46f141450ec9fc43f4e06af51492094ae18`이다.

기준 표본은 74 렌즈/367 대조, source bins 2–4다. 반경별 공통 지지는 60–74개
렌즈 세트다. 세트 동일 가중 random subtraction은 pooled 기준보다 0.2–0.4 Mpc에서
+22.71%, 0.4–0.7 Mpc에서 −19.39%, 0.7–1.0 Mpc에서 −11.97%다. 대조 세트만
동일 가중해도 1.5–2.0 Mpc 값이 −24.01% 달라진다. 이 가중 variant들의 별도
공간 공분산은 산출하지 않았다. DESI `n(z)` 중앙 곡선은 SACC 대비 약 −2–4%이며,
lower/upper는 점별 주변 envelope라 joint posterior가 아니다.

따라서 pooled random profile은 estimator 선택에 민감하다. 동일 aperture 지지,
zero-support set 처리, footprint-aware covariance 및 cross-shear가 정리될 때까지
mass fit과 detection 해석을 보류한다. 영상 PSF/depth/mask/sky 매칭, reduced shear,
magnification, member dilution/boost, selection response 및 injection/recovery도
미해결이다. `science_measurement=false` 유지.

- 평가: `analysis/icl_wgl_membership/hsc_debug/wl_s23b_weight_and_desi_nz_sensitivity_assessment_2026-11-11.md`
- 산출물 (Git 외): `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_weight_and_desi_nz_sensitivity_2026-11-11_03/`
- summary SHA-256: `9c66b281e0920b79a5c7854e40081b873621867cfaca2a9991d5e93e031f1aa5`
- 이전 `_01`과 `_02` 시도는 각각 비교 지원 불일치와 문구 오류 때문에 superseded 처리했다. 최종 기준은 `_03`이다.

## 2026-11-11 S23B 반경별 matched-control support 후속

최신 Astrafek runner는 승인된 80-lens/400-control QA 결과의 이미 집계된 표만
읽었다. Raw shape CSV, 이미지, HSC 서버, 네트워크를 읽지 않았다. source bins 2–4,
전역 4/5 support sample은 74 lens sets/367 controls다. 양쪽 역할이 지지되는
반경별 세트 수는 60/67/69/73/74/74/74다. 반경마다 4개 이상 positive controls가
있는 수는 52/59/64/67/70/72/74, 5개 모두인 수는 38/51/60/64/65/68/71이다.
모든 반경에서 lens rows와 ≥4/5 controls를 요구하는 고정 교집합은 52세트다.

고정 52세트의 redshift 중앙값은 전역 74세트와 유사하다 (0.44805 대 0.44705).
그러나 렌즈 aperture별 유지 source-pair 중앙값은 3,415 대 3,118.5이며 shape-weight
합 중앙값도 19,682 대 17,662로 높다. 따라서 지지 기반 표본은 자료 support가 더
좋은 렌즈를 고른다. 이 고정 표본 `DeltaSigma_t` (`10^13 Msun Mpc^-2`)는 일곱
반경에서 23.554, 6.048, 4.003, 3.521, 2.521, 0.675, 1.255다. 반경별 60–74개
지지 표본과 비교하면 최대 −54.8%/+61.8% 차이가 난다. per-set 산포는 descriptive
only이며 covariance를 계산하지 않았다. `science_measurement=false`; 질량 적합과
검출 주장을 보류한다.

- support 감사 runner: `analysis/icl_wgl_membership/hsc_debug/audit_s23b_wl_matched_radial_support.py` (SHA-256 `1c924dd6173addca14e0eec36e8698d687e535f659425e52b93a5c18e1ad3ba3`)
- fixed profile runner: `analysis/icl_wgl_membership/hsc_debug/analyze_s23b_wl_fixed_common_support_profile.py` (SHA-256 `3640ec5b6313a4bf1ef48cf8270676a8c4c466e20e274bbef30a952b8257524d`)
- support summary (Git 외) SHA-256: `c51a4d9268935a0919d0abd3da2f10b6cd16ee7e35f1e6e23a8f3fe8c9e51234`
- fixed profile summary (Git 외) SHA-256: `4d270b81737e5bcb6188941f9f26cd6255614fd7ae8d1e6f2977cf51a4686930`
- 결과: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_matched_radial_support_2026-11-11_02/` 및 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_fixed_common_support_profile_2026-11-11_02/`

이 기록에서 앞서 저장한 per-center radial/tomographic QA TSV의 해시는 문서 오타를
수정했다. 검증된 값은 `f8d604e91c43705cb24152c7e98553efd8838933034a8657bcc603c2f8c11bed`
(64자)다.

## 2026-11-11 fixed-52 matched-random spatial covariance

승인된 shape-row 20개 CSV와 기존 HSC Y3 SACC만 사용해 고정 52개 lens-set
프로파일의 paired endpoint 공간 공분산을 `lageunha` CPU에서 계산했다. 매 세트에서
lens와 matched random profile을 각각 재보정하고, 세트 안의 random center는 source-pair
lensing weight로 pooled한 뒤 바깥 52세트 평균은 동일 가중했다. source bin은 2–4이며
모든 7개 반경에서 lens 지지와 대조 ≥4/5를 요구한다. 새 CAS 조회·다운로드·이미지 분석은
없다.

2.5° equal-area RA/sin(Dec) grid는 64개 점유 셀, 5° grid는 28개다. 두 공분산은 모두
14×14, rank 14, PSD이며 deletion마다 고정 표본의 최소 49/52 및 47/52 세트가 남았다.
이 gate는 근사 grid 안정성 검사다. grid는 연결된 실제 HSC survey region이 아니며 이
공분산으로 검출이나 질량을 주장하지 않는다. tangential profile은 `10^13 Msun Mpc^-2`
단위로 0.1–0.2 Mpc에서 `23.554 ± 5.286 (2.5°) / 4.410 (5°)`, 2–3 Mpc에서
`1.255 ± 0.556 / 0.535`다. 전체 점추정은 이전 fixed-52 프로파일과 최대
0.03125 Msun Mpc⁻² 차이로 일치한다. `science_measurement=false` 유지.

- 실행기: `analysis/icl_wgl_membership/hsc_debug/analyze_s23b_wl_fixed_common_support_spatial_covariance.py`
- 실행기 SHA-256: `662b17470b342f0f2486a2fd23472b7105ba6a76e5354fa37b1aa679a0010256`
- 상세 평가: `analysis/icl_wgl_membership/hsc_debug/wl_s23b_fixed_common_support_spatial_covariance_assessment_2026-11-11.md`
- 최종 산출물 (Git 외): `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_fixed_common_support_spatial_covariance_2026-11-11_02/`
- summary SHA-256: `acdd4da69e4d7fe3dbb52d8584b3996bd836a41196da59aecfe3fcdf3e1fc6e9`
- deletion-support / replicate / covariance matrix SHA-256:
  `810e80bdf9cd28bad42c456968e885ab61d417f68a5c823c5fad5ca98ba6b790` /
  `158bbaf5518b91589c6477673b259a79de9deeab997e8069abd0995bd46a8f27` /
  `bc718cb4c5734e7a981048b735ad6239486aeb99cbc6ee26155228906ceed9b4`
- 첫 `_01`은 2.5° block-support 행이 기록되지 않아 superseded했고, `_02`에 두 폭의 92개
  deletion을 모두 보존했다. 수치 결과는 동일하다.
- 다음 관문: 실측 survey footprint 기반 region jackknife, target/control 영상의 PSF·depth·mask·sky
  일치, `n(z)` 및 selection-response 오차다. 해결 전 mass fit과 검출 주장은 보류한다.

## 2026-11-11 fixed-52 random-center weighting covariance

고정 52개 lens set에서 random source pairs를 source-pair lensing weight로 합치는 기존
방식과 각 matched random center를 따로 보정한 뒤 중심별 동일 가중 평균하는 방식을
비교했다. 입력은 이미 승인된 shape-row CSV 20개와 pinned HSC Y3 SACC `n(z)`뿐이다.
CAS/network/download/image 접근은 없었다. 바깥 lens-set 동일 가중, source bins 2–4,
7개 반경의 최소 4/5 control 지지 및 동일 paired endpoint deletion을 유지했다.

Pooled profile은 직전 fixed-52 결과와 정확히 일치한다. Equal-center profile은
`10^13 Msun Mpc^-2` 단위로 `24.204, 5.900, 4.085, 3.379, 2.551, 0.773, 1.301`이며
pooled 값과의 차이는 paired block error의 최대 1.35배다. 다만 covariance partition이
근사이므로 estimator 동등성은 결론 내리지 않았다. Cross-component zero-profile nominal
p는 2.5°에서 pooled/equal-center `0.500/0.446`, 5°에서 `0.064/0.050`으로 달라 null
판정이 grid 폭에 민감하다.

2.5°와 5°에서 64/28개 셀이 점유됐고 두 방법 및 estimator 차이 covariance 모두 rank
14/14·PSD였다. 삭제 후 최소 공통 lens-set 수는 49/52 및 47/52다. equal-area 셀은 실제
연결 HSC survey field가 아니므로 과학적 covariance 관문을 통과한 것으로 보지 않는다.
깊이·mask·PSF·sky 매칭, `n(z)` 오차, member dilution/boost, selection response와
PSF-주입/회수는 남아 `science_measurement=false` 및 mass-fit 보류를 유지한다.

- 실행기: `analysis/icl_wgl_membership/hsc_debug/analyze_s23b_wl_fixed_equal_control_spatial_covariance.py`
- 상세 평가: `analysis/icl_wgl_membership/hsc_debug/wl_s23b_fixed_equal_control_spatial_covariance_assessment_2026-11-11.md`
- 실행기 SHA-256: `f31f7f30fc9d244d20dd6efaccaa51a9b5b833cf76bbb8c05b310aecb7039fc7`
- 결과 (Git 외): `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_equal_random_center_spatial_covariance_2026-11-11_02/`
- summary SHA-256: `00c8427ef6810b78f0900da1d515e3fca2ed77233bd012ebbeddbf5bc5e5c033`

## 2026-11-11 fixed-52 HSC tract-region influence

승인된 S19A shape 행 20개 CSV, HSC Y3 SACC `n(z)`, 기존 patch-QA만으로
고정 52 lens-set의 근사 트랙 연결요소 삭제 영향도를 계산했다. 전역 840개 트랙 중
분석 표본은 세 요소에 놓였다. 렌즈 중심 수는 28/21/3이며, source-bin 및 보정 행
선택 뒤 남은 대조 중심은 139/105/15개다. 260개 계획 대조 중 259개만 선택 행을
가졌다. 영역 삭제 후 공통 지지는 각각 24/31/49 lens-set이고 모든 유지 세트에서
반경별 대조 수는 최소 4개였다.

두 큰 요소 중 하나를 삭제하면 1.5–2.0 Mpc의 tangential profile이 전체 표본 값에
비해 pooled 기준 −117.6% 또는 +116.1% 이동한다. Equal-center 방식도 −101.1% 또는
+99.6%로 같은 민감도를 보였다. 표본 교집합이 삭제 후 달라지므로 이 값은 순수한
면적 효과가 아닌 영향도 진단이다. 세 요소로 안정적인 covariance를 추정하지 않았고
`science_measurement=false`를 유지한다. 형상 행에 tract/patch ID가 없어 nearest
patch-QA 중심 분류를 사용했으며 공식 필드 경계가 아니다.

- 실행기: `analysis/icl_wgl_membership/hsc_debug/analyze_s23b_fixed52_tract_region_influence.py`
- 평가: `analysis/icl_wgl_membership/hsc_debug/wl_s23b_fixed52_tract_region_influence_assessment_2026-11-11.md`
- 실행기 SHA-256: `d21057f9a5176f22fdd01d1bf18dbb6302ca26effd6892b6528cc6b3df015911`
- 산출물 (Git 외): `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_fixed52_tract_region_influence_2026-11-11_01/`
- summary SHA-256: `afc151a6b189af182f975b7a70a2e8bc176b85fa62ab09bff7a2b8718cb41326`
- 상세 per-bin shifts SHA-256: `b2f1e5ac95c42354cca8c3e593c8207527fcc23dae92e136c03f6593f5640d69`

## 2026-11-11 공식 footprint 경계 파서 준비 및 정적 점검

근사 중심 배정의 한계를 줄이기 위해
`analysis/icl_wgl_membership/hsc_debug/hsc_pdr3_tract_patch_geometry.py`를 추가했다.
HSC PDR3 Wide 4개 tract/patch 좌표 텍스트를 검증·파싱하고, 구면 좌표를 실제
tract/patch 다각형에 할당한다. 파일을 자동 다운로드하지 않으며 아직 실제 공식 파일로
실행하지 않았다. 다운로드 승인이 오기 전까지 데이터 폴더에 파일을 받지 않는다.

정적 점검에서 새 tract로 넘어갈 때 이전 patch ID가 남을 수 있는 상태 오류를
고쳤다. `field` 모드는 tract polygon으로 field만 지정하고 tract/patch ID는
비워 둔다. 합성 또는 공식 경계 데이터로 실행하지 않았으므로 파싱·다각형의
실제 검증은 아직 남아 있다.

- 공식 자료 용량: 정확히 49,225,826 bytes (4개 Wide 필드 좌표 파일)
- 저장 예정: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/hsc_pdr3_wide_tract_patch_geometry_2026-11-11_01/`
- 다음 단계: 승인 후 바이트 수·스키마·다각형 중심을 검증하고, 승인된 고정 52렌즈와
  source rows를 tract에 배정해 tract-block 영향도 및 공분산 조건을 평가한다.
- 그 전까지 weak-lensing 결과는 근사 영역 영향도이며 `science_measurement=false`다.

### 2026-11-11 Astrafek exact-tract 영향도 실행기 준비

기존 고정 52표본 영역 영향도 실행기에 `--geometry-dir` 선택 모드를 추가했다.
기본 patch-QA 근사 모드는 유지된다. exact 모드는 공식 HSC PDR3 Wide tract polygon으로
렌즈·대조 중심과 고유 source 좌표를 할당한 뒤 tract 하나씩 제외해 calibrated profile의
영향도를 계산한다. 반복 object ID의 좌표 대입은 한 번만 수행한다. 활성 행에 미할당 또는
필드 간 모호한 좌표가 있으면 중단하고, 동일 필드 tract overlap은 가장 가까운 포함 tract
중심으로 정한다. 겹침 및 geometry-file hash를 결과에 기록한다.

tract 단위 삭제는 영향도만을 보여준다. 인접 tract 간 상관이 남으므로 실행기는 이 값으로
공분산을 계산하지 않고 `science_measurement=false`를 유지한다. 네 공식 geometry 파일은
아직 수신하지 않았다. 승인 후 `lageunha`에서 실행할 준비가 됐다.

- 실행기: `analysis/icl_wgl_membership/hsc_debug/analyze_s23b_fixed52_tract_region_influence.py`
- 실행기 SHA-256: `878a5c6471af2acc2cbd0cf478f56d5ab08f825268ad9390acbe228df269ef2e`
- 확인: Python 문법, CLI `--help`, diff whitespace, 2-tract 합성 좌표 대입 검사가 통과했다.
- 실제 HSC geometry 입력이나 CSV 분석은 실행하지 않았다.

## 2026-11-11 HSC exact-tract influence and fixed-52 ICL–satellite–WL audit

The earlier approximate patch-center tract influence record above is superseded
for the fixed-52 sample by an exact-footprint assignment using the four HSC
PDR3 Wide tract/patch geometry tables already present on `/gpfs`. The corrected
parser normalizes unwrapped RA, assigns endpoints through patch polygons to
their parent tract, and records the parser and geometry-file hashes. No new
geometry download, HSC query, image read, or shape-row download was performed.

The exact audit reuses the already approved 20 S19A shape-row CSV batches and
the pinned HSC Y3 SACC `n(z)`. It includes 52 fixed lens sets and 260 planned
matched controls (259 have selected source rows), 138,849 lens pairs and
656,440 control pairs, and 493,402 unique source coordinates. All active
endpoints were assigned; same-field overlaps were resolved by nearest
containing patch center. The analysis covers 102
occupied parent-tract blocks. Every deletion retains the four-control floor
and 48–52 common lens sets. Maximum inner-bin equal-control-center tangential
and cross shifts are 1.727 and 2.705 × 10¹³ M☉ Mpc⁻². These are influence
diagnostics, not independent errors or a survey covariance.

- Runner: `analysis/icl_wgl_membership/hsc_debug/analyze_s23b_fixed52_tract_region_influence.py`
- Geometry parser: `analysis/icl_wgl_membership/hsc_debug/hsc_pdr3_tract_patch_geometry.py`
- Assessment: `analysis/icl_wgl_membership/hsc_debug/wl_s23b_fixed52_exact_tract_influence_assessment_2026-11-11.md`
- Output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_fixed52_exact_tract_influence_2026-11-11_02/`
- Summary SHA-256: `7924dfdbe27d9435099bc0d4990c192e28506cf9a50c9ebf0e8ca584d7859e86`

The saved ICL–satellite profiles were then joined to the fixed-52 WL profiles.
All 52 WL clusters have an ICL profile; `i` and `i2` contain disjoint sets of
25 and 27. Requiring at least 10% common support in both ICL annuli leaves 18
and 22, and at least 20% leaves only 5 and 7. The four direct ICL–satellite
partial-correlation intervals include zero. In the 112 pointwise WL
association tests, eight saved percentile intervals excluded zero. A separate
10,000-replicate restricted wild-cluster max-T analysis of the full 116-test
family found 10 unadjusted p-values below 0.05 and zero familywise-adjusted
values below 0.05; the minimum adjusted p-value is 0.5154. This does not adjust
the earlier percentile intervals. The selected sample, approximate 5-degree
blocks, and unresolved photo-z/shear/PSF/mask/member-dilution/source-reuse
systematics prevent physical detection claims. `science_measurement=false`.

- Overlap summary: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_wl_fixed_support_overlap_2026-11-11_01/` (`8145f969b60fe93075d86bd952c9cd0d0e38f9be71be2ff521252048754cf69e`)
- Relation summary: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_satellite_wl_fixed52_relation_2026-11-11_03/` (`a7a9f30107bb800986d5d1aa4bf4e2012497092170e3bc27bf246587f7840e2f`)
- Max-T summary: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_satellite_wl_fixed52_maxt_2026-11-11_01/` (`e1bbe2282adca533791dd3b4a3a82e048e99d74c4c2e1a82b27a85bb69a807c2`)
- Astrafek reports: `analysis/icl_wgl_membership/hsc_debug/s23b_icl_wl_fixed_support_overlap_assessment_2026-11-11.md`, `s23b_icl_satellite_wl_fixed52_relation_assessment_2026-11-11.md`, and `s23b_icl_satellite_wl_fixed52_max_t_assessment_2026-11-11.md`.

## 2026-11-11 S16A shape coverage count: complete and verified

To assess whether the ICL sample can be expanded beyond the S19A-supported
52-cluster overlap, a local planner prepared a count-only S16A CAS plan for all
1,043 centers in the pinned S23B ICL manifest. It uses 3 proper-Mpc cones,
42 batches of at most 25 centers, and returns only `lens_id` plus one integer
distinct-shape-object count per center (1,043 rows; estimated maximum CSV
size about 31.4 kB). It requests no source IDs, positions, shapes, photo-z,
images, or PSF data. The user approved the exact plan SHA-256 below. The
guarded submitter was run on LagEunha and verified the plan/hash before CAS
access. All 42 batches completed. Independent QA verified each batch result
hash, its expected row count and IDs, and all 1,043 unique IDs against the
pinned target manifest. Counts total 765,197; 177 centers are positive and 866
are zero. For positive centers, the median count is 4,622 (range 4–10,138).
Within the existing fixed 52-lens S19A sample, 21 have positive S16A counts
and 31 have zero. Counts are saved outside Git. No credential is saved.
Source-row export remains a separate scope requiring its own approval.

The follow-up local spatial/property audit found 177 S16A-positive centers
clustered in a few RA/Dec regions. Among the same 427 targets already counted
in S19A, 32 are positive in both releases, 66 are S19A-only, and only 3 are
S16A-only; at the raw-count ≥1,000 screen, 27 are positive in both and 3 are
S16A-only. Another 142 S16A-positive targets (126 with count ≥1,000) are
outside that 427-center S19A audit, so their S19A status is unknown. Within the
fixed 52-lens ICL/WL sample, S16A count support is 19/25 for filter i and 2/27
for i2, a strong spatial/selection imbalance. Counts do not show usable source
support, and no cross-release lensing result was formed.

Among the 126 outside-audit centers with S16A raw counts >=1,000, 116 have an
existing ICL target/control image pair and 10 do not. These controls are unique
(98 i, 18 i2), but they are not the established five-random-centers-per-lens
weak-lensing control sample. A count-only S16A plan for these 116 ICL controls
is prepared but not submitted: 5 batches, 116 integer-count rows, 36,139 SQL
bytes, and a conservative 4,100-byte result bound. Plan SHA-256:
`832f2de18f5d512396cc84673e458932635c1f766ea91cdc04c7b47125738dfa`.
It requires separate approval and would not authorize source-row export.

Among the 126 outside-audit centers with raw S16A counts ≥1,000, 116 have an
existing ICL target/control image pair and 10 do not. Those 116 saved controls
are unique (98 i and 18 i2); their S16A shape counts have not been queried.
They are one-to-one ICL image controls, not the established five-random-centers
per-lens weak-lensing control set.

S16A and S19A are separate processing/calibration releases. Aggregate counts
can identify catalog coverage only; they cannot be interpreted as usable
lensing support and do not authorize mixing S16A rows with the current S19A
sample. The official HSC S16A page documents the `s16a_wide` shape-table joins,
and the FAQ explains its PDR1 CAS schema.

- Planner: `analysis/icl_wgl_membership/hsc_debug/prepare_s23b_s16a_shape_count_plan.py`
- Submitter: `analysis/icl_wgl_membership/hsc_debug/submit_s23b_s16a_shape_counts.py`
- Plan report: `analysis/icl_wgl_membership/hsc_debug/s23b_s16a_shape_count_plan_assessment_2026-11-11.md`
- Spatial/release-selection report: `analysis/icl_wgl_membership/hsc_debug/s23b_s16a_shape_count_selection_assessment_2026-11-11.md`
- Spatial/release-selection runner: `analysis/icl_wgl_membership/hsc_debug/audit_s23b_s16a_shape_count_selection.py`
- ICL-control count planner: `analysis/icl_wgl_membership/hsc_debug/prepare_s23b_s16a_icl_control_count_plan.py`
- Prepared SQL and plan (outside Git): `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_s16a_shape_count_plan_2026-11-11_01/`
- Results (outside Git): `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_s16a_shape_count_plan_2026-11-11_01/results/`
- Combined count TSV SHA-256: `0945329119d1100b783211fc128cd0a329b449f5acf6648313a1484ac55a7b12`
- Count summary JSON SHA-256: `0cf90929c2a69eda0e62f3a9350359fd430b37bf869b2ead53246129ac51bb8b`
- Exact plan SHA-256: `c426129b150d440c22fca97d1cf2c722414a7f708132a6267563d77ccda997e0`
- Spatial/release-selection results: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_s16a_shape_coverage_selection_audit_2026-11-11_04/`
- Selection summary SHA-256: `c08602e1b45492a92389dabd074e5f854050e8e14c2cd55b4683b5ae6fc0b398`
- Candidate ICL control-center TSV SHA-256: `81c702166871e75570918f770ffad44702cebe9a29920c80d3d73d2263c21706`
- Prepared control-count plan (not submitted): `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_s16a_icl_control_count_plan_2026-11-11_01/`
- Control-count plan SHA-256: `832f2de18f5d512396cc84673e458932635c1f766ea91cdc04c7b47125738dfa`
- `science_measurement=false`; aggregate count is a catalog-availability ceiling,
  not a post-cut weak-lensing sample. S16A cannot be merged with S19A without
  separate release/calibration validation.
- Official references: https://hsc-release.mtk.nao.ac.jp/doc/index.php/s16a-shape-catalog-pdr2/ ; https://hsc-release.mtk.nao.ac.jp/doc/index.php/database/summary_table_usage-2/ ; https://hsc-release.mtk.nao.ac.jp/doc/index.php/faq__pdr3/
