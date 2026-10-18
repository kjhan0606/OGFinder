# 사용자 AI 분석도구 제출 및 Import 규약 v1

이 규약은 사용자 AI가 도구를 작성·검토하는 제출 형식이다. 정적 검사 통과는 AI 승인이나 라이선스 보증이 아니다. 사용자는 검토한 도구를 **내 브라우저 도구로 등록**한 뒤 직접 실행할 수 있다. 프로젝트 공식 도구 편입과 개인 도구 등록은 구분하며, 자동 패키지 설치는 하지 않는다.

웹에서는 PIPELINE 아래 **사용자 도구 Import · 로컬 검증**을 선택한다. 전체 서버/로컬 모드, 저장된 모드 선택, 다른 분석의 실행 위치를 바꾸지 않는다. 제출 파일은 서버에 보내지 않는다. AI 전달 ZIP을 받고 수정된 `submission.json` 포함 폴더를 선택하면 브라우저에서 사전검사 후 JSON 보고서를 저장할 수 있다. 제출 폴더는 하위 폴더 한 단계까지, 파일별 2 MiB, 명세 포함 101개·18 MiB 이내다.

검사 후 **내 브라우저 도구로 등록**을 누르면 해당 도구 카드, 입력 파일 선택, 설정 JSON, **실행 · 웹 로컬**, 로그 및 산출물 다운로드가 생성된다. 코드와 최근 실행 결과는 사이트·계정별 IndexedDB에 저장되며 새로고침 후 복원된다. 동일 id 등록은 이전 버전을 대체하므로 먼저 원본 백업을 권장한다. 도구 원본 JSON 백업을 복원할 때도 사전검사를 다시 거친다. 원본 백업에는 실행 입력/최근 결과가 포함되지 않으며 산출물은 별도로 다운로드한다. 개인 도구는 일반 서버 job 또는 공식 plugin registry에 등록하지 않는다.

실행은 `sandbox=allow-scripts` (allow-same-origin 없음) iframe의 Worker 안에서 수행한다. 엔진 파일은 부모 페이지가 고정 경로·해시로 확인하여 전달한다. 격리 문서의 CSP는 HTTP/WebSocket 등 연결을 허용하지 않으며 blob 엔진 자료만 허용한다. 사용자 코드에는 앱 쿠키·IndexedDB·다른 도구·서버 파일 접근을 제공하지 않는다. 단순 AST 검사를 sandbox로 취급하지 않는다.

현재 개인 도구는 Python 표준 라이브러리, 제출 폴더의 자체 Python 모듈, 고정 NumPy/SciPy 런타임을 사용할 수 있다. 다른 wheel/native 실행파일은 설치하지 않으며 실행 시 지원하지 않는 의존성을 명시적으로 거부한다. 입력 파일 합계 128 MiB, 반환 산출물 64 MiB, 로그 1 MiB, 엔진 준비 60초, 실행은 manifest의 seconds(최대 600초)다. 실행 취소/초과시간은 iframe과 Worker를 종료한다. 메모리 목표는 브라우저의 강제 RAM 상한이 아니며 큰 배열은 브라우저 메모리를 소진할 수 있다.

정적 배포는 생성된 `user-tool-sandbox.html`, `user-tool-assets.json`, `user-tool-runtime/`도 포함한다. sandbox 문서는 동일 사이트에서 frame으로 열 수 있어야 하며, 모든 페이지에 일괄 `X-Frame-Options: DENY` 또는 부모와 동일한 CSP를 덮어씌우면 안 된다. 내장 Web 서버는 해당 파일에만 별도 CSP/프레임 정책을 적용한다. 다른 페이지의 기존 정책은 유지한다.

Standalone은 Tools 메뉴의 **User tool: AI submission kit...**, **User tool: local Import review...**를 이용한다. 별도 Web 서버가 필요 없고 Python 3.10 이상으로 로컬에서 검사한다. CLI는 `python3 tools/user_tool_import.py --review <제출폴더>` 및 `--kit <새 ZIP 경로>`다. 메뉴 변경은 Tcl 소스가 반영된 빌드/실행본에서 사용한다. 로컬 검사는 sandbox 실행 시험이 아니며 코드가 실행되지 않는다.

## 사용자 → 사용자 AI → Astrafex 검증 담당 AI

1. `examples/user_tool_submission/` 폴더 전체와 이 문서를 사용자 AI에게 전달한다. 폴더의 `AI_INSTRUCTIONS.md`에 연구 목적·입력·출력·원하는 실행 위치를 적는다.
2. 사용자 AI는 `submission.json`, `adapter.py`, `TEST_PLAN.md`를 수정하고, 필요한 소스·합성 fixture·기대값 파일을 추가한다. 제3자 실행파일이나 wheel은 동봉하지 않고 출처·정확한 버전·배포물 해시를 명세에 기록한다.
3. 모든 제출 파일을 `files`에 기록하고 변경 후 SHA-256을 갱신한다. 검증기는 명시된 파일만 읽는다. 누락된 파일은 등록 대상이 아니며 검토자가 명세 외 코드를 실행해서는 안 된다.
4. 수정된 폴더와 정적 검사 보고서를 Astrafex 검증 담당 AI에게 전달한다. 계정, 비밀번호, 토큰, 실제 관측자료는 넣지 않는다.
5. 담당 AI는 아래 검증 단계를 수행한다. 결과는 `rejected`, `needs_revision`, `browser_verified`, `standalone_verified`로 구분한다. 마지막 두 상태도 검증한 버전·입력 범위만 의미한다. Web 서버에는 사용자 제출 코드를 등록하지 않는다.

## 제출 명세

`submission.json`은 UTF-8 JSON이다. 예제의 필드는 모두 유지한다.

| 필드 | 계약 |
|---|---|
| `schema`, `id`, `version` | `astrafex-user-tool/1`, `user_` 접두 식별자, 도구 버전 |
| `target` | `browser`, `standalone`, `both` 중 요청 위치. both는 브라우저와 standalone. 기존 server 값은 읽되 대상 수정 필요로 보고하며 Web 서버에 설치하지 않음 |
| `entrypoint` | v1은 `adapter.py:run` 고정 |
| `purpose`, `provenance` | 과학 목적, 작성자·코드 출처·사용권·한계. AI 생성도 기재 |
| `inputs`, `parameters`, `outputs` | 이름, 형식, 단위, 필수 여부/기본값, 유효 범위, NaN·결측 처리 등 |
| `dependencies` | 직접·전이 의존성 모두. `name`, 정확한 `version`, `imports`(최상위 모듈명), `source`, `license`, `license_source`, `sha256`, `kind`(`python`/`wasm`/`executable`) |
| `permissions` | 네트워크, subprocess, 입력 외 읽기, 출력 외 쓰기 요구를 각각 명시 |
| `limits` | 시간(초), 메모리(MiB), 출력(MiB)의 시험 목표 |
| `files` | 상대경로 → 파일 바이트 SHA-256. `submission.json` 자신은 제외 |
| `test_plan` | 명세에 포함된 시험 계획 파일 경로 |

URL은 공급자 원문과 해당 버전 배포물의 출처를 가리킨다. `latest`, 버전 범위, 해시 없는 의존성은 허용하지 않는다. Python 표준 라이브러리는 의존성 목록에서 제외하되 검증 보고서에 import 목록을 남긴다. 동일 이름의 로컬 파일로 표준/외부 모듈을 가리지 않는다. 설치된 런타임에 이미 있는 NumPy/SciPy도 사용한다면 정확한 버전과 배포물 출처를 기록한다.

## Python 및 데이터 계약

```python
def run(inputs: dict, params: dict, output_dir: str) -> dict:
    # inputs: 논리 이름 -> 검증 담당 호스트가 제공한 입력 파일 경로
    # params: 명세에 정의한 JSON 값
    # output_dir: 호스트가 새로 만든 해당 실행 전용 출력 폴더
    return {"artifacts": [{"name": "result.json", "kind": "json"}],
            "summary": {"science_validated": False}}
```

- 입력은 읽기 전용이다. 산출물은 `output_dir` 바로 아래 안전한 파일명으로만 쓴다. 절대경로·`..`·숨김파일·심볼릭 링크 산출물은 금지한다. JSON에는 NaN/Infinity 대신 `null`과 결측 사유를 쓴다.
- 모듈 import 시 계산·다운로드·설치·파일 쓰기 등 부작용이 없어야 한다. 실행은 `run` 내부에서만 한다. 실패는 예외로 전달하고, 실패 산출물을 완료 결과로 표시하지 않는다.
- `pip`, `micropip`, shell, `eval`/`exec`, 숨겨진 동적 import, 환경변수의 자격증명 접근을 넣지 않는다. 필요한 외부 프로그램은 의존성과 권한에 선언하고 standalone 대상 검토를 요청한다.
- 단위, 좌표계, 픽셀 원점, 영점, 마스크/variance 의미를 명시한다. 난수 seed를 인자로 받고 작은 합성 기준자료로 기대값·허용 오차를 제시한다. 실제 자료의 과학적 타당성을 합성 시험만으로 주장하지 않는다.
- v1 제출 adapter는 제품의 plugin manifest 자체가 아니다. 검증 후 담당 AI가 현재 `ds10core` CLI/manifest 및 브라우저 allowlist에 연결하는 adapter를 구현하고, 파이프라인·산출물·체크포인트까지 시험한다.

## 현재 브라우저 기준

기준 버전은 `tools/browser_runtime_lock.json`이 권위 원본이다. 현재 Pyodide 0.27.7 / Python 3.12.7 / NumPy 2.0.2 / SciPy 1.14.1이며 CPU Web Worker에서 실행한다. 임의의 Linux 바이너리·C 확장 wheel은 그대로 import할 수 없다. 브라우저용 배포물과 전이 의존성의 빌드·라이선스 검토가 필요하다. Worker 자체를 악성 코드 격리 경계로 취급하지 않는다.

현재 제품 한도는 FITS 128 MiB, TSV 32 MiB, 산출물/체크포인트 256 MiB, 계산 600초다. 여러 영상·부속 맵·영역 마스크는 아직 연결되지 않았다. 사용자 AI가 이 제한을 숨기거나 서버로 자동 전송하는 우회를 넣어서는 안 된다. 완전 오프라인을 가정하지 않는다.

Web 라이선스 판정에는 동봉할 `docs/licensing_rules.md`, `docs/browser_runtime_licenses.md`가 우선한다. Standalone 대상은 별도로 해당 프로젝트의 배포·의존성 라이선스 정책을 검토한다. 실행 위치 변경이 라이선스 예외를 뜻하지 않는다. 알려지지 않은 사용권은 보류하며 allowlist나 보안 정책을 제출 AI가 수정하지 않는다.

## 담당 AI의 Import 검증 절차

1. **정적 사전검사:** 신뢰하는 저장소의 `tools/validate_user_tool.py`로 파일 해시·경로·Python 구문·import 선언을 검사한다. 검사 중 제출 Python을 import하거나 설치하지 않는다. 보고서의 `accepted_for_execution`은 항상 `false`다.
2. **독립 검토:** 명세 밖 의존성, import 시 부작용, 동적 호출, 파일/네트워크 접근, 리소스 사용, 전체 전이 의존성·라이선스 원문·배포물 해시를 확인한다. AST 검사만으로 안전하다고 판정하지 않는다.
3. **격리 시험:** 검토 전 코드를 제품 Worker나 운영 서버에서 실행하지 않는다. 별도 일회용 sandbox에서 자격증명 없이, 네트워크 차단, 읽기 전용 입력, 전용 출력, CPU/메모리/시간 제한을 실제로 적용한다. 프로세스 분리나 Python venv만으로 격리됐다고 주장하지 않는다. 적절한 격리가 없으면 실행 단계는 보류한다. HPC 로그인 노드에서는 계산하지 않는다.
4. **수치·실패 시험:** TEST_PLAN의 기대값, 허용 오차, 잘못된 입력, 취소·초과시간·메모리 실패, 입력 불변성을 검증한다. 브라우저 대상은 실제 Chromium/WASM 실행 및 native 비교를 수행하고 의도하지 않은 API/네트워크 요청이 없는지 확인한다.
5. **제품 연결:** 개인 도구 카드는 입력·설정·실행·로그·결과 저장을 제공한다. 공식 과학 레시피에 편입하려면 검증된 adapter의 앞뒤 단계 데이터 연결과 전체 체크포인트까지 별도로 검증한다. 개인 도구 등록만으로 공식 레시피의 중간 단계나 서버 plugin에 자동 편입하지 않는다.
6. **보고:** 도구/파일/의존성 해시, 검증 환경, 실행 명령·결과, 수치 비교, 미검증 항목, 지원 모드를 기록한다. 코드·의존성·권한 변경 시 재검증한다. 실패를 숨기거나 검증 기준을 완화해 통과시키지 않는다.

검사 예시(저장소 루트에서, 결과는 stdout):

```bash
python3 tools/validate_user_tool.py examples/user_tool_submission > /tmp/user-tool-review.json
```

종료 코드 0은 정적 사전검사 통과, 1은 수정 필요다. 이 명령은 패키지를 다운로드·설치·실행하거나 제품에 등록하지 않는다.
