# 사용자 AI 도구 Import 사전검사

Tools 메뉴에서 `User tool: AI submission kit...`로 사용자 AI에게 전달할 ZIP을
저장한다. 수정한 제출 폴더는 `User tool: local Import review...`로 검사하고
JSON 보고서를 새 파일로 저장한다. 기존 ZIP/보고서를 덮어쓰지 않는다.

Python 3.10 이상 필요. Web 서버나 인접 ds10-web checkout은 필요 없다.
메뉴는 Tcl 소스 변경을 반영한 실행본에서 사용할 수 있으며, 기존 바이너리는
아래 CLI를 바로 사용할 수 있다.

```bash
python3 tools/user_tool_import.py --kit /tmp/my-new-ai-kit.zip
python3 tools/user_tool_import.py --review /path/to/submission > /tmp/review.json
```

사전검사는 제출 코드를 import/실행/설치하지 않는다. 성공해도 실행 허가가
아니며 담당 AI의 의존성·라이선스 검토, 격리 시험, 수치 검증 및 제품 연결이
필요하다. 실제 분석 실행 호스트나 자동 plugin 등록은 이 기능에 포함되지 않는다.
웹에서는 이 기능만 브라우저에서 처리하고 전체 서버/로컬 모드와 다른 분석의 실행 위치를 유지한다. Web 서버에 사용자 코드를 등록하지 않는다.

규약 원본은 `tools/user_tool_kit/docs/user_tool_import.md`, 예제는
`tools/user_tool_kit/examples/user_tool_submission/`이다. 검증기 원본은
ds10-web의 `tools/validate_user_tool.py`이며 여기의 `tools/user_tool_validation.py`는
동일한 first-party 사본이다. 규약/검증기 변경 시 kit와 사본을 함께 갱신하고
`python3 regression/test_user_tool_import.py`를 실행한다. Web 코드는 OGFinder를 import하지 않는다.
