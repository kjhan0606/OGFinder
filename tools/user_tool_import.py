"""Standalone personal tools: preflight, registration, and explicit native execution."""
import argparse
import json
from pathlib import Path
import zipfile
from user_tool_validation import validate

FILES = ('docs/user_tool_import.md', 'docs/licensing_rules.md',
         'docs/browser_runtime_licenses.md', 'tools/browser_runtime_lock.json',
         'tools/validate_user_tool.py',
         'examples/user_tool_submission/AI_INSTRUCTIONS.md',
         'examples/user_tool_submission/submission.json',
         'examples/user_tool_submission/adapter.py',
         'examples/user_tool_submission/TEST_PLAN.md')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--review', metavar='DIRECTORY')
    group.add_argument('--kit', metavar='NEW_ZIP')
    group.add_argument('--register', metavar='DIRECTORY')
    group.add_argument('--run', metavar='TOOL_ID')
    group.add_argument('--list', action='store_true')
    parser.add_argument('--registry')
    parser.add_argument('--digest')
    parser.add_argument('--work')
    parser.add_argument('--inputs-json', default='{}')
    parser.add_argument('--params-json', default='{}')
    parser.add_argument('--report-only', action='store_true', help='Return 0 after producing a review JSON, including revision reports (GUI)')
    args = parser.parse_args()
    if args.register or args.run or args.list:
        from user_tool_registry import registry_path, register, run
        registry = registry_path(args.registry)
        if args.register: result = register(args.register, registry)
        elif args.run:
            if not args.work: parser.error('--run requires --work (new output directory)')
            result = run(args.run, registry, args.work, json.loads(args.inputs_json), json.loads(args.params_json), args.digest)
        else:
            result = [json.loads(p.read_text()) for p in sorted(registry.glob('user_*/plugin.json'))] if registry.exists() else []
        print(json.dumps(result, ensure_ascii=False, indent=2))
        return 1 if isinstance(result,dict) and result.get('status') in ('failed', 'cancelled') else 0
    if args.kit:
        base = Path(__file__).resolve().parent / 'user_tool_kit'
        with zipfile.ZipFile(args.kit, 'x', zipfile.ZIP_DEFLATED) as archive:
            for name in FILES:
                if name == 'examples/user_tool_submission/submission.json':
                    manifest = json.loads((base / name).read_text())
                    manifest['target'] = 'standalone'
                    archive.writestr(name, json.dumps(manifest, ensure_ascii=False, indent=2) + '\n')
                else:
                    archive.write(base / name, name)
        print(str(Path(args.kit).resolve()))
        return 0
    report = validate(args.review)
    # CLI reports issues in JSON; exit 1 for automation. GUI reads the JSON even on 1.
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 1 if report['errors'] and not args.report_only else 0


if __name__ == '__main__':
    raise SystemExit(main())
