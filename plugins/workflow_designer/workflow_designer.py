#!/usr/bin/env python3
"""Small JSON bridge for the guided workflow window; no science is run here."""
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from ogfkit.science_workflow import checkpoint, recipes, restore, validate_catalog


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    commands.add_parser('catalog')
    save = commands.add_parser('save')
    save.add_argument('--request', required=True)
    save.add_argument('--destination', required=True)
    load = commands.add_parser('restore')
    load.add_argument('--manifest', required=True)
    load.add_argument('--work')
    args = parser.parse_args()
    if args.command == 'catalog':
        items = recipes()
        validate_catalog(items, ROOT)
        answer = dict(recipes=items)
    elif args.command == 'save':
        request = json.loads(Path(args.request).read_text(encoding='utf-8'))
        state = request['state']
        definitions = {r['id']: r for r in recipes()}
        state['recipe_definition'] = definitions[state['recipe']]
        manifest = checkpoint(args.destination, state, request['outputs'], request['inputs'])
        answer = dict(manifest=str(manifest))
    else:
        answer = restore(args.manifest, args.work)
    print(json.dumps(answer, ensure_ascii=False))


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, KeyError) as exc:
        print(f'workflow: {exc}', file=sys.stderr)
        sys.exit(1)
