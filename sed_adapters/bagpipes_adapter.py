"""SED-fit adapter. ``command`` runs a caller-supplied program with the JSON protocol in script_adapter.py.

The shipped script calls the user-installed package only when that import works. A missing package exits without a bagpipes-labeled fit. Filter files are still exported when a caller command needs them. The in-tree SED model is analytic."""
import sys

from . import common, script_adapter

CODE = 'bagpipes'
TASKS = ('sed_fit',)


def process(records, params=None, task='sed_fit'):
    return script_adapter.process_script(CODE, 'bagpipes_fit.py', records, params, task, 'Bagpipes (%s)', True)


def check(params=None):
    return script_adapter.check(CODE, 'bagpipes', params)


def run(records, params, context=None):
    return process(records, params, (context or {}).get('task', 'sed_fit'))[0]


def main():
    return common.adapter_main(process)


if __name__ == '__main__':
    sys.exit(main())
