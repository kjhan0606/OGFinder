"""Bagpipes adapter (SED fitting at the record's redshift).  See script_adapter.py; real fits: `pip install bagpipes` in the interpreter `python`
(nautilus sampler; no filter curves are shipped by Bagpipes - they are exported from an EAZY FILTER.RES (`eazy_data`/`filters_res`) or read from `filter_dir/<BAND>.dat`)."""
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
