"""Prospector adapter (SED fitting at the record's redshift; python-fsps + sedpy + dynesty/emcee in the interpreter `python`, SPS_HOME must point to the FSPS data).
See script_adapter.py: the shipped native/prospector_fit.py runs a parametric (delay-tau) model with dynesty; `command` replaces it."""
import sys

from . import common, script_adapter

CODE = 'prospector'
TASKS = ('sed_fit',)


def process(records, params=None, task='sed_fit'):
    return script_adapter.process_script(CODE, 'prospector_fit.py', records, params, task, 'Prospector (%s)', False)


def check(params=None):
    return script_adapter.check(CODE, 'prospect', params)


def run(records, params, context=None):
    return process(records, params, (context or {}).get('task', 'sed_fit'))[0]


def main():
    return common.adapter_main(process)


if __name__ == '__main__':
    sys.exit(main())
