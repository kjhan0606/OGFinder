"""Native local execution with resource limits; not an OS security sandbox."""
import importlib.util
import json
from pathlib import Path
import sys
import resource


def main():
    r = json.loads(Path(sys.argv[1]).read_text());m=r['manifest']
    resource.setrlimit(resource.RLIMIT_CPU, (min(m['limits']['seconds'],600), min(m['limits']['seconds'],600)+1))
    memory=min(max(m['limits']['memory_mib'],128),8192)*1024*1024
    resource.setrlimit(resource.RLIMIT_AS, (memory,memory))
    size=min(m['limits']['output_mib'],64)*1024*1024
    resource.setrlimit(resource.RLIMIT_FSIZE, (size,size))
    for dep in m['dependencies']:
        if dep['kind']!='python': raise ValueError('Native external binaries require a separate adapter')
        from importlib.metadata import version
        if version(dep['name'])!=dep['version']: raise ValueError('Dependency version mismatch')
    sys.path.insert(0,r['package'])
    spec=importlib.util.spec_from_file_location('personal_adapter',str(Path(r['package'])/'adapter.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    params={**{k:v['default'] for k,v in m['parameters'].items() if isinstance(v,dict) and 'default' in v},**r['params']}
    result=module.run(r['inputs'],params,r['output'])
    names=[a['name'] for a in result['artifacts']]
    (Path(r['output'])/'result-index.json').write_text(json.dumps({'artifacts':names}))


if __name__=='__main__': main()
