"""Original contract example: finite JSON values to a mean; no import side effects."""
import json
import math
from pathlib import Path


def run(inputs: dict, params: dict, output_dir: str) -> dict:
    values = json.loads(Path(inputs['values']).read_text(encoding='utf-8'))
    if not isinstance(values, list) or not values or any(
        type(v) not in (int, float) or not math.isfinite(v) for v in values
    ):
        raise ValueError('Expected a nonempty list of finite numbers')
    mean = math.fsum(v / len(values) for v in values)
    Path(output_dir, 'result.json').write_text(
        json.dumps({'count': len(values), 'mean': mean}, allow_nan=False), encoding='utf-8')
    return {'artifacts': [{'name': 'result.json', 'kind': 'json'}],
            'summary': {'science_validated': False}}
