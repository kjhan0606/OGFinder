"""EXAMPLE python_callable service (not a model): module:function(records, params, context) -> list of dicts."""


def score(records, params, context):
    scale = float(params.get('scale', 1.0))
    return [{'id': r['id'], 'REALBOGUS_SCORE': min(1.0, scale * (r['x'] or 0) / 1000.0)} for r in records]
