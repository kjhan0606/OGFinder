"""Atomic, input-bound per-frame checkpoints for local member photometry."""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def atomic_json(path, value):
    fd,tmp=tempfile.mkstemp(prefix=path.name+'.tmp.',dir=path.parent)
    try:
        with os.fdopen(fd,'w') as stream:
            json.dump(value,stream,allow_nan=False)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(tmp,path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


class FrameCheckpoint:
    def __init__(self, directory, config, resume=False):
        self.directory=Path(directory)
        if self.directory.exists() and not resume:
            raise FileExistsError('checkpoint directory exists; use --resume after inspecting the prior run')
        self.directory.mkdir(parents=True,exist_ok=resume)
        self.lock=(self.directory/'run.lock').open('a')
        try:
            fcntl.flock(self.lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            self.config_digest=digest(config)
            path=self.directory/'configuration.json'
            if path.exists():
                saved=json.loads(path.read_text())
                if saved != config:
                    raise ValueError('checkpoint configuration, code or input hashes changed')
            else:
                if resume:
                    raise ValueError('cannot resume a checkpoint directory without its configuration')
                atomic_json(path,config)
        except BaseException:
            self.lock.close()
            raise

    def path(self, target_id, band):
        name=f'{target_id}_{band}'
        if not re.fullmatch(r'[A-Za-z0-9_-]+',name):
            raise ValueError('invalid checkpoint frame ID')
        return self.directory/(name+'.json')

    def load(self,target_id,band):
        path=self.path(target_id,band)
        if not path.exists():
            return None
        wrapper=json.loads(path.read_text())
        payload=wrapper['payload']
        if wrapper.get('configuration_sha256') != self.config_digest or wrapper.get('payload_sha256') != digest(payload):
            raise ValueError('checkpoint payload hash/configuration mismatch')
        frame=payload['frame_record']
        if frame['target_id'] != target_id or frame['filter'] != band:
            raise ValueError('checkpoint frame identity mismatch')
        return payload

    def save(self,target_id,band,payload):
        path=self.path(target_id,band)
        if path.exists():
            raise FileExistsError('refusing to overwrite an existing frame checkpoint')
        atomic_json(path,dict(configuration_sha256=self.config_digest,payload_sha256=digest(payload),payload=payload))
