import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from member_profile_checkpoint import FrameCheckpoint


class CheckpointTests(unittest.TestCase):
    def test_resume_preserves_payload_and_rejects_changed_inputs(self):
        with TemporaryDirectory() as tmp:
            path=Path(tmp)/'frames'
            c=FrameCheckpoint(path,{'input':'hash1'})
            data=dict(frame_record={'target_id':'A','filter':'i'},profiles=[{'flux':2.5}])
            c.save('A','i',data)
            self.assertEqual(c.load('A','i'),data)
            with self.assertRaises(FileExistsError): c.save('A','i',data)
            c.lock.close()
            d=FrameCheckpoint(path,{'input':'hash1'},resume=True)
            self.assertEqual(d.load('A','i'),data)
            d.lock.close()
            with self.assertRaisesRegex(ValueError,'changed'):
                FrameCheckpoint(path,{'input':'different'},resume=True)

    def test_corrupt_payload_is_rejected(self):
        with TemporaryDirectory() as tmp:
            c=FrameCheckpoint(Path(tmp)/'frames',{'input':'hash1'})
            data=dict(frame_record={'target_id':'A','filter':'i'},profiles=[])
            c.save('A','i',data)
            path=c.path('A','i')
            wrapper=json.loads(path.read_text()); wrapper['payload']['profiles']=[999]
            path.write_text(json.dumps(wrapper))
            with self.assertRaisesRegex(ValueError,'mismatch'): c.load('A','i')
            c.lock.close()

    def test_explicit_resume_and_exclusive_writer_are_required(self):
        with TemporaryDirectory() as tmp:
            path=Path(tmp)/'frames'
            c=FrameCheckpoint(path,{'input':'hash1'})
            with self.assertRaises(FileExistsError): FrameCheckpoint(path,{'input':'hash1'})
            with self.assertRaises(BlockingIOError): FrameCheckpoint(path,{'input':'hash1'},resume=True)
            with self.assertRaises(ValueError): c.path('../outside','i')
            c.lock.close()


if __name__=='__main__': unittest.main()
