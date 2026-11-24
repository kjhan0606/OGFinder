import os
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch
from run_s23b_member_coverage_after_diagnostics import runner_identity, wait_runner, run


class RunnerWaitTests(unittest.TestCase):
    def test_live_process_has_start_identity(self):
        self.assertTrue(runner_identity(os.getpid()).isdigit())

    def test_pid_reuse_does_not_extend_wait(self):
        with patch('run_s23b_member_coverage_after_diagnostics.runner_identity',side_effect=['old','new']), patch('run_s23b_member_coverage_after_diagnostics.time.sleep') as sleep:
            wait_runner(100,'old')
            sleep.assert_called_once_with(5)

    def test_missing_process_ends_wait(self):
        with patch('run_s23b_member_coverage_after_diagnostics.runner_identity',return_value=None), patch('run_s23b_member_coverage_after_diagnostics.time.sleep') as sleep:
            wait_runner(100,'old')
            sleep.assert_not_called()

    def test_wrong_runner_rejected_before_creating_work(self):
        with TemporaryDirectory() as tmp:
            root=Path(tmp)
            with self.assertRaisesRegex(ValueError,'does not identify'):
                run(root,os.getpid(),100)
            self.assertFalse((root/'coverage_followup').exists())


if __name__=='__main__': unittest.main()
