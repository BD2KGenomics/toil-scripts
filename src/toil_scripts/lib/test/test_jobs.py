import tempfile

import os
from toil.job import Job


def test_map_job():
    from toil_scripts.lib.jobs import map_job
    work_dir = tempfile.mkdtemp()
    options = Job.Runner.getDefaultOptions(os.path.join(work_dir, 'test_store'))
    options.workDir = work_dir
    samples = [x for x in xrange(200)]
    j = Job.wrapJobFn(map_job, _test_batch, samples, 'a', 'b', 'c', disk='1K')
    Job.Runner.startToil(j, options)


def _test_batch(job, sample, a, b, c):
    assert str(sample).isdigit()
    assert a == 'a'
    assert b == 'b'
    assert c == 'c'
