import os
from toil.job import Job
from toil_scripts.lib import get_work_directory


def test_map_job(tmpdir):
    from toil_scripts.lib.jobs import map_job
    work_dir = get_work_directory()
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
