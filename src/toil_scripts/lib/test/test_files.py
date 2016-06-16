import os
import tarfile
from toil.job import Job


def test_tarball_files(tmpdir):
    from toil_scripts.lib.files import tarball_files
    work_dir = str(tmpdir)
    fpath = os.path.join(work_dir, 'output_file')
    with open(fpath, 'wb') as fout:
        fout.write(os.urandom(1024))
    tarball_files(output_dir=work_dir, tar_name='test.tar', file_paths=[fpath])
    assert os.path.exists(os.path.join(work_dir, 'test.tar'))


def test_move_files(tmpdir):
    from toil_scripts.lib.files import move_files
    work_dir = str(tmpdir)
    os.mkdir(os.path.join(work_dir, 'test'))
    fpath = os.path.join(work_dir, 'output_file')
    with open(fpath, 'wb') as fout:
        fout.write(os.urandom(1024))
    move_files([fpath], os.path.join(work_dir, 'test'))
    assert os.path.exists(os.path.join(work_dir, 'test', 'output_file'))


def test_consolidate_tarballs_job(tmpdir):
    options = Job.Runner.getDefaultOptions(os.path.join(str(tmpdir), 'test_store'))
    Job.Runner.startToil(Job.wrapJobFn(_consolidate_tarball_job_setup), options)


def _consolidate_tarball_job_setup(job):
    from toil_scripts.lib.files import consolidate_tarballs_job
    # Create test file
    work_dir = job.fileStore.getLocalTempDir()
    fpath = os.path.join(work_dir, 'output_file')
    with open(fpath, 'wb') as fout:
        fout.write(os.urandom(1024))
    # Create test tarballs
    fpath1 = os.path.join(work_dir, 'test1.tar.gz')
    fpath2 = os.path.join(work_dir, 'test2.tar.gz')
    with tarfile.open(fpath1, 'w:gz') as f_out:
        f_out.add(fpath)
    with tarfile.open(fpath2, 'w:gz') as f_out:
        f_out.add(fpath)
    id1 = job.fileStore.writeGlobalFile(fpath1)
    id2 = job.fileStore.writeGlobalFile(fpath2)
    job.addChildJobFn(consolidate_tarballs_job, dict(test1=id1, test2=id2))
