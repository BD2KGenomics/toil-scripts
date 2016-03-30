import os
import subprocess
import shutil


def test_rnaseq_cgl(tmpdir):
    work_dir = str(tmpdir)
    create_config(work_dir)
    subdir = '/mnt/ephemeral/toil-scripts/rnaseq/'
    os.makedirs(os.path.join(subdir, 'workDir'))
    # URLs for sample inputs
    star_index = 'https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/rnaseq_cgl/ci/starIndex_chr6.tar.gz'
    rsem_ref = 'https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/rnaseq_cgl/ci/rsem_ref_chr6.tar.gz'
    # Call Pipeline
    try:
        subprocess.check_call(['python', '-m', 'toil_scripts.rnaseq_cgl.rnaseq_cgl_pipeline',
                               os.path.join(subdir, 'jstore'),
                               '--config', os.path.join(work_dir, 'config.txt'),
                               '--retryCount', '1',
                               '--s3_dir', 'cgl-driver-projects/test/releases',
                               '--workDir', os.path.join(subdir, 'workDir'),
                               '--starIndex', star_index,
                               '--rsemRef', rsem_ref,
                               '--ci-test'])
    finally:
        shutil.rmtree(subdir)


def create_config(path):
    """Creates config file for test at path"""
    fpath = os.path.join(path, 'config.txt')
    with open(fpath, 'w') as f:
        f.write('ci_test,'
                'https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/rnaseq_cgl/ci/test_sample_chr6.tar.gz')