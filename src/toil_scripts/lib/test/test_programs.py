import os


def test_which():
    from toil_scripts.lib.programs import which
    assert which('python')


def test_docker_call(tmpdir):
    from toil_scripts.lib.programs import docker_call
    work_dir = str(tmpdir)
    parameter = ['--help']
    tool = 'quay.io/ucsc_cgl/samtools'
    docker_call(work_dir=work_dir, parameters=parameter, tool=tool)
    # Test outfile
    fpath = os.path.join(work_dir, 'test')
    with open(fpath, 'w') as f:
        docker_call(tool='ubuntu', env=dict(foo='bar'), parameters=['printenv', 'foo'], outfile=f)
    assert open(fpath).read() == 'bar\n'
