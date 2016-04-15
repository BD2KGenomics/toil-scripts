import os


def test_flatten():
    from toil_scripts.lib import flatten
    x = [(1, 2), (3, 4, (5, 6))]
    y = (1, (2, (3, 4, (5))))
    assert flatten(x) == [1, 2, 3, 4, 5, 6]
    assert flatten(y) == [1, 2, 3, 4, 5]


def test_partitions():
    from toil_scripts.lib import partitions
    x = [z for z in xrange(100)]
    assert len(list(partitions(x, 10))) == 10
    assert len(list(partitions(x, 20))) == 5
    assert len(list(partitions(x, 100))) == 1
    assert list(partitions([], 10)) == []


def test_copy_to_output_dir(tmpdir):
    from toil_scripts.lib import copy_to_output_dir
    work_dir = str(tmpdir)
    subdir = os.path.join(work_dir, 'subdir')
    os.mkdir(subdir)
    fpath = os.path.join(work_dir, 'test')
    with open(fpath, 'wb') as fout:
        fout.write(os.urandom(1024))
    copy_to_output_dir(output_dir=subdir, fpaths=[fpath], uuid='uuid')
    assert os.path.exists(os.path.join(subdir, 'uuid.test'))
