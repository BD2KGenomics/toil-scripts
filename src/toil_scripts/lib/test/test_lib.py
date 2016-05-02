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
