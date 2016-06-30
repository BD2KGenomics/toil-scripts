import os
import tempfile


def flatten(x):
    """
    Flattens a nested array into a single list

    :param list x: The nested list/tuple to be flattened.
    """
    result = []
    for el in x:
        if hasattr(el, "__iter__") and not isinstance(el, basestring):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result


def partitions(l, partition_size):
    """
    >>> list(partitions([], 10))
    []
    >>> list(partitions([1,2,3,4,5], 1))
    [[1], [2], [3], [4], [5]]
    >>> list(partitions([1,2,3,4,5], 2))
    [[1, 2], [3, 4], [5]]
    >>> list(partitions([1,2,3,4,5], 5))
    [[1, 2, 3, 4, 5]]

    :param list l: List to be partitioned
    :param int partition_size: Size of partitions
    """
    for i in xrange(0, len(l), partition_size):
        yield l[i:i + partition_size]


def get_work_directory():
    """
    Creates a a working directory at a location that is suitable for storing files
    temporarily. The caller is responsible for removing the directory.

    :return: The path to directory
    :rtype: str
    """
    try:
        basedir = os.environ['TOIL_WORKDIR']
    except KeyError:
        basedir = '/mnt/ephemeral'
        if not os.path.exists(basedir):
            basedir = None
    return tempfile.mkdtemp(dir=basedir)


class UserError(Exception):
    pass


def require(expression, message):
    if not expression:
        raise UserError('\n\n' + message + '\n\n')
