# FIXME: replace with bd2k.util.iterables.flatten
import os
import shutil


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


# FIXME: replace with bd2k.util humanization function
def sizeof_fmt(num, suffix='B'):
    """
    Humanize input
    >>> sizeof_fmt(1024)
    '1.0KB'
    >>> sizeof_fmt(1024**2)
    '1.0MB'
    >>> sizeof_fmt(1024**3)
    '1.0GB'
    >>> sizeof_fmt(1024**4)
    '1.0TB'

    :param int num: number to be converted to human-readable string
    :param str suffix: Size suffix; defaults to B for byte
    :return: Human-readable string
    :rtype: str
    """
    for unit in ['','K','M','G','T','P','E','Z']:
        if abs(num) < 1024.0:
            return "%3.1f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, 'Y', suffix)


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


def copy_to_output_dir(output_dir, fpaths, uuid=None):
    """
    A list of files to move from work_dir to output_dir.

    :param str output_dir: Directory to move files to
    :param list fpaths: List of filepaths to move
    :param str uuid: Optional UUID to tag files with
    """
    for fpath in fpaths:
        if uuid is None:
            shutil.copy(fpath, os.path.join(output_dir, os.path.basename(fpath)))
        else:
            shutil.copy(fpath, os.path.join(output_dir, '{}.{}'.format(uuid, os.path.basename(fpath))))
