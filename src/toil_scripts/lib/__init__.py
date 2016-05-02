# FIXME: replace with bd2k.util.iterables.flatten
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
    for unit in ['', 'K', 'M', 'G', 'T', 'P', 'E', 'Z']:
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


class UserError(Exception):
    pass


def require(expression, message):
    if not expression:
        raise UserError('\n\n' + message + '\n\n')
