import argparse
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


class UserError(Exception):
    pass


def require(expression, message):
    if not expression:
        raise UserError('\n\n' + message + '\n\n')


def required_length(nmin, nmax):
    """
    For use with argparse's action argument. Allows setting a range for nargs.
    Example: nargs='+', action=required_length(2, 3)

    :param int nmin: Minimum number of arguments
    :param int nmax: Maximum number of arguments
    :return: RequiredLength object
    """
    class RequiredLength(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if not nmin <= len(values) <= nmax:
                msg = 'argument "{f}" requires between {nmin} and {nmax} arguments'.format(
                    f=self.dest, nmin=nmin, nmax=nmax)
                raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)
    return RequiredLength
