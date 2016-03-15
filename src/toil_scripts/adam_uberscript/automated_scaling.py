import sys
import threading
import time
from collections import OrderedDict
from itertools import count

import boto.sdb
from boto.exception import SDBResponseError
from boto.sdb.domain import Domain
from collections import Iterator


class Model(object):
    """
    Base class for models that piggy-back their persistence to a Toil AWS job store.
    """

    @classmethod
    def _get_toil_jobstore_domain(cls):
        aws, region, domain = sys.argv[1].split(':', 3)
        conn = boto.sdb.connect_to_region(region)
        domain = conn.get_domain(domain + '--files')
        return domain


class Sample(object):
    def __init__(self, sample_id, index, nodes):
        super(Sample, self).__init__()
        self.sample_id = sample_id
        self.index = index
        self.nodes = nodes

    @classmethod
    def from_attribute(cls, name, value):
        index, nodes = map(int, value.split(':', 1))
        return cls(sample_id=str(name), index=index, nodes=nodes)

    def to_attribute(self):
        return self.sample_id, '%i:%i' % (self.index, self.nodes)

    def __repr__(self):
        return 'Sample(sample_id={0.sample_id!r}, index={0.index!r}, nodes={0.nodes!r}'.format(self)


class Samples(Model):
    """
    Tracks the number of Toil worker nodes requested per sample and a prioritization ordering among those samples in
    which they are assigned existing nodes.

    >>> sdb = boto.sdb.connect_to_region('us-west-2')
    >>> domain = sdb.create_domain('SamplesTest')
    >>> samples = Samples.load(domain); samples
    Samples(version=0, value=OrderedDict())
    >>> cluster_size = ClusterSize.load(domain)
    >>> cluster_size.value, cluster_size.version
    (0, 0)

    Make the cluster sufficiently large to prevent increase_nodes() from hanging
    >>> cluster_size.value = 8; cluster_size.save()
    >>> samples.increase_nodes('b', 1, domain)
    >>> samples = Samples.load(domain); samples # doctest: +NORMALIZE_WHITESPACE
    Samples(version=1, value=OrderedDict([('b', Sample(sample_id='b', index=1, nodes=1)]))

    >>> samples.increase_nodes('a', 1, domain)
    >>> samples = Samples.load(domain); samples # doctest: +NORMALIZE_WHITESPACE
    Samples(version=2, value=OrderedDict([('b', Sample(sample_id='b', index=1, nodes=1),
                                          ('a', Sample(sample_id='a', index=2, nodes=1)]))

    >>> samples.increase_nodes('c', 1, domain)
    >>> samples = Samples.load(domain); samples # doctest: +NORMALIZE_WHITESPACE
    Samples(version=3, value=OrderedDict([('b', Sample(sample_id='b', index=1, nodes=1),
                                          ('a', Sample(sample_id='a', index=2, nodes=1),
                                          ('c', Sample(sample_id='c', index=3, nodes=1)]))

    >>> samples.increase_nodes('a', 2, domain)
    >>> samples = Samples.load(domain); samples # doctest: +NORMALIZE_WHITESPACE
    Samples(version=4, value=OrderedDict([('b', Sample(sample_id='b', index=1, nodes=1),
                                          ('c', Sample(sample_id='c', index=3, nodes=1),
                                          ('a', Sample(sample_id='a', index=4, nodes=3)]))

    >>> samples.increase_nodes('b', 3, domain)
    >>> samples = Samples.load(domain); samples # doctest: +NORMALIZE_WHITESPACE
    Samples(version=5, value=OrderedDict([('c', Sample(sample_id='c', index=3, nodes=1),
                                          ('a', Sample(sample_id='a', index=4, nodes=3),
                                          ('b', Sample(sample_id='b', index=5, nodes=4)]))

    >>> samples.decrease_nodes('a', 3, domain)
    >>> samples = Samples.load(domain); samples # doctest: +NORMALIZE_WHITESPACE
    Samples(version=6, value=OrderedDict([('c', Sample(sample_id='c', index=3, nodes=1),
                                          ('a', Sample(sample_id='a', index=4, nodes=0),
                                          ('b', Sample(sample_id='b', index=5, nodes=4)]))

    >>> samples.decrease_nodes('b', 4, domain)
    >>> samples = Samples.load(domain); samples # doctest: +NORMALIZE_WHITESPACE
    Samples(version=7, value=OrderedDict([('c', Sample(sample_id='c', index=3, nodes=1),
                                          ('a', Sample(sample_id='a', index=4, nodes=0),
                                          ('b', Sample(sample_id='b', index=5, nodes=0)]))

    >>> domain.delete()
    True
    """

    def __init__(self, domain, version, value):
        self.domain = domain
        self.version = version
        self.value = value

    @classmethod
    def load(cls, domain):
        attributes = domain.get_attributes('nodes_per_sample', consistent_read=True)
        version = int(attributes.pop('version', '0'))
        samples = [Sample.from_attribute(k, v) for k, v in attributes.iteritems()]
        return cls(domain, version, cls._reindex(samples))

    @classmethod
    def _reindex(cls, samples):
        return OrderedDict((sample.sample_id, sample) for sample in
                           sorted(samples, key=lambda sample: (sample.index, sample.sample_id)))

    def save(self):
        attributes = dict((sample.to_attribute() for sample in self.value.itervalues()),
                          version=str(self.version + 1))
        try:
            self.domain.put_attributes('nodes_per_sample', attributes,
                                       expected_value=('version', str(self.version) if self.version else False))
            self.value = self._reindex(self.value.values())
            self.version += 1
        except SDBResponseError as e:
            if e.error_code == 'ConditionalCheckFailed':
                return False
            else:
                raise
        else:
            return True

    @classmethod
    def increase_nodes(cls, sample_id, n, domain=None):
        """
        Increases the node count for the specified sample and returns once the nodes become available. This method
        can ONLY be called from within a job function or method in a Toil script because it assumes that it has
        access to the worker process' command line arguments.
        """
        domain = domain or cls._get_toil_jobstore_domain()

        while True:
            # load nodes per sample from sdb
            self = cls.load(domain)

            # update or insert new nodes
            if sample_id in self.value:
                self.value[sample_id].index = self.version + 1
                self.value[sample_id].nodes += n
            else:
                self.value[sample_id] = Sample(sample_id=sample_id, index=self.version + 1, nodes=n)

            # attempt to write back to sdb; retry if failed
            if self.save():
                break

        while True:
            cluster_size = ClusterSize.load(domain)
            sum = 0
            # noinspection PyUnboundLocalVariable
            for sample in self.value.itervalues():
                sum += sample.nodes
                if sample.sample_id == sample_id and sum <= cluster_size.value:
                    return
                if sum > cluster_size.value:
                    break
            self = cls.load(domain)
            time.sleep(10)

    @classmethod
    def decrease_nodes(cls, sample_id, n, domain=None):
        """
        Decrease the node count for the given sample
        """
        domain = domain or cls._get_toil_jobstore_domain()
        while True:
            self = cls.load(domain)
            sample = self.value[sample_id]
            sample.nodes = min(sample.nodes - n, 0)
            if self.save():
                break

    def __repr__(self):
        return 'Samples(version={0.version!r}, value={0.value!r})'.format(self)


class ConcurrentModificationException(Exception):
    pass


class NoSuchSingletonException(Exception):
    def __init__(self, item_name):
        super(NoSuchSingletonException, self).__init__('Singleton %s does not exist' % item_name)


class SingletonModel(Model):
    default = None
    """
    Subclasses should set this to prevent an exception during load() when the singleton doesn't exist
    """

    @classmethod
    def create(cls, domain, value):
        return cls(domain, value, version=0)

    @classmethod
    def load(cls, domain):
        item_name = cls._item_name()
        attributes = domain.get_attributes(item_name, consistent_read=True)
        try:
            value = attributes['value']
        except KeyError:
            if cls.default is None:
                raise NoSuchSingletonException(item_name)
            else:
                value = cls.default
        else:
            value = str(value)  # suppress unicode
            value = cls._from_string(value)
        return cls(domain, value, version=int(attributes.get('version', '0')))

    @classmethod
    def load_in_toil(cls):
        """
        Same as load() but to be used from within a job function of a Toil user script.
        """
        return cls.load(cls._get_toil_jobstore_domain())

    def __init__(self, domain, value, version):
        self.domain = domain
        self.value = value
        self.version = version

    def save(self, force=False):
        """
        :param force: If False, only save this singleton if the underlying item was not modified since this instance
        was loaded from it and raise ConcurrentModificationException if it was. If True, don't detect concurrent
        modifications and simply write the current value and version.
        """
        try:
            item_name = self._item_name()
            self.domain.put_attributes(
                item_name=item_name,
                expected_value=None if force else ('version', str(self.version) if self.version else False),
                attributes=dict(value=self._from_string(self.value), version=str(self.version + 1)))
        except SDBResponseError as e:
            if e.error_code == 'ConditionalCheckFailed':
                raise ConcurrentModificationException
            else:
                raise

    @classmethod
    def _from_string(cls, value):
        return value

    @classmethod
    def _to_string(cls, value):
        return str(value)

    @classmethod
    def _item_name(cls):
        return cls.__name__


class SingletonIntegerModel(SingletonModel):
    """
    >>> sdb = boto.sdb.connect_to_region('us-west-2')
    >>> domain = sdb.create_domain('SingletonIntegerModelTest')
    >>> SingletonIntegerModel.load(domain).value
    Traceback (most recent call last):
    ...
    NoSuchSingletonException: Singleton SingletonIntegerModel does not exist
    >>> domain.delete()
    True
    """

    @classmethod
    def _from_string(cls, value):
        return int(value)


class ClusterSize(SingletonIntegerModel):
    """
    Caches the actual cluster size so we can avoid listing instances all the time
    """
    default = 0


class SparkMasterAddress(SingletonModel):
    """
    The current address of the external Spark master.

    >>> sdb = boto.sdb.connect_to_region('us-west-2')
    >>> domain = sdb.create_domain('SparkMasterAddressTest')
    >>> a1=SparkMasterAddress.load(domain)
    >>> a1.value
    ''
    >>> a1.value = 'foo'
    >>> a1.save()
    >>> a2 = a1.load(domain)
    >>> a2.value
    'foo'
    >>> a2.value = 'bar'
    >>> a2.save()
    >>> a2 = SparkMasterAddress.load(domain)
    >>> a2.value
    'bar'
    >>> a1.value
    'foo'
    >>> a1.save()
    Traceback (most recent call last):
    ...
    ConcurrentModificationException
    >>> a1.save(force=True)
    >>> a2 = SparkMasterAddress.load(domain)
    >>> a2.value
    'foo'
    >>> domain.delete()
    True
    """
    default = ''


class Semaphore(Model):
    """
    A semaphore implementation that uses optimistic locking of a single item in a SimpleDB domain. Instances of this 
    class are not thread-safe. 
    
    >>> sdb = boto.sdb.connect_to_region('us-west-2')
    >>> domain = sdb.create_domain('SemaphoreTest')
    >>> sem1 = Semaphore.create(domain,'foo',overwrite=True)
    
    >>> sem1.value
    0
    
    >>> for attempt in sem1.release(1): 
    ...     assert False # shouldn't happen since no one else is using the semaphore 
    >>> sem1.value
    1
    
    >>> for attempt, concurrent in sem1.acquire(1): 
    ...     raise # shouldn't happen since no one else is using the semaphore and we know what its value is
    >>> sem1.value
    0
    
    Create another semaphore instance representing the same SDB item and start a thread filling it up.
    
    >>> sem2 = Semaphore.load(domain,'foo')
    >>> def writer():
    ...     for i in range(10):
    ...         time.sleep(1)
    ...         for attempt in sem2.release(1):
    ...             pass # Retry immediately on concurrent updates.
    >>> thread = threading.Thread(target=writer)
    >>> thread.start()
    >>> for i in range(10):
    ...     for attempt, concurrent in sem1.acquire(1):
    ...         if not concurrent: # Retry immediately on concurrent updates ...
    ...             time.sleep(1)  # ... but sleep when semaphore value is too small.
    
    The writer's instance must reflect the most recent release().
         
    >>> thread.join()
    >>> sem2.value > 0
    True
    
    The reader's instance must reflect the most recent acquire().
    
    >>> sem1.value
    0

    >>> domain.delete()
    True
    """

    @classmethod
    def create(cls, domain, name, value=0, overwrite=False):
        """
        Create a new semaphore of the specified name in the given SDB domain.
         
        :param Domain domain: the SimpleDB domain in which to creating the item to use for the semaphore 
        :param str name: the name of the semaphore 
        :param int value: the initial value of the semaphore
        :param bool overwrite: True if an existing semaphore of the same value should be overwritten
        :return: an instance of this class representing the initial value of the semaphore 
        """
        self = cls(domain, name, value=value, version=0)
        self._save(expected_value=None if overwrite else ['version', False])
        return self

    @classmethod
    def load(cls, domain, name):
        """
        Return an instance of this class representing the current value of the semaphore with the specified name in the 
        given domain. 
        
        :param Domain domain: the SimpleDB domain in which to creating the item to use for the semaphore
        :param str name: the name of the semaphore
        :return: an instance of this class representing the current value of the semaphore
        """
        self = cls(domain, name)
        self._load()
        return self

    @classmethod
    def load_in_toil(cls, name):
        """
        Same as load() but to be used from within a job function of a Toil script.
        """
        return cls.load(cls._get_toil_jobstore_domain(), name)

    def acquire(self, delta=1):
        """
        Repeatedly attempt to decrease the semaphore by the given amount until successful. 
        
        Yield a a tuple ``(attempt, concurrent)`` for every failed attempt where ``attempt`` is a 0-based index of 
        the current attempt and ``concurrent`` is True if the attempt failed because the underlying SDB item was 
        updated concurrently, or False if the attempt failed because the semaphore value isn't large enough to 
        acquire the given delta.

        :param int delta: the amount to decrease the semaphore by 
        :rtype: Iterator[(int,bool)] 

        Using acquire() in a non-blocking fashion:
        
        for (attempt, concurrent) in sem.acquire():
           reurn False
        return True
         
        Delaying retries (recommended):
         
        for (attempt, concurrent) in sem.acquire():
                time.sleep(1 if concurrent else 10)

        """
        for attempt in count():
            if self.value < delta:
                yield attempt, False
                self._load()
            else:
                self.value -= delta
                self.version += 1
                if self._update():
                    break
                else:
                    yield attempt, True

    def release(self, delta=1):
        """
        Repeatedly attempt to increase the semaphore by the given amount until successful. 
        
        For every attempt that fails due to a concurrent update to the underlying SDB item, yield the 0-based index 
        of the current attempt. 

        :param int delta: the amount to decrease the semaphore by 
        :rtype: Iterator[(int,bool)] 
        """
        for attempt in count():
            self.value += delta
            self.version += 1
            if self._update():
                break
            yield attempt

    def __init__(self, domain, name, value=None, version=None):
        self.domain = domain
        self.name = name
        self.value = value
        self.version = version

    def _update(self):
        try:
            self._save(expected_value=['version', self.version - 1])
        except SDBResponseError as e:
            if e.error_code == 'ConditionalCheckFailed':
                self._load()
                return False
            else:
                raise
        else:
            return True

    def _load(self):
        attributes = self.domain.get_attributes(item_name=self._item_name,
                                                consistent_read=True)
        self.version = int(attributes['version'])
        self.value = int(attributes['value'])

    def _save(self, **put_attributes_kwargs):
        attributes = dict(version=str(self.version),
                          value=str(self.value))
        self.domain.put_attributes(item_name=self._item_name,
                                   attributes=attributes,
                                   **put_attributes_kwargs)

    @property
    def _item_name(self):
        return 'semaphore_' + self.name
