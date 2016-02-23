#!/usr/bin/env python
#
#   Copyright (C) 2015 Diamond Light Source, Markus Gerstel
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# A simple cache controller. Caching only one file at a time.

from __future__ import division
import dxtbx.filecache
import os
import threading

class simple_controller():
  '''A simple cache controller. Caching one file at a time.'''

  def __init__(self):
    # Reference to the identifier of the last seen object
    # and the relevant cache object
    self._cache_tag = None
    self._cache = None

    # Lock for concurrent access
    self._lock = threading.Lock()

    # Keep the current PID to detect the use of multiprocessing parallelization
    # which breaks caching assumptions
    self._pid = os.getpid()

  def __del__(self):
    '''Garbage collection.
       Since only a cache controller can run .open() on lazy cache objects,
       the lazy cache object can be told to close as soon as possible.
       File handles can still be open legitimately.'''
    if self._cache is not None and self._pid == os.getpid():
      self._cache.close()

  def check(self, tag, open_method):
    '''The main cache controller access method. Checks if an object with name
       "tag" is cached. If so, returns a (pseudo-, ie. cached) file handle to
       this object.
       Otherwise, create a cache first, using the passed open_method()
       function, which returns a (true) file handle.'''
    with self._lock:
      currentpid = os.getpid()
      if currentpid != self._pid:
#       print "(%5s) Controller %s detected PID change, flushing %s" % (os.getpid(), format(id(self), '#x'), self._cache_tag)
        # Change the PID and drop reference to the previous cache.
        self._cache = None
        self._cache_tag = None
        self._pid = currentpid
        # NB: Explicitly do not close the cache in this case. Would be pointless.
      if tag != self._cache_tag:
        if self._cache is not None:
#         print "(%5s) Controller %s closing previous cache on %s" % (os.getpid(), format(id(self), '#x'), self._cache_tag)
          self._cache.close()
        self._cache_tag = tag
#       print "(%5s) Controller %s opening cache on %s" % (os.getpid(), format(id(self), '#x'), self._cache_tag)
        self._cache = dxtbx.filecache.lazy_file_cache(open_method())
#     else:
#       print "(%5s) Controller %s reporting cache hit on %s" % (os.getpid(), format(id(self), '#x'), self._cache_tag)
      return self._cache.open()
