"""
Job scheduler

Dispatch jobs with each submitted calculation

Common methods:
  results(): return an iterator yielding ( identifier, result ) tuples
  submit(target, args = (), kwargs = {}): submits job, return identifier
  is_empty(): returns whether there are no more jobs or results waiting
  is_full(): returns whether the number of currently processing jobs is at maximum
  shutdown(): no more job is submitted
  resume(): continues to process waiting jobs
  join(): shutdown, and finish processing all currently running jobs
  terminate(): kills all processing

Scheduler methods:
  job_count(): return number of jobs currently running
  process_count(): return number of currently active processes
"""

from __future__ import division

import time
from collections import deque
from Queue import Empty

from libtbx.queuing_system_utils import result
from libtbx.scheduling import identifier


# Capacity modes
class limited(object):
  """
  Limited number of jobs
  """

  def __init__(self, njobs):

    self.njobs = njobs


  def is_full(self, njobs):

    return self.njobs <= njobs


class unlimited(object):
  """
  Unlimited number of jobs (to be used with submission queue

  Note: this is a singleton object
  """

  @staticmethod
  def is_full(njobs):

    return False


def job_cycle(outqueue, jobid, target, args, kwargs):

  try:
    value = target( *args, **kwargs )

  except Exception, e:
    res = result.error( exception = e )

  else:
    res = result.success( value = value )

  outqueue.put( ( jobid, res ) )


class manager(object):
  """
  Job scheduler
  """

  def __init__(self, inqueue, job_factory, capacity, waittime = 0.01):

    self.inqueue = inqueue
    self.job_factory = job_factory
    self.capacity = capacity

    self.waittime = waittime

    self.process_data_for = {}
    self.waiting_results = set()
    self.waiting_jobs = deque()
    self.completed_results = deque()

    self.resume()


  def job_count(self):

    return len( self.process_data_for ) + len( self.waiting_jobs )


  def process_count(self):

    return len( self.process_data_for )


  def is_empty(self):

    return not (
      self.waiting_jobs or self.process_data_for or self.waiting_results
      or self.completed_results
      )


  def is_full(self):

    return self.capacity.is_full( njobs = self.job_count() )


  def results(self):

    while not self.is_empty():
      while not self.completed_results:
        self.wait()
        self.poll()

      yield self.completed_results.popleft()


  def submit(self, target, args = (), kwargs = {}):

    jobid = identifier()
    self.waiting_jobs.append( ( jobid, target, args, kwargs ) )
    return jobid


  def shutdown(self):

    self.active = False


  def resume(self):

    self.active = True


  def join(self):

    self.shutdown()

    while self.process_data_for:
      self.poll()
      self.wait()

    self.poll()


  def terminate(self):

    self.shutdown()

    for process in self.process_data_for.values():
      if process.is_alive():
        if hasattr( process, "terminate" ): # Thread has no terminate
          try:
            process.terminate()

          except Exception:
            pass

    self.join()


  # Internal methods
  def wait(self):

    time.sleep( self.waittime )


  def poll(self):

    # Check existing jobs
    for jobid in self.process_data_for.keys():
      process = self.process_data_for[ jobid ]

      if not process.is_alive():
        self.finish_job( jobid = jobid )


    # Submit new jobs
    while ( not self.capacity.is_full( njobs = self.process_count() )
      and self.waiting_jobs and self.active ):
      ( jobid, target, args, kwargs ) = self.waiting_jobs.popleft()

      process = self.job_factory(
        target = job_cycle,
        args = ( self.inqueue, jobid, target, args, kwargs ),
        )

      process.start()
      self.process_data_for[ jobid ] = process

    # Collect results
    while True:
      try:
        ( jobid, res ) = self.inqueue.get( timeout = self.waittime )

      except Empty:
        break

      if jobid in self.process_data_for:
        self.finish_job( jobid = jobid )

      self.waiting_results.remove( jobid )
      self.completed_results.append( ( jobid, res ) )


  def finish_job(self, jobid):

    process = self.process_data_for[ jobid ]
    process.join()
    exit_code = getattr( process, "exitcode", 0 ) # Thread has no "exitcode" attribute

    if exit_code != 0:
      err = getattr(
        process,
        "err",
        RuntimeError( "exit code = %s" % exit_code ),
        )
      self.completed_results.append( ( jobid, result.error( exception = err ) ) )

    else:
      self.waiting_results.add( jobid )

    del self.process_data_for[ jobid ]


class creator(object):
  """
  Information to create and destroy manager
  """

  def __init__(self, job_factory, queue_factory, capacity, waittime = 0.01):

    self.job_factory = job_factory
    self.queue_factory = queue_factory
    self.capacity = capacity
    self.waittime = waittime


  def create(self):

    return manager(
      inqueue = self.queue_factory.create(),
      job_factory = self.job_factory,
      capacity = self.capacity,
      waittime = self.waittime,
      )


  def destroy(self, manager):

    manager.terminate()
    manager.join()
    self.queue_factory.destroy( manager.inqueue )
