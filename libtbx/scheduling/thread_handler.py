from __future__ import division

import threading
import Queue


class exception_capturing_thread(threading.Thread):
  """
  A Thread-specialization that adds the exitcode attribute

  http://stackoverflow.com/questions/986616/python-thread-exit-code
  """

  def __init__(
    self,
    group=None,
    target=None,
    name=None,
    propagate_error_message = False,
    args=(),
    kwargs={},
    ):

    super( exception_capturing_thread, self ).__init__(group, target, name, args, kwargs)
    self.propagate_error_message = propagate_error_message


  def run(self):

    try:
      super( exception_capturing_thread, self ).run()

    except Exception, e:
      self.exitcode = 1
      self.err = e

      if self.propagate_error_message:
        import traceback
        import sys
        traceback.print_exc( file = sys.stderr )

    else:
      self.exitcode = 0
      self.err = None


class qfactory(object):
  """
  Creator pattern for Queue.Queue, also include destruction

  Note this is a singleton object
  """

  @staticmethod
  def create():

    return Queue.Queue()


  @staticmethod
  def destroy(queue):

    pass
