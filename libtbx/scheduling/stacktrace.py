from __future__ import division

__prev_excepthook = None
__last_exception = ( None, None )


def set_last_exception(exception, printout):

  global __last_exception
  __last_exception = ( exception, printout )


def exc_info():

  return __last_exception


def enable():

  global __prev_excepthook

  if __prev_excepthook is None:
    import sys
    __prev_excepthook = sys.excepthook
    sys.excepthook = stacktrace_excepthook
    set_last_exception( None, None )


def disable():

  global __prev_excepthook

  if __prev_excepthook is not None:
    import sys
    sys.excepthook = __prev_excepthook
    __prev_excepthook = None
    set_last_exception( None, None )


def stacktrace_excepthook(ex_cls, ex, tb):
  """
  Prints traceback printout from module global if exists
  """

  data = exc_info()

  if ex is data[0]:
    import traceback
    import sys
    sys.stderr.write(
      "".join( data[1] + traceback.format_exception_only( ex_cls, ex ) )
      )

  else:
    __prev_excepthook( ex_cls, ex, tb )


# Information handling
def no_crash_info(exception):

  raise exception


class stacktrace_info(object):
  """
  Crash printout
  """

  def __init__(self, printout, header):

    self.header = header
    self.printout = printout


  def __call__(self, exception):

    set_last_exception( exception, [ self.header ] + self.printout )
    raise exception


  @classmethod
  def from_stderr(cls, message):

    return cls( printout = [ message ], header = "Stacktrace:\n" )


  @classmethod
  def from_traceback(cls, lines):

    return cls( printout = lines, header = "Traceback (most recent call last):\n" )


class traceback_info(object):
  """
  Transparent traceback transformation on pickling
  """

  def __init__(self, traceback):

    self.traceback = traceback
    self.raise_handler = self.raise_with_traceback
    self.getstate_handler = self.getstate_with_traceback


  def __call__(self, exception):

    self.raise_handler( exception = exception )


  def __getstate__(self):

    return self.getstate_handler()


  def __setstate__(self, state):

    self.stacktrace = state
    self.raise_handler = self.raise_with_stacktrace
    self.getstate_handler = self.getstate_with_stacktrace


  def raise_with_traceback(self, exception):

    raise exception, None, self.traceback


  def getstate_with_traceback(self):

    import traceback
    return stacktrace_info.from_traceback(
      lines = traceback.format_tb( self.traceback )
      )


  def raise_with_stacktrace(self, exception):

    self.stacktrace( exception = exception )


  def getstate_with_stacktrace(self):

    return self.stacktrace
