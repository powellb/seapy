#!/usr/bin/env python
"""
  Wrap your code with a time limit to prevent something from taking too long
  (getting into an infinite loop, etc.)

  **Examples**

  >>> from timeout import timeout
  >>> with timeout(seconds=3):
  >>>     do something

  Taken and slightly modified from Thomas Ahle at:
  <http://stackoverflow.com/questions/2281850/timeout-function-if-it-takes-too-long-to-finish>

"""


import errno
import os
import signal

class TimeoutError(Exception):
    pass

class timeout:
    def __init__(self, seconds=1, minutes=None, error_message='Timeout'):
        self.seconds = seconds
        if minutes is not None:
            self.seconds = minutes*60
        self.error_message = error_message
    def handle_timeout(self, signum, frame):
        raise TimeoutError(self.error_message)
    def __enter__(self):
        signal.signal(signal.SIGALRM, self.handle_timeout)
        signal.alarm(self.seconds)
    def __exit__(self, type, value, traceback):
        signal.alarm(0)


