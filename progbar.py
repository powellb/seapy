# Sebastian Raschka 2014
#
# Taken from:
# Sebastian Raschka 01/27/2014
# PyPrind - Python Progress Indicator module
#        The PyPrind (Python Progress Indicator) module lets you visualize the
#        progress of a programming task in Python via a progress bar or a
#        percentage indicator.
#        Progress bars are visualized via a `ProgBar()` object, and
#        alternatively, the progress can be shown as an percentage via the
#        `ProgPercent()` object.
#
# Progress Bar class to instantiate a progress bar object
# that is printed to the standard output screen to visualize the
# progress in a iterative Python procedure
#
# EXAMPLE:
#         n = 1000000
#         my_bar = seapy.ProgBar(n, width=40, stream=2)
#         for i in range(n):
#             my_bar.update()

from math import floor
import time
import sys
import os

class Prog():
    def __init__(self, iterations, track_time, stream, title, monitor):
        """ Initializes tracking object. """
        self.cnt = 0
        self.title = title
        self.max_iter = float(iterations) # to support Python 2.x
        self.track = track_time
        self.start = time.time()
        self.end = None
        self.total_time = 0.0
        self.monitor = monitor
        self.stream = stream
        self._stream_out = self._no_stream
        self._stream_flush = self._no_stream
        self._check_stream()
        self._print_title()
        
        if monitor:
            import psutil
            self.process = psutil.Process()

    def _check_stream(self):
        """ Determines which output stream (stdout, stderr, or custom) to use. """
        if self.stream == 1 and os.isatty(sys.stdout.fileno()):
            self._stream_out = sys.stdout.write
            self._stream_flush = sys.stdout.flush
        elif self.stream == 2 and os.isatty(sys.stderr.fileno()):
            self._stream_out = sys.stderr.write
            self._stream_flush = sys.stderr.flush
        elif self.stream is not None and hasattr(self.stream, 'write'):
            self._stream_out = self.stream.write
            self._stream_flush = self.stream.flush
        else:
            print('Warning: No valid output stream.')

    def _elapsed(self):
        """ Returns elapsed time at update. """
        return time.time() - self.start

    def _calc_eta(self):
        """ Calculates estimated time left until completion. """
        elapsed = self._elapsed()
        if self.cnt == 0 or elapsed < 0.001:
            return None
        rate = float(self.cnt) / elapsed
        return (float(self.max_iter) - float(self.cnt)) / rate

    def _calc_percent(self):
        """Calculates the rel. progress in percent with 2 decimal points."""
        return round(self.cnt / self.max_iter * 100, 2)

    def _no_stream(self, text=None):
        """ Called when no valid output stream is available. """
        pass

    def _finish(self):
        """ Determines if maximum number of iterations (seed) is reached. """
        if self.cnt == self.max_iter:
            self.total_time = self._elapsed()
            self.end = time.time()
            if self.track:
                self._stream_out('\nTotal time elapsed: {:.3f} sec'.format(self.total_time))
            self._stream_out('\n')

    def _print_title(self):
        """ Prints tracking title at initialization. """
        if self.title:
            self._stream_out('{}\n'.format(self.title))
            self._stream_flush()

    def __repr__(self):
        str_start = time.strftime('%m/%d/%Y %H:%M:%S', time.localtime(self.start))
        str_end = time.strftime('%m/%d/%Y %H:%M:%S', time.localtime(self.end))
        if not self.monitor:
            return """Title: {}
                      Started: {}
                      Finished: {}
                      Total time elapsed: {:.3f} sec""".format(self.title, str_start, 
                                                               str_end, self.total_time)
        else:
            cpu_total = self.process.cpu_percent()
            mem_total = self.process.memory_percent()
            return """Title: {}
                      Started: {}
                      Finished: {}
                      Total time elapsed: {:.3f} sec
                      CPU %: {:2f}
                      Memory %: {:2f}""".format(self.title, str_start, str_end, self.total_time, cpu_total, mem_total)

    def __str__(self):
        return self.__repr__()

class ProgBar(Prog):
    """
    Initializes a progress bar object that allows visuzalization
    of an iterational computation in the standard output screen. 

    Keyword Arguments:
        iterations (int): number of iterations of the computation
        track_time (bool): default True. Prints elapsed time when loop has finished
        width (int): default 30. Sets the progress bar width in characters.
        stream (int): default 2. Takes 1 for stdout, 2 for stderr, or given stream object
        title (str): default ''. A title for the progress bar
        monitor (bool): default False. Monitors CPU and memory usage if True 
            (requires 'psutil' package).

    """
    def __init__(self, iterations, track_time=True, width=30, stream=2, title='', monitor=False):
        Prog.__init__(self, iterations, track_time, stream, title, monitor)
        self.bar_width = width
        self._adjust_width()
        self.last_progress = 0
        self._print_labels()
        self._print_progress_bar(0)
        if monitor:
            self.process.cpu_percent()
            self.process.memory_percent()

    def _adjust_width(self):
        """Shrinks bar if number of iterations is less than the bar width"""
        if self.bar_width > self.max_iter:
            self.bar_width = int(self.max_iter) 
            # some Python 3.3.3 users specifically
            # on Linux Red Hat 4.4.7-1, GCC v. 4.4.7
            # reported that self.max_iter was converted to
            # float. Thus this fix to prevent float multiplication of chars.

    def _print_labels(self):
        self._stream_out('0% {} 100%\n'.format(' ' * (self.bar_width - 6)))
        self._stream_flush()

    def _print_progress_bar(self, progress):
        remaining = self.bar_width - progress
        self._stream_out('[{}{}]'.format('#' * int(progress), ' ' * int(remaining)))
        # int() fix for Python 2 users
        self._stream_flush()

    def _print_eta(self):
        self._stream_out(' | ETA[sec]: {:.3f} '.format(self._calc_eta()))
        self._stream_flush()

    def _print_bar(self):
        progress = floor(self._calc_percent() / 100 * self.bar_width)
        if progress > self.last_progress:
            self._stream_out('\r')
            self._print_progress_bar(progress)
            if self._calc_eta() and self.track:
                self._print_eta()
        self.last_progress = progress

    def update(self, iterations=1):
        """
        Updates the progress bar in every iteration of the task.

        Keyword arguments:
            iterations (int): default argument can be changed to integer values
                >=1 in order to update the progress indicators more than once 
                per iteration.

        """
        self.cnt += iterations
        self._print_bar()
        self._finish() 
