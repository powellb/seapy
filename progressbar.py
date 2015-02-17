#!/usr/bin/env python
"""
    ASCII Progress Bar that work in IPython notebooks. 
  
    Code take from examples shown at:
  
  <http://nbviewer.ipython.org/github/ipython/ipython/blob/3607712653c66d63e0d7f13f073bde8c0f209ba8/docs/examples/notebooks/Animations_and_Progress.ipynb>
  
    Modified to include a timer and added an iterator that displays a 
    progressbar as it iterates.
"""
from __future__ import print_function
import sys, time
try:
    from IPython.display import clear_output
    have_ipython = True
except ImportError:
    have_ipython = False

class ProgressBar:
    def __init__(self, iterations):
        self.iterations = iterations
        self.prog_bar = '[]'
        self.fill_char = '*'
        self.width = 40
        self.__update_amount(0)
        self.start = time.process_time()
        if have_ipython:
            self.animate = self.animate_ipython
        else:
            self.animate = self.animate_noipython

    def animate_ipython(self, iter):
        print('\r', self, end='')
        sys.stdout.flush()
        self.update_iteration(iter + 1)

    def update_iteration(self, elapsed_iter):
        self.__update_amount((elapsed_iter / float(self.iterations)) * 100.0)
        t=time.process_time()
        delta=t-self.start
        togo=(delta/elapsed_iter)*(self.iterations-elapsed_iter)
        self.prog_bar += '  [%d of %d, %.1f secs elapsed/%.1f secs estimate]' % \
            (elapsed_iter, self.iterations, delta, togo)
        if elapsed_iter > self.iterations:
            print("\r [ COMPLETED %d ITERATIONS IN %.1f SECS ] %s" % \
                (self.iterations, delta, " "*60))


    def __update_amount(self, new_amount):
        percent_done = int(round((new_amount / 100.0) * 100.0))
        all_full = self.width - 2
        num_hashes = int(round((percent_done / 100.0) * all_full))
        self.prog_bar = '[' + self.fill_char * num_hashes + ' ' * (all_full - num_hashes) + ']'
        pct_place = (len(self.prog_bar) // 2) - len(str(percent_done))
        pct_string = '%d%%' % percent_done
        self.prog_bar = self.prog_bar[0:pct_place] + \
            (pct_string + self.prog_bar[pct_place + len(pct_string):])

    def __str__(self):
        return str(self.prog_bar)

class progress:
    """
    Iterator with progress bar. 
    """
    def __init__(self, iterable):
        self.size = len(iterable)
        self.pbar = ProgressBar(self.size+1)
        self.iterator = iter(iterable)
        self.count = 0

    def __getitem__(self, index):
        self.count += 1
        self.pbar.animate(self.count)
        return next(self.iterator)


## Usage
if __name__=="__main__":
    p = ProgressBar(100)
    for i in range(101):
        p.animate(i)
        time.sleep(0.01)

    print("use iterator")
    for i in progress(range(100)):
        time.sleep(0.01)