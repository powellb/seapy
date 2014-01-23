#!/usr/bin/env python
"""
  qserver
  
  This module will use a queuing server to execute a number of tasks
  across multiple threads.
  
  Example
  -------
  
  tasks = (qserver.os_task("list files","ls -1"), \
           qserver.task("my job",my_func,arg1,arg2,arg3))
  qserver.start(tasks)
  
  Written by Brian Powell on 1/17/14
  Copyright (c)2014 University of Hawaii under the BSD-License.
"""
from __future__ import print_function
import Queue
import threading
import subprocess
import sys
import datetime

class task:
    def __init__(self, name, cmd, *args):
        self.cmd = cmd
        self.name = name
        self.args = args
        pass
        
    def run(self):
        if callable(self.cmd):
            self.cmd(*self.args)
        pass
    
    pass
    
class os_task(task):
    def run(self):
        subprocess.call(self.cmd, shell=True)

class process_thread(threading.Thread):
    def __init__(self, queue):
        threading.Thread.__init__(self)
        self.queue = queue
        
    def run(self):
        while True:
            item = self.queue.get()
            print(self.getName()+" running "+ \
                  item.name+" at "+str(datetime.datetime.now())+"\n")
            sys.stdout.flush()
            item.run()
            print(self.getName()+" completed "+item.name+ \
                  " at "+str(datetime.datetime.now())+"\n")
            sys.stdout.flush()
            self.queue.task_done()

def execute(tasks, nthreads=2):
    q = Queue.Queue()
    for i in range(nthreads):
         t = process_thread(q)
         t.daemon = True
         t.start()

    for item in tasks:
        q.put(item)

    q.join()


