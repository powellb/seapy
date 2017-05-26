#!/usr/bin/env python
"""
  This module will execute any number of tasks using a queue of the
  specified number of threads.

  Define a list of qserver.task objects and then tell the server to
  execute the list with a given number of threads.

  **Examples**

  >>> tasks = (qserver.os_task("list files","ls -1"),
  >>>          qserver.task("my job",my_func,arg1,arg2,arg3))
  >>> qserver.execute(tasks, nthreads=2)

"""


import sys
if sys.version_info < (3, 0):
    from Queue import Queue
else:
    from queue import Queue
import threading
import subprocess
import sys
import datetime
from seapy.timeout import timeout, TimeoutError


class task:
    """
    task class simply defines a task for the queue server to process.
    It requires a descriptive name of the task for logging and the command
    to execute.

    Parameters
    ----------

    name : string
        title of the task
    cmd : string
        shell command with arguments to ru

    Returns
    -------
        none

    """

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
    """
    subclass of task to simply call shell commands

    It requires a descriptive name of the task for logging, the method
    to call, and a list of arguments for the method.

    Parameters
    ----------

    name : string
        title of the task
    cmd : method name, function pointer
        method to call
    args : vary [optional]
        arguments to pass to cmd

    Returns
    -------
        none

    """

    def run(self):
        subprocess.call(self.cmd, shell=True)


class process_thread(threading.Thread):
    def __init__(self, queue):
        threading.Thread.__init__(self)
        self.queue = queue

    def run(self):
        while True:
            item = self.queue.get()
            print(self.getName() + " running " +
                  item.name + " at " + str(datetime.datetime.now()) + "\n")
            sys.stdout.flush()
            item.run()
            print(self.getName() + " completed " + item.name +
                  " at " + str(datetime.datetime.now()) + "\n")
            sys.stdout.flush()
            self.queue.task_done()


def execute(tasks, nthreads=2):
    """
    Run the list of tasks in a queued server with the specified number of
    threads.

    Parameters
    ----------
    tasks : list
        list of task classes to execute in the queue server
    nthreads: int
        number of threads to use to process tasks in the queue

    Returns
    -------
    None
    """
    q = Queue()
    for i in range(nthreads):
        t = process_thread(q)
        t.daemon = True
        t.start()

    for item in tasks:
        q.put(item)

    q.join()


