"""
This file implements time handling
including:
    - Stamping, get a time_stamp in string format
    - Timer
"""

import datetime
import time


def get_stamp():
    """
    This method returns a time stamp in string format
    :return: str: time stamp
    """
    now = datetime.datetime.now()
    stamp = f"{now.year}:{int(now.day / 10)}{int(now.day % 10)}:{int(now.month / 10)}{int(now.month % 10)}__"
    stamp += f"{int(now.hour / 10)}{int(now.hour % 10)}:{int(now.minute / 10)}{int(now.minute % 10)}:{int(now.second / 10)}{int(now.second % 10)}__"
    stamp += f"{now.microsecond}"
    return stamp


class Timer:
    def __init__(self):
        self.st = None
        self.delta = None
    
    def start(self):
        self.st = time.time()
    
    def stop(self):
        self.delta = time.time() - self.st




def compare_and_format(ref, comp):
    ratio = int((comp/ref)*100)
    space = 20
    diff = int((ratio - 100) / 10)
    head = ""
    tail = ""
    if diff > 0:
        head = "".join([" " for i in range(space)])
        tail = "".join(["-" for i in range(diff)]) + "b" + "".join([" " for i in range(space-diff - 1)])
    else:
        diff = -diff
        head = "".join([" " for i in range(space-diff -1)]) + "b" + "".join(["-" for i in range(diff)])
        tail = "".join([" " for i in range(space)])
    string = f"\t{head}:a:{tail}b/a ratio: {ratio}%"
    return string
    
