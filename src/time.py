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