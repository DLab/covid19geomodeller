#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from datetime import datetime
from datetime import timedelta

def timeJStoPy(t):
    return datetime.strptime(t[:10],'%Y-%m-%d')


def timetxttoDatetime(t):
    return datetime.strptime(t,'%Y-%m-%d')

