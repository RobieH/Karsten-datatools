import numpy as np
import glob
import mmap
import os
import re
from datetime import datetime, timedelta
class probe:

    def __init__(self, filename):
        self.filename = filename
        self.size = 0
        self.data_name = ''

    def load(self):
        with open(self.filename, 'r') as f:
           self.m = mmap.mmap(f.fileno(),0,prot=mmap.PROT_READ)
        self._line_count()
        self._data_name()
        self._array_data()

    def print_file(self):
        with open(self.filename, 'r') as f:
            m = mmap.mmap(f.fileno(),0,prot=mmap.PROT_READ)
        data = m.readline()
        while data:
            print data
            data = m.readline()

    def _line_count(self):
        lines = 0
        with open(self.filename, 'r') as f:
           m = mmap.mmap(f.fileno(),0,prot=mmap.PROT_READ)
        readline = m.readline
        while readline():
            lines += 1
        self.size = lines - 18

    def _data_name(self):
        index = self.filename.index('_') + 1
        self.data_name = self.filename[index:]

    def _array_data(self):
        with open(self.filename, 'r') as f:
           m = mmap.mmap(f.fileno(),0,prot=mmap.PROT_READ)
        i = 0
        while i<18:
            row = m.readline()
            i+= 1
        #initialize arrays
        data = np.empty((self.size,2))
        readline = m.readline
        for i in xrange(self.size):
            row = m.readline()
            data[i,...] = [float(j) for j in row.split()]
        print (b-a).seconds, (b-a).microseconds
        self.data = data

    def _location(self):
        num = re.search("\d",self.filename).start()
        end = self.filename.index('_')
        self.location = self.filename[num:end]
