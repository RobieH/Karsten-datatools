import numpy as np
import glob
import mmap
import os
import re
import scipy.io as sio
import multiprocessing
class probe:
    """Creates a probe object, allowing the user to efficiently and
    easily load probe data files.  Takes as input the path to the file.
    .. Note: The probe files should be saved with the following
    naming convention: nameNUMBER_DATA.  For example, Location1_ua is
    an appropriate name."""

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
        self._location()

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
           m = mmap.mmap(f.fileno(),0,prot=mmap.PROT_READ)    print indiv_locs

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
        self.data = data

    def _location(self):
        num = re.search("\d",self.filename).start()
        end = self.filename.index('_')
        self.location = self.filename[num:end]

def load_probe(location_files):
    mdict = {}
    location = ''
    for filename in location_files:
        probe_data = probe(filename)
        probe_data.load()
        mdict[probe_data.data_name] = probe_data.data
        location = probe_data.location
    saveName = 'Location_' + str(location)
    sio.savemat(saveName, mdict=mdict, oned_as='column')

if __name__ == '__main__':
    datadir = '/home/robie/Documents/python/probes/'
    #get list of probe files in the directory.
    specifer = 'L*'
    files = glob.glob(datadir + specifer)
    #break the files up into the individual locations
    loc = []
    for File in files:
        num = re.search("\d",File).start()
        end = File.index('_')
        loc.append(int(File[num:end]))
    num_locs = max(loc)
    data_entries = len(files)/num_locs

    loc = np.array(loc)
    data_entries = np.unique(loc)
    indiv_locs = []
    for i in xrange(data_entries):
        curr_loc = np.where(loc == i + 1)[0]
        indiv_locs.append([files[j] for j in curr_loc])

    #use multiprocessing to save the probe files
    pool = multiprocessing.Pool()
    pool.map(load_probe, indiv_locs)
