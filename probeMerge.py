import numpy as np
import glob
import mmap
import os
import re
import scipy.io as sio
import multiprocessing
from datetime import datetime, timedelta

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
        """load all the necessary information about the probe"""

        self._line_count()
        self._data_name()
        self._get_dim()
        self._array_data()
        self._location()
        self._depth_array()

    def print_file(self):
        """Print the data from the file on the terminal"""

        with open(self.filename, 'r') as f:
            m = mmap.mmap(f.fileno(),0,prot=mmap.PROT_READ)
        data = m.readline()
        while data:
            print data
            data = m.readline()

    def _line_count(self):
        """Determine the number of lines in the file, disregarding
        the header."""

        lines = 0
        with open(self.filename, 'r') as f:
           m = mmap.mmap(f.fileno(),0,prot=mmap.PROT_READ)

        readline = m.readline
        while readline():
            lines += 1
        self.size = lines - 18

    def _data_name(self):
        """Determine the name of the data (i.e. ua, va, el)"""
        item = os.path.basename(self.filename)
        index = item.index('_') + 1
        self.data_name = item[index:]

    def _get_dim(self):
        """Determine the dimension of the data"""

        with open(self.filename, 'r') as f:
            m = mmap.mmap(f.fileno(),0,prot=mmap.PROT_READ)
        i = 0
        while i < 18:
            m.readline()
            i += 1
        row = m.readline()
        self.num_cols = len(row.split())
        m.close()

    def _array_data(self):
        """load the data from the file and place it in arrays"""
        with open(self.filename, 'r') as f:
           m = mmap.mmap(f.fileno(),0,prot=mmap.PROT_READ)
        i = 0
        while i<18:
            row = m.readline()
            i+= 1
        #initialize arrays
        data = np.empty((self.size,self.num_cols))
        for i in xrange(self.size):
            row = m.readline()
            data[i,...] = row.split()
        self.time = data[:,0]
        self.data = data[:,1:]

    def _depth_array(self):
        """creates an array with the probe depths.
        currently assumes a linear distribution with the
        first probe on the surface of the ocean."""

        #find the dedyro dance the pain awaypth of the ocean in the file
        with open(self.filename, 'r') as f:
             m = mmap.mmap(f.fileno(),0,prot=mmap.PROT_READ)
        i = 0
        while i < 14:
            m.readline()
            i += 1
        row = m.readline()
        ocean_depth = float(row.split()[-1])
        delta_z = ocean_depth / (self.num_cols - 1)
        z = np.empty(self.num_cols-1)

        for i in xrange(self.num_cols-1):
            z[i] = -i * delta_z
        self.z = z


    def _location(self):
        """Determine the probe location/specifier"""
        item = os.path.basename(self.filename)
        num = re.search("\d", item).start()
        end = item.index('_')
        self.location = item[num:end]

def load_probe(location_files):
    mdict = {}
    location = ''
    for filename in location_files:
        probe_data = probe(filename)
        probe_data.load()
        mdict[probe_data.data_name] = probe_data.data
        mdict['time'] = probe_data.time
        mdict['z'] = probe_data.z
        location = probe_data.location
    saveName = 'Location_' + str(location)
    sio.savemat(saveName, mdict=mdict, oned_as='column')

if __name__ == '__main__':
    datadir = '/home/robie/Documents/python/aidans_a_dick/'
    #get list of probe files in the directory.
    specifer = 'L*'
    files = glob.glob(datadir + specifer)
    #break the files up into the individual locations
    loc = []
    for File in files:
        item = os.path.basename(File)
        num = re.search("\d", item).start()
        end = item.index('_')
        loc.append(int(item[num:end]))
    num_locs = max(loc)
    data_entries = len(files)/num_locs

    loc = np.array(loc)
    data_entries = np.unique(loc)
    indiv_locs = []
    for i in data_entries:
        curr_loc = np.where(loc == i)[0]
        indiv_locs.append([files[j] for j in curr_loc])
    #use multiprocessing to save the probe files
    pool = multiprocessing.Pool()
    pool.map(load_probe, indiv_locs)
