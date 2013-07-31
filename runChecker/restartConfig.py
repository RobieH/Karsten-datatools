from scipy.io import netcdf
import glob

<<<<<<< HEAD
datadir = '/output/'
#get the correct restart file, the second newest one (if there is more than one)
files = glob.glob(datadir + "*.nc")
#find the restart files
restart_files = []
ordinary_files =[]
for i in files:
	if "restart" in i:
		restart_files.append(i)
    else:
        ordinary_files.append(i)
=======
def getRestartTime():
    datadir = '/home/daugue6/capeislerestart/output/'
    #get the correct restart file, the second newest one (if there is more than one)
    files = glob.glob(datadir + "*.nc")
    #find the restart files
    restart_files = []
    for i in files:
        if "restart" in i:
            restart_files.append(i)
>>>>>>> 459f88c604b86bde6b1d45b93c76a3a1df5d0c92

    #get the latest restart file
    #we want to make sure that we are not restarting from a point that is a nan

    #look at the latest ordinary file for the last set of time series
    ord_nums = []
    for  i in ordinary_files:
        ord_nums.append(int(i[-7:-3]))

    newest = ord_nums.index(max(ord_nums))

    last_file = netcdf.netcdf_file(files[newest])
    time_series = last_file.variables['time'].data


    filenums = []
    for i in restart_files:
        filenums.append(int(i[-7:-3]))

    latest = filenums.index(max(filenums))

    #we need the times data from the restart file
    ncid = netcdf.netcdf_file(files[latest],'r')
    Times = ncid.variables['Times'].data

    ind = Times.shape[0] - 1

    #join the elements of the list into a single string

<<<<<<< HEAD
time = "\'"
for i in Times[ind,:]:
	if i == 'T':
		time += '\s'
	else:
		time += i
time += "\'"
=======
    time = "\'"
    for i in Times[ind,:]:
        if i == 'T':
            time += ' '
        else:
            time += i
    time += "\'"
>>>>>>> 459f88c604b86bde6b1d45b93c76a3a1df5d0c92

    name="\'{}\'".format(files[latest])

    return time, name
