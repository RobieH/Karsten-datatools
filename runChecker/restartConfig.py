from scipy.io import netcdf
import glob

datadir = '/home/daugue6/capeislerestart/output/'
#get the correct restart file, the second newest one (if there is more than one)
files = glob.glob(datadir + "*.nc")
#find the restart files
restart_files = []
for i in files:
	if "restart" in i:
		restart_files.append(i)

#get the latest restart file

filenums = []
for i in restart_files:
	filenums.append(int(i[-7:-3]))

latest = filenums.index(max(filenums))

#we need the times data from the restart file
ncid = netcdf.netcdf_file(files[latest],'r')
Times = ncid.variables['Times'].data

ind = Times.shape[0] - 1

#join the elements of the list into a single string

time = "\'"
for i in Times[ind,:]:
	if i == 'T': 
		time += ' '
	else:
		time += i
time += "\'"
print time
