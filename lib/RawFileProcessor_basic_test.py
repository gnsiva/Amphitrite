import os, subprocess, time

rawfolder = '121108_bsa_013.raw'
grain = 2
cpp = 'cppapplication.exe'

#rawfolder = os.path.join(r'E:\workspaces\Amphitrite_2.0\lib',rawfolder)
#cpp = os.path.join(r'E:\workspaces\Amphitrite_2.0\lib',cpp)

print os.getcwd()

#s = '%s %s 0 1 %d 0' %(cpp,rawfolder,grain)
#p = subprocess.Popen(s)

#output = subprocess.check_output(s)
#print output
output = subprocess.check_output(['cppapplication.exe'],shell=True)
print output

#output = subprocess.check_output(['dir'],shell=True)
#print output

