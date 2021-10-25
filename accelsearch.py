import os
import time
import multiprocessing
#os.system("ls *dat | xargs -n 1 accelsearch -zmax 0")

def parallel_func(filelist):
    for file in filelist:
         os.system("accelsearch -ncpus 1 -numharm 16 -zmax 200 " + file)

proc = 28
list = os.popen("ls *.fft | xargs -n 1").read()
list = list.split('\n')
list = list[:-1]

Process = []
for i in range(proc):
    filelist = []
    Process.append(filelist)
for i in range(len(list)):
    Process[i%proc].append(list[i])

for filelist in Process:
    #print filelist
    multiprocessing.Process(target=parallel_func, args=[filelist]).start()

