import os
import time
import multiprocessing
#os.system("ls *dat | xargs -n 1 accelsearch -zmax 0")

def parallel_func(cmdlist):
    for cmd in cmdlist:
         os.system("echo "+cmd+" >> command.txt")
         os.system(cmd)

proc =28
list = []
dedispfile = open("Dedispcmd.txt")
for line in dedispfile.readlines():
     list.append(line)

Process = []
for i in range(proc):
    cmdlist = []
    Process.append(cmdlist)
for i in range(len(list)):
    Process[i%proc].append(list[i])
for cmdlist in Process:
    #print filelist
    multiprocessing.Process(target=parallel_func, args=[cmdlist]).start()

