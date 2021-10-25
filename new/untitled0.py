# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 16:20:49 2021

@author: 24566
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pylab import *
from scipy.optimize import curve_fit
import linecache

##mjd58767
file1=r"C:\\Users\\24566\\Desktop\\DM and flux\\parkes\\mjd58767\\newmjd58767\\mjd58767.calib.64bin.2minpF.txt"
flux = []
error = []
orbital = []
for i in range(9, 79):
    a = linecache.getline(file1, i)
    a1 = a.split()
    flux.append(float(a1[4]))
    error.append(float(a1[5]))
    orbital.append(float(a1[2]))
for i in range(len(orbital)):
    orbital[i] = orbital[i] / 129 + 0.19
new_x_58767 = []
new_y_58767 = []
new_error2_58767 = []
for i in range(len(error)):
    if error[i] != 0:
       new_x_58767.append(orbital[i])
       new_y_58767.append(flux[i])
       new_error2_58767.append(error[i])
y_u=np.mean(new_y_58767)+np.std(new_y_58767)
y_s=[]
for i in range(len(new_y_58767)):
    if new_y_58767[i]>y_u:
        y_s.append(new_y_58767[i])
y_max=np.mean(y_s)
y_norm_58767 = []
error_norm_58767 = []
for i in range(len(new_y_58767)):
    y_norm_58767.append((new_y_58767[i]) / (y_max))
    error_norm_58767.append(new_error2_58767[i] / (y_max))
    #y_norm_58767.append((new_y_58767[i] - np.min(new_y_58767)) / (np.max(new_y_58767) - np.min(new_y_58767)))
    #error_norm_58767.append(new_error2_58767[i] / (np.max(new_y_58767) - np.min(new_y_58767)))
###mjd58796
file1=r"C:\\Users\\24566\\Desktop\\DM and flux\\parkes\\mjd58796\\newmjd58796\\mjd58796.calib.2minpF.txt"
y_58796_1=[]
error_58796_1=[]
x_58796_1=[]
for i in range(9,320):
    a=linecache.getline(file1,i)
    a1=a.split()
    y_58796_1.append(float(a1[4]))
    error_58796_1.append(float(a1[5]))
    x_58796_1.append(float(a1[2]))
for i in range(len(x_58796_1)):
    x_58796_1[i]=x_58796_1[i]/129

new_x_58796=[]
new_y_58796=[]
new_error2_58796=[]
for i in range(len(error_58796_1)):
    #if error2_58796_1[i]<1:
        new_x_58796.append(x_58796_1[i])
        new_y_58796.append(y_58796_1[i])
        new_error2_58796.append(error_58796_1[i])


y_norm_58796=[]
error_norm_58796=[]
for i in range(len(new_y_58796)):
    y_norm_58796.append(new_y_58796[i]/(np.max(new_y_58796)*0.865))
    error_norm_58796.append(new_error2_58796[i]/(np.max(new_y_58796)*0.865))

# y_norm_58796=[]
# error_norm_58796=[]
# for i in range(len(newy)):
#  y_norm_58796.append((newy[i]-np.min(newy))/(np.max(newy)-np.min(newy)))
# error_norm_58796.append(newerror[i]/(np.max(newy)-np.min(newy)))
####分别归一
#y_norm_58796 = y_norm_58796_1 + y_norm_58796_2 + y_norm_58796_3 + y_norm_58796_4 + y_norm_58796_5
#error_norm_58796 = error_norm_58796_1 + error_norm_58796_2 + error_norm_58796_3 + error_norm_58796_4 + error_norm_58796_5
#####
# 58892
file1 = r"C:\\Users\\24566\\Desktop\\DM and flux\\parkes\\mjd58892\\mjd58892.2min.txt"
flux = []
error = []
orbital = []
for i in range(9, 65):
    a = linecache.getline(file1, i)
    a1 = a.split()
    flux.append(float(a1[4]))
    error.append(float(a1[5]))
    orbital.append(float(a1[2]))
for i in range(len(orbital)):
    orbital[i] = orbital[i] / 129 + 0.613
new_x_58892 = []
new_y_58892 = []
new_error2_58892 = []
for i in range(len(error)):
    # if error[i]<flux[i]:
    new_x_58892.append(orbital[i])
    new_y_58892.append(flux[i])
    new_error2_58892.append(error[i])
y_norm_58892 = []
error_norm_58892 = []
for i in range(len(new_y_58892)):
    y_norm_58892.append(new_y_58892[i] / (np.max(new_y_58892) * 0.865))
    error_norm_58892.append(new_error2_58892[i] / (np.max(new_y_58892) * 0.865))
    #y_norm_58892.append((new_y_58892[i] - np.min(new_y_58892)) / (np.max(new_y_58892) - np.min(new_y_58892)))
    #error_norm_58892.append(new_error2_58892[i] / (np.max(new_y_58892) - np.min(new_y_58892)))

########58900
file1 = r"C:\\Users\\24566\\Desktop\\DM and flux\\parkes\\mjd58900\\mjd58900.calib.2min.txt"
x_58900 = []
y_58900 = []
error_58900 = []
for i in range(9, 120):
    a = linecache.getline(file1, i)
    a1 = a.split()
    y_58900.append(float(a1[4]))
    error_58900.append(float(a1[5]))
    x_58900.append(float(a1[2]))
for i in range(len(x_58900)):
    x_58900[i] = x_58900[i] / 129 + 0.63

error2_58900 = []
for i in range(len(error_58900)):
    error2_58900.append(error_58900[i] / y_58900[i])

new_x_58900 = []
new_y_58900 = []
new_error2_58900 = []
for i in range(len(error2_58900)):
    # if error2_58900[i]<1:
    new_x_58900.append(x_58900[i])
    new_y_58900.append(y_58900[i])
    new_error2_58900.append(error_58900[i])

y_norm_58900 = []
error_norm_58900 = []
for i in range(len(new_error2_58900)):
    y_norm_58900.append(new_y_58900[i] / (np.max(new_y_58900) * 0.865))
    error_norm_58900.append(new_error2_58900[i] / (np.max(new_y_58900) * 0.865))
    #y_norm_58900.append((new_y_58900[i] - np.min(new_y_58900)) / (np.max(new_y_58900) - np.min(new_y_58900)))
    #error_norm_58900.append(new_error2_58900[i] / (np.max(new_y_58900) - np.min(new_y_58900)))

##########
# mjd58974
file1 = r"C:\\Users\\24566\\Desktop\\DM and flux\\parkes\\mjd58974\\mjd58974.first.3min.txt"
y_58974 = []
error_58974 = []
x_58974 = []
for i in range(9, 52):
    a = linecache.getline(file1, i)
    a1 = a.split()
    y_58974.append(float(a1[4]))
    error_58974.append(float(a1[5]))
    x_58974.append(float(a1[2]))
for i in range(len(x_58974)):
    x_58974[i] = x_58974[i] / 129
new_x_58974 = []
new_y_58974 = []
new_error2_58974 = []
for i in range(len(error_58974)):
    # if error[i]<flux[i]:
    new_x_58974.append(x_58974[i])
    new_y_58974.append(y_58974[i])
    new_error2_58974.append(error_58974[i])
y_norm_58974 = []
error_norm_58974 = []
for i in range(len(new_y_58974)):
    y_norm_58974.append(new_y_58974[i] / (np.max(new_y_58974) * 0.865))
    error_norm_58974.append(new_error2_58974[i] / (np.max(new_y_58974) * 0.865))
    #y_norm_58974.append((new_y_58974[i] - np.min(new_y_58974)) / (np.max(new_y_58974) - np.min(new_y_58974)))
    #error_norm_58974.append(new_error2_58974[i] / (np.max(new_y_58974) - np.min(new_y_58974)))

################
# mjd59017
file1 = r"C:\\Users\\24566\\Desktop\\DM and flux\\meerkat\\meerkat.64s.txt"
x_59017_2 = []
y_59017_2 = []
error_59017_2 = []
for i in range(9, 122):
    a = linecache.getline(file1, i)
    a1 = a.split()
    y_59017_2.append(float(a1[4]))
    error_59017_2.append(float(a1[5]))
    x_59017_2.append(float(a1[2]))
for i in range(len(x_59017_2)):
    x_59017_2[i] = x_59017_2[i] / 129 + 0.488

error2_59017_2 = []
for i in range(len(error_59017_2)):
    error2_59017_2.append(error_59017_2[i] / y_59017_2[i])
new_x_59017_2 = []
new_y_59017_2 = []
new_error2_59017_2 = []
for i in range(len(error2_59017_2)):
    if error2_59017_2[i] < 1:
        new_x_59017_2.append(x_59017_2[i])
        new_y_59017_2.append(y_59017_2[i])
        new_error2_59017_2.append(error_59017_2[i])
y_norm_59017 = []
error_norm_59017 = []
for i in range(len(new_y_59017_2)):
    y_norm_59017.append(new_y_59017_2[i]/(np.max(new_y_59017_2)*0.865))
    error_norm_59017.append(new_error2_59017_2[i]/(np.max(new_y_59017_2)*0.865))
    #y_norm_59017.append((new_y_59017_2[i] - np.min(new_y_59017_2)) / (np.max(new_y_59017_2) - np.min(new_y_59017_2)))
    #error_norm_59017.append(new_error2_59017_2[i] / (np.max(new_y_59017_2) - np.min(new_y_59017_2)))
# plt.errorbar(new_x_59017_2,new_y_59017_2,yerr=new_error2_59017_2,fmt='-k')

###############################################################################################
# mjd59018
file1 = r"C:\\Users\\24566\\Desktop\\DM and flux\\parkes\\mjd59018\\mjd59018.calib.3min.txt"
y_59018 = []
error_59018 = []
x_59018 = []
for i in range(9, 57):
    a = linecache.getline(file1, i)
    a1 = a.split()
    y_59018.append(float(a1[4]))
    error_59018.append(float(a1[5]))
    x_59018.append(float(a1[2]))
for i in range(len(x_59018)):
    x_59018[i] = x_59018[i] / 129 + 0.43

new_x_59018 = []
new_y_59018 = []
new_error2_59018 = []
for i in range(len(error_59018)):
     if error_59018[i]!=0:
        new_x_59018.append(x_59018[i])
        new_y_59018.append(y_59018[i])
        new_error2_59018.append(error_59018[i])

y_norm_59018 = []
error_norm_59018 = []
for i in range(len(new_y_59018)):
    y_norm_59018.append(new_y_59018[i] / (np.max(new_y_59018) * 0.865))
    error_norm_59018.append(new_error2_59018[i] / (np.max(new_y_59018) * 0.865))
    #y_norm_59018.append((new_y_59018[i] - np.min(new_y_59018)) / (np.max(new_y_59018) - np.min(new_y_59018)))
    #error_norm_59018.append(new_error2_59018[i] / (np.max(new_y_59018) - np.min(new_y_59018)))
###################################################################################################################
# mjd59187
file1 = r"C:\\Users\\24566\\Desktop\\DM and flux\\parkes\\mjd59187\\mjd59187.calib.3min.txt"
y_59187 = []
x_59187 = []
error_59187 = []
for i in range(9, 124):
    a = linecache.getline(file1, i)
    a1 = a.split()
    y_59187.append(float(a1[4]))
    error_59187.append(float(a1[5]))
    x_59187.append(float(a1[2]))
for i in range(len(x_59187)):
    x_59187[i] = x_59187[i] / 129 + 0.512

error2_59187 = []
for i in range(len(error_59187)):
    if error_59187[i] != 0 and x_59187[i] < 2:
        error2_59187.append(error_59187[i] / y_59187[i])
new_x_59187 = []
new_y_59187 = []
new_error2_59187 = []
for i in range(len(error2_59187)):
    # if error2_59187[i]<1:
    new_x_59187.append(x_59187[i])
    new_y_59187.append(y_59187[i])
    new_error2_59187.append(error_59187[i])
y_norm_59187 = []
error_norm_59187 = []
for i in range(len(new_y_59187)):
    y_norm_59187.append((new_y_59187[i] - np.min(new_y_59187)) / (np.max(new_y_59187) - np.min(new_y_59187)))
    error_norm_59187.append(new_error2_59187[i] / (np.max(new_y_59187) - np.min(new_y_59187)))
#############################################
# subplots_adjust(left=0.15,bottom=0.1,top=0.9,right=1.2,hspace=0.2,wspace=0.5)
fig = plt.figure(figsize=(18, 20))
subplots_adjust(left=0.15, bottom=0.1, top=0.9, right=1.2, hspace=0.2, wspace=0.2)

ax1 = subplot2grid((4, 4), (0, 0), colspan=1)
ax1.errorbar(new_x_58767, y_norm_58767, yerr=error_norm_58767, fmt='.k',capsize=3)
plt.axhline(y=0, c="dimgray", ls="--")
plt.title('MJD58767(a)', x=0.1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.ylim(-0.55, 1.4)
##############拟合
####mjd8767
x_58767_1_egress = []
y_58767_1_egress = []
for i in range(len(new_x_58767)):
    if new_x_58767[i] > 0.2 and new_x_58767[i] < 0.5:
        x_58767_1_egress.append(new_x_58767[i])
        y_58767_1_egress.append(y_norm_58767[i])


def func(x, p1, p2):
    return 1 / (np.exp((p1 - x) / p2) + 1)


popt_58767_1_egress, pcov_58767_1_egress = curve_fit(func, x_58767_1_egress, y_58767_1_egress)
x = np.linspace(0.2, 0.6, 30)
plt.plot(x, func(x, *popt_58767_1_egress),"b", linewidth=6, linestyle='--',alpha=1)
# print(popt_58767_1_egress)
###########
x_58767_2_ingress = []
y_58767_2_ingress = []
for i in range(len(new_x_58767)):
    if new_x_58767[i] > 0.6 and new_x_58767[i] < 0.85:
        x_58767_2_ingress.append(new_x_58767[i])
        y_58767_2_ingress.append(y_norm_58767[i])

# def func(x,p1,p2):
# return 1/(np.exp((x-p1)/t1)+1)#+a/(np.exp((p2-x)/t2)+1)

# popt_58767_2_ingress,pcov_58767_2_ingress=curve_fit(func,x_58767_2_ingress,y_58767_2_ingress)

# plt.plot(x,func(x,*popt_58767_2_ingress),"b--",linewidth=2,label="Ingress")
x = np.linspace(0.6, 0.75, 50)
a = 1
p1 = 0.6526
p2 = 0.8507
t1 = 0.01972
t2 = 0.0415
y = a / (np.exp((x - p1) / t1) + 1) + a / (np.exp((p2 - x) / t2) + 1)
plt.plot(x, y, "b", linewidth=6, linestyle='--',alpha=1)
x = np.linspace(0.75, 0.85, 50)
a = 1
p1 = 0.6526
p2 = 0.8507
t1 = 0.01972
t2 = 0.0415
y = a / (np.exp((x - p1) / t1) + 1) + a / (np.exp((p2 - x) / t2) + 1)
plt.plot(x, y, "b", linewidth=6, linestyle='--',alpha=1)
#plt.legend(loc="upper right")
# print(popt_58767_2_egress)
# plt.legend()
###############
ax2 = subplot2grid((4, 4), (1, 0), colspan=4)
plt.xticks([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5])
ax2.errorbar(new_x_58796
             ,y_norm_58796, yerr=error_norm_58796, fmt='.k',capsize=3)
# ax2.errorbar(newx,newy,yerr=newerror,fmt='--k')
plt.title('MJD58796(b)', x=0.03)
plt.xlim(0, 5)
plt.axhline(y=0, c="dimgray", ls="--")
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.ylim(-0.8, 1.4)
#plt.yticks([])
##########拟合，第一段
x_58796_1_in=[]
y_58796_1_in=[]
for i in range(len(new_x_58796)):
      if new_x_58796[i]>0.4 and new_x_58796[i]<1:
            x_58796_1_in.append(new_x_58796[i])
            y_58796_1_in.append(y_norm_58796[i])
def func(x,p1,p2,t1,t2):
    return 1/(np.exp((x-p1)/t1)+1)+1/(np.exp((p2-x)/t2)+1)
popt,pcov=curve_fit(func,x_58796_1_in,y_58796_1_in)
x=np.linspace(0.3,1,50)
plt.plot(x,func(x,*popt),"b", linewidth=6, linestyle='--',alpha=1)
####
x_58796_2_in=[]
y_58796_2_in=[]
for i in range(len(new_x_58796)):
      if new_x_58796[i]>1 and new_x_58796[i]<1.45:
            x_58796_2_in.append(new_x_58796[i])
            y_58796_2_in.append(y_norm_58796[i])
def func(x,p1,p2,t1,t2):
    return 1/(np.exp((x-p1)/t1)+1)+1/(np.exp((p2-x)/t2)+1)

popt,pcov=curve_fit(func,x_58796_2_in,y_58796_2_in)
x=np.linspace(1,1.45,50)
plt.plot(x,func(x,*popt),"b", linewidth=6, linestyle='--',alpha=1)
####
x_58796_2_2_in=[]
y_58796_2_2_in=[]
for i in range(len(new_x_58796)):
      if new_x_58796[i]>1.45 and new_x_58796[i]<2.4:
            x_58796_2_2_in.append(new_x_58796[i])
            y_58796_2_2_in.append(y_norm_58796[i])
def func(x,p1,p2,t1,t2):
    return 1/(np.exp((x-p1)/t1)+1)+1/(np.exp((p2-x)/t2)+1)

popt,pcov=curve_fit(func,x_58796_2_2_in,y_58796_2_2_in)
x=np.linspace(1.45,2.4,20)
plt.plot(x,func(x,*popt),"b", linewidth=6, linestyle='--',alpha=1)
###
x_58796_3_2_in=[]
y_58796_3_2_in=[]
for i in range(len(new_x_58796)):
      if new_x_58796[i]>2.3 and new_x_58796[i]<2.65:
            x_58796_3_2_in.append(new_x_58796[i])
            y_58796_3_2_in.append(y_norm_58796[i])

def func(x,p1,t1):
    return 1/(np.exp((x-p1)/t1)+1)#+0.8/(np.exp((p2-x)/t2)+1)
popt,pcov=curve_fit(func,x_58796_3_2_in,y_58796_3_2_in)
x=np.linspace(2.4,2.55,20)
plt.plot(x,func(x,*popt),"b", linewidth=6, linestyle='--',alpha=1)
###
x_58796_3_2_eg=[]
y_58796_3_2_eg=[]
for i in range(len(new_x_58796)):
      if new_x_58796[i]>2.55 and new_x_58796[i]<2.7:
            x_58796_3_2_eg.append(new_x_58796[i])
            y_58796_3_2_eg.append(y_norm_58796[i])

def func(x,p2,t2):
    return 0.8/(np.exp((p2-x)/t2)+1)
popt,pcov=curve_fit(func,x_58796_3_2_eg,y_58796_3_2_eg)
x=np.linspace(2.55,2.7,20)
plt.plot(x,func(x,*popt),"b", linewidth=6, linestyle='--',alpha=1)
#####拟合第四段
x_58796_4_1_in=[]
y_58796_4_1_in=[]
for i in range(len(new_x_58796)):
      if new_x_58796[i]>3.12 and new_x_58796[i]<3.9:
            x_58796_4_1_in.append(new_x_58796[i])
            y_58796_4_1_in.append(y_norm_58796[i])

def func(x,p1,p2,t1,t2):
    return 0.89/(np.exp((x-p1)/t1)+1)+0.89/(np.exp((p2-x)/t2)+1)

popt,pcov=curve_fit(func,x_58796_4_1_in,y_58796_4_1_in)
x=np.linspace(3.12,3.8,20)
plt.plot(x,func(x,*popt),"b", linewidth=6, linestyle='--',alpha=1)
#x=np.linspace(3.45,3.8,20)
#plt.plot(x,func(x,*popt),"r",linewidth=2,linestyle='--')
##########
x_58796_4_1_in=[]
y_58796_4_1_in=[]
for i in range(len(new_x_58796)):
      if new_x_58796[i]>3.8 and new_x_58796[i]<4:
            x_58796_4_1_in.append(new_x_58796[i])
            y_58796_4_1_in.append(y_norm_58796[i])

a=0.76
def func(x,p1,p2,t1,t2):
    return a/(np.exp((x-p1)/t1)+1)+a/(np.exp((p2-x)/t2)+1)

popt,pcov=curve_fit(func,x_58796_4_1_in,y_58796_4_1_in)
x=np.linspace(3.8,4,20)
plt.plot(x,func(x,*popt),"b", linewidth=6, linestyle='--',alpha=1)
#x=np.linspace(3.95,4,20)
#plt.plot(x,func(x,*popt),"r",linewidth=2,linestyle='--')
print(popt)
print(np.sqrt(np.diag(pcov)))
###########拟合第5段
x_58796_5_1_in=[]
y_58796_5_1_in=[]
for i in range(len(new_x_58796)):
      if new_x_58796[i]>4 and new_x_58796[i]<4.6:
            x_58796_5_1_in.append(new_x_58796[i])
            y_58796_5_1_in.append(y_norm_58796[i])
a=0.8
def func(x,p1,p2,t1,t2):
    return a/(np.exp((x-p1)/t1)+1)+a/(np.exp((p2-x)/t2)+1)

popt,pcov=curve_fit(func,x_58796_5_1_in,y_58796_5_1_in)

x=np.linspace(4,4.6,20)
plt.plot(x,func(x,*popt),"b", linewidth=6, linestyle='--',alpha=1)
#x=np.linspace(4.25,4.65,20)
#plt.plot(x,func(x,*popt),"r",linewidth=2,linestyle='--')
x_58796_5_2_in=[]
y_58796_5_2_in=[]
for i in range(len(new_x_58796)):
      if new_x_58796[i]>4.6 and new_x_58796[i]<4.83:
            x_58796_5_2_in.append(new_x_58796[i])
            y_58796_5_2_in.append(y_norm_58796[i])
def func(x,p1,p2):
    return 0.81/(np.exp((x-p1)/p2)+1)

popt,pcov=curve_fit(func,x_58796_5_2_in,y_58796_5_2_in)
x=np.linspace(4.6,4.83,20)
plt.plot(x,func(x,*popt), "b", linewidth=6, linestyle='--',alpha=1)
#plt.legend(loc="upper right")
# plt.legend()
#
ax3 = subplot2grid((4, 4), (2, 0), colspan=1)
ax3.errorbar(new_x_58892, y_norm_58892, yerr=error_norm_58892, fmt='.k',capsize=3)
plt.title('MJD58892(c)', x=0.1)
# plt.xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
plt.axhline(y=0, c="dimgray", ls="--")
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.ylim(-0.55, 1.4)
plt.ylabel("Normalised flux density", fontsize=25, y=1)
##################拟合
x_58892_in=[]
y_58892_in=[]
for i in range(len(new_x_58892)):
    if new_x_58892[i]>1 and new_x_58892[i]<1.42:
        x_58892_in.append(new_x_58892[i])
        y_58892_in.append(y_norm_58892[i])
def func(x,p1,p2,t1,t2):
    return 1/(np.exp((x-p1)/t1)+1)+1/(np.exp((p2-x)/t2)+1)
popt,pcov=curve_fit(func,x_58892_in,y_58892_in)
x=np.linspace(0.7,1.42,100)
plt.plot(x,func(x,*popt),"b", linewidth=6, linestyle='--',alpha=1)
#x=np.linspace(1.2,1.42,100)
#plt.plot(x,func(x,*popt),"r",linewidth=2,linestyle='--',label="Egress")
#.legend(loc="upper right")

# plt.legend()
#

ax4 = subplot2grid((4, 4), (2, 1), colspan=1)
ax4.errorbar(new_x_58900, y_norm_58900, yerr=error_norm_58900, fmt='.k',capsize=3)
plt.title('MJD58900(d)', x=0.1)
plt.xticks([0.75,1.25,1.75,2.25])
plt.axhline(y=0, c="dimgray", ls="--")
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.ylim(-0.55, 1.4)
#######拟合
x_58900_1_egress = []
y_58900_1_egress = []
for i in range(len(new_x_58900)):
    if new_x_58900[i] > 0.63 and new_x_58900[i] < 0.82:
        x_58900_1_egress.append(new_x_58900[i])
        y_58900_1_egress.append(y_norm_58900[i])


def func_58900(x, p1, p2):
    return 1 / (np.exp((p1 - x) / p2) + 1)


popt_58900, pcov_58900 = curve_fit(func_58900, x_58900_1_egress, y_58900_1_egress)
x = np.linspace(0.63, 1, 50)
plt.plot(x, func_58900(x, *popt_58900),"b", linewidth=6, linestyle='--',alpha=1)
#######
x_58900_2_ingress = []
y_58900_2_ingress = []
for i in range(len(new_x_58900)):
    if new_x_58900[i] > 1 and new_x_58900[i] < 1.5:
        x_58900_2_ingress.append(new_x_58900[i])
        y_58900_2_ingress.append(y_norm_58900[i])
# def func(x,p1,t1):
# return 1/(np.exp((x-p1)/t1)+1)#+1/(np.exp((p2-x)/t2)+1)
# popt_58900,pcov_58900=curve_fit(func,x_58900_2_ingress,y_58900_2_ingress)
p1 = 1.187
p2 = 1.409
t1 = 0.0193
t2 = 0.05518
x = np.linspace(1, 1.25, 50)
y = 1 / (np.exp((x - p1) / t1) + 1) + 1 / (np.exp((p2 - x) / t2) + 1)
plt.plot(x, y, "b", linewidth=6, linestyle='--',alpha=1)
p1 = 1.187
p2 = 1.409
t1 = 0.0193
t2 = 0.05518
x = np.linspace(1.25, 1.65, 50)
y = 1 / (np.exp((x - p1) / t1) + 1) + 1 / (np.exp((p2 - x) / t2) + 1)
plt.plot(x, y, "b", linewidth=6, linestyle='--',alpha=1)
# plt.plot(x,func(x,*popt_58900),"b",linewidth=2,linestyle='--')

########
x_58900_3_ingress = []
y_58900_3_ingress = []
for i in range(len(new_x_58900)):
    if new_x_58900[i] > 1.65 and new_x_58900[i] < 2.16:
        x_58900_3_ingress.append(new_x_58900[i])
        y_58900_3_ingress.append(y_norm_58900[i])
x_58900_3_ingress = [1.6140038759689923, 1.629506976744186, 1.6450116279069769, 1.6605155038759691, 1.6760193798449614,
                     1.6915232558139532, 1.707026356589147, 1.7225302325581393, 1.7380341085271316, 1.8000503875968992,
                     1.815553488372093, 1.8310573643410852, 1.8465612403100775, 1.8620651162790698, 1.877568992248062,
                     1.8930736434108528, 1.9085767441860466, 1.9240813953488374, 1.9395852713178292, 1.955088372093023,
                     1.9705922480620153, 1.986096899224806, 2.0016007751937983, 2.0171046511627906, 2.032608527131783,
                     2.048112403100775, 2.063615503875969, 2.0791201550387597, 2.094624031007752, 2.1101271317829458,
                     2.125631007751938, 2.141135658914729, 2.156639534883721]
y_58900_3_ingress = [1.0620488718622485, 0.8501272564396043, 0.8424136195310564, 1.0132567284229772, 0.9896148760998691,
                     0.9930286203481095, 0.8703316129732127, 0.5028576990750545, 0.06855860792130311,
                     -0.11543033128822487, 0.22062547567859447, 0.04258035966915744, -0.08359601098478553,
                     -0.07397154098963259, 0.2099203105561523, 0.007098476308398994, -0.04565482394666665,
                     -0.19063069065336258, -0.07135550305073433, -0.04973071798059041, 0.11279814230127042,
                     -0.12732791015111508, 0.11774187446818139, -0.005041368275334085, 0.1038128478737102,
                     0.01707960808402276, -0.24158167645460482, 0.10021600083467234, 0.018826361205148023,
                     -0.21237375839897685, 0.010903615727355109, 0.013015442993018947, -0.018698085608500666]


def func_58900_1(x, p1, p2):
    return 1 / (np.exp((x - p1) / p2) + 1)


popt_58900, pcov_58900 = curve_fit(func_58900_1, x_58900_3_ingress, y_58900_3_ingress)
x = np.linspace(1.65, 2.16, 50)
plt.plot(x, func_58900_1(x, *popt_58900), "b", linewidth=6, linestyle='--',alpha=1)

#plt.legend(loc="upper right")
#

ax5 = subplot2grid((4, 4), (2, 2), colspan=1)
ax5.errorbar(new_x_58974, y_norm_58974, yerr=error_norm_58974, fmt='.k',capsize=3)
plt.title('MJD58974(e)', x=0.1)
# plt.xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
plt.axhline(y=0, c="dimgray", ls="--")
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.ylim(-0.55, 1.4)
###########拟合
###MJD58974
x_58974_1_egress=[]
y_58974_1_egress=[]
for i in range(len(new_x_58974)):
    if new_x_58974[i]>0.01 and new_x_58974[i]<0.24:
        x_58974_1_egress.append(new_x_58974[i])
        y_58974_1_egress.append(y_norm_58974[i])
def func_58974(x,p1,p2):
    return 1/(np.exp((p1-x)/p2)+1)
popt_58974,pcov_58974=curve_fit(func_58974,x_58974_1_egress,y_58974_1_egress)
x=np.linspace(0.01,0.36,50)
plt.plot(x,func_58974(x,*popt_58974),"b", linewidth=6, linestyle='--',alpha=1)
####
x_58974_2_ingress=[]
y_58974_2_ingress=[]
for i in range(len(new_x_58974)):
    if new_x_58974[i]>0.36 and new_x_58974[i]<0.6:
        x_58974_2_ingress.append(new_x_58974[i])
        y_58974_2_ingress.append(y_norm_58974[i])

def func_58974(x,p1,p2):
    return 1/(np.exp((x-p1)/p2)+1)
popt_58974,pcov_58974=curve_fit(func_58974,x_58974_2_ingress,y_58974_2_ingress)
x=np.linspace(0.36,0.55,50)
plt.plot(x,func_58974(x,*popt_58974),"b", linewidth=6, linestyle='--',alpha=1)
##
x_58974_2_egress=[]
y_58974_2_egress=[]
for i in range(len(new_x_58974)):
    if new_x_58974[i]>0.55 and new_x_58974[i]<0.81:
        x_58974_2_egress.append(new_x_58974[i])
        y_58974_2_egress.append(y_norm_58974[i])
def func_58974(x,p1,p2):
    return 0.75/(np.exp((p1-x)/p2)+1)
popt_58974,pcov_58974=curve_fit(func_58974,x_58974_2_egress,y_58974_2_egress)
x=np.linspace(0.55,0.81,50)
plt.plot(x,func_58974(x,*popt_58974),"b", linewidth=6, linestyle='--',alpha=1)
#plt.legend(loc="upper right")
#

ax6 = subplot2grid((4, 4), (2, 3), colspan=1)
ax6.errorbar(new_x_59017_2, y_norm_59017, yerr=error_norm_59017, fmt='.k',capsize=3)
plt.title('MJD59017(f)', x=0.1)
# plt.xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
plt.axhline(y=0, c="dimgray", ls="--")
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.ylim(-0.6, 1.4)
########拟合
###mjd59017
x_59017_1_egress=[]
y_59017_1_egress=[]
for i in range(len(new_x_59017_2)):
    if new_x_59017_2[i]>0.53 and new_x_59017_2[i]<1.05:
        x_59017_1_egress.append(new_x_59017_2[i])
        y_59017_1_egress.append(y_norm_59017[i])
def func(x,p1,p2):
    return 1/(np.exp((p1-x)/p2)+1)
popt,pcov=curve_fit(func,x_59017_1_egress,y_59017_1_egress)
x=np.linspace(0.53,1.08,50)
plt.plot(x,func(x,*popt),"b", linewidth=6, linestyle='--',alpha=1)
#########
x_59017_2_ingress=[]
y_59017_2_ingress=[]
for i in range(len(new_x_59017_2)):
    if new_x_59017_2[i]>1.08 and new_x_59017_2[i]<1.41:
        x_59017_2_ingress.append(new_x_59017_2[i])
        y_59017_2_ingress.append(y_norm_59017[i])
a=1
def func(x,p1,p2,t1,t2):
    return a/(np.exp((x-p1)/t1)+1)+a/(np.exp((p2-x)/t2)+1)
popt,pcov=curve_fit(func,x_59017_2_ingress,y_59017_2_ingress)
x=np.linspace(1.08,1.41,60)
plt.plot(x,func(x,*popt),"b", linewidth=6, linestyle='--',alpha=1)
#x=np.linspace(1.31,1.41,50)
#plt.plot(x,func(x,*popt),"r",linewidth=2,linestyle='--',label="Egress")
#plt.legend(loc="upper right")
#


ax7 = subplot2grid((4, 4), (3, 0), colspan=1)
ax7.errorbar(new_x_59018, y_norm_59018, yerr=error_norm_59018, fmt='.k',capsize=3)
plt.title('MJD59018(g)', x=0.1)
# plt.([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
plt.axhline(y=0, c="dimgray", ls="--")
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.ylim(-0.6, 1.4)

#######拟合
###MJD59018

x_59018_2_ingress=[]
y_59018_2_ingress=[]
for i in range(len(new_x_59018)):
    if new_x_59018[i]>1.22 and new_x_59018[i]<1.4:
        x_59018_2_ingress.append(new_x_59018[i])
        y_59018_2_ingress.append(y_norm_59018[i])

#x_59018_2_ingress=[1.23372790700000,1.25698372100000,1.27830077500000,1.30155736400000,1.32520077500000,1.34806976700000,1.36744961200000,1.39651938000000,1.41589845000000]
#y_59018_2_ingress=[1,0.724277068000000,0.865501009000000,0.500672495000000,0.518157364000000,0.248789509000000,0.313046402000000,0.164088769000000,0.249899126000000]
def func_59018(x,p1,p2):
    return 1/(np.exp((x-p1)/p2)+1)

popt_59018,pcov_59018=curve_fit(func_59018,x_59018_2_ingress,y_59018_2_ingress)
x=np.linspace(1.22,1.4,50)
plt.plot(x,func_59018(x,*popt_59018),"b", linewidth=6, linestyle='--',alpha=1)
####
x_59018_1_egress=[]
y_59018_1_egress=[]
for i in range(len(new_x_59018)):
    if new_x_59018[i]>0.5 and new_x_59018[i]<1:
        x_59018_1_egress.append(new_x_59018[i])
        y_59018_1_egress.append(y_norm_59018[i])
def func_59018(x,p1,p2):
    return 1/(np.exp((p1-x)/p2)+1)
popt_59018,pcov_59018=curve_fit(func_59018,x_59018_1_egress,y_59018_1_egress)
x=np.linspace(0.5,1.2,50)
plt.plot(x,func_59018(x,*popt_59018),"b", linewidth=6, linestyle='--',alpha=1)
#########

#plt.legend(loc="upper right")
###########
file1=r"C:\\Users\\24566\\Desktop\\DM and flux\\parkes\\mjd59188\\mjd59188.calib.2minpF.txt"
x_59188_1=[]
y_59188_1=[]
error_59188_1=[]
for i in range(9,63):
    a=linecache.getline(file1,i)
    a1=a.split()
    y_59188_1.append(float(a1[4]))
    error_59188_1.append(float(a1[5]))
    x_59188_1.append(float(a1[2]))
for i in range(len(x_59188_1)):
    x_59188_1[i]=x_59188_1[i]/129+0.569

#error2_59358_1=[]
#for i in range(len(error_59358_1)):
    #error2_59358_1.append(error_59358_1[i]/y_59358_1[i])
new_x_59188_1=[]
new_y_59188_1=[]
new_error_59188_1=[]
for i in range(len(error_59188_1)):
    #if error2_59017_2[i]<1:
        new_x_59188_1.append(x_59188_1[i])
        new_y_59188_1.append(y_59188_1[i])
        new_error_59188_1.append(error_59188_1[i])
y_u=np.mean(new_y_59188_1)+np.std(new_y_59188_1)
y_s=[]
for i in range(len(new_y_59188_1)):
    if new_y_59188_1[i]>y_u:
        y_s.append(new_y_59188_1[i])
y_max=np.mean(y_s)

y_norm_59188=[]
error_norm_59188=[]
for i in range(len(new_y_59188_1)):
    y_norm_59188.append(new_y_59188_1[i]/(y_max))
    error_norm_59188.append(new_error_59188_1[i]/(y_max))

ax8 = subplot2grid((4, 4), (3, 1), colspan=1)
ax8.errorbar(new_x_59188_1,y_norm_59188,yerr=error_norm_59188, fmt='.k',capsize=3)
plt.title('MJD59188(h)', x=0.1)
# plt.([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
plt.axhline(y=0, c="dimgray", ls="--")
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.ylim(-0.6, 1.4)
plt.xlabel("Orbital phase", fontsize=25)

x_59188_1_egress=[]
y_59188_1_egress=[]
for i in range(len(new_x_59188_1)):
    if new_x_59188_1[i]>0.7 and new_x_59188_1[i]<1.4:
        x_59188_1_egress.append(new_x_59188_1[i])
        y_59188_1_egress.append(y_norm_59188[i])

a=1
def func_59188(x,p1,p2,t1,t2,a):
    return a/(np.exp((x-p1)/t1)+1)+a/(np.exp((p2-x)/t2)+1)
popt_59188,pcov_59188=curve_fit(func_59188,x_59188_1_egress,y_59188_1_egress)
x=np.linspace(0.65,1.38,50)
plt.plot(x,func_59188(x,*popt_59188),"b", linewidth=6, linestyle='--',alpha=1)
#print(popt_59188)

#########59358
file1 = r"C:\\Users\\24566\\Desktop\\DM and flux\\GMRT\\mjd59358.calib.2minpF.txt"
x_59358_1 = []
y_59358_1 = []
error_59358_1 = []
for i in range(9, 139):
    a = linecache.getline(file1, i)
    a1 = a.split()
    y_59358_1.append(float(a1[4]))
    error_59358_1.append(float(a1[5]))
    x_59358_1.append(float(a1[2]))
for i in range(len(x_59358_1)):
    x_59358_1[i] = x_59358_1[i] / 129 + 0.19
print()
# error2_59358_1=[]
# for i in range(len(error_59358_1)):
# error2_59358_1.append(error_59358_1[i]/y_59358_1[i])
new_x_59358_1 = []
new_y_59358_1 = []
new_error_59358_1 = []
for i in range(len(error_59358_1)):
    # if error2_59017_2[i]<1:
    new_x_59358_1.append(x_59358_1[i])
    new_y_59358_1.append(y_59358_1[i])
    new_error_59358_1.append(error_59358_1[i])
y_norm_59358 = []
error_norm_59358 = []
for i in range(len(new_y_59358_1)):
    y_norm_59358.append(new_y_59358_1[i] / (np.max(new_y_59358_1) * 0.865))
    error_norm_59358.append(new_error_59358_1[i] / (np.max(new_y_59358_1) * 0.865))

ax9 = subplot2grid((4, 4), (3, 2), colspan=1)
ax9.errorbar(new_x_59358_1,y_norm_59358,yerr=error_norm_59358, fmt='.k',capsize=3)
plt.title('MJD59358(i)', x=0.1)
plt.xticks([0.25,0.75,1.25,1.75,2.25])
plt.axhline(y=0, c="dimgray", ls="--")
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.ylim(-0.6, 1.4)

########拟合
x_59358_1_ingress = []
y_59358_1_ingress = []
for i in range(len(new_x_59358_1)):
    if new_x_59358_1[i] > 0.25 and new_x_59358_1[i] < 0.5:
        x_59358_1_ingress.append(new_x_59358_1[i])
        y_59358_1_ingress.append(y_norm_59358[i])
a = 1


def func(x, p1, t1):
    return a / (np.exp((x - p1) / t1) + 1)  # +a/(np.exp((p2-x)/t2)+1)


popt, pcov = curve_fit(func, x_59358_1_ingress, y_59358_1_ingress)
x = np.linspace(0.25, 0.5, 50)
plt.plot(x, func(x, *popt), "b", linewidth=6, linestyle='--', alpha=1)
print(popt)
x_59358_1_ingress = []
y_59358_1_ingress = []
for i in range(len(new_x_59358_1)):
    if new_x_59358_1[i] > 0.3 and new_x_59358_1[i] < 0.67:
        x_59358_1_ingress.append(new_x_59358_1[i])
        y_59358_1_ingress.append(y_norm_59358[i])
a = 1


def func(x, p1, t1):
    return a / (np.exp((p1 - x) / t1) + 1)  # +a/(np.exp((p2-x)/t2)+1)


popt, pcov = curve_fit(func, x_59358_1_ingress, y_59358_1_ingress)
x = np.linspace(0.45, 0.66, 50)
plt.plot(x, func(x, *popt), "b", linewidth=6, linestyle='--', alpha=1)

########

x = np.linspace(0.665, 1.05, 50)
y = []
for i in range(len(x)):
    y.append(1 / (np.exp((x[i] - 0.6831) / 0.0082) + 1) + 1 / (np.exp((0.8053 - x[i]) / 0.05586) + 1))
plt.plot(x, y, "b", linewidth=6, linestyle='--', alpha=1)
##########
x_59358_1_ingress = []
y_59358_1_ingress = []
for i in range(len(new_x_59358_1)):
    if new_x_59358_1[i] > 1.15 and new_x_59358_1[i] < 1.5:
        x_59358_1_ingress.append(new_x_59358_1[i])
        y_59358_1_ingress.append(y_norm_59358[i])
a = 1


def func(x, p1, t1):
    return a / (np.exp((x - p1) / t1) + 1)  # +a/(np.exp((p2-x)/t2)+1)


popt, pcov = curve_fit(func, x_59358_1_ingress, y_59358_1_ingress)
x = np.linspace(1.05, 1.5, 50)
plt.plot(x, func(x, *popt), "b", linewidth=6, linestyle='--', alpha=1)
#############
x_59358_1_ingress = []
y_59358_1_ingress = []
for i in range(len(new_x_59358_1)):
    if new_x_59358_1[i] > 1.4 and new_x_59358_1[i] < 1.55:
        x_59358_1_ingress.append(new_x_59358_1[i])
        y_59358_1_ingress.append(y_norm_59358[i])
a = 0.75


def func(x, p1, t1):
    return a / (np.exp((p1 - x) / t1) + 1)  # +a/(np.exp((p2-x)/t2)+1)


popt, pcov = curve_fit(func, x_59358_1_ingress, y_59358_1_ingress)
x = np.linspace(1.4, 1.55, 50)
plt.plot(x, func(x, *popt), "b", linewidth=6, linestyle='--', alpha=1)
########
x_59358_1_ingress = []
y_59358_1_ingress = []
for i in range(len(new_x_59358_1)):
    if new_x_59358_1[i] > 1.55 and new_x_59358_1[i] < 1.75:
        x_59358_1_ingress.append(new_x_59358_1[i])
        y_59358_1_ingress.append(y_norm_59358[i])
a = 0.75


def func(x, p1, t1):
    return a / (np.exp((x - p1) / t1) + 1)  # +a/(np.exp((p2-x)/t2)+1)


popt, pcov = curve_fit(func, x_59358_1_ingress, y_59358_1_ingress)
x = np.linspace(1.55, 2.2, 50)
plt.plot(x, func(x, *popt), "b", linewidth=6, linestyle='--', alpha=1)

plt.savefig('C:\\Users\\24566\\Desktop\\Flux density',bbox_inches = 'tight',format='pdf')
plt.show()
