#Made by Hongyu Gong 2019-09-05
#For PRESTO_master edition
#Use "readfile filename" to get informations about the data
#and edit the parameters in "Parameters" corresponding to your need
##########################   Parameters  ###############################
#Your Presto route
#prestor = '$PRESTO'
#Project base name
#basename = '1131415232'
#Data route
#droute = '/ssd/Pulsar/'+basename+'/fits/'
#Data file formate (.fits, .fil)
#lastname = '.fits'
#Range of DMs in searching
loDM, hiDM = 1.0, 250.0
#Choose to use normal search (zmax=0) or accelerate search (zmax=?int)
#zmax=200
#Length of time used in rfifind (float in seconds)
#rfitime=12.0
#Would you like to remove the 0-dispersion signal? (birds=1, Yes; birds=0, No)
#birds=0
#Number of subbands used in de-dispersion
#numsubbands=512
#Number of harmonics in accelsearch (1,2,4,8,16,32)
numharm=16
#Zap channels (zc=0 if no gaps among channels; zc=1, sign zap channels)
#zc=1
zapchan = '0:19,108:127,128:147,236:255,256:275,364:383,384:403,492:511,512:531,620:639,640:659,748:767,768:787,876:895,896:915,1004:1023,1024:1043,1132:1151,1152:1171,1260:1279,1280:1299,1388:1407,1408:1427,1516:1535,1536:1555,1644:1663,1664:1683,1772:1791,1792:1811,1900:1919,1920:1939,2028:2047,2048:2067,2156:2175,2176:2195,2284:2303,2304:2323,2412:2431,2432:2451,2540:2559,2560:2579,2668:2687,2688:2707,2796:2815,2816:2835,2924:2943,2944:2963,3052:3071'
#Fold all candidates (fcand=1) or not (fcand=0)
#fcand=0
#Use single pulse search (singlep=1) or not (singlep=0)
#singlep=1
#ncpus
#ncpus = 28
########################################################################


import os
import glob
import sys
import re
from numpy import *
#from Pgplot import *
import getopt, sys
import astropy.io as fits
import argparse
import math
import numpy as np
import matplotlib.pyplot as plt
import time
import datetime


def stcmd(cmd):
    #print "%s"%cmd
    os.system("echo $(date '+%Y-%m-%d %H:%M:%S')>> command.txt")
    os.system("echo "+cmd+">> command.txt")
    os.system(cmd)

#  DDplan for midle and high frequency telescope####################################

class observation:
    def __init__(self, dt, f_ctr, BW, numchan, cDM):
        # dt in sec, f_ctr and in MHz
        self.dt = dt
        self.f_ctr = f_ctr
        self.BW = BW
        self.numchan = numchan
        self.chanwidth = BW/numchan
        self.cDM = cDM
    def guess_dDM(self, DM):
        """
        guess_dDM(self, DM):
            Choose a reasonable dDM by setting the maximum smearing across the
                'BW' to equal the sampling time 'dt'.
        """
        return self.dt*0.0001205*self.f_ctr**3.0/(0.5*self.BW)
class dedisp_method:
    def __init__(self, obs, downsamp, loDM, hiDM, dDM, numDMs=0,
                 numsub=0, smearfact=2.0):
        self.obs = obs
        self.downsamp = downsamp
        self.loDM = loDM
        self.dDM = dDM
        self.numsub = numsub
        self.BW_smearing = BW_smear(dDM, self.obs.BW, self.obs.f_ctr)
        self.numprepsub = 0
        if (numsub):  # Calculate the maximum subband smearing we can handle
            DMs_per_prepsub = 2
            while(1):
                next_dsubDM = (DMs_per_prepsub+2) * dDM
                next_ss = subband_smear(next_dsubDM, numsub, self.obs.BW, self.obs.f_ctr)
                # The 0.8 is a small fudge factor to make sure that the subband
                # smearing is always the smallest contribution
                if (next_ss > 0.8*min(self.BW_smearing, 1000.0*obs.dt*downsamp)):
                    self.dsubDM = DMs_per_prepsub*dDM
                    self.DMs_per_prepsub = DMs_per_prepsub
                    self.sub_smearing = subband_smear(self.dsubDM, numsub,
                                                      self.obs.BW, self.obs.f_ctr)
                    break
                DMs_per_prepsub += 2
        else:
            self.dsubDM = dDM
            self.sub_smearing = subband_smear(self.dsubDM, numsub, self.obs.BW, self.obs.f_ctr)
        # Calculate the nominal DM to move to the next method
        cross_DM = self.DM_for_smearfact(smearfact)
        if (cross_DM > hiDM):
            cross_DM = hiDM
        if (numDMs==0):
            self.numDMs = int(ceil((cross_DM-loDM)/dDM))
            if (numsub):
                self.numprepsub = int(ceil(self.numDMs*dDM / self.dsubDM))
                self.numDMs = self.numprepsub * DMs_per_prepsub
        else:
            self.numDMs = numDMs
        self.hiDM = loDM + self.numDMs*dDM
        self.DMs = arange(self.numDMs, dtype='d')*dDM + loDM
    def chan_smear(self, DM):
        """
        Return the smearing (in ms) in each channel at the specified DM
        """
        try:
            DM = where(DM-cDM==0.0, cDM+self.dDM/2.0, DM)
        except TypeError:
            if (DM-cDM==0.0): DM = cDM+self.dDM/2.0
        return dm_smear(DM, self.obs.chanwidth, self.obs.f_ctr, self.obs.cDM)
    def total_smear(self, DM):
        """
        Return the total smearing in ms due to the sampling rate,
        the smearing over each channel, the smearing over each subband
        (if numsub > 0) and the smearing over the full BW assuming the
        worst-case DM error.
        """
        return sqrt((1000.0*self.obs.dt)**2.0 +
                    (1000.0*self.obs.dt*self.downsamp)**2.0 +
                    self.BW_smearing**2.0 +
                    self.sub_smearing**2.0 +
                    self.chan_smear(DM)**2.0)
    def DM_for_smearfact(self, smearfact):
        """
        Return the DM where the smearing in a single channel is a factor smearfact
        larger than all the other smearing causes combined.
        """
        other_smear = sqrt((1000.0*self.obs.dt)**2.0 +
                           (1000.0*self.obs.dt*self.downsamp)**2.0 +
                           self.BW_smearing**2.0 +
                           self.sub_smearing**2.0)
        return smearfact*0.001*other_smear/self.obs.chanwidth*0.0001205*self.obs.f_ctr**3.0 + self.obs.cDM
    def DM_for_newparams(self, dDM, downsamp):
        """
        Return the DM where the smearing in a single channel is causes the same smearing
        as the effects of the new dosnsampling rate and dDM.
        """
        other_smear = sqrt((1000.0*self.obs.dt)**2.0 +
                           (1000.0*self.obs.dt*downsamp)**2.0 +
                           BW_smear(dDM, self.obs.BW, self.obs.f_ctr)**2.0 +
                           self.sub_smearing**2.0)
        return 0.001*other_smear/self.obs.chanwidth*0.0001205*self.obs.f_ctr**3.0
    def plot(self, work_fract):
        DMspan = self.DMs[-1]-self.DMs[0]
        loDM  = self.DMs[0]  + DMspan*0.02
        hiDM  = self.DMs[-1] - DMspan*0.02
        midDM = self.DMs[0]  + DMspan*0.5
        dt_ms = 1000.0*self.obs.dt*self.downsamp
        plotxy(log10(self.total_smear(self.DMs)), self.DMs, width=4)
        ppgplot.pgslw(1)
        ppgplot.pgsch(0.8)
        ppgplot.pgptxt(midDM, log10(1.1*self.total_smear(midDM)), 90.0, 0.0,
                       "%d (%.1f%%)" % (self.numDMs, 100.0*work_fract))
        # Sample time
        plotxy(log10(zeros(self.numDMs)+dt_ms), self.DMs, width=1, color='green')
        ppgplot.pgptxt(loDM, log10(0.85*dt_ms), 0.0, 0.0,
                       "%g" % dt_ms)
        # DM stepsize smearing
        plotxy(log10(zeros(self.numDMs)+self.BW_smearing), self.DMs, width=1, color='red')
        ppgplot.pgptxt(hiDM, log10(0.85*self.BW_smearing), 0.0, 1.0,
                       "%g" % self.dDM)
        # channel smearing
        plotxy(log10(self.chan_smear(self.DMs)), self.DMs, width=1, color='blue')
        # subband smearing
        if (self.numsub):
            plotxy(log10(zeros(self.numDMs)+self.sub_smearing),
                   self.DMs, width=1, color='purple')
            ppgplot.pgptxt(midDM, log10(0.85*self.sub_smearing), 0.0, 0.5,
                           "%g (%d)" % (self.dsubDM, self.numprepsub))
        ppgplot.pgsci(1)
        ppgplot.pgsch(1.0)
    def __str__(self):
        if (self.numsub):
            dDMs.append(self.dDM)
            dsubDMs.append(self.dsubDM)
            downsamps.append(self.downsamp)
            subcalls.append(self.numprepsub)
            startDMs.append(self.loDM)
            dmspercalls.append(self.DMs_per_prepsub)

            return "%9.3f  %9.3f  %6.2f    %4d  %6.2f  %6d  %6d  %6d " % \
                   (self.loDM, self.hiDM, self.dDM, self.downsamp, self.dsubDM,
                    self.numDMs, self.DMs_per_prepsub, self.numprepsub)
        else:
            return "%9.3f  %9.3f  %6.2f    %4d  %6d" % \
                   (self.loDM, self.hiDM, self.dDM, self.downsamp, self.numDMs)

def dm_smear(DM, BW, f_ctr, cDM=0.0):
    """
    dm_smear(DM, BW, f_ctr, cDM=0.0):
        Return the smearing in ms caused by a 'DM' over a bandwidth
        of 'BW' MHz centered at 'f_ctr' MHz.
    """
    return 1000.0*fabs(DM-cDM)*BW/(0.0001205*f_ctr**3.0)

def BW_smear(DMstep, BW, f_ctr):
    """
    BW_smear(DMstep, BW, f_ctr):
        Return the smearing in ms caused by a search using a DM stepsize of
        'DMstep' over a bandwidth of 'BW' MHz centered at 'f_ctr' MHz.
    """
    maxDMerror = 0.5*DMstep
    return dm_smear(maxDMerror, BW, f_ctr)

def guess_DMstep(dt, BW, f_ctr):
    """
    guess_DMstep(dt, BW, f_ctr):
        Choose a reasonable DMstep by setting the maximum smearing across the
        'BW' to equal the sampling time 'dt'.
    """
    return dt*0.0001205*f_ctr**3.0/(0.5*BW)

def subband_smear(subDMstep, numsub, BW, f_ctr):
    """
    subband_smear(subDMstep, numsub, BW, f_ctr):
        Return the smearing in ms caused by a search using a subband
        DM stepsize of 'subDMstep' over a total bandwidth of 'BW' MHz
        centered at 'f_ctr' MHz, and having numsub subbands.
    """
    if (numsub==0): return 0.0
    subBW = BW/numsub
    maxsubDMerror = 0.5*subDMstep
    return dm_smear(maxsubDMerror, subBW, f_ctr)

def total_smear(DM, DMstep, dt, f_ctr, BW, numchan, subDMstep, cohdm=0.0, numsub=0):
    """
    total_smear(DM, DMstep, dt, f_ctr, BW, numchan, subDMstep, cohdm=0.0, numsub=0):
        Return the total smearing in ms due to the sampling rate,
        the smearing over each channel, the smearing over each subband
        (if numsub > 0) and the smearing over the full BW assuming the
        worst-case DM error.
    """
    return sqrt(2 * (1000.0*dt)**2.0 +
                dm_smear(DM, BW/numchan, f_ctr, cohdm)**2.0 +
                subband_smear(subDMstep, numsub, BW, f_ctr)**2.0 +
                BW_smear(DMstep, BW, f_ctr)**2.0)

def dm_steps(loDM, hiDM, obs, cohdm=0.0, numsub=0, ok_smearing=0.0, device="/XWIN"):

    """
    dm_steps(loDM, hiDM, obs, cohdm=0.0, numsub=0, ok_smearing=0.0):
        Return the optimal DM stepsizes (and subband DM stepsizes if
        numsub>0) to keep the total smearing below 'ok_smearing' (in ms),
        for the DMs between loDM and hiDM.  If 'ok_smearing'=0.0, then
        use the best values based only on the data.
    """
    # Allowable DM stepsizes
    allow_dDMs = [0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0,
                  2.0, 3.0, 5.0, 10.0, 20.0, 30.0, 50.0, 100.0, 200.0, 300.0]
    # Allowable number of downsampling factors
    allow_downsamps = [1, 2, 4, 8, 16, 32, 64, 128, 256]

    # Initial values
    index_downsamps = index_dDMs = 0
    downsamp = allow_downsamps[index_downsamps]
    dDM = allow_dDMs[index_dDMs]
    dtms = 1000.0*obs.dt

    # Fudge factor that "softens" the boundary defining
    # if 2 time scales are equal or not
    ff = 1.2

    # This is the array that will hold the de-dispersion plans
    methods = []

    # Minimum possible smearing
    min_tot_smearing = total_smear(loDM+0.5*dDM, dDM, obs.dt, obs.f_ctr,
                                   obs.BW, obs.numchan, allow_dDMs[0], cohdm, 0)
    # Minimum channel smearing
    min_chan_smearing = dm_smear(linspace(loDM, hiDM, 10000),
                                 obs.chanwidth, obs.f_ctr, cohdm).min()
    # Minimum smearing across the obs.BW
    min_BW_smearing = BW_smear(dDM, obs.BW, obs.f_ctr)

    print
    print "Minimum total smearing     : %.3g ms" % min_tot_smearing
    print "--------------------------------------------"
    print "Minimum channel smearing   : %.3g ms" % min_chan_smearing
    print "Minimum smearing across BW : %.3g ms" % min_BW_smearing
    print "Minimum sample time        : %.3g ms" % dtms
    print

    ok_smearing = max([ok_smearing, min_chan_smearing, min_BW_smearing, dtms])
    print "Setting the new 'best' resolution to : %.3g ms" % ok_smearing

    # See if the data is too high time resolution for our needs
    if (ff*min_chan_smearing > dtms or
        ok_smearing > dtms):
        if (ok_smearing > ff*min_chan_smearing):
            print "   Note: ok_smearing > dt (i.e. data is higher resolution than needed)"
            okval = ok_smearing
        else:
            print "   Note: min_chan_smearing > dt (i.e. data is higher resolution than needed)"
            okval = ff*min_chan_smearing

        while (dtms*allow_downsamps[index_downsamps+1] < okval):
            index_downsamps += 1
        downsamp = allow_downsamps[index_downsamps]
        print "         New dt is %d x %.12g ms = %.12g ms" % \
              (downsamp, dtms, dtms*downsamp)

    # Calculate the appropriate initial dDM 
    dDM = guess_DMstep(obs.dt*downsamp, obs.BW, obs.f_ctr)
    print "Best guess for optimal initial dDM is %.3f" % dDM
    while (allow_dDMs[index_dDMs+1] < ff*dDM):
        index_dDMs += 1

    # Create the first method
    methods = [dedisp_method(obs, downsamp, loDM, hiDM,
                             allow_dDMs[index_dDMs], numsub=numsub)]
    numDMs = [methods[-1].numDMs]

    # Calculate the next methods
    while(methods[-1].hiDM < hiDM):

        # Determine the new downsample factor
        index_downsamps += 1
        downsamp = allow_downsamps[index_downsamps]
        eff_dt = dtms*downsamp
        # Determine the new DM step
        while (BW_smear(allow_dDMs[index_dDMs+1], obs.BW, obs.f_ctr) < ff*eff_dt):
            index_dDMs += 1
        dDM = allow_dDMs[index_dDMs]

        # Get the next method
        methods.append(dedisp_method(obs, downsamp, methods[-1].hiDM,
                                     hiDM, dDM, numsub=numsub))
        numDMs.append(methods[-1].numDMs)

    # Calculate the DMs to search and the smearing at each
    total_numDMs = sum(numDMs)
    DMs = zeros(total_numDMs, dtype='d')
    total_smears = zeros(total_numDMs, dtype='d')

    # Calculate the DMs and optimal smearing for all the DMs
    for ii, offset in enumerate(add.accumulate([0]+numDMs[:-1])):
        DMs[offset:offset+numDMs[ii]] = methods[ii].DMs
        total_smears[offset:offset+numDMs[ii]] = methods[ii].total_smear(methods[ii].DMs)

    # Calculate the predicted amount of time that will be spent in searching
    # this batch of DMs as a fraction of the total
    work_fracts = [meth.numDMs/float(meth.downsamp) for meth in methods]
    work_fracts = asarray(work_fracts)/sum(work_fracts)

    # The optimal smearing
    tot_smear = total_smear(DMs, allow_dDMs[0], obs.dt, obs.f_ctr,
                            obs.BW, obs.numchan, allow_dDMs[0], cohdm, 0)
    
    if (numsub):
        print "\n  Low DM    High DM     dDM  DownSamp  dsubDM   #DMs  DMs/call  calls  WorkFract"
    else:
        print "\n  Low DM    High DM     dDM  DownSamp   #DMs  WorkFract"


#    ij=0
#    print methods
    for method, fract in zip(methods, work_fracts):
        print method, "  %.4g" % fract

#        method.plot(fract)
    print "\n\n"
#    closeplot()

def usage():
    print """
usage:  DDplan.py [options]
  [-h, --help]                    : Display this help
  [-o outfile, --outfile=outfile] : Output .eps plot file (default is xwin)
  [-l loDM, --loDM=loDM]          : Low DM to search   (default = 0 pc cm-3)
  [-d hiDM, --hiDM=HIDM]          : High DM to search  (default = 1000 pc cm-3)
  [-f fctr, --fctr=fctr]          : Center frequency   (default = 1400MHz)
  [-b BW, --bw=bandwidth]         : Bandwidth in MHz   (default = 300MHz)
  [-n #chan, --numchan=#chan]     : Number of channels (default = 1024)
  [-c cDM, --cohdm=cDM]           : Coherent DM in each chan  (default = 0.0)
  [-t dt, --dt=dt]                : Sample time (s)    (default = 0.000064 s)
  [-s subbands, --subbands=nsub]  : Number of subbands (default = #chan) 
  [-r resolution, --res=res]      : Acceptable time resolution (ms)
  The program generates a good plan for de-dispersing raw data.  It
  trades a small amount of sensitivity in order to save computation costs.

"""
################################################################################
#DDplan for low frequency telescopes############################################
#! /usr/bin/env python3



def dd_plan(centrefreq, bandwidth, nfreqchan, timeres, lowDM, highDM, min_DM_step=0.02):
    """
    Work out the dedisperion plan

    Parameters
    ----------
    centrefreq: float
        The center frequency of the observation in MHz
    bandwidth: float
        The bandwidth of the observation in MHz
    nfreqchan: int
        The number of frequency channels
    timeres: float
        The time resolution of the observation in ms
    lowDM: float
        The lowest dispersion measure
    highDM: float
        The highest dispersion measure
    min_DM_step: float
        Will overwrite the minimum DM step with this value

    Returns
    -------
    DD_plan_array: list list
        dedispersion plan format:
        [[low_DM, high_DM, DM_step, nDM_step, timeres, downsample]]
    """

    DD_plan_array = []
    freqres = bandwidth / float(nfreqchan)
    previous_DM = lowDM

    #number of time samples smeared over before moving to next D_dm
    smear_fact = 3.

    #Loop until you've made a hit your range max
    D_DM = 0.
    downsample = 1
    while D_DM < round(highDM, 2):
        #calculate the DM where the current time resolution equals the
        #dispersion in a frequency channel (a bit of an overkill)

        #Dm smear over a frequency channel
        dm_smear = previous_DM * freqres * 8.3 * 10.**6 / centrefreq**3
        total_smear = math.sqrt(timeres**2 +
                                dm_smear**2)


        D_DM = smear_fact * timeres * centrefreq**3 /\
               (8.3 * 10.**6 * freqres)

        #difference in DM that will double the effective width (eq 6.4 of pulsar handbook)
        #TODO make this more robust
        #DM_step = math.sqrt( (2.*timeres)**2 - timeres**2 )/\
        #          (8.3 * 10**6 * bandwidth / centrefreq**3)
        DM_step = smear_fact * total_smear * centrefreq**3 /\
                  (8.3 * 10.**6 * 0.5 * bandwidth)


        #round to nearest 0.01
        DM_step = round(DM_step, 2)
        if DM_step < min_DM_step:
            #set DM to 0.01 as a zero DM doesn't make sense
            DM_step = min_DM_step


        if D_DM > highDM:
            #last one so range from to max
            D_DM = highDM
        #range from last to new
        D_DM = round(D_DM, 2)
        nDM_step = int((D_DM - previous_DM) / DM_step)
        if D_DM > lowDM:
            DD_plan_array.append([ previous_DM, D_DM, DM_step, nDM_step, timeres, downsample ])
            previous_DM = D_DM

        #Double time res to account for incoherent dedispersion
        timeres *= 2.
        downsample *= 2

    return DD_plan_array

def decimal(number):
    temp = number.split(".")[1]
    decimalFacter = len(temp)
    return decimalFacter-1 

def creatbirds(birds):
  birdsf = open(birds,"w")
  birdsf.write("#Freq\tWidth\t#harm\tgrow?\tbarry?\n")
  birdsf.close()

def writebirds(basename):
    f1 = basename+"_topo_DM0.00_ACCEL_0"
    h1 = open(f1,'r')
    info = h1.readlines()
    usefulInfo = []

    f2 = basename+".birds"
    m1 = open(f2,"a+")
    for temp in info:
        if temp == "\n":
            break
        usefulInfo.append(temp)

    i = 0
    for temp in usefulInfo:
        if i >= 3:  
            tempInfo = temp.split("\n")
            tempInfo2 = tempInfo[0].split(" ")
            finalInfo = []
            for character in tempInfo2:
                if character != "":
                    finalInfo.append(character)
            parFreq = finalInfo[6].split("(")[0]
            deci = decimal(parFreq)
        
            parWidth = "0."+"0"*deci+finalInfo[6].split("(")[1].split(")")[0]
            parHarm = finalInfo[4]
            parGrow = "0"
            parBarry = "0"
            print parFreq +"  "+parWidth+"  "+finalInfo[4]
            m1.write(parFreq+"\t"+parWidth+"\t"+parHarm+"\t"+parGrow+"\t"+parBarry+"\n")       
        i += 1
    h1.close()
    m1.close()

def findInfoInTxt(txtName,keyInfo):
    h2 = open(txtName,'r')
    info = h2.readlines()

    for temp in info:
        tempPar = temp.split("=")
        #print "tempPar:"
        #print tempPar        
        if len(tempPar) != 1 and tempPar[0].strip() == keyInfo:
             tempPar2 = tempPar[1].strip()
             tempPar3 = tempPar2.split("\n")[0]
             #print tempPar3
             return tempPar3
             h2.close()
    print "**********************************"
    print "we can not find the SpectraNumber"
    h2.close()

################################## main #######################################
if __name__ == '__main__':
  with open('parameters.txt', 'r') as pars:
          parall = pars.read().splitlines()
          prestor = parall[0]
          basename = parall[1]
          lastname = parall[2]
          zmax = int(parall[3])
          rfitime = float(parall[4])
          birds = int(parall[5])
          numsubbands = int(parall[6])
          zc = int(parall[7])
          fcand = int(parall[8])
          singlep = int(parall[9])
          ncpus = int(parall[10])
          droute = parall[11]
  stcmd("export PRESTO=%s" % prestor)
  fitsfile = glob.glob(r"%s*%s" % (droute, lastname))
  fitsname = fitsfile[0]
  rawfiles = droute+"*"+lastname
  fftfiles = basename
  targetTxtFile =  basename.split(".")[0] + ".txt"
  stcmd("readfile %s > %s" % (fitsname, targetTxtFile))
  fctr = float(findInfoInTxt(targetTxtFile,"Central freq (MHz)"))
  numchan = int(findInfoInTxt(targetTxtFile,"Number of channels"))
  BW = float(findInfoInTxt(targetTxtFile,"Total Bandwidth (MHz)"))
  dt = float(findInfoInTxt(targetTxtFile,"Sample time (us)"))/(10**6)
  cDM = 0.0
  ok_smearing = 0.0
  numout = int(findInfoInTxt(targetTxtFile,"Spectra per file")) 
  #Judge whether use search in low frequency
  if fctr > 350.0:
      lof=0
  else:
      lof=1
#  lof=1
  # The name of the maskfile to apply (if no mask, use None)
  maskfile = basename+"_rfifind.mask"
  # The name of topo file
  topofile = basename+"*topo*"
  # The name of topo_DM0.00 file 
  topo0 = basename+"_topo_DM0.00"
  topo0dat = basename+"_topo_DM0.00.dat"
  topo0fft = basename+"_topo_DM0.00.fft"
  rfiinf = basename+"_rfifind.inf"
  rfiinf2 = basename+".inf"
  birds = basename+".birds"
  zaplist = basename+".zaplist"
  
  stcmd("rm -rf *DM*")
  stcmd("rm -rf *rfi*")
  stcmd("rm -rf subbands")
  stcmd("rm -rf rawfft")
  stcmd("rm -rf *.txt")
  stcmd("rm -rf figures")
 
  timefile = open('time.txt','w')
  timefile.write("Pipeline begin\n")
  nowtime = time.asctime( time.localtime(time.time()) )
  timefile.write(nowtime+'\n')
  timefile.write("Begin rfifind\n")
  starttime = datetime.datetime.now()
  totalstart = datetime.datetime.now()
  nowtime = time.asctime( time.localtime(time.time()) )
  timefile.write(nowtime+'\n')
  rfifind_thread_num = 28
  if zc == 0:
      stcmd("rfifind -noclip -ncpus %d -time %f -o %s %s" % (rfifind_thread_num, rfitime, basename, rawfiles))
  else:
      stcmd("rfifind -noclip -ncpus %d -time %f -o %s -zapchan %s %s" % (rfifind_thread_num, rfitime, basename, zapchan, rawfiles))
  timefile.write("End rfifind\n")
  endtime = datetime.datetime.now()
  nowtime = time.asctime( time.localtime(time.time()) )
  timefile.write(nowtime+'\n')
  durationt = float((endtime-starttime).seconds)/60
  timefile.write("rfifind duration: %f minutes\n\n" % durationt)
  if birds == 1:
      stcmd("prepdata -nobary -o %s -dm 0.0 -mask %s -numout %d %s" % (topo0, maskfile, numout, rawfiles))
      stcmd("realfft %s" % topo0dat)
      stcmd("accelsearch -ncpus 10 -numharm 32 -zmax 0 %s" % topo0dat);
      stcmd("cp %s %s" % (rfiinf, rfiinf2))
      creatbirds(birds)
      writebirds(basename)
      stcmd("makezaplist.py %s" % birds)
      stcmd("zapbirds -zap -zapfile %s %s" % (zaplist, topo0fft))
  timefile.write("Begin de-dispersion\n")
  nowtime = time.asctime( time.localtime(time.time()) )
  starttime = datetime.datetime.now()
  timefile.write(nowtime+'\n')
  if lof == 0:
    # dDM steps from DDplan.py
    dDMs = []
    # dsubDM steps
    dsubDMs   = []
    # downsample factors
    downsamps = []
    # number of calls per set of subbands
    subcalls  = []
    # The low DM for each set of DMs 
    startDMs  = []
    # DMs/call
    dmspercalls = []
    
    device = "/xwin"

    obs = observation(dt, fctr, BW, numchan, cDM)
      # The following function creates the de-dispersion plan
      # The ok_smearing values is optional and allows you to raise the floor
      # and provide a level of smearing that you are willing to accept (in ms)
    dm_steps(loDM, hiDM, obs, cDM, numsubbands, ok_smearing, device)

      # The following is an instance of an "observation" class
      # Here's one for a "best" resolution GBT search using the SPIGOT
      # Check out how many DMs you need!  Cool.  ;-)
      #                    dt     f_ctr   BW  numchan
      #obs = observation(0.00008192, 350.0, 50.0, 2048)
      #dm_steps(0.0, 500.0, obs, ok_smearing=0.3) # Create an X-window is the default
      #dm_steps(0.0, 500.0, obs, ok_smearing=0.3, device="GBT_350_DD.ps/CPS")

      # Example for Parkes survey
      #obs = observation(0.000250, 1374.0, 288.0, 96)
      #dm_steps(0.0, 500.0, obs) # Create an X-window is the default
      #dm_steps(0.0, 1500.0, obs, device="PKMB_DD.ps/CPS")

    # Number of subbands
    nsub = numsubbands

    # Loop over the DDplan plans
    for dDM, dsubDM, downsamp, subcall, startDM, dmspercall in \
            zip(dDMs, dsubDMs, downsamps, subcalls, startDMs, dmspercalls):
        # Get our downsampling right
        subdownsamp = downsamp/2
        datdownsamp = 2
        if downsamp < 2: subdownsamp = datdownsamp = 1
        # Loop over the number of calls
        for ii in range(subcall):
            subDM = startDM + (ii+0.5)*dsubDM
            # First create the subbands
            if maskfile:
                stcmd("prepsubband -mask %s -sub -subdm %.2f -nsub %d -downsamp %d -o %s %s" %
                          (maskfile, subDM, nsub, subdownsamp, basename, rawfiles))
            else:
                stcmd("prepsubband -sub -subdm %.2f -nsub %d -downsamp %d -o %s %s" %
                          (subDM, nsub, subdownsamp, basename, rawfiles))
            # And now create the time series
            loDM = startDM + ii*dsubDM
            subnames = basename+"_DM%.2f.sub[0-9]*"%subDM
            stcmd("prepsubband -numout %d -lodm %.2f -dmstep %.2f -numdms %d -downsamp %d -o %s %s" %
                      (numout/downsamp, loDM, dDM, dmspercall, datdownsamp, basename, subnames))

  else:
    parser = argparse.ArgumentParser(description="""
    Used to calculate the Dedispersion plan for low-frequency telescopes. Inspired by PRESTO's DDplan.py.
      """)
    parser.add_argument('-f', '--centrefreq', type=float, default=fctr,
                        help='Centre frequency of the observation in MHz.')
    parser.add_argument('-b', '--bandwidth', type=float, default=BW,
                        help='Bandwidth of the observation in MHz.')
    parser.add_argument('-nf', '--nfreqchan', type=int, default=numchan,
                        help='Number of frequency channels.')
    parser.add_argument('-t', '--timeres', type=float, default=0.1,
                        help='Time resolution in ms.')
    parser.add_argument('-ld', '--lowDM', type=float, default=1.,
                        help='Lowest DM of the required range.')
    parser.add_argument('-hd', '--highDM', type=float, default=250.,
                        help='Highest DM of the required range.')
    parser.add_argument('-o', '--obsid', type=int,
                        help='The MWA observation ID of an observation. Using this command will get the require observation parameters.')
    parser.add_argument('-m', '--min_DM_step', type=float, default=0.02,
                        help='The  minimun DM step size, default 0.02')
    parser.add_argument('-p', '--plot', action='store_true',
                        help='Plot the sensitivty of the DM plan')
    parser.add_argument('--time', type=int, default=4800,
                        help='Time in seconds to calculate the sensitivity')
    #parser.add_argument()
    args=parser.parse_args()

    if args.obsid:
        #get the centrefreq from the obsid metadata
        beam_meta_data = get_common_obs_metadata(args.obsid)
        obs, ra, dec, dura, [xdelays, ydelays], centrefreq, channels = beam_meta_data

        args.centrefreq = channels

    DedispList = []
    Process = []

    DD_plan_array = dd_plan( args.centrefreq, args.bandwidth, args.nfreqchan, args.timeres, args.lowDM, args.highDM, min_DM_step=args.min_DM_step)
    print(" low DM | high DM | DeltaDM | Nsteps | Downsamp | Effective time resolution (ms)")
    total_steps = 0
    for d in DD_plan_array:
        print("{0:7.1f} | {1:7.1f} | {2:7.2f} | {3:6d} | {4:8d} | {5:7.3f}".\
               format(d[0], d[1], d[2], d[3], d[5], d[4]))
        total_steps += d[3]
 #   print("Total DM steps required: {}".format(total_steps))
 #         index_downsamps = 0
 #         while (args.timeres*allow_downsamps[index_downsamps+1] <= d[4]*2):
 #             index_downsamps += 1    
 #         downsamp = allow_downsamps[index_downsamps]
        dispstart = d[0]
        while (dispstart+100*d[2] < d[1]):
              #stcmd("prepsubband -ncpus %d -mask *.mask -lodm %f -nsub %d -dmstep %f -numdms 100 -numout %d -downsamp %d -zerodm -o %s %s" %(ncpus, dispstart, numsubbands, d[2], numout, d[5], basename, rawfiles))
              DedispList.append("prepsubband -ncpus %d -mask *.mask -lodm %f -nsub %d -dmstep %f -numdms 100 -numout %d -downsamp %d -zerodm -o %s %s" %(1, dispstart, numsubbands, d[2], numout, d[5], basename, rawfiles))
         
              dispstart = dispstart+100*d[2]
        #stcmd("prepsubband -ncpus %d -mask *.mask -lodm %f -nsub %d -dmstep %f -numdms %d -numout %d -downsamp %d -zerodm -o %s %s" %(ncpus, dispstart, numsubbands, d[2], (d[1]-dispstart)/d[2], numout, d[5], basename, rawfiles))
        DedispList.append("prepsubband -ncpus %d -mask *.mask -lodm %f -nsub %d -dmstep %f -numdms %d -numout %d -downsamp %d -zerodm -o %s %s" %(1, dispstart, numsubbands, d[2], (d[1]-dispstart)/d[2], numout, d[5], basename, rawfiles))
    for cmd in DedispList:
         os.system("echo "+cmd+" >> Dedispcmd.txt")
    os.system("python dedispersion.py")
       
    print("Total DM steps required: {}".format(total_steps))
  timefile.write("End de-dispersion\n")
  nowtime = time.asctime( time.localtime(time.time()) )
  endtime = datetime.datetime.now()
  timefile.write(nowtime+'\n')
  durationt = float((endtime-starttime).seconds)/60
  timefile.write("de-dispersion duration: %f minutes\n\n" % durationt)

  timefile.write("Begin FFT\n")
  nowtime = time.asctime( time.localtime(time.time()) )
  starttime = datetime.datetime.now()
  timefile.write(nowtime+'\n')
  stcmd("rm %s" % topofile)
  stcmd("ls *.dat | xargs -n 1 realfft")
  stcmd("mkdir rawfft")
  timefile.write("End FFT\n")
  nowtime = time.asctime( time.localtime(time.time()) )
  endtime = datetime.datetime.now()
  timefile.write(nowtime+'\n')
  durationt = float((endtime-starttime).seconds)/60
  timefile.write("FFT duration: %f minutes\n\n" % durationt)
  timefile.write("Begin remove rednoise\n")
  starttime = datetime.datetime.now()
  nowtime = time.asctime( time.localtime(time.time()) )
  timefile.write(nowtime+'\n')
  #Remove rednoise
  if lof == 0:
    for dDM, dsubDM, subcall, startDM, dmspercall in \
          zip(dDMs, dsubDMs, subcalls, startDMs, dmspercalls):  
        for ii in range(subcall*dmspercall):
            fftDM = startDM + ii*dDM
            fftfiles2 = fftfiles + '_DM%.2f.fft' % fftDM
            fftfiles3 = fftfiles + '_DM%.2f_red.fft' % fftDM 
            stcmd("rednoise %s" % fftfiles2)
            stcmd("mv %s ./rawfft" % fftfiles2)
            stcmd("mv %s %s" % (fftfiles3, fftfiles2))
  else:
    for d in DD_plan_array:
        for ii in range(int(d[3])):
            fftDM = d[0] + ii*d[2]
            fftfiles2 = fftfiles + '_DM%.2f.fft' % fftDM
            fftfiles3 = fftfiles + '_DM%.2f_red.fft' % fftDM
            stcmd("rednoise %s" % fftfiles2)
            stcmd("mv %s ./rawfft" % fftfiles2)
            stcmd("mv %s %s" % (fftfiles3, fftfiles2))
  timefile.write("End remove rednoise\n")
  nowtime = time.asctime( time.localtime(time.time()) )
  endtime = datetime.datetime.now()
  timefile.write(nowtime+'\n')
  durationt = float((endtime-starttime).seconds)/60
  timefile.write("remove rednoise duration: %f minutes\n\n" % durationt)
  if birds == 1:
      stcmd("prepdata -o tmp %s | grep Average > Average.txt" % basename)
      AveTopocentricVel = findInfoInTxt("Average.txt","Average topocentric velocity (c)")
      stcmd("ls *.fft | xargs -n 1 zapbirds -zap -zapfile %s -baryv %s" % (zaplist,AveTopocentricVel))
  timefile.write("Begin accelsearch\n")
  starttime = datetime.datetime.now()
  nowtime = time.asctime( time.localtime(time.time()) )
  timefile.write(nowtime+'\n')
  #stcmd("ls *.fft | xargs -n 1 accelsearch -ncpus %d -numharm %d -zmax %d" % (ncpus,numharm,zmax))
  stcmd("python accelsearch.py")
  timefile.write("End accelsearch\n")
  nowtime = time.asctime( time.localtime(time.time()) )
  endtime = datetime.datetime.now()
  timefile.write(nowtime+'\n')
  durationt = float((endtime-starttime).seconds)/60
  timefile.write("accelsearch duration: %f minutes \n\n" % durationt)

  stcmd("mkdir subbands")
  stcmd("mv *_red.inf subbands/")
  #stcmd("mv tmp* subbands/")
  stcmd("cp -r $PRESTO/examplescripts/ACCEL_sift.py .")
  stcmd("python ACCEL_sift.py > cands.txt")
  stcmd('''cat cands.txt | grep _ACCEL_ | awk '{print $2"\t"$8}' > candslist.txt''')
  if fcand == 1:
      PrepfoldList = []
      timefile.write("Begin prepfold\n")
      nowtime = time.asctime( time.localtime(time.time()) )
      starttime = datetime.datetime.now()
      timefile.write(nowtime+'\n')
      candsfile = 'candslist.txt'
      DMsearch,Psearch = [],[]
      with open(candsfile, 'r') as cands:
          candslines = cands.readlines()
          for candsline in candslines:
              candsval = [float(cval) for cval in candsline.split()]
              DMsearch.append(candsval[0])
              Psearch.append(candsval[1])
              PrepfoldList.append("prepfold -ncpus 1 -n 128 -noxwin -noclip -o %s_DM%.2f -p %f -dm %f -nosearch  -mask *.mask %s" % (basename, candsval[0], candsval[1]/1000, candsval[0], rawfiles))
      for cmd in PrepfoldList:
         os.system("echo "+cmd+" >> Prepfold.txt")
      os.system("python prepfold.py")
      timefile.write("End prepfold\n")
      nowtime = time.asctime( time.localtime(time.time()) )
      endtime = datetime.datetime.now()
      timefile.write(nowtime+'\n')
      durationt = float((endtime-starttime).seconds)/60
      timefile.write("prepfold duration: %f minutes \n\n" % durationt) 
  if singlep == 1:
      timefile.write("Begin single pulse search\n")
      nowtime = time.asctime( time.localtime(time.time()) )
      starttime = datetime.datetime.now()
      timefile.write(nowtime+'\n')
      stcmd("python /ssd/Pulsar/lxl/single_pulse_search_opt.py -b *.dat")
      timefile.write("End single pulse search\n")
      nowtime = time.asctime( time.localtime(time.time()) )
      endtime = datetime.datetime.now()
      timefile.write(nowtime+'\n')
      durationt = float((endtime-starttime).seconds)/60
      timefile.write("single pulse search duration: %f minutes \n\n" % durationt)
  stcmd("mkdir figures")
  stcmd("mv *.pfd ./figures/")
  stcmd("mv *.bestprof ./figures/")
  stcmd("mv *.ps ./figures/")
  timefile.write("End pipeline\n")
  nowtime = time.asctime( time.localtime(time.time()) )
  totalend = datetime.datetime.now()
  timefile.write(nowtime+'\n')
  durationt = float((totalend-totalstart).seconds)/60
  timefile.write("Total duration: %f minutes\n" % durationt)
  timefile.close()
