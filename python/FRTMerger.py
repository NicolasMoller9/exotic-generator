# -*- coding: utf-8 -*-

# FRTMerger
# Merge Signal DOMLaunchMap with fixed rate trigger data
#
# (c) 2012 Simon Zierke
#
#$Id$
#

from icecube import icetray, dataclasses, dataio, trigger_sim
import random
import glob
import math
import os

class SlowExoHit():
  def __init__(self, OMKey, time, info):
     self.omkey = OMKey
     self.time = time
     self.info = info
  def __str__(self):
    return "%s, Time %s ms" % (self.omkey, self.time / icetray.I3Units.millisecond)
  
class FRTMerger(icetray.I3ConditionalModule):

  def __init__(self, context):
    icetray.I3ConditionalModule.__init__(self,context)

    self.AddParameter("SignalMap","Map/Mask for input", "SignalRawData")

    self.AddParameter("FRTMap","Name for the input", "InIceRawData")
    self.AddParameter("OutputName","Name for the output", "InIceRawData")

    self.AddParameter("StartWindow","", -10 * icetray.I3Units.microsecond)
    self.AddParameter("EndWindow","", 10 * icetray.I3Units.microsecond)

    self.AddParameter("FrtFileList","", glob.glob(os.path.expandvars("$I3_SRC/exotic-generator/resources/important_files/Run_Untriggered_Background_Recombined.i3.zst")))
    self.AddParameter("DeadTimeSLC","", 2.4 * icetray.I3Units.microsecond)
    self.AddParameter("DeadTimeHLC","", 6.45 * icetray.I3Units.microsecond)
    self.AddParameter("ErrorCount","Maximal number of tries", 5)
    self.AddOutBox('OutBox')

    # Check if DOMSet6 exists #TODO Why do we care for DeepCore
    DefaultDOMSets = trigger_sim.GetDefaultDOMSets()
    if not trigger_sim.InDOMSet(icetray.OMKey(80,50), 6, DefaultDOMSets):
         raise Exception(NotImplementedError, "InDOMSet error")


  def Configure(self):
        self.SignalMap = self.GetParameter("SignalMap")
        self.FRTMap = self.GetParameter("FRTMap")
        self.OutputName = self.GetParameter("OutputName")

        self.StartWindow = self.GetParameter("StartWindow")
        self.EndWindow = self.GetParameter("EndWindow")

        self.FrtFileList = self.GetParameter("FrtFileList")

        self.DeadtimeSLC = self.GetParameter("DeadTimeSLC")
        self.DeadtimeHLC = self.GetParameter("DeadTimeHLC")

        self.ErrorCount = self.GetParameter("ErrorCount")

        self.Hits = 0
        self.Rejected = 0
        self.RejectedFrame = 0

        frtFile = self.FrtFileList[0]
        frtI3 = dataio.I3File(frtFile)
        # for i in range(random.randint(4,1000)):
        #     frtFrame = frtI3.pop_frame()
        # frtI3.close()
        self.FRTFrame = frtI3


  def Physics(self, frame):
    self.PushFrame(frame)


  def DAQ(self, frame):
    signalHits = []
    def numeric_compare(x, y):
        return x - y
    if frame.Has(self.SignalMap):
        if len(frame[self.SignalMap])>0:
            # print ( frame[self.SignalMap])
            for omkey, launchVector in frame[self.SignalMap]:
                for launch in launchVector:
                    signalHits += [SlowExoHit(omkey, launch.time, launch)]

            signalHits=sorted(signalHits, key = lambda n: n.time)

            self.Error = False
            self.FRTFile = ""
            domLaunchSeriesMap = self.MixFRT(signalHits)

            self.Name = "FrtMerger"
            frame[self.OutputName] = domLaunchSeriesMap
            frame[self.Name + "_Rejected"] = dataclasses.I3Double(self.RejectedFrame)
            frame[self.Name + "_Error"] = icetray.I3Bool(self.Error)
            frame[self.Name + "_FRTFile"] = dataclasses.I3String(self.FRTFile)

            self.PushFrame(frame)


  def Finish(self):
    #print ("Hits merged: %i davon %i verworfen (%f)" % (self.Hits, self.Rejected, float(self.Rejected)/self.Hits))
    return True


  def MixFRT(self, signalHits, errorcount = 0):
    frtHits = []
    start = signalHits[0].time
    end = signalHits[-1].time
    length = end - start


    #frtFile = self.FrtFileList[random.randint(4,len(self.FrtFileList)-1)]
    #frtI3 = dataio.I3File(frtFile)


    frtFile = self.FrtFileList[0]
    frtI3 = dataio.I3File(frtFile)
    for i in range(random.randint(5,10000)):
        frtFrame = frtI3.pop_frame()
        try:
            var= frtFrame[self.FRTMap]
        except KeyError:
            frtFrame = frtI3.pop_frame()
        # frtFrame = self.FRTFrame.pop_frame()
    frtI3.close()

    #self.FRTFile = frtFile

    # for trigger in frtFrame["I3TriggerHierarchy"]:
    #     if trigger.key.config_id == 23050:
    #         trigger_start = trigger.time - 5 * icetray.I3Units.millisecond # Trigger starttime is set to 5ms
    #         trigger_length = 10 * icetray.I3Units.millisecond
    #         trigger_end = trigger_start + trigger_length

    # if length > trigger_length:
    #     raise  Exception(NotImplementedError, "Eventlength too long for FRT-Data")

    front = 1000 * icetray.I3Units.microsecond
    # offset = random.uniform(front, max(front, trigger_length - length - front) )

    # DOMLaunches innerhalb des Triggers herausschreiben
    for omkey, launchVector in frtFrame[self.FRTMap]:
        for launch in launchVector:
            # if launch.time > trigger_start and launch.time < trigger_end:
            #     launch.time += start - trigger_start - offset + (end-start)/2.

            frtHits += [SlowExoHit(omkey, launch.time, launch)]
            self.Hits += 1

    frtHits=sorted(frtHits, key = lambda n:  n.time)

    # HLC Start- und Endzeit vom Signal
    timeWindow_Signal = self.CheckTriggerCandidates(signalHits)

    # HLC Start- und Endzeit vom FRT
    timeWindow = self.CheckTriggerCandidatesWindow(frtHits, timeWindow_Signal[0], timeWindow_Signal[1])

    # Mindestens ganzes Signal soll in den Daten sein
    timeWindow = (min(timeWindow[0], start), max(timeWindow[1], end))



    self.RejectedFrame = 0

    ## Deadtime erhalten, Signal bevorzugen
    for signal in signalHits:
        for frt in frtHits:
            if frt.omkey == signal.omkey:
                first = min([signal, frt], key=lambda x: x.time)
                if (first.info.lc_bit and abs(frt.time - signal.time) < self.DeadtimeHLC) or \
                (not first.info.lc_bit and abs(frt.time - signal.time) < self.DeadtimeSLC):

                    self.RejectedFrame += 1
                    self.Rejected += 1
                    frtHits.pop(frtHits.index(frt))


    domLaunchSeriesMap = dataclasses.I3DOMLaunchSeriesMap()

    # Ränder müssen größer als 500mus sein
    # if errorcount < self.ErrorCount:
    #     if  timeWindow[0] - (start - offset + (end-start)/2.) < 500 * icetray.I3Units.microsecond:
    #         frtHits = None
    #         return self.MixFRT(signalHits, errorcount + 1)
    #
    #     if (start + trigger_length - offset + (end-start)/2.) - timeWindow[1] < 500 * icetray.I3Units.microsecond:
    #         frtHits = None
    #         return self.MixFRT(signalHits, errorcount + 1)
    # else:
    #     self.Error = True

    # Signale in timeWindow mergen
    while len(signalHits) or len(frtHits):
        hit = self.PopNextHit(signalHits, frtHits)
        if hit.time < timeWindow[0]+ self.StartWindow or hit.time > timeWindow[1] + self.EndWindow:
            continue

        if domLaunchSeriesMap.has_key(hit.omkey):
            domLaunchSeriesMap[hit.omkey].append(hit.info)
        else:
            domLaunchSeries = dataclasses.I3DOMLaunchSeries()
            domLaunchSeries.append(hit.info)
            domLaunchSeriesMap[hit.omkey] = domLaunchSeries

    return domLaunchSeriesMap


  def PopNextHit(self, ListA, ListB):
    ## Abbruch bei leeren Listen
    if len(ListB) == 0:
        return ListA.pop(0)

    if len(ListA) == 0:
        return ListB.pop(0)

    ## Nächsten Hit ausgeben
    if ListA[0].time <= ListB[0].time:
        return ListA.pop(0)

    if ListA[0].time > ListB[0].time:
        return ListB.pop(0)


  def CheckTriggerCandidates(self, hitlist):
    self.debug = False

    self.t_proximity = 2500
    self.t_min = 0
    self.t_max = 500000

    self.OneHitList = []
    self.TwoHitList = []
    self.muon_time_window = -1

    DefaultDOMSets = trigger_sim.GetDefaultDOMSets()
    hitlist=sorted(hitlist, key = lambda hit:  hit.time)
    for hit in hitlist:
        if hit.info.lc_bit == True:
            if trigger_sim.InDOMSet(hit.omkey, 6, DefaultDOMSets): #TODO Why look only hits which were inside DeepCore
                self.CheckOneHitList(hit)
    
    if len(self.TwoHitList)>1:
        self.TwoHitList=sorted(self.TwoHitList, key = lambda n:  n.time)
        return(self.TwoHitList[0].time, self.TwoHitList[-1].time)
    return (hitlist[0].time, hitlist[-1].time)


  def CheckTriggerCandidatesWindow(self, hitlist, start=None, ende=None):
    self.debug = False

    self.t_proximity = 2500
    self.t_min = 0
    self.t_max = 500000

    self.OneHitList = []
    self.TwoHitList = []
    self.muon_time_window = -1

    DefaultDOMSets = trigger_sim.GetDefaultDOMSets()
    for hit in hitlist:
        if hit.info.lc_bit == True:
            if trigger_sim.InDOMSet(hit.omkey, 6, DefaultDOMSets):
                self.CheckOneHitList(hit)

    if len(self.TwoHitList)>1:
        n = 0
        while n < len(self.TwoHitList[:-1]):
            if self.TwoHitList[n].time > start:
                break

            if abs(self.TwoHitList[n].time - min(self.TwoHitList[n+1].time, start) ) > self.t_max:
                for i in range(n+1):
                    #print i , len(self.TwoHitList)
                    self.TwoHitList.pop(0)
                n = 0
            else:
                n += 1

        for n in range(len(self.TwoHitList[:-1])):
            if self.TwoHitList[n].time < ende:
                continue
            if abs(self.TwoHitList[n].time - self.TwoHitList[n+1].time ) > self.t_max:
                self.TwoHitList = self.TwoHitList[:n+1]
                break

        return (min(self.TwoHitList[0].time, start), max(self.TwoHitList[-1].time, ende))

    return (start, ende)

#
# Parts of SLOPTools - TupleTagger by Emanuel Jacobi
# modified for IC86-2011
# removed tuple check
#
  def HLCPairCheck(self, hit1, hit2):
    if hit1.omkey.string == hit2.omkey.string:
        if abs(hit1.omkey.om - hit2.omkey.om) <= 2:
            if self.debug:
                print ("self.debug: HLC found!!! String: %i, DOM: %i, Time: %f" % (hit1.omkey.string, hit1.omkey.om, hit1.time))
                print ("self.debug: The other is String: %i, DOM: %i, Time: %f" % (hit2.omkey.string, hit2.omkey.om, hit2.time))
            return True
    return False


  def CheckOneHitList(self, newHit):
    if len(self.OneHitList) == 0:                        # one hit list is empty. just add current hit to list
        self.OneHitList.append(newHit)
    else:                                                # one hit list contains stuff. compare current hit to the one in the list
        while math.fabs(newHit.time - self.OneHitList[0].time) > 1000.:
            self.OneHitList.pop(0)
            if len(self.OneHitList) == 0: break
        i=0
        while i < len(self.OneHitList):                  # iterate over one hit list to form HLCs
            if self.HLCPairCheck(self.OneHitList[i], newHit):
                payload=self.OneHitList[i]
                self.CheckTwoHitList(payload)
                self.OneHitList.pop(i)                   # delete the used hit from OneHitList after checking TwoHitList
            else:
                i+=1
        self.OneHitList.append(newHit)                   # at the end add the current hitPayload for further comparisons


## Kann vereinfacht werden
  def CheckTwoHitList(self, payload):
    if self.debug: print ("self.debug: ", payload.omkey, payload.time, payload.info, len(self.TwoHitList))
    if len(self.TwoHitList) == 0:                        # hit list is empty
        if self.muon_time_window == -1:
            self.TwoHitList.append(payload)                      # add to list (no time_window set)
            if self.debug: print ("self.debug:  adding hit %f" % payload.time)
        else:
            if payload.time - self.muon_time_window <= self.t_proximity:
                self.muon_time_window = payload.time                 # not adding hit, due to t_proximity, set time window to new hit
                if self.debug: print ("self.debug:  not adding hit, due to t_proximity")
            else:
                self.TwoHitList.append(payload)                  # add to list (time window set)
                self.muon_time_window = -1                         # reset time window
                if self.debug: print ("self.debug:  adding hit %f" % payload.time)
    else:                                                   # hit list is not empty
        if self.muon_time_window == -1:
            if payload.time - self.TwoHitList[-1].time <= self.t_proximity:
                self.muon_time_window = payload.time                 # t_proximity test failed, discard current hit, set time window, remove last hit
                self.TwoHitList.pop()
                if self.debug: print ("self.debug:  not adding hit, due to t_proximity, deleting last hit")
            else:                                         # t_proximity test passed
                if payload.time - self.TwoHitList[-1].time < self.t_max:
                    self.TwoHitList.append(payload)              # add to list
                    if self.debug: print ("self.debug:  adding hit %f" % payload.time)
                else:
                    self.CheckTriggerStatus()     # this hit is too late to make a tuple, but check previous hits in list
                    self.TwoHitList.append(payload)              # now add current hit to list
                    if self.debug: print ("self.debug:  adding hit %f after checking previous hits" % payload.time)
        else:                                            # time window is set
                 if payload.time - self.muon_time_window <= self.t_proximity:
                     self.muon_time_window = payload.time        # t_proximity test failed, disregard current pulse, set time window
                     if self.debug: print ("self.debug:  not adding hit, due to t_proximity")
                 else:                                      # t_proximity test passed
                     self.muon_time_window = -1                    # this is not a muon, set window back to -1
                     if payload.time - self.TwoHitList[-1].time < self.t_max:
                         self.TwoHitList.append(payload)          # add to list
                         if self.debug: print ("self.debug:  adding hit %f" % payload.time)
                     else:
                         self.CheckTriggerStatus() # this hit is too late to make a tuple, but check previous hits in list
                         self.TwoHitList.append(payload)          # now add current hit to list
                         if self.debug: print ("self.debug:  adding hit %f after checking previous hits" % payload.time)

  def CheckTriggerStatus(self):
    if self.debug:
        print ("Checking trigger status")
        #for slowexohit in self.TwoHitList:
            #print slowexohit
    # the c++ version loops here through the trigger_container_vector (which holds the tuples)
    # and issues the trigger if n>=min_tuples. after that the container is cleared
    #self.TwoHitList = []

