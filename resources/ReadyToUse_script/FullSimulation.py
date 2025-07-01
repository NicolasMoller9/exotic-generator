# Usage exemples:
#1) Monopole with proton decay but no luminescence
#python3 FullSimulation.py -pID -41 -b_ 0.001 -lamb 10 -keepP 0 -folder test/ -nevents 10 -num 1
#2) Exotic (ex: shadow charge) with only luminescence
#python3 FullSimulation.py -pID -42 -Q 1 -b_ 0.001 -keepP 1 -folder test/ -nevents 10 -num 1
#3) Shadow charge with beta range
#python3 FullSimulation.py -pID -42 -Q 1 -b_ 0.001 -bplus 0.1 -keepP 1 -folder test/ -nevents 100 -num 1

import os
import math
import time
import argparse
import sys
import icecube
from icecube import icetray, exotic_generator, dataclasses, dataio, trigger_sim, wavedeform, WaveCalibrator, DOMLauncher
from icecube.icetray import I3Tray, load, I3Units
from icecube.dataclasses import I3VectorDoubleDouble, make_pair
from icecube.hdfwriter import I3HDFTableService
from icecube.exotic_generator.FRTMerger import FRTMerger
from datetime import datetime
from os.path import expandvars
load("xppc")
load("ppc")
# Are you runing on gpu?
gpu = False
if gpu: os.putenv("OGPU", "1")  # makes sure only GPUs are used (with OpenCL version)
else: os.putenv("OCPU", "1")  # makes sure only GPUs are used (with OpenCL version)
# Path to your output storage space 
output_folder = './'


# The simulation is split in 3 parts:
# I) Exotic generation and propagation
# II) Light simulation using Photon Propagation Code (PPC)
# III) Detector simulation
# Comment out part II or II + III if not needed

### PARAMETERS FOR PART I #########################################################################################################################################################
gcdFile = '/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_2020.Run134142.Pass2_V0.i3.gz'  # Path to the gcd file containing detector information. Can be also found in: '$I3_SRC/exotic-generator/resources/important_files/' 
DiskD = 1000 # m - Generation-disk distance from IceCube center
DiskR = 800  # m - Generation-disk radius
### PARAMETERS FOR PART II ########################################################################################################################################################
icemodel = "spice_3.2"    # Specify ice model
os.putenv("PPCTABLESDIR", expandvars(f"$I3_SRC/exotic-generator/resources/important_files/ppc_config_files/{icemodel}/"))
os.putenv("ICEMODELDIR", expandvars(f"$I3_SRC/ice-models/resources/models/ICEMODEL/{icemodel}/"))
os.putenv("PPCHOLEICE", expandvars("$I3_SRC/exotic-generator/resources/important_files/ppc_config_files/HoleIceParametrization/as.h2-50cm"))
UnshadowedFraction = 0.94  # Effect of the shadowing of the cables (default clsim value = 0.94)
DOMEfficiency = 1
DOMOversizeFactor = 1      # An oversize of 5 will make DOMs 5 times bigger to catch more photons which are later scaled down
check_method = 1           # 2 options to check wether a photon is inside the detector and thus worth propagating: 1 (faster) or 0 (loop through each DOM)
dnde_tau = [[0.999999684981697, 2.44], [3.1343842585129187e-07, 196.1], [9.551499701160983e-10, 5030.0], [6.247272615115733e-10, 56100.0]] # from [1]
tau_dnde = I3VectorDoubleDouble([ make_pair(dnde_tau[i][1],dnde_tau[i][0]) for i in range(len(dnde_tau)) ])  # Decay time and strength for luminescence
### PARAMETERS FOR PART III #######################################################################################################################################################
FRT = True           # Decide to add background with FRT Merger or not
OnlySLOP = False      # Decide to keep only the SLOP triggers and discard the others
throw_untrig = True # Decide to throw away events that do not trigger
###################################################################################################################################################################################


startTime = datetime.now()
print(f'Simulation starting time: {startTime}')

if __name__ == "__main__":
    # 0) INITIALIZATION
    # 0.1) Define the input parameters
    parser = argparse.ArgumentParser()
    parser.add_argument("-pID",     type = int,   dest="part_type",   required=True,     help="Define which exotic to simulate: -41=Monopole, -42=ShadowCharge, -43=Qball")
    parser.add_argument("-M",       type = float, dest="mass",        required=False,    help="Mass of the particle [GeV]. Allowed range=[1e5,1e17] for Monopoles and [1e17,1e29] for Qballs")
    parser.add_argument("-Q",       type = int,   dest="charge",      default=1,         help="Particle Charge. For monopoles it represents g/gD")
    parser.add_argument("-b_",      type = float, dest="beta_min",    required=True,     help="minimum beta value")
    parser.add_argument("-bplus",   type = float, dest="beta_max",    required=False,    help="maximum beta value. If not set, all exotics will have the same speed (beta_min)")
    parser.add_argument("-lamb",    type = float, dest="lambdaa",     default=-1,        help="mean free path for proton decay catalysis [m]. Use -1 to calculate lambda from the cross section value.")
    parser.add_argument("-sig0",    type = float, dest="sigmaa",      required=False,    help="cross section for proton decay catalysis [m]. Exemple: 1e-31 [m^2]")
    parser.add_argument("-keepP",   type = int,   dest="keepPrimary", default=1,         help="Wether to keep the primary particle track, allowing it to produce light (=1), or make it dark (=0)")
    parser.add_argument("-folder",                dest="folder",      default='test',    help="folder to save the output. Will be created in 'output_folder' if doesn't exist yet")
    parser.add_argument("-nevents", type = int,   dest="evento",      required=True,     help="number of events. The number of exotic generated will be nevents-3")
    parser.add_argument("-num",     type = float, dest="numbe",       required=True,     help="Id of the simulation. Used for the output file name.")
    args = parser.parse_args()

    # 0.2) Catch the input parameters
    pID = int(args.part_type)
    part_name = {-41:'monopoles', -42:'shadow charges', -43:'Qballs'}[pID]
    if args.mass is None:
        if pID==-41:
            Mass = 1e11
            print("Warning: Simulating Monopoles without setting the mass! Set to default value of 1e11 GeV/c^2. Ignore this warning if you are not simulating energy loss.")
        elif pID==-43:
            Mass = 1e17
            print("Warning: Simulating Qballs without setting the mass! Set to default value of 1e17 GeV/c^2. Ignore this warning if you are not simulating energy loss.")
        else: Mass = 0
    else: Mass = args.mass
    q = args.charge
    b_ = args.beta_min
    if args.beta_max is None: bp = b_
    else: bp = args.beta_max
    c = args.lambdaa
    if args.sigmaa is None: sig = float('nan')
    else: sig = args.sigmaa
    keep_prim = [False,True][args.keepPrimary]
    output_path = output_folder+str(args.folder)
    os.makedirs(output_path, exist_ok=True)  # Create the output directory if doesn't exist already
    n_events = args.evento
    numb = int(args.numbe)
    tray = I3Tray()
    
    # 0.3) Build a random seed to be able to reproduce the simualtion
    ti_ran=time.time()*10**6
    ti_ran=int(ti_ran%10000000)
    print(f'Seed of the simulation is: {ti_ran}\n')
    tray.AddService("I3GSLRandomServiceFactory", Seed=ti_ran)  # Configurate a random number generator from IceCube software, based on the GSL library
    # Note: When using "Add.Service", the added module (here the random number generator) is registered as a global service in the pipeline. All the other services of the pipeline can then access it automatically
    
    # 0.4) Generate the 3D grid, axis system and detector geometry from the GCD file
    tray.AddModule("I3InfiniteSource", "infinite", Prefix = gcdFile, Stream=icetray.I3Frame.DAQ)

 
    ###### PART I: PARTICLE GENERATION & PROPAGATION #####################################################################################################################
    part = 1
    
    print ("--------------PART1: PARTICLE GENERATION & PROPAGATION--------------")
    print (f'--> This script will generate {n_events} {part_name} with the following parameters: \n- Mass = {Mass} GeV\n- Speed/c = [{b_} , {bp}]')
    if c==-1:
        if math.isnan(sig):
            print("No mean free path or cross section given. No catalysis process will be simulated.")
        else:
            print(f'- Lambda will be calculated from a fixed value of the cross section ({sig} m^2)')
    else:
        print(f'- MeanFreePath = {c} m')
    print(f'--> The {part_name} are generated inside of a disk of radius {DiskR} m located {DiskD} m away from the detector center.\n')
    
    ###GENERATOR##################################
    tray.AddModule("I3ExoticGenerator",'Generator',
              NEvents=n_events+3,                       # The number of events generated is given by NEvents-3
              PartID=pID,                               # Particle identifier: {Monopole:-41,ShadowCharge:-42,Qball:-43}
              Mass=Mass,                                # Particle mass in GeV
              Charge=q,                                 # Particle charge (for monopole, charge=n=g/gD)
              BetaRange=[b_, bp],                       # Must be inside the range [1e-6,0.99995]
              PowerLawIndex=0.,                         # int that determines the slope of the beta distribution between beta_min and beta_max (useless if beta_min=beta_max). Ex: 5 (more statistics at low speed) or -5 (more statistics at high speed)
              Disk_dist=DiskD*icetray.I3Units.m,        # Distance between generation disk and icecube center, in m
              Disk_rad=DiskR*icetray.I3Units.m,         # Radius of the generation disk, in m
              
              # Optional parameters:
              #Length=2.5*DiskD,                        # Length of the exotic track. Can be `NaN` or any length. The default value is calculated to twice the disk distance. Set to `-1` to indicate that the default value, `2 * Disk_dist`, should be used.
              #ZenithRange=[60.0 * I3Units.deg, 90.0 * I3Units.deg],  # List of lower and upper zenith bound. Use this parameter to restrict the direction of the exotic
              #AzimuthRange=[0.0 * I3Units.deg, 20.0 * I3Units.deg],  # List of lower and upper azimuth bound. Use this parameter to restrict the direction of the exotic
              #Rad_on_disk=10*icetray.I3Units.m,        # Set the radius coordinate of the starting position on the generation disk. Randomized if `NaN`
              #Azi_on_disk=20*I3Units.deg,              # Set the azimuth coordinate of the starting position on the generation disk. Randomized if `NaN`
              #StartTime=20*I3Units.ns,                 # The time, measured from the beginning of the event, the exotic particle should be started
              #ShiftCenter=([10.0 * icetray.I3Units.m, 10.0 * icetray.I3Units.m, 10.0 * icetray.I3Units.m])  # Shifts the exotic. This is useful to explore different geometries. To shift according to the center of DeepCore (IC86-I SLOP trigger only acts on DC)
    )
   

    ###PROPAGATOR##################################
    tray.AddModule("I3ExoticPropagator",'Propagator',

                # Catalysis:
                MeanFreePath = c * I3Units.m,                          # (default=-1) Mean free path of the catalysed process, in m. If -1 --> calculate from the value of cross section 
                #CrossSection = sig * I3Units.m * I3Units.m,            # (default=nan) Cross section of the catalysed process, in m^2
                #ScaleEnergy = True,                                    # (default=False) Whether to set the mean free path (lambda) to 1 meter and scale up the energy by 1/lambda. The overall light output stays comparable, while the number of secondary particles in the `I3MCTree` is reduced. This saves computing resources especially for short mean free paths
                #EnergyScaleFactor = 0.1,                               # (default=1) Factor between 0 and 1 to scale down the cascade energy
                #UseCorrectDecay = True,                                # (default=False) Ex: Proton decay: Whether to simulate back-to-back positrons (460 MeV) and neutral pions (480 MeV) instead of just one positron. The overall light output is similar, but correct decay has twice as much secondary particles in the `I3MCTree`. This option cannot be used together with `ScaleEnergy` or `EnergyScaleFactor`, since the energies are hard coded
                KeepPrimary = keep_prim,                               # (default=True) Wether to keep the primary particle track or not

                # Track:
                #MaxDistanceFromCenter = 1100 * I3Units.m,              # (default=1300m) How far beyond the detector to propagate the exotic. If the start of the exotic is further from the detector than this value, the propagator will IGNORE the parameter and propagate until it reaches the same distance away on the other side of detector
                #CalculateEnergy = True,                                # (default=False) Calculate energy loss of the exotic iteratively and reduce its speed accordingly
                #StepSize = 10 * I3Units.m,                             # (default=nan) Size of the constant speed steps, in m. If set, and if calculating energy loss, this will override `MinLength` and `MaxLength`. Otherwise, `MinLength` and `MaxLength` are used to set the lower and upper bounds on the track segment lengths
                #MinLength = 1.0 * I3Units.m,                           # (default=nan) Assuming stepsize is `NaN`, this represents the smallest segment the propagator will generate
                #MaxLength = 10.0 * I3Units.m,                          # (default=nan) Assuming stepsize is `NaN`, this represents the largest segment the propagator will generate
                #Profiling = True,                                      # (default=False) If `true`, adds a profile (type `I3VectorDouble`) of the monopole speed for each track segment to the frame
                #ParticleCheck = False                                  # (default=True) Produce sanity checks on the added segments
                )

    
    #input_path = 'test/'
    #inp_file = f'{input_path}/part1_charge{q}_beta{b_}lambda{c}_nevents_{n_events}_{numb}_slop.i3.zst'
    #tray.Add("I3Reader", Filename=inp_file)
    ###### PART II ################################################################################################################################################
    part = 2
    print ("--------------PART2: LIGHT PRODUCTION & PROPAGATION WITH PPC--------------")
    print(f"--> Ice model used: {icemodel}")
    print(f"--> UnshadowedFraction: {UnshadowedFraction}")
    print(f"--> DOMEfficiency: {DOMEfficiency}")
    print(f"--> DOMOversizeFactor: {DOMOversizeFactor}")
    mth = ['dom loop','cylinder'][check_method]
    print(f'--> Method used to check if light inside the detector: {mth}')
    # More details about the light propagation can be found in the ppc code: ppc\private\ppc\i3ppc.cxx'
    
    tray.Add('i3ppc','ppc',                                 # Add the module i3ppc to the pipeline 'tray' with the internal name 'ppc' 
                    MCTree = "I3MCTree",                    # Specify that the module will take its input data in the structures named 'I3MCTree' of each frame
                    infoName = "PPCInfoDict",               # Name of the ppc info dictionary. Will not be created if set to empty string.
                    gpu = -1,                               # GPU to use. Write n-1 to use the nième gpu, and -1 (default) to use everything accessible
                    cyl = check_method,                     # PPC checks if the track/cascade is close to the detector, and ignore it if too far away. 
                                                            # To verify this, it can use either a cylinder of 300m radius and 600m height centered on IceCube center, 
                                                            # or loop through all the DOMS and consider a sphere of 300m radius around them.
                                                            # cyl=1 => cylinder (faster, but doesn't contain the edges of IceCube), cyl=0 => DOM loop (longer)
                    keep=False,                             # Keep the events even if they don't produce light inside of the geometry
                    verbose=False,                          # Print information messages during the simulation for easier debuging
                    photons=True,                           # Save photons that hit DOMs
                    pseries=False,                          # If false, save photons that hit DOMs into Photonseries instead of McTree
                    charge=q,                               # Charge of the particle (1 by default)
                    tau_dnde_vec=tau_dnde,                  # Vector of pairs of luminescence decay time and dnde, tau in ns, dNdE in eV^-1
                    efficiency_scaling_factor=UnshadowedFraction*DOMEfficiency,
                    oversize=DOMOversizeFactor)

    
    ###### PART 3: DETECTOR SIMULATION ##################################################################################################################################################
    part = 3
    print("--------------PART3: DOM RESPONSE AND TRIGGER--------------")
    print(f'--> Background simulation: {str(FRT)}')

    # 3.1) SIMULATE PMT RESPONSE
    pmt_config = {'Input':"MCPESeriesMap",                     # The input is the PhotoElectron Series Map, containing information on the photons received at each PMT. (Each photoelectron represents a detected photon)
                  'Output':"I3MCPulseSeriesMap",               # Will create a new key named I3MCPulseSeriesMap with the information on the PMT pulses. The pulses represent the response of the PMTs to the detected PEs. Pulses are modeled to account for the PMT's electronics, such as pulse shapes and timing.
                  'MergeHits':True,                            # Determines whether multiple PEs arriving in close succession should be combined into a single pulse. This setting is useful to reflect the realistic behavior of PMTs, which can produce a single output for closely timed photons
                  'RandomServiceName' : 'I3RandomService'}     # Use a consistent random service to ensure reproducibility.
    tray.AddModule('PMTResponseSimulator',"_pmt",**pmt_config) 
    # Adds the keys: I3MCPulseSeriesMap & I3MCPulseSeriesMapParticleIDMap
    
    # 3.2) SIMULATE DOM RESPONSE
    dom_config = {'Input':'I3MCPulseSeriesMap',                # The input is the simulated PMT responses. Used to model how the DOM hardware processes the pulses.
                  'Output': "SignalRawData",                   # Will create a new key named SignalRawData containing the raw output from the DOM simulation, such as digital signals that include the effects of digitization and electronics noise.
                  'UseTabulatedPT':True,                       # Indicates whether to use pre-tabulated photon timing distributions for efficiency. These tables describe the expected photon arrival times based on the geometry and optical properties of the detector.
                  'RandomServiceName' : 'I3RandomService',
                  }
    tray.AddModule('DOMLauncher','_dommb',**dom_config)
    # Adds the key: SignalRawData

    
    # 3.3) ADD NOISE WITH THE FRT MERGER
    # The FRT merger module adds a key named 'InIceRawData', which is needede for further processing.
    # In case you turned off the background simulation, we thus need to create the InIceRawData frame by copying 'SignalRawData'
    class add_key(icetray.I3Module):
        def DAQ(self, frames):
            frames['InIceRawData'] = frames['SignalRawData']
            self.PushFrame(frames)
    if FRT:                        # FRTMerger takes SignalRawData as input (the DOM signal) and returns a new key named "InIceRawData"
        tray.AddModule(FRTMerger)  # It also adds the keys: FrtMerger_Error, FrtMerger_FRTFile, FrtMerger_Rejected
    else:                          # If not using FRT, we need to create the 'InIceRawData' key, which we create from 'SignalRawData'
        tray.Add(add_key)
    # It adds the keys: InIceRawData, (FrtMerger_Error, FrtMerger_FRTFile, FrtMerger_Rejected)
    
    # 3.4) ICECUBE TRIGGER HANDLING
    # 3.4.1) Create a copy of the 'InIceRawData' key before TriggerSim is executed, because TriggerSim modifies InIceRawData to only keep what's inside the trigger
    class CopyInIceRawData(icetray.I3Module):
        def DAQ(self, frame):
            if 'InIceRawData' in frame:
                frame['SignalNoiseRawData'] = frame['InIceRawData']
            else:
                icetray.logging.log_warn("InIceRawData key not found in frame, no backup created.")
            self.PushFrame(frame)
    tray.AddModule(CopyInIceRawData, 'copy_inice_rawdata')
    # It adds the key: SignalNoiseRawData

    # 3.4.2) Class to restrict trigger_sim to only SLOP trigger
    class OnlySLOPTrigger(icetray.I3Module):
        def DetectorStatus(self, frame):
            old_TriggerMap = frame["I3DetectorStatus"].trigger_status
            new_triggerMap = dataclasses.I3TriggerStatusMap()
            for trigger_key, trigger_config in old_TriggerMap.items():
                if trigger_key.config_id == 24002:
                    new_triggerMap[trigger_key] = trigger_config
            frame["I3DetectorStatus"].trigger_status = new_triggerMap
            self.PushFrame(frame)
    if OnlySLOP: tray.AddModule(OnlySLOPTrigger, "only_slop_triggers")
    
    # 3.4.3) add TriggerSim segment
    tray.AddSegment(trigger_sim.TriggerSim,       # Adds the TriggerSim module to the IceTray framework. A "segment" in IceTray is essentially a higher-level processing step, composed of multiple low-level modules or operations.
                '_triggersim',                    # This is the label or instance name given to the TriggerSim module in the processing pipeline. It allows referencing this module later in the pipeline if needed.
                gcd_file=dataio.I3File(gcdFile),  # Detector geometry. This is required for trigger auto-configuration, as it ensures the triggers are applied correctly based on the detector’s state.
                run_id = 1234,
                prune = True,
                time_shift = True,                # Enables time-shifting of events. In IceCube, simulated events often need to be synchronized to match real-time conditions
                filter_mode = throw_untrig)       # If True, will throw away all events that do not pass any trigger conditions
    # It adds the keys: I3EventHeader, I3TriggerHierarchy, I3Triggers, TimeShift
    
    
    # 3.5) Calibrate raw waveform data to produce physics-useful "pulses" that signify particle hits in the detector
    tray.AddModule('I3WaveCalibrator', '_wavecal',              # Converts raw waveforms (low-level signal outputs from DOMs) into calibrated waveforms that account for detector-specific corrections. This step ensures that the waveforms are in a form suitable for further processing, such as extracting pulses.
                Launches='SignalRawData',                       # Specifies the input data source containing the raw waveform data (DOM launches). These are the signals recorded by the DOMs when they detect light, stored in the frame under the key SignalRawData
                Waveforms='CalibratedWaveforms',                # Specifies the output key where the calibrated waveforms will be stored in the frame. Calibration involves correcting for DOM-specific effects like time delays, amplitude scaling, and baseline shifts.
                Errata='CalibrationErrata',                     # Stores any errors or anomalies encountered during calibration, such as DOMs with bad calibration data or unexpected waveform properties. This data can be used to filter out or flag problematic DOMs.
                WaveformRange='CalibratedWaveformRange')        # Specifies the time range of the calibrated waveforms. This range is determined during calibration and ensures downstream modules only process valid portions of the waveform.
    # Adds the keys: CalibratedWaveformRange, CalibratedWaveforms

    tray.AddModule('I3Wavedeform', '_wavedeform',               # processes the calibrated waveforms and extracts pulses, which are discrete time points representing particle interactions in the detector. This step is crucial for reconstructing particle events.
                Waveforms='CalibratedWaveforms',                # Specifies the input key where the calibrated waveforms (output from I3WaveCalibrator) are read. These are the cleaned and corrected waveforms that are ready for pulse extraction.
                WaveformTimeRange='CalibratedWaveformRange',    # Specifies the valid time range of the calibrated waveforms, ensuring only valid portions are analyzed for pulse extraction.
                Output='InIcePulses')                           # Specifies the output key where the extracted pulses are stored in the frame. These pulses represent the light hits detected by DOMs and are the primary data used for reconstructing the energy, direction, and type of particle interactions in IceCube.
    # Adds the keys: InIcePulses, InIcePulsesTimeRange


    ###### 4) WRITE THE OUTPUT
    file_name = f"part{part}_charge{q}_beta{b_}_lambda{c}_nevents{n_events}_{numb}_slop.i3.zst"
    tray.Add("I3Writer", Filename=os.path.join(output_path, file_name),
            Streams = [icetray.I3Frame.DAQ, icetray.I3Frame.TrayInfo, icetray.I3Frame.Geometry, icetray.I3Frame.DetectorStatus, icetray.I3Frame.Calibration])
    tray.Execute(n_events+3)
    
print(f'\nThe output has been written. Open it with: dataio-pyshovel {output_path}/{file_name}\n')
endTime = datetime.now()
print(f'Simulation ending time: {endTime}')
print(f'The running time was {endTime-startTime}')


# [1]: https://wiki.icecube.wisc.edu/index.php/Luminescence_measurement_in_SpiceCore