#!/usr/bin/env python3
#
#  Exotic Generator
#
#  $Id$
#
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import collections
import numpy as np
from icecube import icetray, dataclasses, dataio, exotic_generator
from scipy.optimize import curve_fit


# config
nevents = int(1e5)  
part_id = -41       # monopole
mctreename = "I3MCTree"
# betaRange=[0.001]
betaRange = [0.001, 0.995]
powerLawIndex = 0.5
# powerLawIndex=float('nan')
mass = (1e11) * icetray.I3Units.GeV
disk_dist = 1000 * icetray.I3Units.m
disk_rad = 850 * icetray.I3Units.m
flux = 1e-16

# tray
tray = icetray.I3Tray()

# Add the random generator service
tray.AddService("I3GSLRandomServiceFactory", "random")

# Make some empty frames
tray.AddModule("I3InfiniteSource", "infinite", Stream=icetray.I3Frame.DAQ)

tray.AddModule("I3ExoticGenerator", "generator")(
    ("NEvents", nevents),
    ("PartID", part_id),
    ("TreeName", mctreename),
    ("Disk_dist", disk_dist),
    ("Disk_rad", disk_rad),
    ("InfoName", "EXOInfoDict"),
    ("BetaRange", betaRange),
    ("powerLawIndex", powerLawIndex),
    ("Mass", mass),
    # Those parameters are optional. Here's some choices though
    #   ("ZenithRange", [0*icetray.I3Units.degree,180*icetray.I3Units.degree]),
    #   ("AzimuthRange", [0*icetray.I3Units.degree,360*icetray.I3Units.degree]),

    #   ("Rad_on_disk",0*icetray.I3Units.m),
    #   ("Azi_on_disk",0*icetray.I3Units.radian),
    #   ("Length",1000*icetray.I3Units.m),
)

# tray.AddModule("I3Writer","writer")(
#    ("filename", "exo-gen.i3")
#    )

resultsdict = collections.defaultdict(list)  # type: collections.defaultdict[str,list]
resultsdict["Counter"] = 0  # type: ignore[assignment]


def getResultsFromMCTree(frame):
    if frame.Has(mctreename):
        tree = frame.Get(mctreename)
        particleType = dataclasses.I3Particle.ParticleType.Monopole  # Monopole
        prim = dataclasses.get_most_energetic(tree, particleType)
        if prim is None:
            print(tree)
            print("WARNING: Couldn't find monopole in tree!")
            return False
    else:
        raise RuntimeError("Fatal Error: No %s in Frame", mctreename)

    if frame.Has("EXOInfoDict"):
        exo = frame.Get("EXOInfoDict")
    else:
        raise RuntimeError("Something really went wrong")

    resultsdict["Counter"] += 1  # type: ignore[arg-type]
    resultsdict["zenith"].append(prim.dir.zenith)
    resultsdict["azimuth"].append(prim.dir.azimuth)
    resultsdict["energy"].append(prim.energy)
    resultsdict["fit_status"].append(prim.fit_status)
    resultsdict["location_type"].append(prim.location_type)
    resultsdict["shape"].append(prim.shape)
    resultsdict["speed"].append(exo["Beta"])
    resultsdict["time"].append(prim.time)
    resultsdict["type"].append(prim.type)
    resultsdict["weight"].append(exo["Weight"])
    resultsdict["radius"].append(prim.pos.r)
    resultsdict["phi"].append(prim.pos.phi)
    resultsdict["theta"].append(prim.pos.theta)
    resultsdict["x"].append(prim.pos.x)
    resultsdict["y"].append(prim.pos.y)
    resultsdict["z"].append(prim.pos.z)


tray.AddModule(getResultsFromMCTree, "getResults", Streams=[icetray.I3Frame.DAQ])
tray.Execute()
tray.Finish()

# ---------------------------------------------------------------------------------------


if len(betaRange) == 1:
    infos = "beta=%g, diskDist=%g m, distRad=%g m" % (betaRange[0], disk_dist, disk_rad)
else:
    infos = "betamin=%g,betamax=%g, diskDist=%g m, distRad=%g m" % (
    betaRange[0], betaRange[1], disk_dist, disk_rad)

print("\nValidation for %s: " % infos)
print("Simulated %d events and found %s events in frame" % (nevents, resultsdict["Counter"]))
if nevents != resultsdict["Counter"]:
    print("ERROR: number of simulated events doesn't match!")


def findDeviations(arr):
    dict = {}
    for a in arr:
        if str(a) not in dict:
            dict[str(a)] = 1
        else:
            dict[str(a)] += 1
    #for key, val in dict.iteritems():
    #    print("%s\tfound %d times" % (key, val))
    if len(dict) > 1:
        print("ERROR: More than 1 value found for this variable!")


if len(betaRange) == 1:
    print("\nTest energy:")
    findDeviations(resultsdict["energy"])

print("\nTest location_type:")
findDeviations(resultsdict["location_type"])

print("\nTest shape:")
findDeviations(resultsdict["shape"])

if len(betaRange) == 1:
    print("\nTest speed:")
    findDeviations(resultsdict["speed"])

print("\nTest time:")
findDeviations(resultsdict["time"])

print("\nTest fit_status:")
findDeviations(resultsdict["fit_status"])

print("\nTest type:")
findDeviations(resultsdict["type"])


# ----------------------



def power_law(x, a):
    lam = powerLawIndex
    if lam > 0:
        return a * x ** -lam
    elif lam < 0:
        return a * (sum(betaRange)-x) ** lam
    else:
        print('ERROR: Lambda = 0')

def plotRandomDistributedVariables(filename, title, arr, infos="", weights=None, fit=False, ylim=False):
    print("\nCreate and save a plot for random distributed variable: %s" % title)
    mi = np.min(arr)
    ma = np.max(arr)
    borders = [mi, ma]
    binNumber = 50
    binning = [i * float(borders[1] - borders[0]) / binNumber + borders[0] for i in range(binNumber + 1)]

    if weights is not None:
        if len(weights) != len(arr):
            print(f"WARNING: weights size {len(weights)} does not match arr size {len(arr)}. Ignoring weights.")
            weights = None
        else:
            discRadius = float(disk_rad * 100)  # m -> cm
            discArea = math.pi * discRadius ** 2
            w = flux * discArea * 4 * math.pi
            weights = np.array(weights) / sum(weights)
            weights *= w
        #print(weights)

    hist, bins = np.histogram(arr, bins=binning, weights=weights)
    center = (bins[:-1] + bins[1:]) / 2
    width = 0.7 * (bins[1] - bins[0])

    fig, ax = plt.subplots()
    ax.bar(center, hist, align='center', width=width, color="b")

    if fit:
        # Filter for non-zero bins for fitting
        nonzero = hist > 0
        x_fit = center[nonzero]
        y_fit = hist[nonzero]
        try:
            popt, pcov = curve_fit(power_law, x_fit, y_fit)
            a_fit = popt[0]
        except RuntimeError:
            print("Curve fitting failed.")
            a_fit, b_fit = None, None
        if a_fit is not None:
            x_plot = np.linspace(np.min(x_fit), np.max(x_fit), 500)
            y_plot = power_law(x_plot, *popt)
            if powerLawIndex>0: label=fr"Fit: ${a_fit:.2f}x^{{-\lambda}}$"
            else: label=fr"Fit: ${a_fit:.2f}(x_{{max}}+x_{{min}}-x)^\lambda$"
            ax.plot(x_plot, y_plot, 'r--', label=label, linewidth=2)
            plt.legend()

    ax.set_title(f"{len(arr)} events\n{infos}", fontsize=12, pad=15)
    ax.set_xlabel(title)
    ax.set_ylabel("a.u.")
    if ylim:
        mean = (np.min(hist) + np.max(hist))/2
        diff = np.max(hist) - np.min(hist)
        ymin = max(0., mean - 10*diff)
        ymax = mean + 10*diff
        ax.set_ylim(ymin, ymax)
    ax.tick_params(axis='both', labelsize=12)
    fig.tight_layout()
    fig.savefig(filename)
    print(f"Saved: {filename}")


plotRandomDistributedVariables("cos_zenith.png", "cos(zenith)", np.cos(resultsdict["zenith"]), infos=infos, weights=resultsdict["weight"], ylim=True)

plotRandomDistributedVariables("azimuth.png", "azimuth", np.array(resultsdict["azimuth"]), infos=infos,weights=resultsdict["weight"], ylim=True)

plotRandomDistributedVariables("x.png", "x", np.array(resultsdict["x"]),infos=infos, weights=resultsdict["weight"])

plotRandomDistributedVariables("y.png", "y", np.array(resultsdict["y"]),infos=infos, weights=resultsdict["weight"])

plotRandomDistributedVariables("z.png", "z", np.array(resultsdict["z"]),infos=infos, weights=resultsdict["weight"])

plotRandomDistributedVariables("radius.png", "radius", np.array(resultsdict["radius"]),infos=infos, weights=resultsdict["weight"])

plotRandomDistributedVariables("phi.png", "phi", np.array(resultsdict["phi"]),infos=infos, weights=resultsdict["weight"], ylim=True)

plotRandomDistributedVariables("cos_theta.png", "cos(theta)", np.cos(np.array(resultsdict["theta"])),infos=infos, weights=resultsdict["weight"], ylim=True)

plotRandomDistributedVariables("speed.png", "speed", np.array(resultsdict["speed"]), infos=infos, fit=True)

plotRandomDistributedVariables("speed_weighted.png", "speed (weighted)", np.array(resultsdict["speed"]),infos=infos, weights=resultsdict["weight"], ylim=True)

# With negative power law index:
#plotRandomDistributedVariables("speed_rev.png", "speed", np.array(resultsdict["speed"]), infos=infos, fit=True)

#plotRandomDistributedVariables("speed_rev_weighted.png", "speed (weighted)", np.array(resultsdict["speed"]),infos=infos, weights=resultsdict["weight"], ylim=True)