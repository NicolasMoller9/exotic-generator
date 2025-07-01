#!/usr/bin/env python3
import unittest
from icecube import icetray, dataclasses, dataio
from icecube import exotic_generator
from icecube import icetray
from icecube.icetray import I3Units
import math
import numpy
import scipy.stats


class TestExoticGeneratorKnownBugTests(unittest.TestCase):
    def test_max_length_of_particles_is_correct(self):
        tray = icetray.I3Tray()

        def frame_tester(fr, unit_tester=self):
            for particle in fr["I3MCTree"]:
                if particle.type == "monopole":
                    unit_tester.assertEqual(particle.length, 10 * I3Units.m)

        tray.AddModule("I3InfiniteSource", "infinite", Stream=icetray.I3Frame.DAQ)
        tray.AddService("I3MTRandomServiceFactory", "random")
        tray.AddModule("I3ExoticGenerator", "generator", Mass=1e12*I3Units.GeV, Nevents=1, BetaRange=[0.5])
        tray.AddModule("I3ExoticPropagator",
                       MaxLength             = 10 * I3Units.m,
                       MinLength             = 10 * I3Units.m,
                       StepSize              = float("NAN"),
                       )
        tray.AddModule(frame_tester, Streams=[icetray.I3Frame.DAQ])
        tray.Execute()


class TestExoticPropagator_working_usecases(unittest.TestCase):
    def test_max_length_of_particles_is_correct(self):
        tray = icetray.I3Tray()
        tray.AddModule("I3InfiniteSource", "infinite", Stream=icetray.I3Frame.DAQ)
        tray.AddService("I3MTRandomServiceFactory", "random")
        tray.AddModule("I3ExoticGenerator", "generator", Mass=1e11*I3Units.GeV, Nevents=1, BetaRange=[0.5], Disk_dist=1000 * I3Units.m, Disk_rad = 0)
        tray.AddModule("I3ExoticPropagator",
                       MaxLength             = 10 * I3Units.m,
                       MinLength             = 1 * I3Units.m,
                       StepSize              = float("NAN"),
                       )
        tray.Execute()

    def test_check_that_Stepsize_gets_properly_overwritten(self):
        tray = icetray.I3Tray()
        def frame_tester(fr, unit_tester=self):
            for particle in fr["I3MCTree"]:
                if particle.type == "monopole":
                    unit_tester.assertEqual(particle.length, 10 * I3Units.m)

        tray.AddModule("I3InfiniteSource", "infinite", Stream=icetray.I3Frame.DAQ)
        tray.AddService("I3MTRandomServiceFactory", "random")
        tray.AddModule("I3ExoticGenerator", "generator", Mass=1e12*I3Units.GeV, Nevents=1, BetaRange=[0.5])
        tray.AddModule("I3ExoticPropagator",
                       MaxLength             = 10 * I3Units.m,
                       MinLength             = 10 * I3Units.m,
                       StepSize              = 100 * I3Units.m,
                       )
        tray.AddModule(frame_tester, Streams=[icetray.I3Frame.DAQ])
        tray.Execute()


class TestMonopoleLogFatalShouldWork(unittest.TestCase):
    def test_max_length_of_particles_is_correct(self):
        tray = icetray.I3Tray()
        tray.AddModule("I3InfiniteSource", "infinite", Stream=icetray.I3Frame.DAQ)
        tray.AddService("I3MTRandomServiceFactory", "random")
        tray.AddModule("I3ExoticGenerator", "generator", Mass=1e11*I3Units.GeV, Nevents=1, BetaRange=[0.5], Disk_dist=1000 * I3Units.m, Disk_rad = 0)
        tray.AddModule("I3ExoticPropagator",
                           MaxLength             = 1 * I3Units.m,
                           MinLength             = 10 * I3Units.m,
                           StepSize              = float("NAN"),
                           )
        with self.assertRaisesRegex(Exception, "Oops, MaxLength<MinLength."):
            tray.Execute()

#TODO
"""
* This test is designed to comprehensively check the interface of the Propagate module through various ways
* (1) That it only propagates exotics (done by adding non-exotic primary/secondaries to test generation tree)
* (2) That the number of children remains as expected - constant if no proton decay, and the correct number of
*     segments if they are monopoles (the num_child tests)
* (3) That setting the exotic to beyond the max length will cause it to propagate beyond the detector the same amount
*     as it started (tree1 above)
* (4) That, otherwise, it will propagate to the set max distance from the detector center (tree3 above)
* (5) That direction is preserved
* (6) That when CalculateNextLength gives something larger than user-defined max length,
*     length is set to this instead (tree 1)
*/
"""


if __name__ == '__main__':
    unittest.main()
