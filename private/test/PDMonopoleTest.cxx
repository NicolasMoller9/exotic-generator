/**
 * Monopole Secondaries Test Suite
 * Copyright (c) 2004 IceCube Collaboration
 * Version $Id: PDMonopoleTest.cxx 125118 2014-10-27 17:22:09Z jacobi $
 *
 * @file MonopoleSecondariesTest.cxx
 * @date $Date: 2008-01-18 00:43:41 (Fri, 18 Jan 2008)$
 * @author Alex Olivas <olivas@icecube.umd.edu>
 * @author Brian Christy <bchristy@icecube.umd.edu>
 * @author Emanuel Jacobi <emanuel.jacobi@desy.de>
 * @brief Tests for the Monopoles catalyzing Proton Decay
 * @todo Test tree name interface
 * @todo Test no random service produces log_fatal
 */

#include <I3Test.h>

#include <exotic-generator/I3ExoticPropagatorUtils.h>
#include <exotic-generator/I3ExoticPropagator.h>
#include <exotic-generator/I3ExoticGenerator.h>


#include <icetray/I3Frame.h>
#include <icetray/I3Tray.h>
#include <icetray/I3Units.h>
#include <dataclasses/I3Double.h>
#include <dataclasses/physics/I3Particle.h>
#include <dataclasses/physics/I3MCTree.h>
#include <dataclasses/physics/I3MCTreeUtils.h>
#include <dataclasses/I3Constants.h>
#include <phys-services/I3MTRandomService.h>

TEST_GROUP(slow - monopole);

namespace exo_sec_io_test {
    /**
     *Takes a mctree from the frame that contains a
     * propagated SLOP and checks that everything
     * has expected behavior
     * Used by GeneralSecondariesTest below
     */
    class TestSecondaries : public I3Module {
    public:
        TestSecondaries(const I3Context &ctx) : I3Module(ctx) { AddOutBox("OutBox"); }

        void DAQ(I3FramePtr frame) {
          const I3MCTree &testTree = frame->Get<I3MCTree>();
          //I3MCTree::iterator tree_iter;
          int num_sec_child = 0;

          const std::vector<I3Particle> primaries = I3MCTreeUtils::GetPrimaries(testTree);
          ENSURE(primaries.size() == 1);
          ENSURE(primaries[0].GetType() == I3Particle::Monopole);
          ENSURE_EQUAL(primaries[0].GetLength(), 200 * I3Units::m);
          for(auto const &tree_iter: testTree){
            if (I3MCTreeUtils::HasParent(testTree, tree_iter)) {
              if (tree_iter.GetType() == I3Particle::EPlus) {
                ++num_sec_child;
                const I3Particle &parent = I3MCTreeUtils::GetParent(testTree, tree_iter);
                ENSURE(parent.GetType() == I3Particle::Monopole);
                double dist_from_parent = sqrt(pow(parent.GetX() - tree_iter.GetX(), 2)
                                               + pow(parent.GetY() - tree_iter.GetY(), 2)
                                               + pow(parent.GetZ() - tree_iter.GetZ(), 2));
                ENSURE(dist_from_parent <= 200 * I3Units::m);
              } else {
                FAIL("Was only expecting delta electrons...has propagator been extended?  If so, need to update this test");
              }
            }
          }//End For loop over tree
          int Num_Exp = static_cast<int>(200 / 0.715115);
          double sigma = sqrt(Num_Exp);
          ENSURE_DISTANCE(num_sec_child, Num_Exp, 3 * sigma);
          PushFrame(frame, "OutBox");
        }
    };

    /**
     * Used to make sure that scaling energy option working correctly
     */
    class CheckSecEnergy : public I3Module {
    public :
        CheckSecEnergy(const I3Context &ctx) :
                I3Module(ctx),
                expEnergy_(NAN) {
          AddParameter("ExpectedEnergy", "Expected Energy", expEnergy_);
          AddOutBox("OutBox");
        }

        void Configure() {
          GetParameter("ExpectedEnergy", expEnergy_);
          if (isnan(expEnergy_)) {
            FAIL("oops - test module didn't know what to expect");
          }
        }

        void DAQ(I3FramePtr frame) {
          for (auto const &tree_iter : frame->Get<I3MCTree>()) {
            if (tree_iter.GetType() != I3Particle::Monopole) {
              //ENSURE_EQUAL(tree_iter->GetEnergy(),expEnergy_);
              ENSURE_DISTANCE(tree_iter.GetEnergy(), expEnergy_,
                              1e-8);  // If mfp is too low, ENSURE_EQUAL fails due to rounding error
            }
          }
          PushFrame(frame, "OutBox");
        }

    private:
        double expEnergy_;
    };
}
I3_MODULE(exo_sec_io_test::TestSecondaries);
I3_MODULE(exo_sec_io_test::CheckSecEnergy);

void TestingMFPConstraints(double mfpTest, bool shouldPass) {
  const int NFRAMES(3);
  bool pass(true);

  std::vector<double> betaRange{0.001};//10^-3

  try {
    I3Tray tray;
    tray.AddService("I3MTRandomServiceFactory", "random");

    tray.AddModule("I3InfiniteSource", "infinite");
    tray.AddModule("I3ExoticGenerator", "mono-gen")
            ("NEvents", 1)
            ("BetaRange", betaRange);
    tray.AddModule("I3ExoticPropagator", "exo-prop")
            ("MeanFreePath", mfpTest * I3Units::m);

    tray.AddModule("TrashCan", "trash");

    tray.Execute(NFRAMES + 4);
    tray.Finish();
  }
  catch (...) {
    pass = false;
  }

  if (shouldPass != pass) { FAIL("MFP Test Passed when it should have failed or vice versa"); }

}

void TestingScaleEnergy(double mfpTest, bool scaleTest, double energyExp) {
  const int NFRAMES(3);

  std::vector<double> betaRange{0.001};//10^-3

  try {
    I3Tray tray;

    tray.AddService("I3MTRandomServiceFactory", "random");

    tray.AddModule("I3InfiniteSource", "infinite");
    tray.AddModule("I3ExoticGenerator", "mono-gen")
            ("NEvents", 1)
            ("BetaRange", betaRange);
    tray.AddModule("I3ExoticPropagator", "exo-prop")
            ("MeanFreePath", mfpTest * I3Units::meter)
            ("ScaleEnergy", scaleTest);
    tray.AddModule("exo_sec_io_test::CheckSecEnergy", "sec-test")
            ("ExpectedEnergy", energyExp);
    tray.AddModule("TrashCan", "trash");
    tray.Execute(NFRAMES + 4);
    tray.Finish();
  } catch (...) {
    FAIL("Scale Energy Test has something wrong with script");
  }
}

// # My guess is that short mean free paths leads to either long running times where
// # ctest times out (mildly OK) or there is a region of phase space that causes
// # things to go horribly pear-shaped (very bad).
// # In two tests TestMFP and TestScaling I set the mpfs to something reasonably
// # long while keeping the test conditions intact.

TEST(TestMFP) {
  TestingMFPConstraints(-2, false);
  TestingMFPConstraints(0, false);
  TestingMFPConstraints(NAN, false);
  I3MTRandomServicePtr random(new I3MTRandomService(time(NULL)));
  double mfpTest{1};
  TestingMFPConstraints(mfpTest, true);
}

TEST(TestScaling) {
  I3MTRandomServicePtr random(new I3MTRandomService(time(NULL)));
  double mfpTest{1};
  TestingScaleEnergy(mfpTest, false, 0.94 * I3Units::GeV);
  TestingScaleEnergy(mfpTest, true, (0.94 * I3Units::GeV) / (mfpTest));
  //Check that it doesn't scale energy when mfp>1
  double mfpTest2{1.1};
  TestingScaleEnergy(mfpTest2, false, 0.94 * I3Units::GeV);
  TestingScaleEnergy(mfpTest2, true, 0.94 * I3Units::GeV);
}

TEST(GeneralSecondariesTest) {
  const int NFRAMES(3);

  std::vector<double> betaRange{0.001};//10^-3

  I3Tray tray;

  tray.AddService("I3MTRandomServiceFactory", "random");

  tray.AddModule("I3InfiniteSource", "infinite");
  tray.AddModule("I3ExoticGenerator", "mono-gen")
          ("NEvents", 100)
          ("Disk_dist", 395 * I3Units::m)
          ("BetaRange", betaRange) //Corresponds to beta of 0.01
          ("Length", 200 * I3Units::m);
  tray.AddModule("I3ExoticPropagator", "exo-prop")
          ("MeanFreePath", 0.715115 * I3Units::meter);

  tray.AddModule("exo_sec_io_test::TestSecondaries", "test");
  tray.AddModule("TrashCan", "trash");

  tray.Execute(NFRAMES + 4);
  tray.Finish();

}

TEST(DefaultSecondariesTest) {
  const int NFRAMES(3);

  std::vector<double> betaRange{0.001};//10^-3

  I3Tray tray;
  bool pass(true);
  try {

    tray.AddService("I3MTRandomServiceFactory", "random");

    tray.AddModule("I3InfiniteSource", "infinite");
    tray.AddModule("I3ExoticGenerator", "mono-gen")
            ("NEvents", 100)
            ("BetaRange", betaRange);
    tray.AddModule("I3ExoticPropagator", "exo-prop")
            ("MeanFreePath", 0.5 * I3Units::meter);
    tray.AddModule("TrashCan", "trash");

    tray.Execute(NFRAMES + 4);
    tray.Finish();
  } catch (...) { pass = false; }
  if (!pass) { FAIL("Something went wrong with just a default set of modules!"); }
}
