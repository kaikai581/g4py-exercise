#!/usr/bin/env python
'''
This script simulates monoenergetic 6 MeV protons on beryllium target to generate epithermal neutrons.
The goal is to study neutron kinematics and proton-beryllium interactions.

Shih-Kai Lin, Feb 2021. shihkailin78@gmail.com
'''

# python3 style print
from __future__ import print_function

# Physical Units
from Geant4 import cm, deg, m, MeV

import argparse
import Geant4
import sys

# ------------------------------------------------------------------
# Detector Construction
class MyDetectorConstruction(Geant4.G4VUserDetectorConstruction):
    '''
    My Detector Construction
    '''

    def __init__(self, target_thickness=1000e-6*m):
        Geant4.G4VUserDetectorConstruction.__init__(self)
        self.target_thickness = target_thickness
        # a member variable to store physical volumes is necessary to avoid this error:
        # ReferenceError: Attempt to return dangling pointer to object of type: G4VPhysicalVolume
        # An alternative way used in other examples is to make the volumes global.
        self.vols = []
        
    def Construct(self):
        '''
        Construct geometry
        '''
        # World (cylinder of air)
        world_s = Geant4.G4Tubs('World', 0*cm, 5*m, 5*m, 0*deg, 360*deg)
        world_l = Geant4.G4LogicalVolume(world_s, Geant4.gNistManager.FindOrBuildMaterial('G4_Air'),
                                         'World')
        world_p = Geant4.G4PVPlacement(None,                   #no rotation
                                       Geant4.G4ThreeVector(), #at (0,0,0)
                                       world_l,                #its logical volume
                                       'World',                #its name
                                       None,                   #its mother  volume
                                       False,                  #no boolean operation
                                       0,                      #copy number
                                       True)                   #check overlap
        self.vols += [world_s, world_l, world_p]

        # target (disk of beryllium)
        target_halflen = self.target_thickness/2
        target_s = Geant4.G4Tubs('Target', 0*cm, 6*cm, target_halflen, 0*deg, 360*deg)
        target_l = Geant4.G4LogicalVolume(target_s, Geant4.gNistManager.FindOrBuildMaterial('G4_Be'),
                                          'Target')
        Geant4.G4PVPlacement(None,
                             Geant4.G4ThreeVector(0,0,target_halflen),
                             target_l,
                             'Target',
                             world_l,
                             False,
                             0,
                             True)

        # Return world volume
        return world_p

# ------------------------------------------------------------------
class MyPrimaryGeneratorAction(Geant4.G4VUserPrimaryGeneratorAction):
    '''
    My Primary Generator Action.
    '''

    def __init__(self, zstart=0*m, ke=6*MeV):
        Geant4.G4VUserPrimaryGeneratorAction.__init__(self)
        self.particleGun = Geant4.G4ParticleGun(1)
        self.particleGun.SetParticleByName('proton')
        self.particleGun.SetParticleEnergy(ke)
        self.particleGun.SetParticleMomentumDirection(Geant4.G4ThreeVector(0., 0., 1.))
        self.particleGun.SetParticlePosition(Geant4.G4ThreeVector(0.,0.,zstart))

    def GeneratePrimaries(self, event):
        self.particleGun.GeneratePrimaryVertex(event)

# ------------------------------------------------------------------
class MySimulation:
    '''
    My Simple Simulation class.
    '''
    def __init__(self):
        '''
        The class constructor takes care of run manager construction
        happening in the main function of a C++ counterpart.

        Note: Mandatory classes have to be stored in member variables.
              Otherwise "TypeError: 'NoneType' object is not callable" will result from Geant4.gRunManager.
        '''
        # configure the simulation
        self.args = self.configure_simulation()

        # optional: specify the random seed used by this run
        if self.args.random_seed:
            Geant4.HepRandom.setTheSeed(self.args.random_seed)

        # construct geometry and assign to run manager
        self.detector = MyDetectorConstruction(target_thickness=self.args.target_thickness*m)
        Geant4.gRunManager.SetUserInitialization(self.detector)

        # specify physics list by name
        # factory = Geant4.G4PhysListFactory()
        # physlist_name = self.args.physlist
        # if not factory.IsReferencePhysList(physlist_name):
        #     print('Provided reference physics list name {} does not exist.'.format(physlist_name))
        #     print('Please use a valid name and run the simulation again.')
        #     sys.exit(-1)
        # self.physlist = factory.GetReferencePhysList(physlist_name)
        self.physlist = Geant4.QBBC()
        Geant4.gRunManager.SetUserInitialization(self.physlist)
    
        # set event generator actions
        self.pga = MyPrimaryGeneratorAction()
        Geant4.gRunManager.SetUserAction(self.pga)

        # study secondary paricles produced by proton-beryllium interactions
        self.stackingaction = MyStackingAction()
        Geant4.gRunManager.SetUserAction(self.stackingaction)

        # initialize
        Geant4.gRunManager.Initialize()
    
    def configure_simulation(self):
        '''
        Configure this simulation with command line arguments.
        '''
        parser = argparse.ArgumentParser(description='Particle beam on fixed target simulation.')
        parser.add_argument('-n', '--nevents', help='Number of beam particles to simulate', type=int, default=10)
        parser.add_argument('-pl', '--physlist', help='Name of the reference physics list', type=str, default='QBBC')
        parser.add_argument('-rs', '--random_seed', help='Random seed', type=int)
        parser.add_argument('-tt', '--target_thickness', help='Thickness of the beryllium target, unit in meter', type=float, default=1000e-6)

        return parser.parse_args()
    
    def run(self):
        '''
        Simulate!
        '''
        # run!
        Geant4.gRunManager.BeamOn(self.args.nevents)

        # close
        Geant4.gTerminate()

# ------------------------------------------------------------------
class MyStackingAction(Geant4.G4UserStackingAction):
    '''
    My stacking action for recording neutrons generated by
    proton-beryllium interactions.
    '''
    def __init__(self):
        '''
        Constructor.
        '''
        Geant4.G4UserStackingAction.__init__(self)
        self.nprim = 1
    
    def ClassifyNewTrack(self, track):
        # if it's a primary particle, continue
        if track.GetParentID() == 0:
            # print(track.GetDefinition().GetParticleName(), self.nprim)
            self.nprim += 1
            return Geant4.G4ClassificationOfNewTrack.fUrgent
        
        # print secondary information
        eid = Geant4.gRunManager.GetRunManager().GetCurrentEvent().GetEventID()
        print('Secondary particle found!', track.GetDefinition().GetParticleName())
        print('\tEvent ID: {}\tParent ID: {}\tTrack ID: {}'.format(eid, track.GetParentID(), track.GetTrackID()))
        
        return Geant4.G4ClassificationOfNewTrack.fKill

if __name__ == '__main__':
    '''
    Run the simulation.
    '''
    mysim = MySimulation()
    mysim.run()
    # mysim.finalize()