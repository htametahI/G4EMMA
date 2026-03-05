//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .        *
// ********************************************************************
//

#ifndef EMMAS3Detector_h
#define EMMAS3Detector_h 1

#include "EMMAS3Hit.hh"          
#include "G4VSensitiveDetector.hh" 

class G4HCofThisEvent;           
class G4Step;                    
class G4TouchableHistory;        

// Sensitive detector for the upstream S3 annular silicon detector.
class EMMAS3Detector : public G4VSensitiveDetector
{
  public:
    EMMAS3Detector(                    
      const G4String& name,           
      const G4String& hitsCollectionName, 
      G4int nofRings);                 
    virtual ~EMMAS3Detector();       

    virtual void   Initialize(G4HCofThisEvent* hce); 
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history); 
    virtual void   EndOfEvent(G4HCofThisEvent* hce); 

  private:
    EMMAS3HitsCollection* fHitsCollection; 
    G4int fNofRings;                       

};

#endif
