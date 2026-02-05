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

#include "EMMAS3Hit.hh"          // Upstream S3 hit definition.
#include "G4VSensitiveDetector.hh" // Geant4 sensitive detector base class.

class G4HCofThisEvent;           // Forward declaration for Geant4 event hit collection container.
class G4Step;                    // Forward declaration for Geant4 step type.
class G4TouchableHistory;        // Forward declaration for Geant4 touchable history.

// Sensitive detector for the upstream S3 annular silicon detector.
class EMMAS3Detector : public G4VSensitiveDetector
{
  public:
    EMMAS3Detector(                    // Constructor.
      const G4String& name,            // Sensitive detector name.
      const G4String& hitsCollectionName, // Output hits collection name.
      G4int nofRings);                 // Number of S3 rings.
    virtual ~EMMAS3Detector();         // Destructor.

    virtual void   Initialize(G4HCofThisEvent* hce); // Geant4 per-event hit collection initialization.
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history); // Geant4 step callback.
    virtual void   EndOfEvent(G4HCofThisEvent* hce); // Geant4 end-of-event callback.

  private:
    EMMAS3HitsCollection* fHitsCollection; // Per-event S3 hits collection.
    G4int fNofRings;                       // Number of physical S3 rings.
};

#endif
