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

#include "EMMAS3Detector.hh"  

#include "G4HCofThisEvent.hh"     
#include "G4ParticleDefinition.hh" 
#include "G4SDManager.hh"         
#include "G4Step.hh"              
#include "G4StepPoint.hh"        
#include "G4Track.hh"             
#include "G4ios.hh"               

EMMAS3Detector::EMMAS3Detector(
  const G4String& name,
  const G4String& hitsCollectionName,
  G4int nofRings)
  : G4VSensitiveDetector(name),  
    fHitsCollection(0),          
    fNofRings(nofRings)         
{
  collectionName.insert(hitsCollectionName); 
}

EMMAS3Detector::~EMMAS3Detector()
{
}

void EMMAS3Detector::Initialize(G4HCofThisEvent* hce)
{
  fHitsCollection = new EMMAS3HitsCollection( 
    SensitiveDetectorName,
    collectionName[0]);

  const G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); 
  hce->AddHitsCollection(hcID, fHitsCollection); 

  for (G4int i = 0; i < fNofRings + 1; ++i) { 
    fHitsCollection->insert(new EMMAS3Hit()); 
  }
}

G4bool EMMAS3Detector::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  const G4Track* track = step->GetTrack(); 
  const G4ParticleDefinition* particle = track->GetDefinition(); 
  if (track->GetParentID() != 0) { // Keep only the generated primary triton track to align with target-exit observables.
    return false; // Ignore secondary tritons to avoid mixing track-level quantities in one output row, uncomment this if we need secondary effects
  }

  const G4int z = particle->GetAtomicNumber(); 
  const G4int a = particle->GetAtomicMass();   
  if (z != 1 || a != 3) { // Keep only tritons (Z=1,A=3), regardless of particle-type string.
    return false; 
  }

  const G4double edep = step->GetTotalEnergyDeposit(); // Get this-step deposited energy in S3 silicon.
  if (edep < 0.0) { 
    return false; 
  }

  const G4StepPoint* preStepPoint = step->GetPreStepPoint(); // Access pre-step state.
  const G4int ringNumber = preStepPoint->GetTouchableHandle()->GetCopyNumber(); // Current ring copy number.

  if (ringNumber < 0 || ringNumber >= fNofRings) { // Guard against unexpected geometry copy numbers.
    G4cerr << "EMMAS3Detector: invalid ring copy number " << ringNumber << G4endl; 
    return false; // Ignore malformed hit indexing.
  }

  EMMAS3Hit* ringHit = (*fHitsCollection)[ringNumber]; // Resolve per-ring hit object.
  EMMAS3Hit* totalHit = (*fHitsCollection)[fHitsCollection->entries() - 1]; // Resolve total accumulator hit.

  const G4double kineticEnergy = preStepPoint->GetKineticEnergy(); 
  const G4ThreeVector direction = preStepPoint->GetMomentumDirection(); 
  const G4double theta = direction.theta();
  const G4double phi = direction.phi();     

  ringHit->AddEdep(edep); // Accumulate deposited energy in this ring.
  ringHit->SetKinematicsIfUnset(kineticEnergy, theta, phi); // Record first-hit kinematics for this ring.

  totalHit->AddEdep(edep); // Accumulate total deposited energy over all rings.
  totalHit->SetKinematicsIfUnset(kineticEnergy, theta, phi); // Keep first-hit kinematics for total slot too.

  return true; // Signal that this step produced a valid S3 hit update.
}

void EMMAS3Detector::EndOfEvent(G4HCofThisEvent*)
{
  if (verboseLevel > 1) { 
    const G4int nofHits = fHitsCollection->entries();
    G4cout << "\n-------->S3 Hits Collection has " << nofHits << " entries." << G4endl; 
    for (G4int i = 0; i < nofHits; ++i) { 
      (*fHitsCollection)[i]->Print(); 
    }
  }
}
