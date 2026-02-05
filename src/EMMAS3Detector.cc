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

#include "EMMAS3Detector.hh"      // Local S3 SD declaration.

#include "G4HCofThisEvent.hh"     // Event hit collection container.
#include "G4ParticleDefinition.hh" // Particle definition access.
#include "G4SDManager.hh"         // SD manager interface.
#include "G4Step.hh"              // Step information.
#include "G4StepPoint.hh"         // Step-point information.
#include "G4Track.hh"             // Track information.
#include "G4ios.hh"               // Geant4 I/O stream definitions.

EMMAS3Detector::EMMAS3Detector(
  const G4String& name,
  const G4String& hitsCollectionName,
  G4int nofRings)
  : G4VSensitiveDetector(name),  // Build base SD.
    fHitsCollection(0),          // Initialize hit collection pointer.
    fNofRings(nofRings)          // Cache ring count for bounds checks and per-event allocation.
{
  collectionName.insert(hitsCollectionName); // Register one hits collection name for this SD.
}

EMMAS3Detector::~EMMAS3Detector()
{
}

void EMMAS3Detector::Initialize(G4HCofThisEvent* hce)
{
  fHitsCollection = new EMMAS3HitsCollection( // Allocate fresh per-event S3 hit collection.
    SensitiveDetectorName,
    collectionName[0]);

  const G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); // Resolve collection id.
  hce->AddHitsCollection(hcID, fHitsCollection); // Register collection in current event container.

  for (G4int i = 0; i < fNofRings + 1; ++i) { // Create one hit slot per ring plus one total slot.
    fHitsCollection->insert(new EMMAS3Hit()); // Insert empty hit object.
  }
}

G4bool EMMAS3Detector::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  const G4Track* track = step->GetTrack(); // Access current track.
  const G4ParticleDefinition* particle = track->GetDefinition(); // Access particle identity.
  if (track->GetParentID() != 0) { // Keep only the generated primary triton track to align with target-exit observables.
    return false; // Ignore secondary tritons to avoid mixing track-level quantities in one output row.
  }

  const G4int z = particle->GetAtomicNumber(); // Particle atomic number (works for ions/light ions).
  const G4int a = particle->GetAtomicMass();   // Particle atomic mass number.
  if (z != 1 || a != 3) { // Keep only tritons (Z=1,A=3), regardless of particle-type string.
    return false; // Ignore all nuclei except tritons.
  }

  const G4double edep = step->GetTotalEnergyDeposit(); // Get this-step deposited energy in S3 silicon.
  if (edep < 0.0) { // Guard only against unphysical negative deposits.
    return false; // Ignore invalid negative-energy steps.
  }

  const G4StepPoint* preStepPoint = step->GetPreStepPoint(); // Access pre-step state.
  const G4int ringNumber = preStepPoint->GetTouchableHandle()->GetCopyNumber(); // Current ring copy number.

  if (ringNumber < 0 || ringNumber >= fNofRings) { // Guard against unexpected geometry copy numbers.
    G4cerr << "EMMAS3Detector: invalid ring copy number " << ringNumber << G4endl; // Debug warning.
    return false; // Ignore malformed hit indexing.
  }

  EMMAS3Hit* ringHit = (*fHitsCollection)[ringNumber]; // Resolve per-ring hit object.
  EMMAS3Hit* totalHit = (*fHitsCollection)[fHitsCollection->entries() - 1]; // Resolve total accumulator hit.

  const G4double kineticEnergy = preStepPoint->GetKineticEnergy(); // Capture triton kinetic energy at entry.
  const G4ThreeVector direction = preStepPoint->GetMomentumDirection(); // Capture triton direction at entry.
  const G4double theta = direction.theta(); // Compute polar angle in radians.
  const G4double phi = direction.phi();     // Compute azimuthal angle in radians.

  ringHit->AddEdep(edep); // Accumulate deposited energy in this ring.
  ringHit->SetKinematicsIfUnset(kineticEnergy, theta, phi); // Record first-hit kinematics for this ring.

  totalHit->AddEdep(edep); // Accumulate total deposited energy over all rings.
  totalHit->SetKinematicsIfUnset(kineticEnergy, theta, phi); // Keep first-hit kinematics for total slot too.

  return true; // Signal that this step produced a valid S3 hit update.
}

void EMMAS3Detector::EndOfEvent(G4HCofThisEvent*)
{
  if (verboseLevel > 1) { // Optional verbose printout for debugging.
    const G4int nofHits = fHitsCollection->entries(); // Number of allocated hit slots this event.
    G4cout << "\n-------->S3 Hits Collection has " << nofHits << " entries." << G4endl; // Header line.
    for (G4int i = 0; i < nofHits; ++i) { // Loop over all hit slots.
      (*fHitsCollection)[i]->Print(); // Print each hit slot.
    }
  }
}
