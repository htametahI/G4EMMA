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

#include "EMMAS3Hit.hh"      // Local hit declaration.

#include "G4UnitsTable.hh"   // G4BestUnit helper.
#include "G4ios.hh"          // Geant4 stream definitions.

#include <iomanip>           // std::setw formatting.

G4Allocator<EMMAS3Hit> EMMAS3HitAllocator; // Define Geant4 allocator instance.

EMMAS3Hit::EMMAS3Hit()
  : G4VHit(),               // Construct base hit.
    fEdep(0.0),             // Initialize deposited energy accumulator.
    fKineticEnergy(0.0),    // Initialize kinetic energy storage.
    fTheta(0.0),            // Initialize theta storage.
    fPhi(0.0),              // Initialize phi storage.
    fHasKinematics(false)   // Initialize kinematics flag.
{
}

EMMAS3Hit::~EMMAS3Hit()
{
}

EMMAS3Hit::EMMAS3Hit(const EMMAS3Hit& right)
  : G4VHit()                  // Construct base hit.
{
  fEdep = right.fEdep;        // Copy accumulated deposited energy.
  fKineticEnergy = right.fKineticEnergy; // Copy stored kinetic energy.
  fTheta = right.fTheta;      // Copy stored theta.
  fPhi = right.fPhi;          // Copy stored phi.
  fHasKinematics = right.fHasKinematics; // Copy "kinematics recorded" flag.
}

const EMMAS3Hit& EMMAS3Hit::operator=(const EMMAS3Hit& right)
{
  if (this != &right) {       // Guard against self-assignment.
    fEdep = right.fEdep;      // Copy accumulated deposited energy.
    fKineticEnergy = right.fKineticEnergy; // Copy stored kinetic energy.
    fTheta = right.fTheta;    // Copy stored theta.
    fPhi = right.fPhi;        // Copy stored phi.
    fHasKinematics = right.fHasKinematics; // Copy "kinematics recorded" flag.
  }
  return *this;               // Return assigned object reference.
}

G4int EMMAS3Hit::operator==(const EMMAS3Hit& right) const
{
  return (this == &right) ? 1 : 0; // Pointer equality semantics, consistent with existing hits.
}

void EMMAS3Hit::Print()
{
  G4cout
    << "S3 hit: Edep=" << std::setw(7) << G4BestUnit(fEdep, "Energy") // Print deposited energy.
    << " Ekin=" << std::setw(7) << G4BestUnit(fKineticEnergy, "Energy") // Print entry kinetic energy.
    << " theta(rad)=" << fTheta   // Print stored theta in radians.
    << " phi(rad)=" << fPhi       // Print stored phi in radians.
    << G4endl;
}
