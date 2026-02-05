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

#ifndef EMMAS3Hit_h
#define EMMAS3Hit_h 1

#include "G4Allocator.hh"       // Geant4 hit allocator support.
#include "G4THitsCollection.hh" // Geant4 typed hits collection.
#include "G4VHit.hh"            // Geant4 base hit class.
#include "globals.hh"           // Geant4 basic typedefs.

// This hit stores per-ring triton observables for the upstream S3 detector.
class EMMAS3Hit : public G4VHit
{
  public:
    EMMAS3Hit();                               // Default constructor.
    EMMAS3Hit(const EMMAS3Hit& right);         // Copy constructor.
    virtual ~EMMAS3Hit();                      // Destructor.

    const EMMAS3Hit& operator=(const EMMAS3Hit& right); // Assignment operator.
    G4int operator==(const EMMAS3Hit& right) const;     // Equality operator.

    inline void* operator new(size_t);         // Custom allocator new.
    inline void  operator delete(void* hit);   // Custom allocator delete.

    virtual void Draw() {}                     // No custom drawing for this hit.
    virtual void Print();                      // Print helper for debug output.

    void AddEdep(G4double edep);               // Accumulate deposited energy.
    void SetKinematicsIfUnset(                 // Store entry kinematics once.
      G4double kineticEnergy,
      G4double theta,
      G4double phi);

    G4double GetEdep() const;                  // Getter for accumulated energy deposit.
    G4double GetKineticEnergy() const;         // Getter for entry kinetic energy.
    G4double GetTheta() const;                 // Getter for entry theta.
    G4double GetPhi() const;                   // Getter for entry phi.
    G4bool HasKinematics() const;              // Getter for "kinematics recorded" flag.

  private:
    G4double fEdep;                            // Total energy deposited in this ring.
    G4double fKineticEnergy;                   // Triton kinetic energy at first hit in this ring.
    G4double fTheta;                           // Triton polar angle at first hit in this ring.
    G4double fPhi;                             // Triton azimuthal angle at first hit in this ring.
    G4bool fHasKinematics;                     // True after first-hit kinematics is stored.
};

typedef G4THitsCollection<EMMAS3Hit> EMMAS3HitsCollection; // Event-level collection type.

extern G4Allocator<EMMAS3Hit> EMMAS3HitAllocator; // Geant4 allocator instance declaration.

inline void* EMMAS3Hit::operator new(size_t)
{
  void* hit = (void*)EMMAS3HitAllocator.MallocSingle(); // Allocate one hit from Geant4 allocator.
  return hit;                                            // Return allocated pointer.
}

inline void EMMAS3Hit::operator delete(void* hit)
{
  EMMAS3HitAllocator.FreeSingle((EMMAS3Hit*)hit); // Release one hit to Geant4 allocator.
}

inline void EMMAS3Hit::AddEdep(G4double edep)
{
  fEdep += edep; // Add this step contribution to the total deposited energy.
}

inline void EMMAS3Hit::SetKinematicsIfUnset(
  G4double kineticEnergy,
  G4double theta,
  G4double phi)
{
  if (!fHasKinematics) {     // Only record entry kinematics the first time this ring is hit.
    fKineticEnergy = kineticEnergy; // Save kinetic energy at ring entry.
    fTheta = theta;                 // Save theta at ring entry.
    fPhi = phi;                     // Save phi at ring entry.
    fHasKinematics = true;          // Mark kinematics as initialized.
  }
}

inline G4double EMMAS3Hit::GetEdep() const
{
  return fEdep; // Return accumulated deposited energy.
}

inline G4double EMMAS3Hit::GetKineticEnergy() const
{
  return fKineticEnergy; // Return first-hit kinetic energy.
}

inline G4double EMMAS3Hit::GetTheta() const
{
  return fTheta; // Return first-hit theta.
}

inline G4double EMMAS3Hit::GetPhi() const
{
  return fPhi; // Return first-hit phi.
}

inline G4bool EMMAS3Hit::HasKinematics() const
{
  return fHasKinematics; // Return whether first-hit kinematics has been recorded.
}

#endif
