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

#include "G4Allocator.hh"       
#include "G4THitsCollection.hh" 
#include "G4VHit.hh"            
#include "globals.hh"           

// This hit stores per-ring triton observables for the upstream S3 detector.
class EMMAS3Hit : public G4VHit
{
  public:
    EMMAS3Hit();                               
    EMMAS3Hit(const EMMAS3Hit& right);        
    virtual ~EMMAS3Hit();                      

    const EMMAS3Hit& operator=(const EMMAS3Hit& right);
    G4int operator==(const EMMAS3Hit& right) const;     

    inline void* operator new(size_t);        
    inline void  operator delete(void* hit);   

    virtual void Draw() {}                     
    virtual void Print();                      

    void AddEdep(G4double edep);               
    void SetKinematicsIfUnset(                
      G4double kineticEnergy,
      G4double theta,
      G4double phi);

    G4double GetEdep() const;                  
    G4double GetKineticEnergy() const;         
    G4double GetTheta() const;                 
    G4double GetPhi() const;                   
    G4bool HasKinematics() const;              

  private:
    G4double fEdep;                           
    G4double fKineticEnergy;                 
    G4double fTheta;                           
    G4double fPhi;                             
    G4bool fHasKinematics;                   
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
    fKineticEnergy = kineticEnergy; 
    fTheta = theta;                 
    fPhi = phi;                     
    fHasKinematics = true;          
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
