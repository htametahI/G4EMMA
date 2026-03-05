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

#include "EMMAS3Hit.hh"      

#include "G4UnitsTable.hh" 
#include "G4ios.hh"         

#include <iomanip>          

G4Allocator<EMMAS3Hit> EMMAS3HitAllocator; 

EMMAS3Hit::EMMAS3Hit()
  : G4VHit(),              
    fEdep(0.0),             
    fKineticEnergy(0.0),    
    fTheta(0.0),          
    fPhi(0.0),              
    fHasKinematics(false)  
{
}

EMMAS3Hit::~EMMAS3Hit()
{
}

EMMAS3Hit::EMMAS3Hit(const EMMAS3Hit& right)
  : G4VHit()                  
{
  fEdep = right.fEdep;       
  fKineticEnergy = right.fKineticEnergy; 
  fTheta = right.fTheta;      
  fPhi = right.fPhi;         
  fHasKinematics = right.fHasKinematics; 
}

const EMMAS3Hit& EMMAS3Hit::operator=(const EMMAS3Hit& right)
{
  if (this != &right) {      
    fEdep = right.fEdep;      
    fKineticEnergy = right.fKineticEnergy; 
    fTheta = right.fTheta;    
    fPhi = right.fPhi;       
    fHasKinematics = right.fHasKinematics; 
  }
  return *this;               
}

G4int EMMAS3Hit::operator==(const EMMAS3Hit& right) const
{
  return (this == &right) ? 1 : 0; // Compare hit pointers 
}

void EMMAS3Hit::Print()
{
  G4cout
    << "S3 hit: Edep=" << std::setw(7) << G4BestUnit(fEdep, "Energy") 
    << " Ekin=" << std::setw(7) << G4BestUnit(fKineticEnergy, "Energy") 
    << " theta(rad)=" << fTheta   
    << " phi(rad)=" << fPhi      
    << G4endl;
}
