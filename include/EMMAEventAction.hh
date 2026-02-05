//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: EMMAEventAction.hh,v 1.6 2006-06-29 16:31:04 gunter Exp $
// --------------------------------------------------------------
//
#ifndef EMMAEventAction_h
#define EMMAEventAction_h 1


#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#ifdef G4ANALYSIS_USE

#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TFolder.h>
#include <TTree.h>
#include <TROOT.h>
#include <TAxis.h>

#endif // G4ANALYSIS_USE

class EMMAEventActionMessenger;

class EMMAEventAction : public G4UserEventAction
{
  public:
    EMMAEventAction();
    virtual ~EMMAEventAction();

  public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

  private:
    G4int DHC2ID;
    G4int s3HCID;
    G4ThreeVector localPos;
    G4double theta;
    G4double fp_pos[2],fp_theta;
    static const G4int kMaxS3Rings = 64;
    G4double fS3EnergySigmaMeV;
    G4long fS3RowsWritten;
    G4long fS3EventsWithHits;
    G4long fS3MultiHitEvents;
    G4long fS3ExtraRowsFromMultiHits;

    EMMAEventActionMessenger* messenger;
    G4int verboseLevel;

#ifdef G4ANALYSIS_USE
		TFile* rootfile;
    TTree* fp_tree;
		TH2F* fp_hitpos;
		TH1F* fp_hitangle;
    G4int s3_nhit;
    G4int s3_ring[kMaxS3Rings];
    G4double s3_energy_true[kMaxS3Rings];
    G4double s3_energy_smeared[kMaxS3Rings];
    G4double s3_edep[kMaxS3Rings];
    G4double s3_theta[kMaxS3Rings];
    G4double s3_phi[kMaxS3Rings];
    G4double s3_total_edep;
#endif // G4ANALYSIS_USE

  private:
    void ResetS3EventBuffers();

  public:
    inline void SetVerbose(G4int val) { verboseLevel = val; }
    inline G4int GetVerbose() const { return verboseLevel; }
};

#endif
