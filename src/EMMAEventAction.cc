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
// $Id: EMMAEventAction.cc,v 1.10 2007-05-17 09:55:14 duns Exp $
// --------------------------------------------------------------
//

#include "EMMAEventAction.hh"          // Event action declaration.
#include "EMMAEventActionMessenger.hh" // Event-action messenger declaration.
#ifdef G4ANALYSIS_USE
#include "EMMAAnalysisManager.hh"      // ROOT output manager.
#endif // G4ANALYSIS_USE

#include "EMMADriftChamberHit.hh"      // Focal-plane hit definition.
#include "EMMAS3Hit.hh"                // Upstream S3 hit definition.

#include "G4Event.hh"                  // Event type.
#include "G4EventManager.hh"           // Event manager interface.
#include "G4HCofThisEvent.hh"          // Event hit container.
#include "G4PrimaryParticle.hh"        // Primary particle interface.
#include "G4RunManager.hh"             // AbortRun interface for accepted-row target stopping.
#include "G4SDManager.hh"              // Sensitive detector manager.
#include "G4SystemOfUnits.hh"          // CLHEP unit constants.
#include "G4Trajectory.hh"             // Retained from original file for compatibility.
#include "G4TrajectoryContainer.hh"    // Retained from original file for compatibility.
#include "G4UImanager.hh"              // Retained from original file for compatibility.
#include "G4VHitsCollection.hh"        // Retained from original file for compatibility.
#include "G4VVisManager.hh"            // Retained from original file for compatibility.
#include "G4ios.hh"                    // Geant4 stream definitions.
#include "Randomize.hh"                // G4RandGauss smearing helper.

#include <algorithm>                   // std::max, std::min.
#include <fstream>                     // std::ofstream for plain-text S3 output.
#include <cmath>                       // std::sqrt, std::acos, std::atan2.

using namespace std;                     // Keep original namespace style.

namespace {
// ------------------------------------------------------------------
// Hardcoded output gate on triton angle right before S3 entry.
// Set these values directly in code (no UI macro control).
// ------------------------------------------------------------------
const G4bool kApplyS3PreHitThetaOutputGate = true; // Enable/disable plain-text row filter on pre-S3 theta.
const G4bool kWriteOnlySingleHitS3Events = true;   // If true, plain-text S3 output keeps only events with exactly one S3 ring hit.
const G4double kS3PreHitThetaMinDeg = 131.531771;       // Lower theta bound (deg) for writing S3 rows.
const G4double kS3PreHitThetaMaxDeg = 132.357455;       // Upper theta bound (deg) for writing S3 rows.
// const G4double kS3PreHitThetaMinDeg = 130.0;       // Lower theta bound (deg) for writing S3 rows.
// const G4double kS3PreHitThetaMaxDeg = 160.0;       // Upper theta bound (deg) for writing S3 rows.
// const G4double kManualGateMinDeg = 158.83874018; // Manual fallback min theta if ring gate is off.
// const G4double kManualGateMaxDeg = 160.46334506;

bool IsS3PreHitThetaAccepted(G4double thetaDeg)
{
  if (!std::isfinite(thetaDeg)) {
    return false;
  }
  if (!kApplyS3PreHitThetaOutputGate) {
    return true;
  }
  const G4double thetaMin = std::min(kS3PreHitThetaMinDeg, kS3PreHitThetaMaxDeg);
  const G4double thetaMax = std::max(kS3PreHitThetaMinDeg, kS3PreHitThetaMaxDeg);
  return thetaDeg >= thetaMin && thetaDeg <= thetaMax;
}
} // namespace

//global variable
extern G4String UserDir;                // Keep original external declaration.
extern G4double beamEnergyAtTarget;     // Beam energy at reaction depth from primary generator.
extern G4int tritonLabGateRingNumber;   // 0-based S3 ring copy number selected for rejection sampling.
extern G4String s3TritonOutputFileName; // Plain-text output file path set by reaction initialization.
extern G4double ejectileEnergy;         // Triton (product #4) kinetic energy at reaction generation in LAB.
extern G4double ejectileDirX;           // Triton generation-direction x-component in LAB.
extern G4double ejectileDirY;           // Triton generation-direction y-component in LAB.
extern G4double ejectileDirZ;           // Triton generation-direction z-component in LAB.
extern G4double tritonExitTargetEnergy; // True triton kinetic energy immediately after target exit.
extern G4double tritonExitTargetTheta;  // Triton theta immediately after target exit.
extern G4double tritonExitTargetPhi;    // Triton phi immediately after target exit.
extern G4bool enforceS3AcceptedRowTarget; // Enable/disable stopping run once accepted output rows reach target.
extern G4long s3AcceptedRowTarget;        // Requested accepted output-row target.
extern G4long s3AcceptedRowCount;         // Running accepted output-row count.
extern G4long s3GlobalEventSerial;        // Monotonic event serial written into output eventID column.

EMMAEventAction::EMMAEventAction()
  : DHC2ID(-1),                          // Initialize focal-plane HC id.
    s3HCID(-1),                          // Initialize S3 HC id as unresolved.
    localPos(),                          // Initialize focal-plane position buffer.
    theta(0.0),                          // Initialize focal-plane theta buffer.
    fp_pos{0.0, 0.0},                    // Initialize focal-plane x/y output.
    fp_theta(0.0),                       // Initialize focal-plane theta output.
    fS3EnergySigmaMeV(0.07),              // Hardcoded Gaussian sigma = 100 keV = 0.1 MeV.
    fS3RowsWritten(0),                   // Count rows written to the plain-text S3 file.
    fS3EventsWithHits(0),                // Count events with at least one S3 ring hit.
    fS3MultiHitEvents(0),                // Count events with more than one S3 ring hit.
    fS3ExtraRowsFromMultiHits(0),        // Count extra rows produced by multi-hit events.
    messenger(0),                        // Initialize messenger pointer.
    verboseLevel(1)                      // Keep original default verbose behavior.
#ifdef G4ANALYSIS_USE
  , rootfile(0)                          // Initialize ROOT file pointer.
  , fp_tree(0)                           // Initialize ROOT tree pointer.
  , fp_hitpos(0)                         // Initialize focal-plane position histogram pointer.
  , fp_hitangle(0)                       // Initialize focal-plane angle histogram pointer.
  , s3_nhit(0)                           // Initialize S3 hit count output.
  , s3_ring{}                            // Zero-initialize S3 ring-id array.
  , s3_energy_true{}                     // Zero-initialize S3 true-energy array.
  , s3_energy_smeared{}                  // Zero-initialize S3 smeared-energy array.
  , s3_edep{}                            // Zero-initialize S3 deposited-energy array.
  , s3_theta{}                           // Zero-initialize S3 theta array.
  , s3_phi{}                             // Zero-initialize S3 phi array.
  , s3_total_edep(0.0)                   // Initialize S3 total Edep output.
#endif // G4ANALYSIS_USE
{
  // Resolve the focal-plane drift chamber collection id once.
  G4String colName;                                              // Temporary collection-name buffer.
  G4SDManager* SDman = G4SDManager::GetSDMpointer();             // Access global SD manager.
  DHC2ID = SDman->GetCollectionID(colName = "detector/detectorColl"); // Cache focal-plane HC id.

  // Build the event-action messenger used by /mydet/verbose.
  messenger = new EMMAEventActionMessenger(this);                // Allocate messenger instance.

  // Reset event output buffers to clean defaults.
  ResetS3EventBuffers();                                         // Initialize S3 per-event arrays.
  fp_pos[0] = -9999.0;                                           // Initialize fp x sentinel.
  fp_pos[1] = -9999.0;                                           // Initialize fp y sentinel.
  fp_theta = -9999.0;                                            // Initialize fp theta sentinel.

#ifdef G4ANALYSIS_USE
  // Get analysis manager singleton and root file handle.
  EMMAAnalysisManager* analysisManager = EMMAAnalysisManager::getInstance(); // Analysis manager singleton.
  rootfile = analysisManager->getRootfile();                                 // Output ROOT file pointer.

  // Create original focal-plane histograms.
  fp_hitpos = new TH2F("hitpos", "Focal plane hit position", 6000, -30, 30, 6000, -30, 30); // FP x-y histogram.
  fp_hitpos->GetXaxis()->SetTitle("X position (mm)");                                          // FP x-axis title.
  fp_hitpos->GetYaxis()->SetTitle("Y position (mm)");                                          // FP y-axis title.
  fp_hitangle = new TH1F("hitangle", "Focal plane hit angle", 100, 0, 10);                  // FP theta histogram.
  fp_hitangle->GetXaxis()->SetTitle("Theta (deg)");                                            // FP theta x-axis title.
  fp_hitangle->GetYaxis()->SetTitle("Counts");                                                  // FP theta y-axis title.

  // Get event tree and create focal-plane + S3 branches.
  fp_tree = analysisManager->getfpRoottree();                                                    // Shared event tree pointer.
  fp_tree->Branch("fp_pos", &fp_pos, "fp_pos[2]/D");                                           // FP x,y branch.
  fp_tree->Branch("fp_theta", &fp_theta, "fp_theta/D");                                        // FP theta branch.
  fp_tree->Branch("s3_nhit", &s3_nhit, "s3_nhit/I");                                           // Number of S3 ring hits.
  fp_tree->Branch("s3_ring", s3_ring, "s3_ring[s3_nhit]/I");                                   // 1-based ring id per S3 hit.
  fp_tree->Branch("s3_energy_true", s3_energy_true, "s3_energy_true[s3_nhit]/D");             // True triton kinetic energy (MeV).
  fp_tree->Branch("s3_energy_smeared", s3_energy_smeared, "s3_energy_smeared[s3_nhit]/D");    // Smeared triton kinetic energy (MeV).
  fp_tree->Branch("s3_edep", s3_edep, "s3_edep[s3_nhit]/D");                                   // Deposited energy per ring (MeV).
  fp_tree->Branch("s3_theta", s3_theta, "s3_theta[s3_nhit]/D");                                // Theta per ring hit (deg).
  fp_tree->Branch("s3_phi", s3_phi, "s3_phi[s3_nhit]/D");                                      // Phi per ring hit (deg).
  fp_tree->Branch("s3_total_edep", &s3_total_edep, "s3_total_edep/D");                         // Event total S3 deposited energy (MeV).
  fp_tree->Branch("s3_sigma_mev", &fS3EnergySigmaMeV, "s3_sigma_mev/D");                        // Hardcoded Gaussian resolution sigma (MeV).
#endif // G4ANALYSIS_USE
}

EMMAEventAction::~EMMAEventAction()
{
#ifdef G4ANALYSIS_USE
  // Persist and close ROOT output.
  EMMAAnalysisManager::dispose();      // Dispose analysis manager singleton.
  rootfile->Write();                   // Write histograms and trees to file.
  rootfile->Close();                   // Close ROOT file handle.
  delete rootfile;                     // Release ROOT file pointer.
#endif // G4ANALYSIS_USE
  if (fS3EventsWithHits > 0 || fS3RowsWritten > 0) { // Print S3 multi-hit diagnostics once at run end.
    G4cout << "S3 output summary: rows=" << fS3RowsWritten
           << ", eventsWithHits=" << fS3EventsWithHits
           << ", multiHitEvents=" << fS3MultiHitEvents
           << ", extraRowsFromMultiHits=" << fS3ExtraRowsFromMultiHits
           << G4endl;
  }
  delete messenger;                    // Release messenger pointer.
}

void EMMAEventAction::BeginOfEventAction(const G4Event* evt)
{
  ++s3GlobalEventSerial;               // Advance monotonic serial used for output eventID tracking.

  // Reset S3 outputs at the start of each event.
  ResetS3EventBuffers();               // Clear S3 output arrays.

  // Reset focal-plane outputs at the start of each event.
  fp_pos[0] = -9999.0;                 // FP x sentinel for no-hit event.
  fp_pos[1] = -9999.0;                 // FP y sentinel for no-hit event.
  fp_theta = -9999.0;                  // FP theta sentinel for no-hit event.

  // Keep original event-start printout.
  G4cout << "Start event: " << evt->GetEventID() << G4endl; // Event header print.
}

void EMMAEventAction::EndOfEventAction(const G4Event* evt)
{
  // Get event-level hit container.
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent(); // Geant4 event hit container.

  // Resolve focal-plane drift chamber hits pointer.
  EMMADriftChamberHitsCollection* DHC2 = 0;       // FP collection pointer.
  if (HCE) {                                      // Only query if event has hit container.
    DHC2 = (EMMADriftChamberHitsCollection*)(HCE->GetHC(DHC2ID)); // Retrieve FP collection.
  }

  // Resolve S3 collection id lazily in case SD registration order changes.
  if (s3HCID < 0) {                                                        // Resolve once on first event.
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();                 // Access SD manager.
    s3HCID = sdManager->GetCollectionID("S3Detector/S3HitsCollection");    // Try fully-qualified collection name first.
    if (s3HCID < 0) {                                                      // Fallback for unqualified name setups.
      s3HCID = sdManager->GetCollectionID("S3HitsCollection");             // Try bare collection name.
    }
  }

  // Resolve S3 hits pointer.
  EMMAS3HitsCollection* s3Hits = 0;                                       // S3 collection pointer.
  if (HCE && s3HCID >= 0) {                                               // Only query if container and id are valid.
    s3Hits = (EMMAS3HitsCollection*)(HCE->GetHC(s3HCID));                 // Retrieve S3 collection.
  }

  // Collect focal-plane data.
  G4bool hasFocalPlaneData = false;                                       // Track whether this event has FP output data.
  if (DHC2) {                                                             // Guard against missing FP collection.
    const G4int n_hit = DHC2->entries();                                  // Number of FP hits.
    for (G4int i2 = 0; i2 < 5 && !hasFocalPlaneData; ++i2) {              // Keep original layer scan priority.
      for (G4int i1 = 0; i1 < n_hit; ++i1) {                              // Scan available FP hits.
        EMMADriftChamberHit* aHit = (*DHC2)[i1];                          // Current FP hit pointer.
        if (!aHit) continue;                                              // Skip null hit slots.
        if (aHit->GetLayerID() != i2) continue;                           // Keep original layer filter logic.

        hasFocalPlaneData = true;                                         // Mark FP payload as available.
#ifdef G4ANALYSIS_USE
        localPos = aHit->GetLocalPos();                                   // Local focal-plane position.
        theta = aHit->GetTheta()/CLHEP::deg;                              // Focal-plane theta in degrees.
        if (fp_hitpos) fp_hitpos->Fill(localPos.x()/CLHEP::mm, localPos.y()/CLHEP::mm); // Fill FP x-y histogram.
        if (fp_hitangle) fp_hitangle->Fill(theta);                        // Fill FP theta histogram.
        fp_pos[0] = localPos.x()/CLHEP::mm;                               // Store FP x output.
        fp_pos[1] = localPos.y()/CLHEP::mm;                               // Store FP y output.
        fp_theta = theta;                                                 // Store FP theta output.
#endif // G4ANALYSIS_USE
        break;                                                            // Keep first matching hit only.
      }
    }
  }

  // Collect S3 per-ring triton observables.
  if (s3Hits) {                                                           // Guard against missing S3 collection.
    const G4int nEntries = s3Hits->entries();                             // Number of allocated S3 hit entries.
    const G4int nRingEntries = (nEntries > 0) ? (nEntries - 1) : 0;      // Last entry is total accumulator slot.

    for (G4int ringIndex = 0; ringIndex < nRingEntries; ++ringIndex) {    // Loop over physical ring slots.
      EMMAS3Hit* ringHit = (*s3Hits)[ringIndex];                          // Current ring hit pointer.
      if (!ringHit) continue;                                             // Skip null hit slots.
      if (!ringHit->HasKinematics()) continue;                            // Require stored first-hit kinematics.
      if (s3_nhit >= kMaxS3Rings) continue;                               // Protect fixed-size output arrays.

      const G4double trueEnergyMeV = ringHit->GetKineticEnergy() / CLHEP::MeV; // True triton kinetic energy in MeV.
      const G4double sigmaMeV = fS3EnergySigmaMeV;                               // Fixed Gaussian sigma in MeV.
      const G4double smearedRawMeV = G4RandGauss::shoot(trueEnergyMeV, sigmaMeV); // Gaussian-smeared energy.
      const G4double smearedEnergyMeV = (smearedRawMeV > 0.0) ? smearedRawMeV : 0.0; // Clamp negative smears.

      s3_ring[s3_nhit] = ringIndex + 1;                                   // Store 1-based ring id.
      s3_energy_true[s3_nhit] = trueEnergyMeV;                            // Store true kinetic energy.
      s3_energy_smeared[s3_nhit] = smearedEnergyMeV;                      // Store smeared kinetic energy.
      s3_edep[s3_nhit] = ringHit->GetEdep() / CLHEP::MeV;                 // Store deposited energy.
      s3_theta[s3_nhit] = ringHit->GetTheta() / CLHEP::deg;               // Store theta in degrees.
      s3_phi[s3_nhit] = ringHit->GetPhi() / CLHEP::deg;                   // Store phi in degrees.
      ++s3_nhit;                                                           // Increment S3 ring-hit count.
    }

    if (nEntries > 0) {                                                    // Guard total slot access.
      EMMAS3Hit* totalHit = (*s3Hits)[nEntries - 1];                      // Resolve total accumulator slot.
      if (totalHit) {                                                      // Guard null total hit pointer.
        s3_total_edep = totalHit->GetEdep() / CLHEP::MeV;                 // Store event total deposited energy.
      }
    }
  }

  const G4bool passesS3MultiplicityWriteFilter =                          // Optional event-level multiplicity filter for plain-text S3 rows.
    (!kWriteOnlySingleHitS3Events) || (s3_nhit == 1);                    // Keep all events, or only exactly-one-ring-hit events.

  G4int s3OutputRowsThisEvent = 0;                                        // Number of rows that pass output filters.
  if (passesS3MultiplicityWriteFilter) {                                  // Apply optional single-hit requirement before row-level theta gate.
    for (G4int i = 0; i < s3_nhit; ++i) {                                 // Scan stored S3 hits to count output-eligible rows.
      if (IsS3PreHitThetaAccepted(s3_theta[i])) {                         // Apply output gate using pre-S3 theta observable.
        ++s3OutputRowsThisEvent;                                          // Count rows that pass angular output filtering.
      }
    }
  }

  G4int s3RowsToWriteThisEvent = s3OutputRowsThisEvent;                   // Rows actually written after strict target-cap truncation.
  if (enforceS3AcceptedRowTarget && s3AcceptedRowTarget > 0) {            // Apply option-2 accepted-row target cap.
    const G4long remainingRows = s3AcceptedRowTarget - s3AcceptedRowCount; // Accepted rows still needed to hit target.
    if (remainingRows <= 0) {                                             // Target already satisfied before this event.
      s3RowsToWriteThisEvent = 0;                                         // Suppress further output rows.
    } else if (remainingRows < s3RowsToWriteThisEvent) {                  // This event would overshoot target.
      s3RowsToWriteThisEvent = static_cast<G4int>(remainingRows);         // Truncate output rows to land exactly on target.
    }
  }

  if (s3RowsToWriteThisEvent > 0) {                                       // Update summary counters for actual output rows.
    ++fS3EventsWithHits;                                                  // Count events that produced at least one output row.
    fS3RowsWritten += s3RowsToWriteThisEvent;                             // Add rows that passed output filtering and target cap.
    if (s3RowsToWriteThisEvent > 1) {                                     // Multi-hit in output file: more than one written row per event.
      ++fS3MultiHitEvents;                                                // Count multi-hit event.
      fS3ExtraRowsFromMultiHits += (s3RowsToWriteThisEvent - 1);         // Rows above one-per-event baseline.
    }
  }
  s3AcceptedRowCount += s3RowsToWriteThisEvent;                           // Update global accepted-row count used by option-2 stop logic.

#ifdef G4ANALYSIS_USE
  // Append plain-text file rows with requested triton observables.
  if (!s3TritonOutputFileName.empty() && s3RowsToWriteThisEvent > 0) {    // Only write when run mode enabled it and output-gate accepted rows.
    std::ofstream s3Outfile(s3TritonOutputFileName, std::ios::app);        // Append one line per accepted S3 ring hit.
    if (s3Outfile.is_open()) {                                             // Guard against file-open failures.
      s3Outfile.precision(17);                                             // Keep full precision for analysis.
      const G4double tritonGenDirNorm = std::sqrt(                         // Triton generation-direction magnitude.
        ejectileDirX * ejectileDirX +
        ejectileDirY * ejectileDirY +
        ejectileDirZ * ejectileDirZ);
      G4double tritonGenThetaDeg = -9999.0;                                // Triton-generation theta sentinel if direction invalid.
      G4double tritonGenPhiDeg = -9999.0;                                  // Triton-generation phi sentinel if direction invalid.
      if (tritonGenDirNorm > 0.0) {                                        // Convert valid generated triton direction to lab angles.
        const G4double tritonGenCosTheta = std::max(-1.0, std::min(1.0, ejectileDirZ / tritonGenDirNorm)); // Clamp acos argument.
        tritonGenThetaDeg = std::acos(tritonGenCosTheta) / CLHEP::deg;     // Generated-triton polar angle in degrees.
        tritonGenPhiDeg = std::atan2(ejectileDirY, ejectileDirX) / CLHEP::deg; // Generated-triton azimuthal angle in degrees.
      }
      G4int rowsWrittenThisEvent = 0;                                      // Count rows emitted in this event after target-cap truncation.
      for (G4int i = 0; i < s3_nhit; ++i) {                                // Loop over ring hits recorded this event.
        if (rowsWrittenThisEvent >= s3RowsToWriteThisEvent) break;         // Stop once strict target cap for this event is reached.
        if (!IsS3PreHitThetaAccepted(s3_theta[i])) continue;               // Keep only requested pre-S3 theta interval in output file.
        const G4double tritonTargetEnergyLossMeV = (ejectileEnergy - tritonExitTargetEnergy) / CLHEP::MeV; // Energy loss in target.
        s3Outfile << beamEnergyAtTarget / CLHEP::MeV << ","                // Beam energy at target (MeV).
                  << tritonExitTargetEnergy / CLHEP::MeV << ","            // Triton kinetic energy immediately at target exit (MeV).
                  << tritonExitTargetTheta / CLHEP::deg << ","             // Triton theta immediately at target exit (deg).
                  << tritonExitTargetPhi / CLHEP::deg << ","               // Triton phi immediately at target exit (deg).
                  << s3_edep[i] << ","                                     // S3 deposited energy (MeV).
                  << s3_theta[i] << ","                                    // Triton theta just before entering S3 (deg).
                  << s3_phi[i] << ","                                      // Triton phi just before entering S3 (deg).
                  << ejectileEnergy / CLHEP::MeV << ","                    // Triton kinetic energy at reaction generation (MeV).
                  << tritonGenThetaDeg << ","                              // Triton theta at reaction generation (deg).
                  << tritonGenPhiDeg << ","                                 // Triton phi at reaction generation (deg).
                  << s3_energy_smeared[i] << ","                           // S3 energy with Gaussian resolution (MeV).
                  << tritonTargetEnergyLossMeV                              // Triton energy loss in target (MeV).
                  << G4endl;
        ++rowsWrittenThisEvent;                                            // Track emitted rows so we can enforce event-level truncation.
      }
    }
  }
#endif // G4ANALYSIS_USE

  if (enforceS3AcceptedRowTarget && s3AcceptedRowTarget > 0 &&             // Stop generation once exact accepted-row target has been reached.
      s3AcceptedRowCount >= s3AcceptedRowTarget) {
    G4RunManager::GetRunManager()->AbortRun(true);                         // Abort current BeamOn loop after this event.
  }

  // Verbose printout for selected events only.
  const G4bool doVerbosePrint = !(verboseLevel == 0 || evt->GetEventID() % verboseLevel != 0); // Verbose selector.
  if (doVerbosePrint) {                                                   // Apply original verbose gating.
    if (evt->GetPrimaryVertex(0)) {                                       // Guard primary vertex access.
      G4PrimaryParticle* primary = evt->GetPrimaryVertex(0)->GetPrimary(0); // Primary particle pointer.
      if (primary) {                                                       // Guard primary particle pointer.
        G4cout << G4endl                                                   // Start verbose block.
               << "\t>>> Event " << evt->GetEventID() << " >>> Simulation truth : "
               << primary->GetG4code()->GetParticleName()
               << " " << primary->GetMomentum() << G4endl;
      }
    }
    if (DHC2) {                                                            // Print FP hit count if collection exists.
      G4cout << "Stopping Block has " << DHC2->entries() << " hits." << G4endl;
    }
    if (s3Hits) {                                                          // Print S3 hit summary if collection exists.
      G4cout << "S3 has " << s3_nhit << " triton ring hits." << G4endl;
    }
  }

#ifdef G4ANALYSIS_USE
  // Write one tree entry for events with either focal-plane or S3 data.
  if (fp_tree && (hasFocalPlaneData || s3_nhit > 0)) {                    // Keep tree compact and relevant.
    fp_tree->Fill();                                                       // Commit event payload.
  }
#endif // G4ANALYSIS_USE
}

void EMMAEventAction::ResetS3EventBuffers()
{
#ifdef G4ANALYSIS_USE
  s3_nhit = 0;                                                             // Reset S3 hit count.
  s3_total_edep = 0.0;                                                     // Reset event total deposited energy.
  for (G4int i = 0; i < kMaxS3Rings; ++i) {                                // Loop over fixed output slots.
    s3_ring[i] = -1;                                                       // Reset ring id sentinel.
    s3_energy_true[i] = 0.0;                                               // Reset true energy slot.
    s3_energy_smeared[i] = 0.0;                                            // Reset smeared energy slot.
    s3_edep[i] = 0.0;                                                      // Reset deposited energy slot.
    s3_theta[i] = -9999.0;                                                 // Reset theta sentinel.
    s3_phi[i] = -9999.0;                                                   // Reset phi sentinel.
  }
#endif // G4ANALYSIS_USE
}
