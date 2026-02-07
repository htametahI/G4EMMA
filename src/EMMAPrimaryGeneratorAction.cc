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
// $Id: EMMAPrimaryGeneratorAction.cc,v 1.5 2006-06-29 16:33:05 gunter Exp $
// --------------------------------------------------------------
//

#include "EMMAPrimaryGeneratorAction.hh"
#include "EMMAPrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"

#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4IonTable.hh"
#include "G4NucleiProperties.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4Box.hh"

#include <string> //words and sentences
#include <fstream> //Stream class to both read and write from/to files
#include <sstream>
#include <algorithm>
#include <cmath>
#include <limits>
using namespace std;

namespace {
bool IsFinite(G4double value)
{
  return std::isfinite(value);
}

bool IsFiniteVector(const G4ThreeVector& vec)
{
  return std::isfinite(vec.x()) && std::isfinite(vec.y()) && std::isfinite(vec.z());
}

// ------------------------------------------------------------------
// Hardcoded S3 ring-gate controls (set values here directly in code).
// ------------------------------------------------------------------
const G4bool kApplyTritonLabAngleGate = true;    // Enable reaction-point pre-sampling gate to reduce wasted transport.
const G4bool kUseS3RingGate = false;             // Ignored when reaction-point gate is disabled.
const G4int kSelectedS3Ring = 0;                 // 0-based ring copy number to keep (requested: ring 0).
const G4bool kGenerateTritonToS3 = true;         // If true, generate product #4 (triton) as tracked primary.
const G4int kS3RingCount = 24;                   // Must match S3 geometry ring count.
const G4double kS3InnerRadius = 11.0 * mm;       // Must match S3 geometry inner radius.
const G4double kS3OuterRadius = 35.0 * mm;       // Must match S3 geometry outer radius.
const G4double kS3DistanceFromTarget = 31. * mm; // Must match upstream S3 target distance.
const G4double kManualGateMinDeg = 130.5; // Manual fallback min theta if ring gate is off.
const G4double kManualGateMaxDeg = 133.0; // Manual fallback max theta if ring gate is off.
const G4int kS3AcceptedRowSamplingBudgetFactor = 200; // Max generated events = factor * requested accepted output rows.
const G4int kReactionBeamOnChunkEvents = 5000; // Reaction is run in chunks; reduces upfront beam pre-generation.
const G4int kBeamTopUpChunkEvents = 5000; // Extra beam samples are generated in small chunks on demand.
// const G4double kManualGateMinDeg = 130.0; // Manual fallback min theta if ring gate is off.
// const G4double kManualGateMaxDeg = 160.0; // Manual fallback max theta if ring gate is off.

bool ComputeS3RingLabAngleRangeDeg(
  G4int ringNumberZeroBased,
  G4int ringCount,
  G4double innerRadius,
  G4double outerRadius,
  G4double distanceFromTarget,
  G4double& thetaMinDeg,
  G4double& thetaMaxDeg)
{
  if (ringCount <= 0 || ringNumberZeroBased < 0 || ringNumberZeroBased >= ringCount) {
    return false; // Ring index is outside valid 0..ringCount-1.
  }
  if (distanceFromTarget <= 0.0 || outerRadius <= innerRadius) {
    return false; // Geometry values are not physically valid.
  }

  const G4double ringWidth = (outerRadius - innerRadius) / ringCount; // Uniform radial pitch.
  const G4double ringInnerRadius = innerRadius + ringNumberZeroBased * ringWidth; // Selected-ring inner edge.
  const G4double ringOuterRadius = ringInnerRadius + ringWidth; // Selected-ring outer edge.

  // Upstream detector is at negative z, so theta = 180 deg - atan(r/|z|).
  const G4double thetaAtInnerEdgeDeg = 180.0 - std::atan2(ringInnerRadius, distanceFromTarget) / deg; // Inner-edge theta.
  const G4double thetaAtOuterEdgeDeg = 180.0 - std::atan2(ringOuterRadius, distanceFromTarget) / deg; // Outer-edge theta.

  thetaMinDeg = std::min(thetaAtInnerEdgeDeg, thetaAtOuterEdgeDeg); // Lower angular boundary.
  thetaMaxDeg = std::max(thetaAtInnerEdgeDeg, thetaAtOuterEdgeDeg); // Upper angular boundary.
  return true; // Valid range computed successfully.
}
}

bool IsTritonLabAngleAccepted(G4double tritonDirZ, G4double minDeg, G4double maxDeg)
{
  if (!IsFinite(tritonDirZ)) {
    return false;
  }
  tritonDirZ = std::max(-1.0, std::min(1.0, tritonDirZ));
  G4double thetaLabDeg = std::acos(tritonDirZ) * 180.0 / CLHEP::pi;
  if (!IsFinite(thetaLabDeg)) {
    return false;
  }
  return thetaLabDeg >= minDeg && thetaLabDeg <= maxDeg;
}


// global variables 
G4bool prepareBeam = true;
G4String inTargetFileName;
G4String postTargetFileName;
G4String focalPlaneFileName;
G4double userCharge = 54.; // default value
G4String postDegrader1FileName;
G4double depth;
G4double beamEnergyAtTarget = 0.;
// Ejectiles 
G4String postTargetEjectileFileName;
G4String postTargetRecoilSpectrometerFileName;
G4double ejectileEnergy = 0.;
G4double ejectileDirX = 0.;
G4double ejectileDirY = 0.;
G4double ejectileDirZ = 0.;
G4double tritonExitTargetEnergy = -1.0 * MeV;
G4double tritonExitTargetTheta = -9999.0 * deg;
G4double tritonExitTargetPhi = -9999.0 * deg;
G4double recoilLabEnergy = 0.;
G4double recoilLabDirX = 0.;
G4double recoilLabDirY = 0.;
G4double recoilLabDirZ = 0.;
G4int ejectileZ = 0;
G4int ejectileA = 0;
G4double recoilThetaCM = 0.;
G4double ejectileThetaCM = 0.;
G4String s3TritonOutputFileName;

// Triton gate controls used in simulateTwoBodyReaction rejection sampling.
G4bool applyTritonLabAngleGate = kApplyTritonLabAngleGate; // Enable/disable rejection gate here in code.
G4int tritonLabGateRingNumber = kSelectedS3Ring; // 0-based ring selector used to derive angular gate.
G4double tritonLabAngleMinDeg = kManualGateMinDeg; // Active lower bound (deg), set below.
G4double tritonLabAngleMaxDeg = kManualGateMaxDeg; // Active upper bound (deg), set below.

G4long gateTrialTotal = 0;
G4long gateAcceptedEvents = 0;
G4bool enforceS3AcceptedRowTarget = false; // If true, abort reaction run once accepted S3 output rows hit target.
G4long s3AcceptedRowTarget = 0;            // Target accepted S3 output rows (set from nEvents in reaction mode).
G4long s3AcceptedRowCount = 0;             // Running count of accepted S3 output rows written to plain-text file.
G4long s3GlobalEventSerial = -1;           // Monotonic event serial used as output eventID across possible run restarts.
G4bool s3AcceptedRowAbortActive = false;   // True only during /mydet/doReaction BeamOn; keeps manual runs from aborting.
G4int beamInputRowsLoaded = 0;             // Number of prepared beam rows currently loaded in memory.
G4long reactionGeneratedEventCount = 0;    // Monotonic generated-event index across chunked reaction BeamOn calls.

EMMAPrimaryGeneratorAction::EMMAPrimaryGeneratorAction()  // constructor
{
	sigmaEnergy = 0.*MeV;
	transEmittance = 0.*mm*mrad;
	energyBeam = nullptr;
	posxBeam = nullptr;
	posyBeam = nullptr;
	poszBeam = nullptr;
	dirxBeam = nullptr;
	diryBeam = nullptr;
	dirzBeam = nullptr;

  // Configure triton rejection gate from chosen S3 ring geometry.
  if (applyTritonLabAngleGate && kUseS3RingGate) {
    const G4bool hasValidRange = ComputeS3RingLabAngleRangeDeg(
      tritonLabGateRingNumber,
      kS3RingCount,
      kS3InnerRadius,
      kS3OuterRadius,
      kS3DistanceFromTarget,
      tritonLabAngleMinDeg,
      tritonLabAngleMaxDeg);
    if (!hasValidRange) {
      std::ostringstream msg;
      msg << "Invalid S3 ring gate settings. ring=" << tritonLabGateRingNumber
          << " ringCount=" << kS3RingCount
          << " inner(mm)=" << (kS3InnerRadius / mm)
          << " outer(mm)=" << (kS3OuterRadius / mm)
          << " distance(mm)=" << (kS3DistanceFromTarget / mm);
      G4Exception("EMMAPrimaryGeneratorAction::EMMAPrimaryGeneratorAction", "EMMA0009",
                  FatalException, msg.str().c_str());
    }
    G4cout << "Configured triton gate for S3 ring(copyNo) " << tritonLabGateRingNumber
           << " with theta_lab in [" << tritonLabAngleMinDeg
           << ", " << tritonLabAngleMaxDeg << "] deg" << G4endl;
  } else if (applyTritonLabAngleGate) {
    tritonLabGateRingNumber = -1; // Mark that this run uses manual theta limits, not a ring-derived gate.
    G4cout << "Configured manual triton gate with theta_lab in ["
           << tritonLabAngleMinDeg << ", " << tritonLabAngleMaxDeg
           << "] deg" << G4endl;
  } else {
    tritonLabGateRingNumber = -1; // Mark that no ring-specific gate is active.
    G4cout << "Triton lab-angle gate disabled." << G4endl;
  }
  if (kGenerateTritonToS3) {
    G4cout << "Primary reaction product mode: tracking ejectile/triton (product #4)." << G4endl;
  } else {
    G4cout << "Primary reaction product mode: tracking recoil (product #3)." << G4endl;
  }

	G4int n_particle = 1;
	//G4ParticleGun class generates primary particle(s) with a given momentum and position
	particleGun  = new G4ParticleGun(n_particle);

	//create a messenger for this class
	gunMessenger = new EMMAPrimaryGeneratorMessenger(this);

	// default particle kinematics
	particleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.*m));
	particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));

	G4Step aStep;
	
	energy = aStep.GetTotalEnergyDeposit(); // =0
	Angle = 0;

	// default values
	nEvents = 1e6;
	fZ1=0.,fA1=0., fZ2=0.,fA2=0.;
	fZ3=0.,fA3=0., fZ4=0.,fA4=0.;
	fqmin = 180.*deg;
	fqmax = 0.*deg;
	fCharge3 = 0.;

	inTargetFileName = UserDir + "/ExcitationEnergy/beam.dat";  //Used in EMMASteppingAction

	// alpha-source input file
	useAlphaSource = false;
	G4String text, line;
	ifstream inputfil;
	G4String filename = UserDir + "/UserInput/alphaSource.dat";
	inputfil.open ( filename, ios::in );
	if ( inputfil.is_open() ) {
	  int n=0;
	  while ( inputfil.good() ) {
	    inputfil >> text;
	    if (text=="#") { // skip comments
	      getline (inputfil,line);
	    }
	    else {
	      n = n+1;
	      if (n==1) {if (text=="YES") useAlphaSource = true;}
	      if (n==2) {energyAlphaSource = atof(text.c_str());} 
	      if (n==3) {maxAngleAlphaSource = atof(text.c_str());} 
	    }
	  }
	  inputfil.close();
	}
	else G4cout << "Unable to open " << filename << G4endl; 

}

EMMAPrimaryGeneratorAction::~EMMAPrimaryGeneratorAction()
{
  delete[] energyBeam;
  delete[] posxBeam;
  delete[] posyBeam;
  delete[] poszBeam;
  delete[] dirxBeam;
  delete[] diryBeam;
  delete[] dirzBeam;
  delete particleGun;	//must delete G4ParticleGun
  delete gunMessenger;
}

void EMMAPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ // This method is invoked at the beginning of each event

  G4double randNumb;
  if (prepareBeam) randNumb = G4UniformRand(); //used to chose a random reaction depth
  else randNumb = 1.0;
  
  
  // this is to update target thickness
  // <><><><><><><><><><><><><><><><><><><><> //
  depth = targetThickness * randNumb;
  G4VSolid* targetSolid = new G4Box("target",5.*cm,5.*cm,depth/2.0);
  G4LogicalVolumeStore* logVolStore = G4LogicalVolumeStore::GetInstance();	
  G4LogicalVolume* target = logVolStore->GetVolume("targetLogical",true); //"targetLogical" declared in DetectorConstruction
  target->SetSolid(targetSolid);
  G4RunManager::GetRunManager()->ReOptimize( target );
  // <><><><><><><><><><><><><><><><><><><><> //
		
  G4double Ekin;
  G4ParticleDefinition* particleDef;
  tritonExitTargetEnergy = -1.0 * MeV; // Reset per-event; stepping action fills it when triton exits the target.
  tritonExitTargetTheta = -9999.0 * deg; // Reset per-event triton target-exit theta.
  tritonExitTargetPhi = -9999.0 * deg;   // Reset per-event triton target-exit phi.


  // to simulate just an isotropic alpha source
  if (useAlphaSource) {	//boolean determined from alphaSource.dat input file
    simulateReaction = false;
    
    // Ion (values read in from alphaSource.dat in constructor)
    /* 
       NB: The simulation terminates with a bus error if we pick an alpha particle, so instead
       we "cheat" and use a 6Li(3+) ion instead, while scaling the kinetic energy by a factor of 1.5
    */
    particleDef = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(3,6,0.0); //FindParticle("alpha");
    particleGun->SetParticleDefinition(particleDef);
    // charge
    userCharge = 3; // userCharge is a global variable! (used e.g. in BGFields1-7.cc)
    particleGun->SetParticleCharge(userCharge);    
    // Energy
    Ekin = energyAlphaSource * 1.5;
    particleGun->SetParticleEnergy(Ekin *MeV);
    // Position
    G4double xBeam=0.*m, yBeam=0.*m;
    particleGun->SetParticlePosition(G4ThreeVector(xBeam,yBeam,0.0*m));
    // Sample angle
    G4double x, y, z;
    G4double costh = 1.0 - G4UniformRand() * (1-std::cos(maxAngleAlphaSource*CLHEP::pi/180));
    G4double theta = std::acos(costh);
    G4double phi = G4UniformRand()*CLHEP::twopi;
    x = std::sin(theta) * std::cos(phi);
    y = std::sin(theta) * std::sin(phi);
    z = std::cos(theta);
    particleGun->SetParticleMomentumDirection(G4ThreeVector(x,y,z));

  }


  // BEAM
  if (simulateReaction==false && !useAlphaSource) {
    // Ion
    particleDef = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(beamZ,beamA,0.0);
    particleGun->SetParticleDefinition(particleDef);
    particleGun->SetParticleCharge(userCharge);    
    // Sample energy
    Ekin = energy;
    if (sigmaEnergy>0.) { 
      G4double mean = energy;
      G4double FWHM = sigmaEnergy/100.*energy;
      G4double std = FWHM/2.35; 
      Ekin = CLHEP::RandGauss::shoot(mean,std);
    }
//---------------------------------------------------------------------------------------//
    //energy including spread
    particleGun->SetParticleEnergy(Ekin *MeV); 
    //fixed energy
    // particleGun->SetParticleEnergy(energy *MeV);
//---------------------------------------------------------------------------------------//
    // Sample position
    G4double xBeam=0.*m, yBeam=0.*m;
    G4double r=0;
    G4double rmax = beamSpotDiameter / 2.0;
    if (beamSpotDiameter>0.) {
      r = G4UniformRand() * rmax;
      G4double phi = G4UniformRand()*CLHEP::twopi;
      xBeam = r*std::cos(phi);
      yBeam = r*std::sin(phi);
    }
    // beam particle z emission location is set to 10 Angstroms, i.e. immediately, in front of the target
    G4double zemit = targetZoffset - targetThickness/2 - 10*angstrom;
//---------------------------------------------------------------------------------------//
    //random emittance off optical axis
    particleGun->SetParticlePosition(G4ThreeVector(xBeam,yBeam,zemit));
    //fix emittance location to optical axis
    // particleGun->SetParticlePosition(G4ThreeVector(0,0,zemit));
//---------------------------------------------------------------------------------------//
    // Determine max angle from normalized transverse emittance
    G4double mass = particleDef->GetPDGMass();
    G4double Etot = Ekin + mass;
    G4double gamma = Etot/mass;
    G4double beta = sqrt(1.0-1.0/(gamma*gamma));
    G4double MaxAngle = 0*deg;
    if (rmax>0) {
      MaxAngle = transEmittance/(gamma*beta*rmax);
      MaxAngle = MaxAngle * 180./CLHEP::pi/1000. * deg; // mrad to deg conversion
    }
    //    G4cout << "MAX RADIUS: " << rmax/mm << " mm" << G4endl;
    //    G4cout << "MAX ANGLE: " << MaxAngle/deg << " deg" << G4endl;
    // Sample angle
    G4double x=0., y=0., z=1.;
//---------------------------------------------------------------------------------------//
    //random angles
    /*if (MaxAngle>0.) {
      G4double theta = G4UniformRand() * MaxAngle * sqrt(1.0-(r/rmax)*(r/rmax));
      G4double phi = G4UniformRand()*CLHEP::twopi;
      x = sin(theta) * cos(phi);
      y = sin(theta) * sin(phi);
      z = cos(theta);
    }*/
    //fixed angles
    // Edit: MQ random angles: 
    G4double theta = 0*deg;
    G4double THETA = 0*deg;
    G4double theta0 = Angle;
    G4double phi = 0*deg;
    x = sin(theta) * cos(phi);
    y = sin(theta) * sin(phi);
    z = cos(theta);
    if (MaxAngle>0. && rmax>0.) {
      theta = G4UniformRand() * MaxAngle * sqrt(1.0-(r/rmax)*(r/rmax));
      phi = G4UniformRand()*CLHEP::twopi;
      THETA = theta0 + theta*cos(phi);
      x = sin(THETA);
      y = sin(theta) * sin(phi);
      z = cos(THETA);
    } else {
      x = sin(theta0) * cos(phi);
      y = sin(theta0) * sin(phi);
      z = cos(theta0);
    }
    particleGun->SetParticleMomentumDirection(G4ThreeVector(x,y,z));
    
    G4cout<<"Prim.Gen.Action output "<<"Energy(MeV)= "<<energy <<" z emission location (mm) "
          <<zemit/mm<<" theta(deg)= "<<theta/deg<<" phi(deg)= "<<phi/deg<<G4endl;
//---------------------------------------------------------------------------------------//
  }


  // REACTION (values read in from beam.dat in EMMAapp)
  else if (simulateReaction) {
    const G4long idLong = reactionGeneratedEventCount++; // Monotonic sample index across chunked BeamOn calls.
    if (idLong < 0 || idLong > std::numeric_limits<G4int>::max()) {
      std::ostringstream msg;
      msg << "Reaction generated-event index overflow: " << idLong;
      G4Exception("EMMAPrimaryGeneratorAction::GeneratePrimaries", "EMMA0014",
                  FatalException, msg.str().c_str());
    }
    const G4int id = static_cast<G4int>(idLong);
    if (id >= beamInputRowsLoaded) {
      std::ostringstream msg;
      msg << "Reaction event index " << id
          << " exceeds loaded beam rows (" << beamInputRowsLoaded
          << "). Generate/load more beam rows before this reaction chunk.";
      G4Exception("EMMAPrimaryGeneratorAction::GeneratePrimaries", "EMMA0010",
                  FatalException, msg.str().c_str());
    }
    Ekin = energyBeam[id]; //from initializeReactionSimulation()
    beamEnergyAtTarget = Ekin;
    G4ThreeVector dir(dirxBeam[id],diryBeam[id],dirzBeam[id]);	  
    if (fZ1==0.) {
      G4cout << "ERROR: Two-body reaction not defined" << G4endl;
      exit (EXIT_FAILURE);
    }
    simulateTwoBodyReaction( Ekin, dir );

    // Select which reaction product is tracked as the Geant4 primary.
    G4int generatedZ = fZ3;                   // Default: recoil (product #3), original behavior.
    G4int generatedA = fA3;                   // Default: recoil mass number.
    G4double generatedEx = fExcitationEnergy3; // Default: recoil excitation energy.
    G4double generatedCharge = userCharge;    // Default: recoil charge state from input file.
    if (kGenerateTritonToS3) {
      generatedZ = fZ4;                       // Switch to ejectile/triton (product #4).
      generatedA = fA4;                       // Switch to ejectile mass number.
      generatedEx = 0.0;                      // Triton is generated in ground state.
      generatedCharge = static_cast<G4double>(generatedZ); // Use physical triton charge (+1).
    }

    particleDef = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(generatedZ, generatedA, generatedEx);  // Create selected ion.
    particleGun->SetParticleDefinition(particleDef);

    particleGun->SetParticleEnergy(Ekin*MeV);
    particleGun->SetParticleMomentumDirection(dir);
    particleGun->SetParticleCharge(generatedCharge);
    G4double x=posxBeam[id]*mm, y=posyBeam[id]*mm, z=poszBeam[id]*mm;
    G4double dz=depth/2.;
    z = z-dz; // correction needed because target placement refers to center of target ...
    particleGun->SetParticlePosition(G4ThreeVector( x, y, z ));
    
  }
  
  // Edit by MQ: guarding against exceptions
  particleGun->GeneratePrimaryVertex(anEvent);
  G4ThreeVector direction = particleGun->GetParticleMomentumDirection();

  // Print info:
  G4bool printInfo=false;
  G4double mass = particleDef->GetPDGMass();
  G4double pp = std::sqrt((Ekin + mass)*(Ekin + mass)-(mass*mass));
  G4double charge = particleGun->GetParticleCharge();	
  mass = mass / 931.494061;  // convert to amu
  if (printInfo) {
    G4cout << "\n Momentum " << pp << ", Mass " << mass << " amu, Ekin " 
	   << Ekin << " MeV, Charge " << charge << G4endl;
  }

}





void EMMAPrimaryGeneratorAction::initializeReactionSimulation() // called using /mydet/doReaction
{
  prepareBeam = false;
  simulateReaction = true;
  enforceS3AcceptedRowTarget = true;   // Option 2: keep generating until we collect requested accepted S3 rows.
  s3AcceptedRowTarget = nEvents;       // Requested number of accepted output rows comes from user-set nEvents.
  s3AcceptedRowCount = 0;              // Reset accepted-row counter at run start.
  s3GlobalEventSerial = -1;            // Reset monotonic global event serial written to output file.
  reactionGeneratedEventCount = 0;     // Reset generated-event index for this reaction run.
  beamInputRowsLoaded = 0;             // Reset loaded beam-row count; filled below.
  userCharge = fCharge3; //read in from reaction.dat in EMMAapp
  std::ofstream outfile; 
  focalPlaneFileName = UserDir;
  focalPlaneFileName.append("/ExcitationEnergy/fp_reaction.dat"); //Used in EMMADriftChamberHit
  outfile.open (focalPlaneFileName);
  outfile.close();
  postTargetFileName = UserDir;
  postTargetFileName.append("/ExcitationEnergy/postTarget_reaction.dat"); //Used in EMMASteppingAction
  outfile.open (postTargetFileName);
  outfile.close();
  // Ejectiles
  postTargetEjectileFileName = UserDir;

  std::ostringstream excitationTag;
  excitationTag.setf(std::ios::fixed);
  excitationTag.precision(3);
  excitationTag << fExcitationEnergy3 / MeV;
  G4String excitationSuffix = excitationTag.str();
  excitationSuffix = excitationSuffix.replace(excitationSuffix.find("."), 1, "p");
  std::ostringstream angleTag;
  angleTag.setf(std::ios::fixed);
  angleTag.precision(2);
  angleTag << tritonLabAngleMinDeg << "to" << tritonLabAngleMaxDeg;
  G4String angleSuffix = angleTag.str();
  angleSuffix = angleSuffix.replace(angleSuffix.find("."), 1, "p");
  angleSuffix = angleSuffix.replace(angleSuffix.find("."), 1, "p");
  postTargetEjectileFileName.append("/ExcitationEnergy/Ex");
  postTargetEjectileFileName.append(excitationSuffix);
  postTargetEjectileFileName.append("MeV_");
  postTargetEjectileFileName.append(angleSuffix);
  postTargetEjectileFileName.append("deg.dat");
  outfile.open(postTargetEjectileFileName);
  outfile.close();

  postTargetRecoilSpectrometerFileName = UserDir;
  postTargetRecoilSpectrometerFileName.append("/ExcitationEnergy/SpecEx");
  postTargetRecoilSpectrometerFileName.append(excitationSuffix);
  postTargetRecoilSpectrometerFileName.append("MeV_");
  postTargetRecoilSpectrometerFileName.append(angleSuffix);
  postTargetRecoilSpectrometerFileName.append("deg_angledeg.dat");
  outfile.open(postTargetRecoilSpectrometerFileName);
  outfile.close();

  postDegrader1FileName = UserDir;
  postDegrader1FileName.append("/ExcitationEnergy/postDegrader1_reaction.dat"); //Used in EMMASteppingAction
  outfile.open (postDegrader1FileName);
  outfile.close();

  // Create/clear plain-text S3 triton observables file for this run.
  s3TritonOutputFileName = UserDir;
  s3TritonOutputFileName.append("/ExcitationEnergy/S3_triton_ring_observables.dat");
  std::ofstream s3Outfile(s3TritonOutputFileName, std::ios::trunc);
  s3Outfile << "# beamEnergyAtTarget_MeV,tritonTargetExitEnergy_MeV,tritonTargetExitTheta_deg,tritonTargetExitPhi_deg,s3EnergyLoss_MeV,tritonThetaBeforeS3_deg,tritonPhiBeforeS3_deg,tritonEnergyAtReactionPoint_MeV,tritonThetaAtReactionPoint_deg,tritonPhiAtReactionPoint_deg,s3EnergySmeared_MeV,tritonTargetEnergyLoss_MeV" << G4endl;
  s3Outfile << "# ringGateCopyNo=" << tritonLabGateRingNumber
            << " thetaMinDeg=" << tritonLabAngleMinDeg
            << " thetaMaxDeg=" << tritonLabAngleMaxDeg << G4endl;
  s3Outfile.close();

  // Option 2: run with an oversampling budget and stop early when accepted S3 rows reach target.
  G4long budgetLong = static_cast<G4long>(nEvents) * static_cast<G4long>(kS3AcceptedRowSamplingBudgetFactor); // Oversampling budget.
  if (budgetLong < nEvents) budgetLong = nEvents; // Guard against overflow wrap on multiplication.
  if (budgetLong > std::numeric_limits<G4int>::max()) budgetLong = std::numeric_limits<G4int>::max(); // BeamOn takes int.
  const G4int generatedEventBudget = static_cast<G4int>(budgetLong); // Final generated-event cap for this run.

  // Utility: count how many prepared beam rows are currently available in beam.dat.
  auto countPreparedBeamRows = [&]() -> G4long {
    std::ifstream beamCountFile(inTargetFileName, std::ios::in);
    if (!beamCountFile.is_open()) return 0;
    G4long count = 0;
    G4double e = 0., px = 0., py = 0., pz = 0., dx = 0., dy = 0., dz = 0.;
    while (beamCountFile >> e >> px >> py >> pz >> dx >> dy >> dz) {
      ++count;
    }
    return count;
  };

  // Utility: append additional prepared beam rows until at least requiredRows exist.
  auto ensurePreparedRows = [&](G4long requiredRows) -> bool {
    G4long currentRows = countPreparedBeamRows();
    if (currentRows >= requiredRows) return true;

    G4cout << "Reaction generator: prepared beam rows " << currentRows
           << " < required " << requiredRows
           << ". Generating more beam samples on demand." << G4endl;

    const G4bool savedPrepareBeam = prepareBeam;
    const G4bool savedSimulateReaction = simulateReaction;
    const G4bool savedEnforceTarget = enforceS3AcceptedRowTarget;
    const G4bool savedAbortActive = s3AcceptedRowAbortActive;
    const G4double savedUserCharge = userCharge;

    prepareBeam = true;                // Reuse existing beam-preparation stepping/writing path.
    simulateReaction = false;
    enforceS3AcceptedRowTarget = false;
    s3AcceptedRowAbortActive = false;
    userCharge = beamCharge;           // Beam charge for preparation.

    G4long previousRows = currentRows;
    G4int stagnantPasses = 0;
    while (currentRows < requiredRows && stagnantPasses < 4) {
      G4long missingRows = requiredRows - currentRows;
      G4long batchLong = std::min<G4long>(kBeamTopUpChunkEvents, missingRows);
      if (batchLong < 1) batchLong = 1;
      if (batchLong > std::numeric_limits<G4int>::max()) {
        batchLong = std::numeric_limits<G4int>::max();
      }
      G4RunManager::GetRunManager()->BeamOn(static_cast<G4int>(batchLong));
      currentRows = countPreparedBeamRows();
      if (currentRows <= previousRows) {
        ++stagnantPasses; // Guard against infinite loops when preparation produces no new rows.
      } else {
        stagnantPasses = 0;
      }
      previousRows = currentRows;
    }

    prepareBeam = savedPrepareBeam;
    simulateReaction = savedSimulateReaction;
    enforceS3AcceptedRowTarget = savedEnforceTarget;
    s3AcceptedRowAbortActive = savedAbortActive;
    userCharge = savedUserCharge;

    return currentRows >= requiredRows;
  };

  // Utility: load first requiredRows prepared rows into reaction memory arrays.
  auto loadPreparedRows = [&](G4int requiredRows) -> bool {
    if (requiredRows <= 0) return false;
    delete[] energyBeam; energyBeam = nullptr;
    delete[] posxBeam; posxBeam = nullptr;
    delete[] posyBeam; posyBeam = nullptr;
    delete[] poszBeam; poszBeam = nullptr;
    delete[] dirxBeam; dirxBeam = nullptr;
    delete[] diryBeam; diryBeam = nullptr;
    delete[] dirzBeam; dirzBeam = nullptr;

    energyBeam = new G4double[requiredRows];
    posxBeam = new G4double[requiredRows];
    posyBeam = new G4double[requiredRows];
    poszBeam = new G4double[requiredRows];
    dirxBeam = new G4double[requiredRows];
    diryBeam = new G4double[requiredRows];
    dirzBeam = new G4double[requiredRows];

    std::ifstream beamFile(inTargetFileName, std::ios::in);
    if (!beamFile.is_open()) return false;
    for (G4int i = 0; i < requiredRows; ++i) {
      if (!(beamFile >> energyBeam[i]
                    >> posxBeam[i]
                    >> posyBeam[i]
                    >> poszBeam[i]
                    >> dirxBeam[i]
                    >> diryBeam[i]
                    >> dirzBeam[i])) {
        return false;
      }
    }
    beamInputRowsLoaded = requiredRows;
    return true;
  };

  G4cout << "S3 accepted-row mode enabled: targetRows=" << s3AcceptedRowTarget
         << ", generatedEventBudget=" << generatedEventBudget
         << ", reactionChunkSize=" << kReactionBeamOnChunkEvents
         << ", beamTopUpChunkSize=" << kBeamTopUpChunkEvents << G4endl;

  G4int generatedEventsRemaining = generatedEventBudget;
  while (generatedEventsRemaining > 0 && s3AcceptedRowCount < s3AcceptedRowTarget) {
    const G4int reactionChunk = std::min(generatedEventsRemaining, kReactionBeamOnChunkEvents);
    const G4long requiredRows = reactionGeneratedEventCount + reactionChunk;
    if (requiredRows > std::numeric_limits<G4int>::max()) {
      std::ostringstream msg;
      msg << "Requested prepared-row load exceeds G4int capacity: " << requiredRows;
      G4Exception("EMMAPrimaryGeneratorAction::initializeReactionSimulation", "EMMA0015",
                  FatalException, msg.str().c_str());
    }

    if (!ensurePreparedRows(requiredRows)) {
      std::ostringstream msg;
      msg << "Unable to prepare enough beam rows on demand. requiredRows=" << requiredRows;
      G4Exception("EMMAPrimaryGeneratorAction::initializeReactionSimulation", "EMMA0012",
                  FatalException, msg.str().c_str());
    }
    if (requiredRows > beamInputRowsLoaded) {
      if (!loadPreparedRows(static_cast<G4int>(requiredRows))) {
        std::ostringstream msg;
        msg << "Failed to load prepared beam rows from " << inTargetFileName
            << " up to row " << requiredRows;
        G4Exception("EMMAPrimaryGeneratorAction::initializeReactionSimulation", "EMMA0013",
                    FatalException, msg.str().c_str());
      }
    }

    const G4long generatedBefore = reactionGeneratedEventCount;
    s3AcceptedRowAbortActive = true; // Enable abort-on-target only for this internal reaction run.
    G4RunManager::GetRunManager()->BeamOn(reactionChunk);
    s3AcceptedRowAbortActive = false; // Disable abort-on-target between chunks and after run completion.
    const G4long generatedThisChunk = reactionGeneratedEventCount - generatedBefore;
    if (generatedThisChunk <= 0) {
      G4cout << "Warning: reaction chunk produced zero generated events; stopping early." << G4endl;
      break;
    }
    if (generatedThisChunk >= generatedEventsRemaining) {
      generatedEventsRemaining = 0;
    } else {
      generatedEventsRemaining -= static_cast<G4int>(generatedThisChunk);
    }
  }

  if (s3AcceptedRowCount < s3AcceptedRowTarget) { // Notify user if budget exhausted before target accepted rows.
    G4cout << "Warning: accepted S3 rows (" << s3AcceptedRowCount
           << ") below target (" << s3AcceptedRowTarget
           << "). Increase kS3AcceptedRowSamplingBudgetFactor if needed." << G4endl;
  } else {
    G4cout << "Reached S3 accepted-row target: " << s3AcceptedRowCount << G4endl;
  }
}


void EMMAPrimaryGeneratorAction::initializeBeamSimulation() // called using /mydet/doBeam
{
  prepareBeam = false;
  simulateReaction = false;
  delete[] energyBeam; energyBeam = nullptr;
  delete[] posxBeam; posxBeam = nullptr;
  delete[] posyBeam; posyBeam = nullptr;
  delete[] poszBeam; poszBeam = nullptr;
  delete[] dirxBeam; dirxBeam = nullptr;
  delete[] diryBeam; diryBeam = nullptr;
  delete[] dirzBeam; dirzBeam = nullptr;
  beamInputRowsLoaded = 0;
  reactionGeneratedEventCount = 0;
  s3AcceptedRowAbortActive = false; // Ensure manual runs don't inherit abort-on-target behavior.
  enforceS3AcceptedRowTarget = false; // Disable accepted-row stopping logic outside reaction mode.
  s3AcceptedRowTarget = 0;            // Clear target rows when not in reaction mode.
  s3AcceptedRowCount = 0;             // Clear running accepted-row counter.
  s3GlobalEventSerial = -1;           // Reset global event serial for non-reaction runs.
  s3TritonOutputFileName = ""; // Disable S3 triton text output in pure beam mode.
  userCharge = beamCharge; //read in from beam.dat in EMMAapp
  std::ofstream outfile;
  focalPlaneFileName = UserDir;
  focalPlaneFileName.append("/ExcitationEnergy/fp_beam.dat"); //Used in EMMADriftChamberHit
  outfile.open (focalPlaneFileName);
  outfile.close();	  
  postTargetEjectileFileName = "";
  postTargetRecoilSpectrometerFileName = "";
  postTargetFileName = UserDir;
  postTargetFileName.append("/ExcitationEnergy/postTarget_beam.dat"); //Used in EMMASteppingAction
  outfile.open (postTargetFileName);
  outfile.close();
  postDegrader1FileName = UserDir;
  postDegrader1FileName.append("/ExcitationEnergy/postDegrader1_beam.dat"); //Used in EMMASteppingAction
  outfile.open (postDegrader1FileName);
  outfile.close();

  // simulated nEvents
  G4RunManager::GetRunManager()->BeamOn(nEvents);
}



void EMMAPrimaryGeneratorAction::initializeBeamPreparation() // called using /mydet/doPrepare
{
  prepareBeam = true;
  simulateReaction = false;
  delete[] energyBeam; energyBeam = nullptr;
  delete[] posxBeam; posxBeam = nullptr;
  delete[] posyBeam; posyBeam = nullptr;
  delete[] poszBeam; poszBeam = nullptr;
  delete[] dirxBeam; dirxBeam = nullptr;
  delete[] diryBeam; diryBeam = nullptr;
  delete[] dirzBeam; dirzBeam = nullptr;
  beamInputRowsLoaded = 0;
  reactionGeneratedEventCount = 0;
  s3AcceptedRowAbortActive = false; // Ensure manual runs don't inherit abort-on-target behavior.
  enforceS3AcceptedRowTarget = false; // Disable accepted-row stopping logic outside reaction mode.
  s3AcceptedRowTarget = 0;            // Clear target rows when not in reaction mode.
  s3AcceptedRowCount = 0;             // Clear running accepted-row counter.
  s3GlobalEventSerial = -1;           // Reset global event serial for non-reaction runs.
  s3TritonOutputFileName = ""; // Disable S3 triton text output in beam-preparation mode.
  userCharge = beamCharge; //read in from beam.dat in EMMAapp
  std::ofstream outfile; 
  outfile.open (inTargetFileName); //declared in constructor
  outfile.close();
  postTargetEjectileFileName = "";

  // simulated nEvents
  G4RunManager::GetRunManager()->BeamOn(nEvents);
}


void EMMAPrimaryGeneratorAction::simulateTwoBodyReaction( G4double &Ebeam, G4ThreeVector &dir )
{
  G4int Z1 = fZ1;
  G4int A1 = fA1;
  G4int Z2 = fZ2;
  G4int A2 = fA2;

  // masses of 1+2
  G4double m1 = G4NucleiProperties::GetNuclearMass(A1, Z1);
  G4double m2 = G4NucleiProperties::GetNuclearMass(A2, Z2);

  // Z and A of reaction products (3+4):
  G4int Z3 = fZ3;
  G4int A3 = fA3;
  G4int Z4 = fZ4;
  G4int A4 = fA4;

  // masses of 3+4
  G4double m3 = G4NucleiProperties::GetNuclearMass(A3, Z3);
  G4double m4 = G4NucleiProperties::GetNuclearMass(A4, Z4);

  // take into account excitation energy of fragment 3
  m3 = m3 + fExcitationEnergy3;

  // 4-momentum of beam
  G4double p1n = sqrt( (m1+Ebeam)*(m1+Ebeam) - m1*m1 );
  G4double dirn = sqrt( dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2] );
  G4ThreeVector p1 = dir/dirn*p1n;

  // determine velocity of CM frame relative to LAB frame
  G4LorentzVector lv1(p1,m1+Ebeam);
  G4LorentzVector lv2(0.0,0.0,0.0,m2);   
  G4LorentzVector lv = lv1 + lv2;
  G4ThreeVector bst = lv.boostVector(); // divides the spatial component (p) by the time component (E)

  // transform 4-momenta to CM frame
  lv1.boost(-bst);
  lv2.boost(-bst);
  G4double etot = lv1[3] + lv2[3]; // total energy in CM

  // if energy is insufficient, nothing happens
  if (etot<m3+m4) {
    // G4cout << "NOT ENOUGH ENERGY FOR REACTION" << G4endl;
    // edit: 
    std::ostringstream msg;
    msg << "Not enough energy for reaction. etot=" << etot
        << " m3+m4=" << (m3 + m4);
    G4Exception("EMMAPrimaryGeneratorAction::simulateTwoBodyReaction", "EMMA0005",
                EventMustBeAborted, msg.str().c_str());
    return;
  }

  // Compute c.m. energies and momentum of reaction products (3+4)
  G4double e3 = ( etot*etot + m3*m3 - m4*m4 ) / (2*etot); 
  G4double e4 = etot - e3; 
  G4double pcm = sqrt( e3*e3 - m3*m3 ); 

  // Max and min angles
  G4double fqrmax = (180-fqmin/deg)*deg; //compute recoil c.m. angles from ejectile c.m. angles
  G4double fqrmin = (180-fqmax/deg)*deg;
  G4double thetaCMmin = fqrmin;
  G4double t1 = std::cos(thetaCMmin/rad);
  G4double thetaCMmax = fqrmax;
  G4double t2 = std::cos(thetaCMmax/rad);
  
  // // Sampling of directions in CM system
  // G4double t    = G4UniformRand();
  // G4double phi  = G4UniformRand()*CLHEP::twopi;
  // G4double cost = t1 - (t1-t2)*t;
  // G4double sint = std::sqrt((1.0-cost)*(1.0+cost));
  
  // // Lorentz vectors of reaction products (3+4)
  // G4ThreeVector v3(sint*std::cos(phi),sint*std::sin(phi),cost);
  // v3 = v3 * pcm;
  // G4ThreeVector v4 = -v3;
  // G4LorentzVector lv3(v3.x(),v3.y(),v3.z(),e3);
  // G4LorentzVector lv4(v4.x(),v4.y(),v4.z(),e4);

  
  // // Transform to LAB frame
  // lv3.boost(bst);
  // lv4.boost(bst);

  // Edit (MQ): Rejection sampling on tritons on S3

  const G4int maxGateTrials = 10000;
  G4bool gateAccepted = false;
  G4double cost = 0.0;
  G4int acceptedTrial = -1;
  for (G4int trial = 0; trial < maxGateTrials; ++trial) {
    // Sampling of directions in CM system
    G4double t    = G4UniformRand();
    G4double phi  = G4UniformRand()*CLHEP::twopi;
    cost = t1 - (t1-t2)*t;
    G4double sint = std::sqrt((1.0-cost)*(1.0+cost));
    
    // Lorentz vectors of reaction products (3+4)
    G4ThreeVector v3(sint*std::cos(phi),sint*std::sin(phi),cost);
    v3 = v3 * pcm;
    G4ThreeVector v4 = -v3;
    G4LorentzVector lv3(v3.x(),v3.y(),v3.z(),e3);
    G4LorentzVector lv4(v4.x(),v4.y(),v4.z(),e4);

    // Transform to LAB frame
    lv3.boost(bst);
    lv4.boost(bst);
  
    
    // Return selected product kinematics to the primary generator.
    if (kGenerateTritonToS3) {
      Ebeam = lv4[3] - m4; // Use triton/ejectile (product #4) kinetic energy.
      dir[0] = lv4[0];     // Use triton/ejectile momentum x component.
      dir[1] = lv4[1];     // Use triton/ejectile momentum y component.
      dir[2] = lv4[2];     // Use triton/ejectile momentum z component.
    } else {
      Ebeam = lv3[3] - m3; // Use recoil (product #3) kinetic energy, original behavior.
      dir[0] = lv3[0];     // Use recoil momentum x component.
      dir[1] = lv3[1];     // Use recoil momentum y component.
      dir[2] = lv3[2];     // Use recoil momentum z component.
    }

    // Store ejectile (product #4) kinematics for logging at target exit.
    ejectileEnergy = lv4[3] - m4;
    G4double ejectileP = std::sqrt(lv4[0]*lv4[0] + lv4[1]*lv4[1] + lv4[2]*lv4[2]);
    if (IsFinite(ejectileP) && ejectileP > 0.) {
      ejectileDirX = lv4[0] / ejectileP;
      ejectileDirY = lv4[1] / ejectileP;
      ejectileDirZ = lv4[2] / ejectileP;

    } else {
      ejectileDirX = 0.;
      ejectileDirY = 0.;
      ejectileDirZ = 0.;
    }

    // Store recoil (product #3) lab kinematics for per-event S3 output rows.
    recoilLabEnergy = lv3[3] - m3;
    G4double recoilLabP = std::sqrt(lv3[0]*lv3[0] + lv3[1]*lv3[1] + lv3[2]*lv3[2]);
    if (IsFinite(recoilLabP) && recoilLabP > 0.) {
      recoilLabDirX = lv3[0] / recoilLabP;
      recoilLabDirY = lv3[1] / recoilLabP;
      recoilLabDirZ = lv3[2] / recoilLabP;
    } else {
      recoilLabDirX = 0.;
      recoilLabDirY = 0.;
      recoilLabDirZ = 0.;
    }

    if (!applyTritonLabAngleGate ||
        IsTritonLabAngleAccepted(ejectileDirZ, tritonLabAngleMinDeg, tritonLabAngleMaxDeg)) {
      gateAccepted = true;
      acceptedTrial = trial;
      break;
    }
  }
  // Rejection + efficency tests
  if (!gateAccepted) {
        std::ostringstream msg;
        msg << "Failed to generate event within triton lab-angle gate after "
        << maxGateTrials << " trials.";
    G4Exception("EMMAPrimaryGeneratorAction::simulateTwoBodyReaction", "EMMA0008", EventMustBeAborted, msg.str().c_str());
    return;
  }
  if (applyTritonLabAngleGate) {
    gateTrialTotal += (acceptedTrial + 1);
    ++gateAcceptedEvents;
    if (gateAcceptedEvents % 1000 == 0) {
      const G4double efficiency = static_cast<G4double>(gateAcceptedEvents)
                                  / static_cast<G4double>(gateTrialTotal);
      const G4double avgTrials = static_cast<G4double>(gateTrialTotal)
                                 / static_cast<G4double>(gateAcceptedEvents);
      G4cout << "Triton gate efficiency: " << efficiency
             << " (avg trials/event=" << avgTrials << ")" << G4endl;
    }
  }
  ejectileZ = Z4;
  ejectileA = A4;
  recoilThetaCM = std::acos(cost);
  ejectileThetaCM = CLHEP::pi - recoilThetaCM;

}
