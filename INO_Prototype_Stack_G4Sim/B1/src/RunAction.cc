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
//
/// \file B1/src/RunAction.cc
/// \brief Implementation of the B1::RunAction class

#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
// #include "Run.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4RootAnalysisManager.hh"
#include "G4CsvAnalysisManager.hh"
#include "G4AnalysisManager.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RunAction::RunAction()
 : G4UserRunAction(),
     fEdep(0.),
     fEdep2(0.)
{
 
  auto analysisManager = G4AnalysisManager::Instance();

  // Create directories
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);
  analysisManager->OpenFile("100GeV.root");
  analysisManager->CreateNtuple("mu_e_scatt", "3D hits and time");
  analysisManager->CreateNtupleDColumn("id");
  analysisManager->CreateNtupleDColumn("x");
  analysisManager->CreateNtupleDColumn("y");
  analysisManager->CreateNtupleDColumn("z");
  analysisManager->CreateNtupleDColumn("t");
  analysisManager->CreateNtupleDColumn("Px");
  analysisManager->CreateNtupleDColumn("Py");
  analysisManager->CreateNtupleDColumn("Pz");
  analysisManager->CreateNtupleDColumn("particle_id");
  analysisManager->CreateNtupleDColumn("energy");
  

  // analysisManager->CreateNtupleDColumn("E");
  // analysisManager->CreateNtupleDColumn("theta");
  // analysisManager->CreateNtupleDColumn("phi");
  analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{
  // inform the runManager to save random number seed
  // G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  // // reset accumulables to their initial values
  // G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  // accumulableManager->Reset();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* run)
{
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::AddEdep(G4double edep)
{
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}





















// RunAction::RunAction()
//  : G4UserRunAction(),
//       fEdep(0.),
//       fEdep2(0.)
// {
//   // add new units for dose
//   //
//   const G4double milligray = 1.e-3*gray;
//   const G4double microgray = 1.e-6*gray;
//   const G4double nanogray  = 1.e-9*gray;
//   const G4double picogray  = 1.e-12*gray;

//   new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
//   new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
//   new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
//   new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray);

//   // Register accumulable to the accumulable manager
//   G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
//   accumulableManager->RegisterAccumulable(fEdep);
//   accumulableManager->RegisterAccumulable(fEdep2);
// }


// RunAction::~RunAction()
// {
// }

// void RunAction::BeginOfRunAction(const G4Run*)
// {
//   // inform the runManager to save random number seed
//   G4RunManager::GetRunManager()->SetRandomNumberStore(false);

//   // reset accumulables to their initial values
//   G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
//   accumulableManager->Reset();
//   // G4CsvAnalysisManager *man = G4CsvAnalysisManager::Instance();
//   G4RootAnalysisManager *man = G4RootAnalysisManager::Instance();
//   man->SetNtupleMerging(true); // very important

//   // Open an output file

//   //man->SetNtupleMerging(true); // very important; also enable /run/numberOfThreads 4 in mac file 
//   man->OpenFile("Al_mue.root");
//   // man->OpenFile("mu_e.csv");
//   // Create ntuple
//   man->CreateNtuple("mu_e_scatt", "Position of mu and e");
//   man->CreateNtupleDColumn("event_id");
//   man->CreateNtupleDColumn("mu_x");
//   man->CreateNtupleDColumn("mu_y");
//   man->CreateNtupleDColumn("mu_z");
//   man->CreateNtupleDColumn("e_x");
//   man->CreateNtupleDColumn("e_y");
//   man->CreateNtupleDColumn("e_z");
//   man->CreateNtupleDColumn("mu_E");
//   man->CreateNtupleDColumn("e_E");
//   man->CreateNtupleDColumn("mu_p_x");
//   man->CreateNtupleDColumn("mu_p_y");
//   man->CreateNtupleDColumn("mu_p_z");
//   man->CreateNtupleDColumn("e_p_x");
//   man->CreateNtupleDColumn("e_p_y");
//   man->CreateNtupleDColumn("e_p_z");
//   man->CreateNtupleDColumn("scat_type"); 
//   G4cout<<"............................"<<G4endl;
//   G4cout<<"Added electron Energy branch"<<G4endl;
//   // G4cout<<"Added electron Energy branch"<<G4endl;
 
//   man->FinishNtuple();
// }

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......





// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// void RunAction::EndOfRunAction(const G4Run* run)
// {

//    G4cout<<"Closing file"<<G4endl;
//   G4RootAnalysisManager *man = G4RootAnalysisManager::Instance();
//   // G4CsvAnalysisManager *man = G4CsvAnalysisManager::Instance();
//   man->Write();
//   man->CloseFile();


//   // G4int nofEvents = run->GetNumberOfEvent();
//   // if (nofEvents == 0) return;

//   // // Merge accumulables
//   // G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
//   // accumulableManager->Merge();

//   // // Compute dose = total energy deposit in a run and its variance
//   // //
//   // G4double edep  = fEdep.GetValue();
//   // G4double edep2 = fEdep2.GetValue();

//   // G4double rms = edep2 - edep*edep/nofEvents;
//   // if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;

//   // const auto detConstruction = static_cast<const DetectorConstruction*>(
//   //   G4RunManager::GetRunManager()->GetUserDetectorConstruction());
//   // G4double mass = detConstruction->GetScoringVolume()->GetMass();
//   // G4double dose = edep/mass;
//   // G4double rmsDose = rms/mass;

//   // // Run conditions
//   // //  note: There is no primary generator action object for "master"
//   // //        run manager for multi-threaded mode.
//   // const auto generatorAction = static_cast<const PrimaryGeneratorAction*>(
//   //   G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
//   // G4String runCondition;
//   // if (generatorAction)
//   // {
//   //   const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
//   //   runCondition += particleGun->GetParticleDefinition()->GetParticleName();
//   //   runCondition += " of ";
//   //   G4double particleEnergy = particleGun->GetParticleEnergy();
//   //   runCondition += G4BestUnit(particleEnergy,"Energy");
//   }

// //   // Print
// //   //
// //   if (IsMaster()) {
// //     G4cout
// //      << G4endl
// //      << "--------------------End of Global Run-----------------------";
// //   }
// //   else {
// //     G4cout
// //      << G4endl
// //      << "--------------------End of Local Run------------------------";
// //   }

// //   G4cout
// //      << G4endl
// //      << " The run consists of " << nofEvents << " "<< runCondition
// //      << G4endl
// //      << " Cumulated dose per run, in scoring volume : "
// //      << G4BestUnit(dose,"Dose") << " rms = " << G4BestUnit(rmsDose,"Dose")
// //      << G4endl
// //      << "------------------------------------------------------------"
// //      << G4endl
// //      << G4endl;
// // }

// // //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// void RunAction::AddEdep(G4double edep)
// {
//   fEdep  += edep;
//   fEdep2 += edep*edep;
// }

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// }
