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
/// \file B1/src/SteppingAction.cc
/// \brief Implementation of the B1::SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "G4RootAnalysisManager.hh"
#include "G4SystemOfUnits.hh"


#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4EventManager.hh"
#include <G4VPhysicalVolume.hh>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace B1
{

SteppingAction::SteppingAction(EventAction *eventAction)
    : G4UserSteppingAction(),
      fEventAction(eventAction),
      fScoringVolume(0)
{
  analysisManager = G4AnalysisManager::Instance();
}


SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  const DetectorConstruction* detConstruction= static_cast<const DetectorConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  // if (!fScoringVolume) {
  //   const DetectorConstruction* detConstruction
  //     = static_cast<const DetectorConstruction*>
  //       (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  //   fScoringVolume = detConstruction->GetScoringVolume();
  // }

  // // get volume of the current step
  // G4LogicalVolume* volume= step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  // G4cout<<"I am in: "<<volume->GetName()<<" at "<< step->GetPreStepPoint()->GetPosition().x()<< step->GetPreStepPoint()->GetPosition().y()<<step->GetPreStepPoint()->GetPosition().z()<<G4endl;
    // DetectorConstruction(G4VPhysicalVolume *setWorld = 0)
    
    // fWorld = setWorld;
    // G4cout<<fWorld->GetName()<<G4endl;
    // G4cout<<"Found "<<fWorld->GetLogicalVolume()->GetNoDaughters()<<" daughters in the mother volume ..."<<G4endl; 

  // int nDaughters = 21;
  // for (int i=0; i<nDaughters;i++){
     G4LogicalVolume* volume= step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
     G4double edepStep = step->GetTotalEnergyDeposit();
    //  if(volume->GetDaughter(i)->GetName().contains("GLASS"))
     if(volume->GetName().contains("pickupx"))
     {
      if(step->GetTrack()->GetDefinition()->GetPDGEncoding()== 13 || 11){
    //     const G4VTouchable *touchable = step->GetPreStepPoint()->GetTouchable();
    //    G4cout<<"I am in: "<<volume->GetName()<<" at " <<"x="<< step->GetPreStepPoint()->GetPosition().x()<< "y="<<step->GetPreStepPoint()->GetPosition().y()<<"z="<<step->GetPreStepPoint()->GetPosition().z()<<
    //  " energy: "<<step->GetTotalEnergyDeposit()<<G4endl;

      analysisManager->FillNtupleDColumn(0, G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID());
      analysisManager->FillNtupleDColumn(1, step->GetPreStepPoint()->GetPosition().x());
      analysisManager->FillNtupleDColumn(2, step->GetPreStepPoint()->GetPosition().y());
      analysisManager->FillNtupleDColumn(3, step->GetPreStepPoint()->GetPosition().z());
      analysisManager->FillNtupleDColumn(4, step->GetTrack()->GetGlobalTime());
      analysisManager->FillNtupleDColumn(5, step->GetTrack()->GetMomentum().x());
      analysisManager->FillNtupleDColumn(6, step->GetTrack()->GetMomentum().y());
      analysisManager->FillNtupleDColumn(7, step->GetTrack()->GetMomentum().z());
      analysisManager->FillNtupleDColumn(8, step->GetTrack()->GetDefinition()->GetPDGEncoding()); // added to get particleID
      analysisManager->FillNtupleDColumn(9, step->GetTrack()->GetKineticEnergy());
      analysisManager->AddNtupleRow();


      }
//         if(G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID() % 1000 ==0 )
//         G4cout<<G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID()<<" events done ..."<<G4endl;
//         steppingAction->Reset();
      // G4cout<<fWorld->GetLogicalVolume()->GetDaughter(i)->GetName()<<G4endl;
      // if(fWorld->GetLogicalVolume()->GetDaughter(i)->GetName().contains("GLASS")){
      // G4cout<<"The hit points are "<< step->GetPreStepPoint()->GetPosition().x()<< step->GetPreStepPoint()->GetPosition().y()<<step->GetPreStepPoint()->GetPosition().z()<<G4endl;
      }


}
}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// SteppingAction::~SteppingAction()
// {
// }

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// void SteppingAction::UserSteppingAction(const G4Step *step)
// {
//   if (!fScoringVolume) {
//   const DetectorConstruction* detectorConstruction
//     = static_cast<const DetectorConstruction*>
//       (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
//     fScoringVolume = detectorConstruction->GetScoringVolume();
//   }

//  if(step->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="glass_sens")
//   {

  
//     // G4cout << "The particle ID is " << step->GetTrack()->GetDefinition()->GetPDGEncoding() << G4endl;
//     if (step->GetTrack()->GetDefinition()->GetPDGEncoding() == 13)
//       SetMuPosition(step->GetTrack()->GetPosition(),step->GetTrack()->GetKineticEnergy(),step->GetTrack()->GetMomentum());//.x(),step->GetTrack()->GetMomentum().y());

//       // SetMuPosition(step->GetTrack()->GetPosition(),step->GetTrack()->GetKineticEnergy(),step->GetTrack()->GetMomentumDirection());//.x(),step->GetTrack()->GetMomentum().y());
//     else if (step->GetTrack()->GetDefinition()->GetPDGEncoding() == 11)
//       SetEPosition(step->GetTrack()->GetPosition(),step->GetTrack()->GetKineticEnergy(),step->GetTrack()->GetMomentum());//.x(),step->GetTrack()->GetMomentum().y());

//       // SetEPosition(step->GetTrack()->GetPosition(),step->GetTrack()->GetKineticEnergy(),step->GetTrack()->GetMomentumDirection());//.x(),step->GetTrack()->GetMomentum().y());
//     if(step->GetTrack()->GetDefinition()->GetPDGEncoding() != 13 && step->GetTrack()->GetDefinition()->GetPDGEncoding() != 11)
//       SetScatteringType(step->GetTrack()->GetDefinition()->GetPDGEncoding());
  
//   // //
//     // if(step->GetTrack()->GetDefinition()->GetPDGEncoding()==13){


//       G4RootAnalysisManager *man = G4RootAnalysisManager::Instance();

//     man->FillNtupleDColumn(0,G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID());
//     man->FillNtupleDColumn(1, step->GetPreStepPoint()->GetMuPosition().x());
//     man->FillNtupleDColumn(2, step->GetPreStepPoint()->GetMuPosition().y());
//     man->FillNtupleDColumn(3, step->GetPreStepPoint()->GetMuPosition().z());

//     man->FillNtupleDColumn(4, step->GetPreStepPoint()->GetEPosition().x());
//     man->FillNtupleDColumn(5, step->GetPreStepPoint()->GetEPosition().y());
//     man->FillNtupleDColumn(6, step->GetPreStepPoint()->GetEPosition().z());

//     man->FillNtupleDColumn(7, step->GetPreStepPoint()->GetMuEnergy());

//     man->FillNtupleDColumn(8, step->GetPreStepPoint()->GetEEnergy());

//     man->FillNtupleDColumn(9, step->GetPreStepPoint()->GetMuMomentum().x());
//     man->FillNtupleDColumn(10, step->GetPreStepPoint()->GetMuMomentum().y());
//     man->FillNtupleDColumn(11, step->GetPreStepPoint()->GetMuMomentum().z());

//     man->FillNtupleDColumn(12, step->GetPreStepPoint()->GetEMomentum().x());
//     man->FillNtupleDColumn(13, step->GetPreStepPoint()->GetEMomentum().y());
//     man->FillNtupleDColumn(14, step->GetPreStepPoint()->GetEMomentum().z());

//     man->FillNtupleDColumn(15, step->GetPreStepPoint()->GetScatteringType());
//     // man->FillNtupleDColumn(15, step->GetPreStepPoint()->GetMuMomentumDirection().x());
//     // man->FillNtupleDColumn(16, step->GetPreStepPoint()->GetMuMomentumDirection().y());
//     // man->FillNtupleDColumn(17, step->GetPreStepPoint()->GetMuMomentumDirection().z());

//     // man->FillNtupleDColumn(18, step->GetPreStepPoint()->GetEMomentumDirection().x());
//     // man->FillNtupleDColumn(19, step->GetPreStepPoint()->GetEMomentumDirection().y());
//     // man->FillNtupleDColumn(20, step->GetPreStepPoint()->GetEMomentumDirection().z());
//      man->AddNtupleRow(); 
// //   }
// //       // G4cout<<step->GetPreStepPoint()->GetPosition().x()<<step->GetPreStepPoint()->GetPosition().y()<<step->GetPreStepPoint()->GetPosition().z()<<G4endl;
// //       // G4cout<<" "<<step->GetPreStepPoint()->GetPosition().x() << " "<<step->GetPreStepPoint()->GetPosition().y() << " "<<step->GetPreStepPoint()->GetPosition().z()<<G4endl; 
// //       man->FillNtupleDColumn(0, G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID());
// //       man->FillNtupleDColumn(1, step->GetPreStepPoint()->GetPosition().x());
// //       man->FillNtupleDColumn(2, step->GetPreStepPoint()->GetPosition().y());
// //       man->FillNtupleDColumn(3, step->GetPreStepPoint()->GetPosition().z());
// //       man->FillNtupleDColumn(4, step->GetPreStepPoint()->GetProperTime());
// //       man->FillNtupleDColumn(5, step->GetTrack()->GetDefinition()->GetPDGEncoding());
// //       // man->FillNtupleDColumn(6,step->GetTrack()->GetMomentum().x());
// //       // man->FillNtupleDColumn(7,step->GetTrack()->GetMomentum().y());
// //       // man->FillNtupleDColumn(8,step->GetTrack()->GetMomentum().z());
// //       man->AddNtupleRow();
// //     // }
//   }

//   // // get volume of the current step
//   G4VPhysicalVolume *phys_vol = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
//   G4LogicalVolume *volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

//   G4cout<<"physical volume copy: "<<phys_vol->GetCopyNo()<<G4endl;

// }




//   // if (step->GetPreStepPoint() != NULL)
//   // {
//   //   if (step->GetPreStepPoint()->GetProcessDefinedStep() != NULL)
//   //   {
//   //     G4cout << "ProcName: " << step->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName() << G4endl;
//   //   }
//   // }

//   // if (volume->GetName() == "glass_sens")
//   // {
//   //   // G4cout << "The particle ID is " << step->GetTrack()->GetDefinition()->GetPDGEncoding() << G4endl;
//   //   if (step->GetTrack()->GetDefinition()->GetPDGEncoding() == 13)
//   //     SetMuPosition(step->GetTrack()->GetPosition(),step->GetTrack()->GetKineticEnergy(),step->GetTrack()->GetMomentum());//.x(),step->GetTrack()->GetMomentum().y());

//   //     // SetMuPosition(step->GetTrack()->GetPosition(),step->GetTrack()->GetKineticEnergy(),step->GetTrack()->GetMomentumDirection());//.x(),step->GetTrack()->GetMomentum().y());
//   //   else if (step->GetTrack()->GetDefinition()->GetPDGEncoding() == 11)
//   //     SetEPosition(step->GetTrack()->GetPosition(),step->GetTrack()->GetKineticEnergy(),step->GetTrack()->GetMomentum());//.x(),step->GetTrack()->GetMomentum().y());

//   //     // SetEPosition(step->GetTrack()->GetPosition(),step->GetTrack()->GetKineticEnergy(),step->GetTrack()->GetMomentumDirection());//.x(),step->GetTrack()->GetMomentum().y());
//   //   if(step->GetTrack()->GetDefinition()->GetPDGEncoding() != 13 && step->GetTrack()->GetDefinition()->GetPDGEncoding() != 11)
//   //     SetScatteringType(step->GetTrack()->GetDefinition()->GetPDGEncoding());
//   // }
//   // // // check if we are in scoring volume
//   // if (volume != fScoringVolume) return;

//   // // collect energy deposited in this step
//   // G4double edepStep = step->GetTotalEnergyDeposit();
//   // fEventAction->AddEdep(edepStep);
// }

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......















