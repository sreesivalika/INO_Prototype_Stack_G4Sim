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
/// \file B1/src/DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"
#include <G4SubtractionSolid.hh>
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4GDMLParser.hh"



// DetectorConstruction::DetectorConstruction()
// : G4VUserDetectorConstruction(),
//   fScoringVolume(0)
// { }

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// DetectorConstruction::~DetectorConstruction()
// { }


namespace B1
{
DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // //Mylar
  // G4double density1 = 1.39*g/cm3;
  // G4int nel1 = 3;
  G4Element* elC = G4NistManager::Instance()->FindOrBuildElement("C");
  G4Element* elH = G4NistManager::Instance()->FindOrBuildElement("H");
  G4Element* elO = G4NistManager::Instance()->FindOrBuildElement("O");
  // G4Material* Mylar = new G4Material("Mylar", density1, nel1);
  // Mylar->AddElement(elC, 10);
  // Mylar->AddElement(elH, 8);
  // Mylar->AddElement(elO, 4);

  //R134a
  G4double density2 = 4.25*kg/m3;
  G4int nel2 = 3;
  G4Element* elF = G4NistManager::Instance()->FindOrBuildElement("F");

  // G4double fm1[3] = {0.2613, 0.0517, 0.6870};

  G4Material* R134a = new G4Material("R134a", density2, nel2);
  R134a->AddElement(elC, 2);
  R134a->AddElement(elH, 2);
  R134a->AddElement(elF, 4);
  // R134a->GetIonisation()->SetMeanExcitationEnergy(16.7*eV);

  // Isobutane
  G4double density3 = 2.51*kg/m3;
  G4int nel3 = 2;
  // G4Element* elC = G4NistManager::Instance()->FindOrBuildElement("C");
  // G4Element* elH = G4NistManager::Instance()->FindOrBuildElement("H");
  // G4double fm2[2] = {0.81, 0.19};
  G4Material* isobutane = new G4Material("Isobutane", density3, nel3);
  isobutane->AddElement(elC, 4);
  isobutane->AddElement(elH, 10);
  // isobutane->GetIonisation()->SetMeanExcitationEnergy(19.2*eV);

  //SF6
  G4double density4 = 6.164*kg/m3;
  G4int nel4 = 2;
  G4Element* elS = G4NistManager::Instance()->FindOrBuildElement("S");
  // G4double fm3[2] = {0.33, 0.67};
  G4Material* sf6 = new G4Material("SF6", density4, nel4);
  sf6->AddElement(elS, 1);
  sf6->AddElement(elF, 6);
  // sf6->GetIonisation()->SetMeanExcitationEnergy(19.2*eV);
  
  //Gas mix
  G4double d1 = 4.25*kg/m3;
  G4double d2 = 2.51*kg/m3;
  G4double d3 = 6.164*kg/m3;
  G4double fmass[3] = {0.9515, 0.0451, 0.0034}; //fraction masses
  G4double density = fmass[0]*d1 + fmass[1]*d2 + fmass[2]*d3;
  G4int nel = 3;
  G4Material* gasmix = new G4Material("gasmix", density, nel);
  gasmix->AddMaterial(G4Material::GetMaterial("R134a"), fmass[0]);
  gasmix->AddMaterial(G4Material::GetMaterial("Isobutane"), fmass[1]);
  gasmix->AddMaterial(G4Material::GetMaterial("SF6"), fmass[2]);



  G4int a = 178.0005;

  // Envelope parameters
  //
  G4double env_sizeXY = 2*m, env_sizeZ = 4*m;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_AIR");

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //
  // World
  //
  G4double world_sizeXY = 1.25*env_sizeXY;
  G4double world_sizeZ  = 1.25*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

  G4Box* solidWorld =    
    new G4Box("World",                       //its name
      //  0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*0.5*world_sizeZ);     //its size
        0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size

  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);       //overlaps checking
  
  // //
  // // Envelope
  // //
  // G4Box* solidEnv = new G4Box("Envelope",                    // its name
  //   0.5 * env_sizeXY, 0.5 * env_sizeXY, 0.5 * env_sizeZ);  // its size

  // G4LogicalVolume* logicEnv new G4LogicalVolume(solidEnv,  // its solid
  //   env_mat,                                     // its material
  //   "Envelope");                                 // its name

  // G4VPhysicalVolume* physEnv = new G4PVPlacement(0,  // no rotation
  //   G4ThreeVector(),          // at (0,0,0)
  //   0                 // its logical volume
  //   "Envelope",               // its name
  //   0,               // its mother  volume
  //   false,                    // no boolean operation
  //   0,                        // copy number
  //   checkOverlaps);           // overlaps checking

  


for(G4int i=0; i<=11; i++)
// G4int i=0;
{
  // Box shape
  
  //
  // GAS GAP
  //
  G4Material* GAS_gap_mat = nist->FindOrBuildMaterial("gasmix");

  G4Box* solidGAS_gap =
    new G4Box("GAS_gap",                       //its name
       0.4*0.5*world_sizeXY, 0.4*0.5*world_sizeXY, 0.0004*0.5*world_sizeZ);     //its size

  G4LogicalVolume* logicGAS_gap
   =
    new G4LogicalVolume(solidGAS_gap,         //its solid
                        GAS_gap_mat,          //its material
                        "GAS_gap");           //its name

  G4VPhysicalVolume* physGAS_gap = 
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0, 0, -(0+(i*a))*mm),                    //at position
                    logicGAS_gap,             //its logical volume
                    "GAS_gap",                //its name
                    logicWorld ,             //its mother  volume
                    false,                   //no boolean operation
                    i,                       //copy number
                    checkOverlaps);          //overlaps checking

  G4VisAttributes* Gas_gapVisAttributes = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0)); //cyan
  logicGAS_gap->SetVisAttributes(Gas_gapVisAttributes); 


  //
  // Glass Plates
  //

  // Plate-1:
  G4Material* glass_plate1_mat = nist->FindOrBuildMaterial("G4_GLASS_PLATE");

  // Box shape
  G4Box* solidglass_plate1 =
    new G4Box("glass_plate1",                       //its name
       0.4*0.5*world_sizeXY, 0.4*0.5*world_sizeXY, 0.0006*0.5*world_sizeZ);     //its size

  G4LogicalVolume* logicglass_plate1
   =
    new G4LogicalVolume(solidglass_plate1,         //its solid
                        glass_plate1_mat,          //its material
                        "glass_plate1");           //its name

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0, 0,  -(2.5+(a*i))*mm),                    //at position
                    logicglass_plate1,             //its logical volume
                    "glass_plate1",                //its name
                    logicWorld ,           //its mother  volume
                    false,                   //no boolean operation
                    i,                       //copy number
                    checkOverlaps);          //overlaps checking
  G4VisAttributes* glass_plate1VisAttributes = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));//blue
  logicglass_plate1->SetVisAttributes(glass_plate1VisAttributes);

  // Plate-2:
  G4Material* glass_plate2_mat = nist->FindOrBuildMaterial("G4_GLASS_PLATE");

  // Box shape
  G4Box* solidglass_plate2 =
    new G4Box("glass_plate2",                       //its name
       0.4*0.5*world_sizeXY, 0.4*0.5*world_sizeXY, 0.0006*0.5*world_sizeZ);     //its size

  G4LogicalVolume* logicglass_plate2
   =
    new G4LogicalVolume(solidglass_plate2,         //its solid
                        glass_plate2_mat,          //its material
                        "glass_plate2");           //its name

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0, 0,  -(-2.5+(a*i))*mm),                    //at position
                    logicglass_plate2,             //its logical volume
                    "glass_plate2",                //its name
                    logicWorld ,           //its mother  volume
                    false,                   //no boolean operation
                    i,                       //copy number
                    checkOverlaps);          //overlaps checking
  
  G4VisAttributes* glass_plate2VisAttributes = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));//blue
  logicglass_plate2->SetVisAttributes(glass_plate2VisAttributes);       
  //
  // Graphite Layers
  //

  // Layer-1:
  G4Material* graphite1_mat = nist->FindOrBuildMaterial("G4_GRAPHITE");

  // Box shape
  G4Box* solidgraphite1 =
    new G4Box("graphite1",                       //its name
       0.4*0.5*world_sizeXY, 0.4*0.5*world_sizeXY, 0.00000001*0.5*world_sizeZ);     //its size

  G4LogicalVolume* logicgraphite1
   =
    new G4LogicalVolume(solidgraphite1,         //its solid
                        graphite1_mat,          //its material
                        "graphite1");           //its name

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0, 0,  -(4.000025+(a*i))*mm),                    //at position
                    logicgraphite1,             //its logical volume
                    "graphite1",                //its name
                    logicWorld ,           //its mother  volume
                    false,                   //no boolean operation
                    i,                       //copy number
                    checkOverlaps);          //overlaps checking
  G4VisAttributes* graphite1VisAttributes = new G4VisAttributes(G4Colour(0.0, 0.0, 0.0));//black
  logicgraphite1->SetVisAttributes(graphite1VisAttributes);

  // Layer-2:
  G4Material* graphite2_mat = nist->FindOrBuildMaterial("G4_GRAPHITE");

  // Box shape
  G4Box* solidgraphite2 =
    new G4Box("graphite2",                       //its name
       0.4*0.5*world_sizeXY, 0.4*0.5*world_sizeXY, 0.00000001*0.5*world_sizeZ);     //its size

  G4LogicalVolume* logicgraphite2
   =
    new G4LogicalVolume(solidgraphite2,         //its solid
                        graphite2_mat,          //its material
                        "graphite2");           //its name

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0, 0,  -(-4.000025+(a*i))*mm),                    //at position
                    logicgraphite2,             //its logical volume
                    "graphite2",                //its name
                    logicWorld ,           //its mother  volume
                    false,                   //no boolean operation
                    i,                       //copy number
                    checkOverlaps);          //overlaps checking

  G4VisAttributes* graphite2VisAttributes = new G4VisAttributes(G4Colour(0.0, 0.0, 0.0));//black
  logicgraphite2->SetVisAttributes(graphite2VisAttributes);                  
  //
  // Mylar Layers
  //

  // Layer-1:
  G4Material* mylar1_mat = nist->FindOrBuildMaterial("G4_MYLAR");

  // Box shape
  G4Box* solidmylar1 =
    new G4Box("mylar1",                       //its name
       0.4*0.5*world_sizeXY, 0.4*0.5*world_sizeXY, 0.00000001*0.5*world_sizeZ);     //its size

  G4LogicalVolume* logicmylar1
   =
    new G4LogicalVolume(solidmylar1,         //its solid
                        mylar1_mat,          //its material
                        "mylar1");           //its name

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0, 0,  -(4.000075+(a*i))*mm),                    //at position
                    logicmylar1,             //its logical volume
                    "mylar1",                //its name
                    logicWorld ,           //its mother  volume
                    false,                   //no boolean operation
                    i,                       //copy number
                    checkOverlaps);          //overlaps checking
  G4VisAttributes* mylar1VisAttributes = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0));//magenta
  logicmylar1->SetVisAttributes(mylar1VisAttributes);   

  // Layer-2:
  G4Material* mylar2_mat = nist->FindOrBuildMaterial("G4_MYLAR");

  // Box shape
  G4Box* solidmylar2 =
    new G4Box("mylar2",                       //its name
       0.4*0.5*world_sizeXY, 0.4*0.5*world_sizeXY, 0.00000001*0.5*world_sizeZ);     //its size

  G4LogicalVolume* logicmylar2
   =
    new G4LogicalVolume(solidmylar2,         //its solid
                        mylar2_mat,          //its material
                        "mylar2");           //its name

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0, 0,  -(-4.000075+(a*i))*mm),                    //at position
                    logicmylar2,             //its logical volume
                    "mylar2",                //its name
                    logicWorld ,           //its mother  volume
                    false,                   //no boolean operation
                    i,                       //copy number
                    checkOverlaps);          //overlaps checking
  G4VisAttributes* mylar2VisAttributes = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0));//magenta
  logicmylar2->SetVisAttributes(mylar2VisAttributes);   

  //
  // Honeycomb Layers
  //

  // Layer-1:
  G4Material* hc1_mat = nist->FindOrBuildMaterial("G4_POLYCARBONATE");

  // Box shape
  G4Box* solidhc1 =
    new G4Box("hc1",                       //its name
       0.4*0.5*world_sizeXY, 0.4*0.5*world_sizeXY, 0.001*0.5*world_sizeZ);     //its size

  G4LogicalVolume* logichc1
   =
    new G4LogicalVolume(solidhc1,         //its solid
                        hc1_mat,          //its material
                        "hc1");           //its name

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0, 0,  -(6.50015+(a*i))*mm),                    //at position
                    logichc1,             //its logical volume
                    "hc1",                //its name
                    logicWorld ,           //its mother  volume
                    false,                   //no boolean operation
                    i,                       //copy number
                    checkOverlaps);          //overlaps checking
  G4VisAttributes* hc1VisAttributes = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));//yellow
  logichc1->SetVisAttributes(hc1VisAttributes);   

  // Layer-2:
  G4Material* hc2_mat = nist->FindOrBuildMaterial("G4_POLYCARBONATE");

  // Box shape
  G4Box* solidhc2 =
    new G4Box("hc2",                       //its name
       0.4*0.5*world_sizeXY, 0.4*0.5*world_sizeXY, 0.001*0.5*world_sizeZ);     //its size

  G4LogicalVolume* logichc2
   =
    new G4LogicalVolume(solidhc2,         //its solid
                        hc2_mat,          //its material
                        "hc2");           //its name

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0, 0,  -(-6.50015+(a*i))*mm),                    //at position
                    logichc2,             //its logical volume
                    "hc2",                //its name
                    logicWorld ,           //its mother  volume
                    false,                   //no boolean operation
                    i,                       //copy number
                    checkOverlaps);          //overlaps checking
  G4VisAttributes* hc2VisAttributes = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));//yellow
  logichc2->SetVisAttributes(hc2VisAttributes);   

  //
  // Al Layers
  //

  // Ground-1:
  G4Material* ground1_mat = nist->FindOrBuildMaterial("G4_Al");

  // Box shape
  G4Box* solidground1 =
    new G4Box("ground1",                       //its name
       0.4*0.5*world_sizeXY, 0.4*0.5*world_sizeXY, 0.00000001*0.5*world_sizeZ);     //its size

  G4LogicalVolume* logicground1
   =
    new G4LogicalVolume(solidground1,         //its solid
                        ground1_mat,          //its material
                        "ground1");           //its name

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0, 0,  -(9.0002+(a*i))*mm),                    //at position
                    logicground1,             //its logical volume
                    "ground1",                //its name
                    logicWorld ,           //its mother  volume
                    false,                   //no boolean operation
                    i,                       //copy number
                    checkOverlaps);          //overlaps checking
  G4VisAttributes* ground1VisAttributes = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));//grey
  logicground1->SetVisAttributes(ground1VisAttributes);   

  // Ground-2:
  G4Material* ground2_mat = nist->FindOrBuildMaterial("G4_Al");

  // Box shape
  G4Box* solidground2 =
    new G4Box("ground2",                       //its name
       0.4*0.5*world_sizeXY, 0.4*0.5*world_sizeXY, 0.00000001*0.5*world_sizeZ);     //its size

  G4LogicalVolume* logicground2
   =
    new G4LogicalVolume(solidground2,         //its solid
                        ground2_mat,          //its material
                        "ground2");           //its name

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0, 0,  -(-9.0002+(a*i))*mm),                    //at position
                    logicground2,             //its logical volume
                    "ground2",                //its name
                    logicWorld ,           //its mother  volume
                    false,                   //no boolean operation
                    i,                       //copy number
                    checkOverlaps);          //overlaps checking
  G4VisAttributes* ground2VisAttributes = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));//grey
  logicground2->SetVisAttributes(ground2VisAttributes);   


  //
  // Cu Pickupx Panel
  // Xpick
  G4Material* pickupx_mat = nist->FindOrBuildMaterial("G4_Cu");

  // Box shape
  G4Box* solidpickupx =
    new G4Box("pickupx",                       //its name
       0.4*0.5*world_sizeXY, 0.4*0.5*world_sizeXY, 0.00000001*0.5*world_sizeZ);     //its size

  G4LogicalVolume* logicpickupx
   =
    new G4LogicalVolume(solidpickupx,         //its solid
                        pickupx_mat,          //its material
                        "pickupx");           //its name

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0, 0,  -(-4.000125+(a*i))*mm),                    //at position
                    logicpickupx,             //its logical volume
                    "pickupx",                //its name
                    logicWorld ,           //its mother  volume
                    false,                   //no boolean operation
                    i,                       //copy number
                    checkOverlaps);          //overlaps checking
  G4VisAttributes* pickupxVisAttributes = new G4VisAttributes(G4Colour(0.647, 0.165, 0.165));//brown
  logicpickupx->SetVisAttributes(pickupxVisAttributes); 
  // pickupy
  G4Material* pickupy_mat = nist->FindOrBuildMaterial("G4_Cu");

  // Box shape
  G4Box* solidpickupy =
    new G4Box("pickupy",                       //its name
       0.4*0.5*world_sizeXY, 0.4*0.5*world_sizeXY, 0.00000001*0.5*world_sizeZ);     //its size

  G4LogicalVolume* logicpickupy
   =
    new G4LogicalVolume(solidpickupy,         //its solid
                        pickupy_mat,          //its material
                        "pickupy");           //its name

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0, 0,  -(4.000125+(a*i))*mm),                    //at position
                    logicpickupy,             //its logical volume
                    "pickupy",                //its name
                    logicWorld ,           //its mother  volume
                    false,                   //no boolean operation
                    i,                       //copy number
                    checkOverlaps);          //overlaps checking
  G4VisAttributes* pickupyVisAttributes = new G4VisAttributes(G4Colour(0.647, 0.165, 0.165));//brown
  logicpickupy->SetVisAttributes(pickupyVisAttributes); 

}
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
