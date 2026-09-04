#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "RunActionMessenger.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4GeneralParticleSource.hh"
#include "G4AnalysisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(): G4UserRunAction(), fOutputFileName("./Data") //, fPenCount(0)
{
    //G4AccumulableManager::Instance()->RegisterAccumulable(fPenCount);
    fRunMessenger = new RunActionMessenger(this);
    G4RunManager::GetRunManager()->SetPrintProgress(100000);
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    G4cout << "Using " << analysisManager->GetType() << G4endl;
    analysisManager->SetNtupleMerging(true);
    //Create ntuple
    analysisManager->CreateNtuple("Primary_Hit", "Primary_Hit");
    analysisManager->CreateNtupleIColumn(0,"eventID");
    analysisManager->CreateNtupleDColumn(0,"PhotonDetectedVIS");
    analysisManager->CreateNtupleDColumn(0,"PhotonDetectedUV");
    analysisManager->CreateNtupleDColumn(0,"X");
    analysisManager->CreateNtupleDColumn(0,"Y");
    analysisManager->CreateNtupleDColumn(0,"Z");
    analysisManager->CreateNtupleDColumn(0,"detectorID");

    analysisManager->FinishNtuple(0);

    analysisManager->CreateNtuple("Optical_Properties","Optical_Properties");
    analysisManager->CreateNtupleDColumn(1,"emmited_wabelenght_PEN");
    analysisManager->CreateNtupleDColumn(1,"absorbed_wavelenght_PEN");
    // ADD EVENT ID

    analysisManager->FinishNtuple(1);


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
    delete G4AnalysisManager::Instance();
    delete fRunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{
    G4String fileName   = fOutputFileName;
    G4String fileOutput = fileName;
    fileOutput += ".root";
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->OpenFile(fileOutput);
    //G4AccumulableManager::Instance()->Reset();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{
    auto analysisManager = G4AnalysisManager::Instance();
    // save histograms & ntuple
    G4cout<<"[INFO]: Write the output file "<<fOutputFileName<<".root"<<G4endl;
    analysisManager->Write();
    analysisManager->CloseFile();

 /*    G4AccumulableManager::Instance()->Merge();

    if (IsMaster())
    {
        G4cout << "[INFO]: Write the output file ./Data.root" << G4endl;
        G4cout << "Total de optical photons detectados no PEN: "
               << fPenCount.GetValue() << G4endl;
    } */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
