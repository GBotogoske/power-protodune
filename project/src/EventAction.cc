#include "EventAction.hh"
#include "RunAction.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4AnalysisManager.hh"
#include "TrackingAction.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
: G4UserEventAction(),
  fHCID(-1)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
HitCollection*
EventAction::GetHitsCollection(G4int HCID, const G4Event* event) const
{
  auto hitsCollection = static_cast<HitCollection*>(event->GetHCofThisEvent()->GetHC(HCID));

  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << HCID;
    G4Exception("B4cEventAction::GetHitsCollection()",
      "MyCode0003", FatalException, msg);
  }
  return hitsCollection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::BeginOfEventAction(const G4Event* event)
{
    trackMap.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
   G4int nEvt = event->GetEventID();
    auto analysisManager = G4AnalysisManager::Instance();
    if ( fHCID == -1 ) 
    {
        fHCID = G4SDManager::GetSDMpointer()->GetCollectionID("HitCollection");
    }
    // ***Primary Hit Collection*** //
    auto HC = GetHitsCollection(fHCID, event);
    if(HC->entries()>0)
    {
            G4int PhotonDetected_UV = 0;
            G4int PhotonDetected_VIS = 0;
            double X=0.0;
            double Y=0.0;
            double Z=0.0;
            int iD= (*HC)[0]->GetDetectorID();

            PhotonDetected_UV = (*HC)[0]->GetPhotonCounter_UV();
            PhotonDetected_VIS = (*HC)[0]->GetPhotonCounter_VIS();
            X = (*HC)[0]->GetX();
            Y = (*HC)[0]->GetY();
            Z = (*HC)[0]->GetZ();

            int origin = (*HC)[0]->GetOrigin();
            
            analysisManager->FillNtupleIColumn(0,0,  nEvt);
            analysisManager->FillNtupleDColumn(0,1,  PhotonDetected_VIS);
            analysisManager->FillNtupleDColumn(0,2,  PhotonDetected_UV);
            analysisManager->FillNtupleDColumn(0,3,  X);
            analysisManager->FillNtupleDColumn(0,4,  Y);
            analysisManager->FillNtupleDColumn(0,5,  Z);
            analysisManager->FillNtupleDColumn(0,6,  iD);
            analysisManager->FillNtupleDColumn(0,7,  origin);
            analysisManager->AddNtupleRow(0);
    }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
