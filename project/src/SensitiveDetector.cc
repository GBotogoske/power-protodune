#include "SensitiveDetector.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4TouchableHistoryHandle.hh"
#include "G4TouchableHistory.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4Tubs.hh"
#include "G4String.hh"
#include "G4OpticalPhoton.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4GeometryTolerance.hh"
#include "TrackingAction.hh"
#include <math.h>

using namespace std;


static const G4double GeometricalTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
static const G4double h = CLHEP::h_Planck/((CLHEP::eV)*(CLHEP::ns));
static const G4double c = CLHEP::c_light/((CLHEP::nm)/(CLHEP::ns));

SensitiveDetector::SensitiveDetector(const G4String &SDname,const G4String &HitCollectionName)
  : G4VSensitiveDetector(SDname),
    fHitCollection(NULL),
    fP_uv(0),
    fP_vis(0)
{
  G4cout<<"Creating SD with name: "<<SDname<<G4endl;
  collectionName.insert(HitCollectionName);
}

SensitiveDetector::~SensitiveDetector(){}

void SensitiveDetector::Initialize(G4HCofThisEvent *HCE)
{
    fHitCollection = new HitCollection(GetName(),collectionName[0]);
    G4int HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); //<<-- this is to get an ID for the colletionName[0]
    //G4cout<<"*** "<<fHitCollection->GetName()<<" initialized [ID = "<<HCID<<"]"<<G4endl;
    HCE->AddHitsCollection(HCID, fHitCollection);
}


G4bool SensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *)
{
    G4double eDep   = aStep->GetTotalEnergyDeposit();
    G4int counter = 0;

    G4String thisParticle = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();
    G4String thisVolume   = aStep->GetTrack()->GetVolume()->GetName();

    //G4Track* track = aStep->GetTrack();
    //G4TrackStatus status = track->GetTrackStatus();

    G4int pDetected = 0;
    bool isVIS=false;
    bool isUV=false;
    G4double p;

    if(thisParticle=="opticalphoton")
    {
        G4String procName = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
        G4bool checkAbsorption = G4StrUtil::contains(procName,"Absorption");
        //G4cout << "Process:" << procName << std::endl;

        const auto& sipm_spectrum = SiPMSpectrum::get();
        std::vector<G4double> eff;
        std::vector<G4double> E;
        double reflectance = sipm_spectrum.getReflectivy();
        int n;
        if(checkAbsorption==true and G4StrUtil::contains(thisVolume,"sipmVIS"))
        {   
            //G4cout << "VIS" << std::endl;
            isVIS=true;
            eff = sipm_spectrum.get_effVIS();
            E = sipm_spectrum.get_EVIS();
            n = sipm_spectrum.getNVIS();
        }
        else if(checkAbsorption==true and G4StrUtil::contains(thisVolume,"sipmUV"))
        {   
            //G4cout << "UV" << std::endl;
            isUV=true;
            eff = sipm_spectrum.get_effUV();
            E = sipm_spectrum.get_EUV();
            n = sipm_spectrum.getNUV();
        }
        else
        {
            return true;
        }
       
        auto Ephoton = aStep->GetTrack()->GetTotalEnergy();
        if(Ephoton<=E[0])
        {
            p=eff[0];
        }
        else if(Ephoton>=E[n-1])
        {
            p=eff[n-1]; 
        }
        else
        {
            for (int i = 0; i < n - 1; ++i) 
            {
                if (Ephoton >= E[i] && Ephoton <= E[i + 1])
                {
                    G4double x0 = E[i];
                    G4double x1 = E[i + 1];
                    G4double y0 = eff[i];
                    G4double y1 = eff[i + 1];
                    p = y0 + (y1 - y0) * (Ephoton - x0) / (x1 - x0);
                }
            }
        }
        G4double r = G4UniformRand();  

        p = p/(1-reflectance);

        if(r<p)
        {
            pDetected +=1;
        }
    }
    if(pDetected>=1)
    {
        G4Track* track = aStep->GetTrack();
        G4int ancestorID = track->GetTrackID();  

        // sobe até a primária
        while(trackMap[ancestorID].parentID != 0) 
        {
            ancestorID = trackMap[ancestorID].parentID;
        }

        G4ThreeVector motherPos = trackMap[ancestorID].vertex;

        double X = motherPos.getX()/m;
        double Y = motherPos.getY()/m;
        double Z = motherPos.getZ()/m;
        
        //std::cout << X << " " << Y <<  " " << Z  << std::endl;

        int idDetector = aStep->GetTrack()->GetVolume()->GetCopyNo();

        if(isVIS)
        {
            OneHit *aHit = new OneHit();
            SensitiveDetector::SetCounterStatus_VIS(pDetected);
            G4int TotP_VIS = SensitiveDetector::GetCounterStatus_VIS();
            aHit->SetPhotonCounter_VIS(1);
            aHit->SetPhotonCounter_UV(0);
            aHit->SetX(X);
            aHit->SetY(Y);
            aHit->SetZ(Z);
            aHit->SetDetectorID(idDetector);
            fHitCollection->insert(aHit);
            SensitiveDetector::PrintSDMemoryStatus();
            SensitiveDetector::CleanSDMemory();

        }
        if(isUV)
        {
            OneHit *aHit = new OneHit();
            SensitiveDetector::SetCounterStatus_UV(pDetected);
            G4int TotP_UV = SensitiveDetector::GetCounterStatus_UV();
            aHit->SetPhotonCounter_VIS(0);
            aHit->SetPhotonCounter_UV(1);
            aHit->SetX(X);
            aHit->SetY(Y);
            aHit->SetZ(Z);
            aHit->SetDetectorID(idDetector);
            fHitCollection->insert(aHit);
            SensitiveDetector::PrintSDMemoryStatus();
            SensitiveDetector::CleanSDMemory();
        }
    }
    
    
    return true;
}

void SensitiveDetector::EndOfEvent(G4HCofThisEvent*)
{
   /*  //Fill the hits
    OneHit *aHit = new OneHit();
    G4int TotP_VIS = SensitiveDetector::GetCounterStatus_VIS();
    G4int TotP_UV = SensitiveDetector::GetCounterStatus_UV();

    aHit->SetPhotonCounter_VIS(TotP_VIS);
    aHit->SetPhotonCounter_UV(TotP_UV);

    fHitCollection->insert(aHit);

    SensitiveDetector::PrintSDMemoryStatus();
    SensitiveDetector::CleanSDMemory(); */
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
G4bool SensitiveDetector::IsAnOpticalPhoton(G4Step* aStep)
{
    G4bool flag = false;
    G4Track *aTrack = aStep->GetTrack();
    const G4ParticleDefinition* aDef = aTrack->GetParticleDefinition();
    G4String aName = aDef->GetParticleName();
    if(aName=="opticalphoton")
    {
        flag=true;
    }
    return flag;
}

void SensitiveDetector::CleanSDMemory()
{
    SensitiveDetector::ResetCounterStatus_UV();
    SensitiveDetector::ResetCounterStatus_VIS();
}

void SensitiveDetector::PrintSDMemoryStatus()
{
    
/*     G4int    TotP = SensitiveDetector::GetCounterStatus();
    G4cout<<"O.o.O.o.O.o.O.o.O.o.O.o.O.o.O.o.O.o.O.o.O"<<G4endl;
    G4cout<<" Total Energy deposited = "<<G4BestUnit(TotE,"Energy")<<G4endl;
    G4cout<<" Total Photon detected  = "<<TotP<<G4endl;
    G4cout<<"O.o.O.o.O.o.O.o.O.o.O.o.O.o.O.o.O.o.O.o.O"<<G4endl; */
}
