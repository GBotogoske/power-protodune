#include "SteppingAction.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4LogicalVolume.hh"
#include "G4EventManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"
#include "RunAction.hh"
#include "G4AnalysisManager.hh"

SteppingAction::SteppingAction() {}
SteppingAction::~SteppingAction() {}

G4ThreadLocal int cont_PEN=0;

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    auto track = step->GetTrack();
    G4String particle = track->GetDefinition()->GetParticleName();

    // PreStep volume info
    G4StepPoint *aPrePoint = step->GetPreStepPoint();
    G4VPhysicalVolume *aPrePV = aPrePoint->GetPhysicalVolume();
    G4String PreVolName = "";
    if (aPrePV)
        PreVolName = aPrePV->GetName();

    // PostStep volume info
    G4StepPoint *aPostPoint = step->GetPostStepPoint();
    G4VPhysicalVolume *aPostPV = aPostPoint->GetPhysicalVolume();
    G4String PostVolName = "";
    if (aPostPV)
        PostVolName = aPostPV->GetName();

    if (particle == "opticalphoton")
    {
        if(G4StrUtil::contains(PreVolName,"argon") && G4StrUtil::contains(PostVolName,"PEN"))
        {
            auto Energy = aPostPoint->GetTotalEnergy();
            if(Energy>9*eV && Energy<12*eV)
            {
                auto runAction = const_cast<RunAction*>(static_cast<const RunAction*>(G4RunManager::GetRunManager()->GetUserRunAction()));
                runAction->AddPenCount(1);
                //G4cout << "thaaank you meestre" << G4endl;
            }
        }
        if((G4StrUtil::contains(PreVolName,"Vikuiti")  && G4StrUtil::contains(PostVolName,"argon")) || (G4StrUtil::contains(PreVolName,"Reflector")  && G4StrUtil::contains(PostVolName,"PEN")))
        {
            // G4cout << "oba" << G4endl;
        }
    }

    /* auto track = step->GetTrack();
    auto pos = track->GetPosition();
    auto preVolume = step->GetPreStepPoint()->GetPhysicalVolume(); // <-- seguro
    auto eventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
    auto Energy = step->GetPreStepPoint()->GetTotalEnergy();

    if(preVolume) {
        G4cout << "EventID: " << eventID
               << " TrackID: " << track->GetTrackID()
               << " Volume: " << preVolume->GetName()
               << " Position: " << pos 
               << " Energy: " << Energy << G4endl;
    } else {
        G4cout << "EventID: " << eventID
               << " TrackID: " << track->GetTrackID()
               << " Volume: [NULL - saiu do mundo]"
               << " Position: " << pos << G4endl;
    } */

    
    // This section is responsible for collecting PEN emission and absorbption data

    G4bool get_pen_data = true; // set to false if you don't want to collect PEN data (this will result in much smaller .root files)
    G4AnalysisManager *man = G4AnalysisManager::Instance(); // is a singleton btw
    if(get_pen_data){
        G4bool wls_absorbption = false;
        G4bool wls_emission = false;
        

        if (particle == "opticalphoton"){
            const G4VProcess* post_step_process = step->GetPostStepPoint()->GetProcessDefinedStep();
            if(track->GetTrackStatus()==fStopAndKill && post_step_process){
                if(post_step_process->GetProcessName() == "OpWLS"){
                    G4double absorbed_pen_energy = track->GetKineticEnergy();
                    G4double absorbed_pen_wavelenght = (CLHEP::h_Planck * CLHEP::c_light)/absorbed_pen_energy;
                    G4double absorbed_pen_wavelenght_nm = absorbed_pen_wavelenght/nm;

                    man->FillNtupleDColumn(1,1,absorbed_pen_wavelenght_nm);
                    //G4cout << "absorcao wls pen" << G4endl;
                    wls_absorbption = true;
                }
            }
        }
        

        const std::vector<const G4Track*> *secondaries = step->GetSecondaryInCurrentStep();
        if (secondaries && secondaries->size()>0){
            for(auto secondarie_track: *secondaries){
                G4String secondarie_particle_name = secondarie_track->GetDefinition()->GetParticleName();
                if(secondarie_particle_name == "opticalphoton"){
                    const G4VProcess* creator_process = secondarie_track->GetCreatorProcess();
                    if(creator_process && creator_process->GetProcessName() == "OpWLS"){
                        G4double emmited_pen_energy = secondarie_track->GetKineticEnergy();
                        G4double emmited_pen_wavelenght = (CLHEP::h_Planck * CLHEP::c_light)/emmited_pen_energy;
                        G4double emmited_pen_wavelenght_nm = emmited_pen_wavelenght/nm;

                        man->FillNtupleDColumn(1,0,emmited_pen_wavelenght_nm);
                        wls_emission = true;
                        //G4cout << "emissao wls pen" << G4endl;
                    }
                }
            }
        }
        if(wls_emission && wls_absorbption){
            man->AddNtupleRow(1);
        }
    }
}
