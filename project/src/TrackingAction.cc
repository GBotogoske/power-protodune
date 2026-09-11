#include "TrackingAction.hh"
#include "G4EventManager.hh"

G4ThreadLocal std::map<G4int, TrackInfo> trackMap;

void MyTrackingAction::PreUserTrackingAction(const G4Track* track) 
{
    G4int origin = -1;

    // Volume onde esse track nasceu
    if(track->GetVolume())
    {
        G4String volumeName = track->GetVolume()->GetName();

        if(volumeName.find("inside_argon") != std::string::npos ||
        volumeName.find("Cathode_Hole_argon") != std::string::npos)
        {
            origin = 1;
        }
        else if(volumeName.find("World_argon") != std::string::npos)
        {
            origin = 0;
        }
    }

    trackMap[track->GetTrackID()] = {
        track->GetVertexPosition(),
        track->GetParentID(),
        origin
    };
}