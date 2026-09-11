#include "G4UserTrackingAction.hh"
#include "G4Track.hh"
#include <map>
#include "G4ThreeVector.hh"
#include "globals.hh"

struct TrackInfo
{
    G4ThreeVector vertex;
    G4int parentID;
    G4int volOrigin; //1 inside, 0 outside
};

G4ThreadLocal extern std::map<G4int, TrackInfo> trackMap;

class MyTrackingAction : public G4UserTrackingAction 
{
public:
    virtual void PreUserTrackingAction(const G4Track* track) override;
};
