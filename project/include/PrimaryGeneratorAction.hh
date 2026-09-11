#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

class G4GeneralParticleSource;
class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction();    
    virtual ~PrimaryGeneratorAction();

    // method from the base class
    virtual void GeneratePrimaries(G4Event*);         
  
    // method to access particle gun
    const G4GeneralParticleSource* GetGPS() const { return fGPS; }
  
  private:
    G4GeneralParticleSource*  fGPS;
    G4ParticleGun* fParticleGun;
    G4double cryostat_sizeX;
    G4double cryostat_sizeY;
    G4double cryostat_sizeZ;

    G4double FC_sizeY;

    G4int Nout;
    G4int Nout_edge;
    G4int Nin;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
