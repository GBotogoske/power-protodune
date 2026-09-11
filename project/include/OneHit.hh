//GammaBari G4Week - By Corrado Altomare & Davide Serini :D
//Ex2-Physics By Corrado
//

#ifndef OneHit_h
#define OneHit_h 1

//Include Native G4 Classes
#include "G4VHit.hh"
#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4UnitsTable.hh"
#include "G4Track.hh"

class OneHit : public G4VHit{
  public:

    OneHit();
    virtual ~OneHit();

    // methods from base class

    void SetPhotonCounter_VIS (G4int p)                       {fPhotonCounter_VIS = p;};
    G4int GetPhotonCounter_VIS()                   const { return fPhotonCounter_VIS;};

    void SetPhotonCounter_UV (G4int p)                       {fPhotonCounter_UV = p;};
    G4int GetPhotonCounter_UV()                   const { return fPhotonCounter_UV;};

    void SetDetectorID(G4int id) {fDetectorID = id;};
    G4int GetDetectorID() const { return fDetectorID;};

    void SetX (G4double x)                       {fX=x;};
    G4double GetX()                   const { return fX;};
    void SetY (G4double x)                       {fY=x;};
    G4double GetY()                   const { return fY;};
    void SetZ (G4double x)                       {fZ=x;};
    G4double GetZ()                   const { return fZ;};

    void SetOrigin(G4int origin) {fOrigin = origin;};
    G4int GetOrigin() const {return fOrigin;};

  private:
 
    G4int         fPhotonCounter_VIS;
    G4int         fPhotonCounter_UV;
    G4double fX;
    G4double fY;
    G4double fZ;
    G4int fDetectorID;
    G4int fOrigin;

  };

typedef G4THitsCollection<OneHit> HitCollection;


#endif
