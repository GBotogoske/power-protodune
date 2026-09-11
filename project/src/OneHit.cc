//GammaBari G4Week - By Corrado Altomare & Davide Serini :D
//Ex2-Physics By Corrado
//

//Include Custom Classes
#include "OneHit.hh"

//Include Native G4 Classes
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"

#include <iomanip>

OneHit::OneHit():fPhotonCounter_UV(0),fPhotonCounter_VIS(0),fX(0.0),fY(0.0),fZ(0.0),fDetectorID(0), fOrigin(-1) {}

OneHit::~OneHit() {}
