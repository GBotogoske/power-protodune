#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4VPhysicalVolume.hh"

#include <fstream>
#include <nlohmann/json.hpp>
using json = nlohmann::json;



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fGPS(0),fParticleGun(0)
{
    //fGPS  = new G4GeneralParticleSource();
    json config;
    fParticleGun = new G4ParticleGun(1);
    std::ifstream f("../configuration/detector.json");
    f >> config;
    auto config_cryostat = config["cryostat"];
    this->cryostat_sizeX = config_cryostat["size"][0].get<double>()*cm;
    this->cryostat_sizeY = config_cryostat["size"][1].get<double>()*cm;
    this->cryostat_sizeZ = config_cryostat["size"][2].get<double>()*cm;

    auto config_FC = config["FC"];
    
    this->FC_sizeY = config_FC["vertical_bar"][1].get<double>()*cm;

    this->Nin = 0;
    this->Nout = 0;
    this->Nout_edge = 0;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    //delete fGPS;
    delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    
    /*
    fGPS->SetParticlePosition(G4ThreeVector(0.,0.,0.));
    fGPS->GeneratePrimaryVertex(anEvent);
    */


    //Contas: 
    // Dentro: 14 x 6 voxels  --> 84 x 25000 = 2100000
    // Fora: 14 x 18 - 14 x 6.52 + 2x18 --> 196.72 x 40000 = 7869800

    G4double scale_factor_N = 1.0;

    G4double NinTotal = 2100000*scale_factor_N;
    G4double NoutTotal = 7869800*scale_factor_N;
    G4double NoutTotalEdges = 2*18*40000*scale_factor_N;
    G4double NoutTotalNonEdges =  NoutTotal - NoutTotalEdges;
    G4double NTotal = NinTotal + NoutTotal;
   
    G4ParticleDefinition* photon = G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton");

    float this_probability = 0.3; // quantity of argon light

    for(int i = 0; i < 1 ; i++) 
    {

        fParticleGun->SetParticleDefinition(photon);
        if(G4UniformRand()<=this_probability)
        {
            fParticleGun->SetParticleEnergy(9.69*eV);  //127nm  
        }
        else
        {
            fParticleGun->SetParticleEnergy(7.08*eV);  //175nm
        }

        G4double x_pos = 0*m;//(-0.25 + 0.5*G4UniformRand())*m;
        G4double y_pos = 0*m;// (-1+2*G4UniformRand())*cryostat_sizeY/2;
        G4double z_pos = 0*m;//(-1+2*G4UniformRand())*cryostat_sizeZ/2;
        G4String targetVolumeNameIn = "inside_argon";
        G4String targetVolumeNameIn2 = "Cathode_Hole_argon";
        G4String targetVolumeNameOut = "World_argon";
        G4ThreeVector pos;
        pos = G4ThreeVector(x_pos, y_pos, z_pos);
        G4VPhysicalVolume* volume = nullptr;
    
        bool redo = true;

        do 
        {
            // Sorteia posição dentro do criostato
            x_pos = (-0.25 + 0.5*G4UniformRand())*m;

            y_pos = (-1 + 2*G4UniformRand())*cryostat_sizeY/2;
            
            z_pos = (-1 + 2*G4UniformRand())*cryostat_sizeZ/2;

            pos = G4ThreeVector(x_pos, y_pos, z_pos);

            // Descobre qual volume contém essa posição
            volume = G4TransportationManager::GetTransportationManager()
                        ->GetNavigatorForTracking()
                        ->LocateGlobalPointAndSetup(pos);

            if (!volume)
                continue;

            std::string volumeName = volume->GetName();

            if(volumeName.find(targetVolumeNameIn) != std::string::npos || volumeName.find(targetVolumeNameIn2) != std::string::npos)
            {
                if(this->Nin < NinTotal) 
                {
                    this->Nin++;
                    redo = false;
                }
            }
            else if(volumeName.find(targetVolumeNameOut) != std::string::npos)
            {
                if(std::abs(y_pos) >= this->FC_sizeY/2)
                {
                    // EDGE
                    if(this->Nout_edge < NoutTotalEdges)
                    {
                        this->Nout++;
                        this->Nout_edge++;
                        redo = false;
                    }
                }
                else
                {
                    // NON-EDGE
                    G4double Nout_nonedge = this->Nout - this->Nout_edge;

                    if(Nout_nonedge < NoutTotalNonEdges)
                    {
                        this->Nout++;
                        redo = false;
                    }
                }
            }

        } while (redo);

        fParticleGun->SetParticlePosition(G4ThreeVector(x_pos,y_pos,z_pos));

        G4double theta = std::acos(1 - 2*G4UniformRand()); // de 0 a pi
        G4double phi   = 2*M_PI*G4UniformRand();           // de 0 a 2pi
        G4double x = std::sin(theta)*std::cos(phi);
        G4double y = std::sin(theta)*std::sin(phi);
        G4double z = std::cos(theta);
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(x,y,z));

        G4ThreeVector k = fParticleGun->GetParticleMomentumDirection();

        // escolhe um vetor qualquer não paralelo a k
        G4ThreeVector temp(1,0,0);
        if (std::abs(k.dot(temp)) > 0.99) temp = G4ThreeVector(0,1,0);

        // gera vetor perpendicular
        // gera base ortonormal: u, v no plano perpendicular a k
        G4ThreeVector u = (temp - k*(k.dot(temp))).unit();
        G4ThreeVector v = k.cross(u);

        // sorteia ângulo entre 0 e 2π
        G4double alpha = 2*M_PI*G4UniformRand();
        // combinação linear que continua ortogonal a k
        G4ThreeVector pol = std::cos(alpha)*u + std::sin(alpha)*v;

        fParticleGun->SetParticlePolarization(pol);

        fParticleGun->GeneratePrimaryVertex(anEvent);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
