#include "DetectorConstruction.hh"
#include "OneHit.hh"
#include "OptFileManager.hh"

#include "SensitiveDetector.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"  

#include "G4Material.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4EllipticalTube.hh"
#include "G4ExtrudedSolid.hh"
#include "G4AssemblyVolume.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4NistManager.hh"
#include <math.h>

//#include "G4GDMLParser.hh"

//nota de coordadenas (x,y,z) --> (direcao da tela, direcao do campo, direcao transversal da tela)

// Option to switch on/off checking of volumes overlaps
G4bool checkOverlaps = true;

DetectorConstruction::DetectorConstruction():G4VUserDetectorConstruction()
{
    std::ifstream f("../configuration/detector.json");
    f >> this->config;
}

DetectorConstruction::~DetectorConstruction(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* DetectorConstruction::Construct()
{
    G4NistManager *fNistManager = G4NistManager::Instance();
    fNistManager->SetVerbose(0);
  
    // ----  World -----

    auto config_world = config["world"];
    G4double world_sizeX = config_world["size"][0].get<double>()*cm;
    G4double world_sizeY = config_world["size"][1].get<double>()*cm;
    G4double world_sizeZ = config_world["size"][2].get<double>()*cm;

    G4cout << "LAr World size: " << world_sizeX/m << " .m  X " << world_sizeY/m << " .m  X " << world_sizeZ/m << " .m " << std::endl;

    // -----  Liquid Argon Filling --------
    auto config_LAr = config["liquid_argon"];
    auto config_Reflector = config["Reflector"];
    auto config_FC = config["FC"];

    auto config_cryostat = config["cryostat"];
    G4Material* lar_mat = fNistManager->FindOrBuildMaterial(config_cryostat["filling"].get<string>());
    std::string file_lar_index = config_LAr["refraction_index"].get<string>();
    std::ifstream f_lar_index(std::string("../configuration/")+file_lar_index);
    json data = json::parse(f_lar_index);
    size_t n_rindex = data.size();
    std::vector<G4double> energies_r(n_rindex);
    std::vector<G4double> rindex(n_rindex);

    G4cout << "Liquid Argon properties: " << std::endl;
    G4cout << "Refraction Index: " << std::endl;
    for (size_t i = 0; i < n_rindex; i++) 
    {
        energies_r[i] = data[i]["E"].get<double>()*eV;  
        rindex[i]   = data[i]["r"].get<double>();
        G4cout << energies_r[i]/eV << " eV -- " <<   rindex[i] << std::endl;
    }
       
    std::string file_lar_abslength = config_LAr["abs_length"].get<string>();
    std::ifstream f_lar_file_lar_abslength(std::string("../configuration/")+file_lar_abslength);
    json abs_lar = json::parse(f_lar_file_lar_abslength);
    size_t n_abs = abs_lar.size();
    std::vector<G4double> absE_LAr(n_abs);
    std::vector<G4double> absLen_LAr(n_abs);
    G4cout << "Absorption Length: " << std::endl;
    for (size_t i=0; i<n_abs; i++) 
    {
        absE_LAr[i] = abs_lar[i]["E"].get<double>()*eV;    
        absLen_LAr[i] = abs_lar[i]["l"].get<double>()*cm; 
        G4cout <<   absE_LAr[i]/eV << " eV -- " << absLen_LAr[i]/m << " m"<< std::endl;
    }
    
    auto scint_lar = config_LAr["scintilation"];
    std::string file_lar_scint = scint_lar["scintilation_point"].get<string>();
    std::ifstream f_lar_file_lar_scint(std::string("../configuration/")+file_lar_scint);
    json scint_lar_file = json::parse(f_lar_file_lar_scint);
    size_t n_scint = scint_lar_file.size();    
    std::vector<G4double> scintE_LAr(n_scint);
    std::vector<G4double> scint(n_scint);
    G4cout << "Scintilation Spectrum: " << std::endl;
    for (size_t i=0; i<n_scint; i++) 
    {
        scintE_LAr[i] = scint_lar_file[i]["E"].get<double>()*eV;    
        scint[i] = scint_lar_file[i]["s"].get<double>(); 
        G4cout <<   scintE_LAr[i]/eV << " eV -- " << scint[i] << std::endl;  
    }

    std::string file_lar_rayleigh = config_LAr["rayleigh_scattering"].get<string>();
    std::ifstream f_lar_rayleigh(std::string("../configuration/")+file_lar_rayleigh);
    json data_ray = json::parse(f_lar_rayleigh);
    size_t n_rayleigh = data_ray.size();
    std::vector<G4double> energies_ray(n_rayleigh);
    std::vector<G4double> ray_lar(n_rayleigh);
    G4cout << "Rayleigh Scattering: " << std::endl;
    for (size_t i = 0; i < n_rayleigh; i++) 
    {
        energies_ray[i] = data_ray[i]["E"].get<double>()*eV;  
        ray_lar[i]   = data_ray[i]["l"].get<double>()*cm;
        G4cout <<   energies_ray[i] << " eV -- " << ray_lar[i]/m << " m"<< std::endl;
    }

    G4cout << "LightYield - 0 Field: " << scint_lar["LY_without_field"].get<double>()/MeV << " photons/MeV" << std::endl;
    G4cout << "LightYield - Nominal: " << scint_lar["LY"].get<double>()/MeV << " photons/MeV" << std::endl;
    G4cout << "Slow Component - 0 Field: " << scint_lar["Slow_Comp_without_field"].get<double>()<< " ns" << std::endl;
    G4cout << "Slow Component - Nominal: " << scint_lar["Slow_Comp"].get<double>()<< " ns" << std::endl;
    G4cout << "Slow Amplitude: " << scint_lar["Slow_Amp"].get<double>() << std::endl;
    G4cout << "Fast Amplitude: " << scint_lar["Fast_Amp"].get<double>() << std::endl;
    G4cout << "Resolution: " << scint_lar["Resolution"].get<double>() << std::endl;
    
    auto mpt_LAr = new G4MaterialPropertiesTable();
    mpt_LAr->AddProperty("RINDEX", energies_r, rindex, n_rindex);
    mpt_LAr->AddProperty("ABSLENGTH", absE_LAr, absLen_LAr, n_abs);
    mpt_LAr->AddProperty("SCINTILLATIONCOMPONENT1",scintE_LAr,scint,n_scint);
    mpt_LAr->AddProperty("SCINTILLATIONCOMPONENT2",scintE_LAr,scint,n_scint);
    mpt_LAr->AddConstProperty("SCINTILLATIONYIELD",scint_lar["LY_without_field"].get<double>()/MeV);
    mpt_LAr->AddConstProperty("SCINTILLATIONTIMECONSTANT1",scint_lar["Slow_Comp_without_field"].get<double>()*ns);
    mpt_LAr->AddConstProperty("SCINTILLATIONTIMECONSTANT2",scint_lar["Fast_Comp"].get<double>()*ns);
    mpt_LAr->AddConstProperty("SCINTILLATIONYIELD1", scint_lar["Slow_Amp"].get<double>());
    mpt_LAr->AddConstProperty("SCINTILLATIONYIELD2", scint_lar["Fast_Amp"].get<double>());
    mpt_LAr->AddConstProperty("RESOLUTIONSCALE",scint_lar["Resolution"].get<double>());
    mpt_LAr->AddProperty("RAYLEIGH", energies_ray, ray_lar, n_rayleigh);
 
    //lar_mat->GetIonisation()->SetBirksConstant(scint_lar["Birks"].get<double>()*cm/MeV); // se comentar essa linha, o valor eh o default do argonio liquido
    lar_mat->SetMaterialPropertiesTable(mpt_LAr);

    G4Box* solidWorld = new G4Box("World_argon",0.5*world_sizeX, 0.5*world_sizeY, 0.5*world_sizeZ);
    G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld,lar_mat,"World_argon");
    G4VPhysicalVolume* physicalWorld = new G4PVPlacement(0,G4ThreeVector(),logicWorld,"World_argon",0,false,0,checkOverlaps);

    // ----- Cryostat ----
    G4double cryostat_sizeX = config_cryostat["size"][0].get<double>()*cm;
    G4double cryostat_sizeY = config_cryostat["size"][1].get<double>()*cm;
    G4double cryostat_sizeZ = config_cryostat["size"][2].get<double>()*cm;
    G4double cryostatThickness = config_cryostat["thickness"].get<double>()*cm;

    G4cout << "Cryostat Internal size: " << cryostat_sizeX/m << " m  X " << cryostat_sizeY/m << " m  X " << cryostat_sizeZ/m << " m " << std::endl;
    G4cout << "Cryostat ThicknesG4LogicalVolumes: " << cryostatThickness/m << " m " << std::endl;
    G4cout << "Cryostat Refraction Index: " << config_cryostat["refraction_index"].get<double>() << std::endl;
    G4cout << "Cryostat Absorption Length: " << config_cryostat["abs_length"].get<double>()*cm/m << "m" << std::endl;

    G4Material* cryostat_mat = fNistManager->FindOrBuildMaterial(config_cryostat["material"].get<string>());
    auto mpt_Cryo = new G4MaterialPropertiesTable();
    mpt_Cryo->AddProperty("RINDEX", {1*eV,15*eV}, {config_cryostat["refraction_index"].get<double>(),config_cryostat["refraction_index"].get<double>()}, 2);
    mpt_Cryo->AddProperty("ABSLENGTH", {1*eV,15*eV}, {config_cryostat["abs_length"].get<double>()*cm,config_cryostat["abs_length"].get<double>()*cm}, 2);
    cryostat_mat->SetMaterialPropertiesTable(mpt_Cryo);

    G4Box* solidExternalCryostat = new G4Box("Cryostat_ext",0.5*cryostat_sizeX+cryostatThickness, 0.5*cryostat_sizeY+cryostatThickness, 0.5*cryostat_sizeZ+cryostatThickness);
    G4Box* solidInternalCryostat = new G4Box("Cryostat_int",0.5*cryostat_sizeX, 0.5*cryostat_sizeY, 0.5*cryostat_sizeZ);
    auto solidCryostat = new G4SubtractionSolid("Cryostat", solidExternalCryostat,solidInternalCryostat,nullptr,G4ThreeVector());
    
    G4LogicalVolume* logicCryostat = new G4LogicalVolume(solidCryostat,cryostat_mat,"Cryostat");
    G4VPhysicalVolume* physicalCryostat = new G4PVPlacement(0,G4ThreeVector(),logicCryostat,"Cryostat",logicWorld,false,0,checkOverlaps);
    
    // -----  Liquid Argon && Cryostat Interface Boundary --------

    G4OpticalSurface* surface_cryo_lar = new G4OpticalSurface("surface_cryo_lar");
    surface_cryo_lar->SetModel(unified);
    surface_cryo_lar->SetType(dielectric_metal);
    surface_cryo_lar->SetFinish(polished);
    G4MaterialPropertiesTable* mpt_CryoLar_Surface = new G4MaterialPropertiesTable();

    mpt_CryoLar_Surface->AddProperty("REFLECTIVITY",{0.1*eV,15*eV}, {config_FC["reflectivity"].get<double>(),config_FC["reflectivity"].get<double>()},2);
    surface_cryo_lar->SetMaterialPropertiesTable(mpt_CryoLar_Surface);
    new G4LogicalBorderSurface("LiquidArgon-->Cryo", physicalWorld, physicalCryostat , surface_cryo_lar );

    // --------------  ANODE ----------------------------
    G4double AnodeTickness = 0.000000001*cm;
    G4Box* solidAnodeTemplate = new G4Box("Anode",0.5*cryostat_sizeX, 0.5*AnodeTickness, 0.5*cryostat_sizeZ);
    G4LogicalVolume* logicAnode = new G4LogicalVolume(solidAnodeTemplate,cryostat_mat,"Anode");
    G4VPhysicalVolume* physicalTopAnode = new G4PVPlacement(0,G4ThreeVector(0,0.5*cryostat_sizeY-0.5*AnodeTickness,0),logicAnode,"AnodeUp",logicWorld,false,0,checkOverlaps);
    G4VPhysicalVolume* physicalBottomAnode = new G4PVPlacement(0,G4ThreeVector(0,-(0.5*cryostat_sizeY-0.5*AnodeTickness),0),logicAnode,"AnodeDown",logicWorld,false,0,checkOverlaps);

    // -------------- Vikuiti/Reflector Walls -------------

    std::string file_vikuiti_reflectance = config_Reflector["reflectivity"].get<string>();
    std::ifstream f_vikuiti_reflectance(std::string("../configuration/")+file_vikuiti_reflectance);
    json data_reflectance_vikuiti = json::parse(f_vikuiti_reflectance);
    
    G4cout << "Vikuiti Reflectivity: " << std::endl;

    size_t n_vikuiti = data_reflectance_vikuiti.size();
    std::vector<G4double> energies_vikuiti(n_vikuiti);
    std::vector<G4double> r_vikuiti(n_vikuiti);
    for (size_t i = 0; i < n_vikuiti; i++) 
    {
        energies_vikuiti[i] = data_reflectance_vikuiti[i]["E"].get<double>()*eV;  
        r_vikuiti[i]   = data_reflectance_vikuiti[i]["r"].get<double>();
        G4cout <<   energies_vikuiti[i]/eV << " eV -- " << r_vikuiti[i] << std::endl;  
    }

    G4cout << "Vikuiti Spectrum: " << std::endl;
    std::string file_vikuiti_spectrum = config_Reflector["spectrum"].get<string>();
    std::ifstream f_vikuiti_spectrum(std::string("../configuration/")+file_vikuiti_spectrum);
    json data_spectrum_vikuiti = json::parse(f_vikuiti_spectrum);    
    size_t n_vikuiti_spec = data_spectrum_vikuiti.size();
    std::vector<G4double> energies_vikuiti_spec(n_vikuiti_spec);
    std::vector<G4double> spec_vikuiti(n_vikuiti_spec);
    for (size_t i = 0; i < n_vikuiti_spec; i++) 
    {
        energies_vikuiti_spec[i] = data_spectrum_vikuiti[i]["E"].get<double>()*eV;  
        spec_vikuiti[i]   = data_spectrum_vikuiti[i]["r"].get<double>();
        G4cout <<   energies_vikuiti_spec[i]/eV << " eV -- " << spec_vikuiti[i] << std::endl;  
    }

    G4cout << "Vikuiti Absorption Length: " << std::endl;
    std::string file_vikuiti_abslength = config_Reflector["absorption"].get<string>();
    std::ifstream f_vikuiti_abslength(std::string("../configuration/")+file_vikuiti_abslength);
    json data_abslength_vikuiti = json::parse(f_vikuiti_abslength);    
    size_t n_vikuiti_abslength = data_abslength_vikuiti.size();
    std::vector<G4double> energies_vikuiti_abslength(n_vikuiti_abslength);
    std::vector<G4double> abslength_vikuiti(n_vikuiti_abslength);
    std::cout << n_vikuiti_abslength << std::endl;
    for (size_t i = 0; i < n_vikuiti_abslength; i++) 
    {
        energies_vikuiti_abslength[i] = data_abslength_vikuiti[i]["E"].get<double>()*eV;  
        abslength_vikuiti[i]   = data_abslength_vikuiti[i]["l"].get<double>() * mm;
        G4cout <<   energies_vikuiti_abslength[i]/eV << " eV -- " << abslength_vikuiti[i]/mm << " *mm" << std::endl;  
    }

    auto Vikuiti_refractionindex = config_Reflector["refraction_index"].get<double>();
    auto Reflector_thickness = config_Reflector["thickness"].get<double>()*cm;
    auto Vikuiti_efficiency = config_Reflector["efficiency"].get<double>();
    G4cout << "Vikuiti Refraction Index: " << Vikuiti_refractionindex << std::endl;
    G4cout << "Vikuiti Thickness: " << Reflector_thickness/m << " m" << std::endl;
    G4cout << "Vikuiti Efficiency: " << Vikuiti_efficiency << std::endl;
    G4Material* Reflector_mat = fNistManager->FindOrBuildMaterial(config_Reflector["material"].get<string>());

    auto Vikuiti_isWLS = config_Reflector["isWLS"].get<int>();
    auto mpt_Reflector = new G4MaterialPropertiesTable();

    G4cout << "VIkuiti is a " << Vikuiti_isWLS << " a WLS\n";

    mpt_Reflector->AddProperty("RINDEX", {0.1*eV,15*eV} , {Vikuiti_refractionindex,Vikuiti_refractionindex},2);
    if(Vikuiti_isWLS==0)
    {
        mpt_Reflector->AddProperty("ABSLENGTH",energies_vikuiti_abslength,abslength_vikuiti , n_vikuiti_abslength);
    }
    else
    {
        mpt_Reflector->AddProperty("WLSABSLENGTH",energies_vikuiti_abslength,abslength_vikuiti , n_vikuiti_abslength);
        mpt_Reflector->AddProperty("WLSCOMPONENT", energies_vikuiti_spec, spec_vikuiti, n_vikuiti_spec);
        mpt_Reflector->AddConstProperty("WLSTIMECONSTANT",0.2*ns);
        mpt_Reflector->AddConstProperty("WLSMEANNUMBERPHOTONS",Vikuiti_efficiency);
    }
  
    Reflector_mat->SetMaterialPropertiesTable(mpt_Reflector);

    G4int only_half = config["only_half"].get<int>();

    G4double vikuiti_sizeY = (only_half == 1) ?  (0.5*cryostat_sizeY-AnodeTickness) :  (cryostat_sizeY-2*AnodeTickness);
    G4double vikuitiY = (only_half == 1) ? -vikuiti_sizeY/2. : 0.;

    G4Box* solidVikuitiLong = new G4Box("VikuitLong",0.5*cryostat_sizeX, 0.5*vikuiti_sizeY, 0.5*Reflector_thickness);
    G4Box* solidVikuitiEnd = new G4Box("VikuitEnd", 0.5*Reflector_thickness, 0.5*vikuiti_sizeY, 0.5*cryostat_sizeZ);
    G4LogicalVolume* logicVikuitiLong = new G4LogicalVolume(solidVikuitiLong,Reflector_mat,"VikuitiLong");
    G4LogicalVolume* logicVikuitiEnd = new G4LogicalVolume(solidVikuitiEnd,Reflector_mat,"VikuitiEnd");
    G4VPhysicalVolume* physicalVikuitiLong1 = new G4PVPlacement(0,G4ThreeVector(0,vikuitiY,-(cryostat_sizeZ/2-0.5*Reflector_thickness)),logicVikuitiLong,"LogicVikuiteLong",logicWorld,true,1,checkOverlaps);
    G4VPhysicalVolume* physicalVikuitiLong2 = new G4PVPlacement(0,G4ThreeVector(0,vikuitiY,(cryostat_sizeZ/2-0.5*Reflector_thickness)),logicVikuitiLong,"LogicVikuiteLong",logicWorld,true,2,checkOverlaps);
    G4VPhysicalVolume* physicalVikuitiEnd1 = new G4PVPlacement(0,G4ThreeVector(-(cryostat_sizeX/2-0.5*Reflector_thickness),vikuitiY,0),logicVikuitiEnd,"LogicVikuiteEnd",logicWorld,true,1,checkOverlaps);
    G4VPhysicalVolume* physicalVikuitiEnd2 = new G4PVPlacement(0,G4ThreeVector((cryostat_sizeX/2-0.5*Reflector_thickness),vikuitiY,0),logicVikuitiEnd,"LogicVikuiteEnd",logicWorld,true,2,checkOverlaps);

    auto visVikuiti = new G4VisAttributes();
    visVikuiti->SetForceSolid(true);
    visVikuiti->SetColor(G4Colour(1.,1.,1.));
    logicVikuitiEnd->SetVisAttributes(visVikuiti);
    logicVikuitiLong->SetVisAttributes(visVikuiti);

    // -----  Liquid Argon && ANODE && Vikuiti Interface Boundary --------

    G4OpticalSurface* surface_anode_lar = new G4OpticalSurface("surface_anode_lar");
    surface_anode_lar->SetModel(unified);
    surface_anode_lar->SetType(dielectric_metal);
    surface_anode_lar->SetFinish(polished);
    G4MaterialPropertiesTable* mpt_AnodeLar_Surface = new G4MaterialPropertiesTable();

    G4cout << "Anode Reflectivity: " << config_LAr["reflectivity_anode"].get<double>() << std::endl;
    mpt_AnodeLar_Surface->AddProperty("REFLECTIVITY",{1.0*eV,15*eV},{config_LAr["reflectivity_anode"].get<double>(),config_LAr["reflectivity_anode"].get<double>()},2);
    surface_anode_lar->SetMaterialPropertiesTable(mpt_AnodeLar_Surface);

    G4OpticalSurface* surface_vikuiti_lar = new G4OpticalSurface("surface_vikuiti_lar");
    surface_vikuiti_lar->SetModel(unified);
    surface_vikuiti_lar->SetType(dielectric_metal);
    surface_vikuiti_lar->SetFinish(polished);
    G4MaterialPropertiesTable* mpt_VikuitiLar_Surface = new G4MaterialPropertiesTable();
    mpt_VikuitiLar_Surface->AddProperty("REFLECTIVITY",energies_vikuiti,r_vikuiti,n_vikuiti);
    if(Vikuiti_isWLS==1)
    {
        mpt_VikuitiLar_Surface->AddProperty("TRANSMITTANCE",{1.0*eV,15*eV},{1,1},2);
    }
    surface_vikuiti_lar->SetMaterialPropertiesTable(mpt_VikuitiLar_Surface);

    G4OpticalSurface* surface_vikuiti_cryo = new G4OpticalSurface("surface_vikuiti_cryo");
    surface_vikuiti_cryo->SetModel(unified);
    surface_vikuiti_cryo->SetType(dielectric_metal);
    surface_vikuiti_cryo->SetFinish(polished);
    G4MaterialPropertiesTable* mpt_VikuitiCryo_Surface = new G4MaterialPropertiesTable();
    mpt_VikuitiCryo_Surface->AddProperty("REFLECTIVITY",{1.0*eV,15*eV},{0,0},2);
    surface_vikuiti_cryo->SetMaterialPropertiesTable(mpt_VikuitiCryo_Surface); 

    new G4LogicalBorderSurface("LiquidArgon-->AnodeTop", physicalWorld, physicalTopAnode , surface_anode_lar );
    new G4LogicalBorderSurface("LiquidArgon-->AnodeBottom", physicalWorld, physicalBottomAnode , surface_anode_lar );
    new G4LogicalBorderSurface("LiquidArgon-->VikuitiLong1", physicalWorld, physicalVikuitiLong1 , surface_vikuiti_lar );
    new G4LogicalBorderSurface("LiquidArgon-->VikuitiLong2", physicalWorld, physicalVikuitiLong2 , surface_vikuiti_lar );
    new G4LogicalBorderSurface("LiquidArgon-->VikuitiEnd1", physicalWorld, physicalVikuitiEnd1 , surface_vikuiti_lar );
    new G4LogicalBorderSurface("LiquidArgon-->VikuitiEnd2", physicalWorld, physicalVikuitiEnd2 , surface_vikuiti_lar );
    new G4LogicalBorderSurface("VikuitiLong1-->Cryo", physicalVikuitiLong1, physicalCryostat , surface_vikuiti_cryo );
    new G4LogicalBorderSurface("VikuitiLong2-->Cryo", physicalVikuitiLong2, physicalCryostat , surface_vikuiti_cryo );
    new G4LogicalBorderSurface("VikuitiEnd1-->Cryo", physicalVikuitiEnd1, physicalCryostat , surface_vikuiti_cryo );
    new G4LogicalBorderSurface("VikuitiEnd2-->Cryo", physicalVikuitiEnd2, physicalCryostat , surface_vikuiti_cryo );

    // ----------------------------------------------------------------------------
    //Assembling the field cage
    G4Material* FC_mat = fNistManager->FindOrBuildMaterial(config_FC["material"].get<string>());

    G4cout << "FC Refraction Index: " << config_FC["refraction_index"].get<double>() << std::endl;
    G4cout << "FC Absorption Length: " << config_FC["abs_length"].get<double>()*cm/m << "m" << std::endl;

    auto mpt_FC = new G4MaterialPropertiesTable();
    mpt_FC->AddProperty("RINDEX", {0.1*eV,15*eV}, {config_FC["refraction_index"].get<double>(),config_FC["refraction_index"].get<double>()}, 2);
    mpt_FC->AddProperty("ABSLENGTH", {0.1*eV,15*eV}, {config_FC["abs_length"].get<double>()*cm,config_FC["abs_length"].get<double>()*cm}, 2);
    FC_mat->SetMaterialPropertiesTable(mpt_FC);

    G4double d_cryoX =
    config_FC["distance_to_wall_x"].get<double>()*cm;

    G4double d_cryoZ =
        config_FC["distance_to_wall_z"].get<double>()*cm;
    //Starting with the vertical bars
    G4double FCv_sizeX = config_FC["vertical_bar"][0].get<double>()*cm;
    G4double FCv_sizeY = config_FC["vertical_bar"][1].get<double>()*cm;
    G4double FCv_sizeZ = config_FC["vertical_bar"][2].get<double>()*cm;
    G4double FCv_distance = config_FC["vertical_bar_distance"].get<double>()*cm;

    G4cout << "FC Distance from Membrane X: "
       << d_cryoX/m << " m" << std::endl;

    G4cout << "FC Distance from Membrane Z: "
       << d_cryoZ/m << " m" << std::endl;   
    G4cout << "FC Vertical bar dimensions: " << FCv_sizeX/m << " m X " <<  FCv_sizeY/m << " m X " << FCv_sizeZ/m << " m" << std::endl;
    G4cout << "FC Vertical bar distance: " << FCv_distance/m << "m" << std::endl;

    G4Box* solidBARv = new G4Box("vertical_bar",0.5*FCv_sizeX, 0.5*FCv_sizeY, 0.5*FCv_sizeZ);
    G4Box* solidBARv2 = new G4Box("vertical_bar2",0.5*FCv_sizeZ, 0.5*FCv_sizeY, 0.5*FCv_sizeX);
  
    int n_barv_long =
    (cryostat_sizeX-2*d_cryoX)/FCv_distance;
    G4UnionSolid* FCv_united =  new G4UnionSolid("FCv_longwall_template", solidBARv, solidBARv, 0, G4ThreeVector(FCv_distance,0,0));
    for(int i=1;i<n_barv_long/2;i++)
    {
        FCv_united = new G4UnionSolid("FCv_longwall_template", FCv_united, solidBARv, 0, G4ThreeVector((i+1)*FCv_distance,0,0));
        FCv_united = new G4UnionSolid("FCv_longwall_template", FCv_united, solidBARv, 0, G4ThreeVector((-i)*FCv_distance,0,0));     
    } 

    int n_barv_end =
    (cryostat_sizeZ-2*d_cryoZ)/FCv_distance;
    G4UnionSolid* FCv_united_end =  new G4UnionSolid("FCv_endwall_template", solidBARv2, solidBARv2, 0, G4ThreeVector(0,0,FCv_distance));
    for(int i=1;i<n_barv_end/2;i++)
    {
        FCv_united_end = new G4UnionSolid("FCv_endwall_template", FCv_united_end, solidBARv2, 0, G4ThreeVector(0,0,(i+1)*FCv_distance));
        FCv_united_end = new G4UnionSolid("FCv_endwall_template", FCv_united_end, solidBARv2, 0, G4ThreeVector(0,0,(-i)*FCv_distance));     
    } 

    double barH_eixoX1 = config_FC["horizontal_bar_1_majorAxis"].get<double>()*cm;
    double barH_eixoY1 = config_FC["horizontal_bar_1_minorAxis"].get<double>()*cm;
    double barH_eixoX2 = config_FC["horizontal_bar_2_majorAxis"].get<double>()*cm;
    double barH_eixoY2 = config_FC["horizontal_bar_2_minorAxis"].get<double>()*cm;
 
    double cut_value1 = config_FC["cut_horizontal_bar1"].get<double>()*cm;
    double cut_value2 = config_FC["cut_horizontal_bar2"].get<double>()*cm;

    G4cout << "FC Horizontal Bar 1 elipse major axis: " << barH_eixoX1/m  << " m" << std::endl;
    G4cout << "FC Horizontal Bar 1 elipse minor axis: " << barH_eixoY1/m  << " m" << std::endl;
    G4cout << "FC Horizontal Bar 1 length/cut: " << cut_value1/m  << " m" << std::endl;
    G4cout << "FC Horizontal Bar 2 elipse major axis: " << barH_eixoX2/m  << " m" << std::endl;
    G4cout << "FC Horizontal Bar 2 elipse minor axis: " << barH_eixoY2/m  << " m" << std::endl;
    G4cout << "FC Horizontal Bar 2 length/cut: " << cut_value2/m  << " m" << std::endl;

    double barH_length =
        cryostat_sizeX-2*d_cryoX-FCv_sizeX;

    double barH_length_end =
        cryostat_sizeZ-2*d_cryoZ-FCv_sizeZ-2*cut_value1;

    std::vector<G4TwoVector> ellipse1;
    double lastY=1e30; // inicializa com valor impossível
    int npoints = 120;
    for (int i = 0; i < npoints; i++) 
    {
        double phi = 2*CLHEP::pi*i/npoints;
        double x = (barH_eixoX1/2.0) * cos(phi);
        double y = (barH_eixoY1/2.0) * sin(phi);
        y = std::max(y, -(cut_value1 - barH_eixoY1/2));
        G4TwoVector point(x, y);
        // só adiciona se o ponto for diferente do anterior
        if (lastY != y) {
            ellipse1.push_back(point);
            lastY = y;
        }
    }

    std::vector<G4TwoVector> ellipse2;
    lastY = 1e30;
    for (int i = 0; i < npoints; i++) 
    {
        double phi = 2*CLHEP::pi*i/npoints;
        double x = (barH_eixoX2/2.0) * cos(phi);
        double y = (barH_eixoY2/2.0) * sin(phi);
        y = std::max(y, -(cut_value2 - barH_eixoY2/2));

        G4TwoVector point(x, y);
        if (lastY != y) {
            ellipse2.push_back(point);
            lastY = y;
        }
    }

    auto solidBARh1 = new G4ExtrudedSolid("solidEllipse1", ellipse1,  barH_length/2, G4TwoVector(0,0), 1.0,   G4TwoVector(0,0), 1.0);  
    auto solidBARh2 = new G4ExtrudedSolid("solidEllipse2", ellipse2,  barH_length/2, G4TwoVector(0,0), 1.0,   G4TwoVector(0,0), 1.0);  
    G4double FCh_distance = config_FC["horizontal_bar_distance"].get<double>()*cm;
    G4cout << "FC Horizontal bar distance: " << FCh_distance/m << " m" << std::endl;
    int n_barh_long = (FCv_sizeY)/FCh_distance;    

    G4UnionSolid* FCh_united =  new G4UnionSolid("FCh_longwall_template", solidBARh2, solidBARh2, 0, G4ThreeVector(FCh_distance,0,0));
    double sep_bar_type = config_FC["sep_bar_type"].get<double>()*cm;
    G4cout << "FC Horizontal bar type separation: " << sep_bar_type/m << " m" << std::endl;

    auto solidBARh1_end = new G4ExtrudedSolid("solidEllipse1_end", ellipse1,  barH_length_end/2, G4TwoVector(0,0), 1.0,   G4TwoVector(0,0), 1.0);  
    auto solidBARh2_end = new G4ExtrudedSolid("solidEllipse2_end", ellipse2,  barH_length_end/2, G4TwoVector(0,0), 1.0,   G4TwoVector(0,0), 1.0); 
    G4UnionSolid* FCh_united_end =  new G4UnionSolid("FCh_endwall_template", solidBARh2_end, solidBARh2_end, 0, G4ThreeVector(FCh_distance,0,0));

    for(int i=1;i<n_barh_long/2;i++)
    {
        if(((i+1)*FCh_distance - FCh_distance/2) >= sep_bar_type)
        {
            FCh_united = new G4UnionSolid("FCh_longwall_template", FCh_united, solidBARh1, 0, G4ThreeVector((i+1)*FCh_distance,0,0));
            FCh_united = new G4UnionSolid("FCh_longwall_template", FCh_united, solidBARh1, 0, G4ThreeVector((-i)*FCh_distance,0,0));
            FCh_united_end = new G4UnionSolid("FCh_endwall_template", FCh_united_end, solidBARh1_end, 0, G4ThreeVector((i+1)*FCh_distance,0,0));
            FCh_united_end = new G4UnionSolid("FCh_endwall_template", FCh_united_end, solidBARh1_end, 0, G4ThreeVector((-i)*FCh_distance,0,0));        
        }
        else
        {
            FCh_united = new G4UnionSolid("FCh_longwall_template", FCh_united, solidBARh2, 0, G4ThreeVector((i+1)*FCh_distance,0,0));
            FCh_united = new G4UnionSolid("FCh_longwall_template", FCh_united, solidBARh2, 0, G4ThreeVector((-i)*FCh_distance,0,0)); 
            FCh_united_end = new G4UnionSolid("FCh_endwall_template", FCh_united_end, solidBARh2_end, 0, G4ThreeVector((i+1)*FCh_distance,0,0));
            FCh_united_end = new G4UnionSolid("FCh_endwall_template", FCh_united_end, solidBARh2_end, 0, G4ThreeVector((-i)*FCh_distance,0,0));     
        }     
    } 
    G4RotationMatrix* rotTopY = new G4RotationMatrix();
    rotTopY->rotateY(90*deg);
    G4RotationMatrix* rotTopX = new G4RotationMatrix();
    rotTopX->rotateX(90*deg);
    G4RotationMatrix* minus_rotTopX = new G4RotationMatrix();
    minus_rotTopX->rotateX(-90*deg);
    G4RotationMatrix* rotTopZ = new G4RotationMatrix();
    rotTopZ->rotateZ(90*deg);
    G4RotationMatrix* minus_rotTopZ = new G4RotationMatrix();
    minus_rotTopZ->rotateZ(-90*deg);
    G4RotationMatrix* rotFinal1  = new G4RotationMatrix(); 
    *rotFinal1= (*rotTopY) * (*rotTopX);
    G4RotationMatrix* rotFinal2  = new G4RotationMatrix(); 
    *rotFinal2= (*rotTopY) * (*minus_rotTopX);
   
    G4LogicalVolume* LogicalFCv_wall = new G4LogicalVolume(FCv_united,FC_mat,"FC");
    G4LogicalVolume* LogicalFCv_end = new G4LogicalVolume(FCv_united_end,FC_mat,"FC");
    G4LogicalVolume* LogicalFCh_wall = new G4LogicalVolume(FCh_united,FC_mat,"FC");
    G4LogicalVolume* LogicalFCh_end = new G4LogicalVolume(FCh_united_end,FC_mat,"FC");
    G4VisAttributes* visFC = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5)); 
    visFC->SetForceSolid(true); 
    LogicalFCv_wall->SetVisAttributes(visFC);
    LogicalFCv_end->SetVisAttributes(visFC);
    LogicalFCh_wall->SetVisAttributes(visFC);
    LogicalFCh_end->SetVisAttributes(visFC);

    G4VPhysicalVolume* physicalFC1_v = new G4PVPlacement(0, G4ThreeVector(-FCv_distance/2, 0, cryostat_sizeZ/2 - d_cryoZ), LogicalFCv_wall, "FCv_long_wall", logicWorld, true, 1, checkOverlaps);
    G4VPhysicalVolume* physicalFC2_v = new G4PVPlacement(0, G4ThreeVector(-FCv_distance/2, 0, -cryostat_sizeZ/2 + d_cryoZ), LogicalFCv_wall, "FCv_long_wall", logicWorld, true, 2, checkOverlaps);
    G4VPhysicalVolume* physicalFC3_v = new G4PVPlacement(0, G4ThreeVector(cryostat_sizeX/2 - d_cryoX, 0, -FCv_distance/2), LogicalFCv_end, "FCv_end_wall", logicWorld, true, 1, checkOverlaps);
    G4VPhysicalVolume* physicalFC4_v = new G4PVPlacement(0, G4ThreeVector(-cryostat_sizeX/2 + d_cryoX, 0, -FCv_distance/2), LogicalFCv_end, "FCv_end_wall", logicWorld, true, 2, checkOverlaps);

    G4VPhysicalVolume* physicalFC1_h = new G4PVPlacement(rotFinal1, G4ThreeVector(0, -FCh_distance/2, cryostat_sizeZ/2 - d_cryoZ - FCv_sizeZ/2 - (cut_value1 - barH_eixoY1/2)), LogicalFCh_wall, "FCh_long_wall", logicWorld, true, 1, checkOverlaps);
    G4VPhysicalVolume* physicalFC2_h = new G4PVPlacement(rotFinal2, G4ThreeVector(0, +FCh_distance/2, -(cryostat_sizeZ/2 - d_cryoZ - FCv_sizeZ/2 - (cut_value1 - barH_eixoY1/2))), LogicalFCh_wall, "FCh_long_wall", logicWorld, true, 2, checkOverlaps);
    G4VPhysicalVolume* physicalFC3_h = new G4PVPlacement(rotTopZ, G4ThreeVector(-(cryostat_sizeX/2 - d_cryoX - FCv_sizeX/2 - (cut_value1 - barH_eixoY1/2)), +FCh_distance/2, 0), LogicalFCh_end, "FCh_end_wall", logicWorld, true, 1, checkOverlaps);
    G4VPhysicalVolume* physicalFC4_h = new G4PVPlacement(minus_rotTopZ, G4ThreeVector(cryostat_sizeX/2 - d_cryoX - FCv_sizeX/2 - (cut_value1 - barH_eixoY1/2), -FCh_distance/2, 0), LogicalFCh_end, "FCh_end_wall", logicWorld, true, 2, checkOverlaps);        

    // -----  Liquid Argon && FC Interface Boundary --------

    G4MaterialPropertiesTable* mpt_FCLar_Surface = new G4MaterialPropertiesTable();
    G4cout << "FC Reflectivity: " << config_FC["reflectivity"].get<double>() << std::endl;

    new G4LogicalBorderSurface("LiquidArgon-->FC1v", physicalWorld, physicalFC1_v , surface_cryo_lar );
    new G4LogicalBorderSurface("LiquidArgon-->FC2v", physicalWorld, physicalFC2_v , surface_cryo_lar );
    new G4LogicalBorderSurface("LiquidArgon-->FC3v", physicalWorld, physicalFC3_v , surface_cryo_lar );
    new G4LogicalBorderSurface("LiquidArgon-->FC4v", physicalWorld, physicalFC4_v , surface_cryo_lar );
    new G4LogicalBorderSurface("LiquidArgon-->FC1h", physicalWorld, physicalFC1_h , surface_cryo_lar );
    new G4LogicalBorderSurface("LiquidArgon-->FC2h", physicalWorld, physicalFC2_h , surface_cryo_lar );
    new G4LogicalBorderSurface("LiquidArgon-->FC3h", physicalWorld, physicalFC3_h , surface_cryo_lar );
    new G4LogicalBorderSurface("LiquidArgon-->FC4h", physicalWorld, physicalFC4_h , surface_cryo_lar );

    //Inserting the Acrylic --------------------------------------
    auto config_Acrylic = config["Acrylic"];

    auto density_acry = config_Acrylic["density"].get<double>()*g/cm3;
    G4cout << "Acrylic Density: " << density_acry*cm3/g << " g/cm3" << std::endl;
    std::vector<G4int>    natoms_acrylic = config_Acrylic["n_elements"].get<std::vector<G4int>>();
    std::vector<G4String> elements_acrylic = config_Acrylic["elements"].get<std::vector<G4String>>();

    auto mpt_Acry = new G4MaterialPropertiesTable();
    std::string file_acrylic_index = config_Acrylic["refraction_index"].get<string>();
    std::ifstream f_acrylic_index(std::string("../configuration/")+file_acrylic_index);
    
    G4cout << "Acrylic Refraction Index: " << std::endl;
    json data_acry = json::parse(f_acrylic_index);
    size_t n_rindex_acry = data_acry.size();
    std::vector<G4double> energies_r_acry(n_rindex_acry);
    std::vector<G4double> rindex_acry(n_rindex_acry); 
    for (size_t i = 0; i < n_rindex_acry; i++) 
    {
        energies_r_acry[i] = data_acry[i]["E"].get<double>()*eV;  
        rindex_acry[i]   = data_acry[i]["r"].get<double>();
        G4cout <<   energies_r_acry[i]/eV << " eV -- " << rindex_acry[i] << std::endl;  
    }
    G4cout << "Acrylic Absorption Length: " << std::endl;
    std::string file_acrylic_length = config_Acrylic["abs_length"].get<string>();
    std::ifstream f_acrylic_length(std::string("../configuration/")+file_acrylic_length);
    json data_acry_length = json::parse(f_acrylic_length);
    size_t n_abs_acry_len = data_acry_length.size() ;
    std::vector<G4double> absE_acry(n_abs_acry_len);
    std::vector<G4double> absLen_acry(n_abs_acry_len);
    for (size_t i=0; i<n_abs_acry_len; i++) 
    {
        absE_acry[i] = data_acry_length[i]["E"].get<double>()*eV;    
        absLen_acry[i] = data_acry_length[i]["l"].get<double>()*m;
        G4cout <<   absE_acry[i]/eV << " eV -- " << absLen_acry[i]/m << " m" <<std::endl;   
       
    }
    int n_abs_elem = 3;
    G4cout << "Acrylic Elements: " << std::endl;
    G4Material* Acry_mat = new G4Material("My_Acrylic",density_acry,n_abs_elem);
    for (size_t i=0; i<n_abs_elem; i++) 
    {
        G4Element* elem = fNistManager->FindOrBuildElement(elements_acrylic[i]);
        Acry_mat->AddElement(elem,natoms_acrylic[i]);
        G4cout <<  elements_acrylic[i] << " -- " << natoms_acrylic[i] << std::endl;   
    }

    mpt_Acry->AddProperty("RINDEX", energies_r_acry, rindex_acry, n_rindex_acry);
    mpt_Acry->AddProperty("ABSLENGTH", absE_acry, absLen_acry, n_abs_acry_len);
   
    Acry_mat->SetMaterialPropertiesTable(mpt_Acry);

    double acrylic_X =
    cryostat_sizeX - 2*d_cryoX - FCv_sizeX - 2*cut_value1;

    double acrylic_Y =
        (only_half == 1) ? FCv_sizeY/2 : FCv_sizeY;

    double acrylic_Z =
        cryostat_sizeZ - 2*d_cryoZ - FCv_sizeZ - 2*cut_value1;

    double acrylic_Ypos = (only_half == 1) ? -FCv_sizeY/4 : 0;

    double acry_thickness = config_Acrylic["thickness"].get<double>()*cm;
    G4cout << "Acrylic Thickness: " << acry_thickness/m << " m" << std::endl;
    auto acrylic_wall_long = new G4Box("Acrylic_long", 0.5*acrylic_X, 0.5*acrylic_Y, 0.5*acry_thickness);
    auto acrylic_wall_end = new G4Box("Acrylic_end",  0.5*acry_thickness, 0.5*acrylic_Y, 0.5*acrylic_Z);
    auto logical_Acrylic_wall_long = new G4LogicalVolume(acrylic_wall_long, Acry_mat , "Acrylic_long");
    auto logical_Acrylic_wall_end = new G4LogicalVolume(acrylic_wall_end, Acry_mat , "Acrylic_end");

    G4VPhysicalVolume* physicalAcrylical1 = new G4PVPlacement(
    0,
    G4ThreeVector(
        0,
        acrylic_Ypos,
        cryostat_sizeZ/2
        - d_cryoZ
        - 0.5*FCv_sizeZ
        - cut_value1
        - 0.5*acry_thickness
    ),
    logical_Acrylic_wall_long,
    "Acrylical_long",
    logicWorld,
    true,
    1,
    checkOverlaps
    );

    G4VPhysicalVolume* physicalAcrylical2 = new G4PVPlacement(
        0,
        G4ThreeVector(
            0,
            acrylic_Ypos,
            -(cryostat_sizeZ/2
            - d_cryoZ
            - 0.5*FCv_sizeZ
            - cut_value1
            - 0.5*acry_thickness)
        ),
        logical_Acrylic_wall_long,
        "Acrylical_long",
        logicWorld,
        true,
        2,
        checkOverlaps
    );

    G4VPhysicalVolume* physicalAcrylical3 = new G4PVPlacement(
        0,
        G4ThreeVector(
            cryostat_sizeX/2
            - d_cryoX
            - 0.5*FCv_sizeX
            - cut_value1
            - 0.5*acry_thickness,
            acrylic_Ypos,
            0
        ),
        logical_Acrylic_wall_end,
        "Acrylical_end",
        logicWorld,
        true,
        1,
        checkOverlaps
    );

    G4VPhysicalVolume* physicalAcrylical4 = new G4PVPlacement(
        0,
        G4ThreeVector(
            -(cryostat_sizeX/2
            - d_cryoX
            - 0.5*FCv_sizeX
            - cut_value1
            - 0.5*acry_thickness),
            acrylic_Ypos,
            0
        ),
        logical_Acrylic_wall_end,
        "Acrylical_end",
        logicWorld,
        true,
        2,
        checkOverlaps
    );
        
    G4VisAttributes* visAcry = new G4VisAttributes(G4Colour(0.8, 0.1, 0.1)); 
    visAcry->SetForceSolid(true); 
    logical_Acrylic_wall_long->SetVisAttributes(visAcry);
    logical_Acrylic_wall_end->SetVisAttributes(visAcry);

    // interface between acrylic, metal and argon
    G4OpticalSurface* surface_acrylic_lar = new G4OpticalSurface("surface_acry_lar");
    surface_acrylic_lar->SetModel(unified);
    surface_acrylic_lar->SetType(dielectric_dielectric);
    surface_acrylic_lar->SetFinish(polished);
    G4MaterialPropertiesTable* mpt_AcryLar_Surface = new G4MaterialPropertiesTable();
    mpt_AcryLar_Surface->AddProperty("REFLECTIVITY",{0.1*eV,15*eV}, {1,1},2);
    surface_acrylic_lar->SetMaterialPropertiesTable(mpt_AcryLar_Surface);

    new G4LogicalBorderSurface("Acrylic1-->Argon", physicalAcrylical1, physicalWorld , surface_acrylic_lar );
    new G4LogicalBorderSurface("Acrylic2-->Argon", physicalAcrylical2, physicalWorld , surface_acrylic_lar );
    new G4LogicalBorderSurface("Argon-->Acrylic1", physicalWorld ,  physicalAcrylical1,surface_acrylic_lar );
    new G4LogicalBorderSurface("Argon-->Acrylic2", physicalWorld , physicalAcrylical2, surface_acrylic_lar ); 
    new G4LogicalBorderSurface("Acrylic3-->Argon", physicalAcrylical3, physicalWorld , surface_acrylic_lar );
    new G4LogicalBorderSurface("Acrylic4-->Argon", physicalAcrylical4, physicalWorld , surface_acrylic_lar );
    new G4LogicalBorderSurface("Argon-->Acrylic3", physicalWorld ,  physicalAcrylical3,surface_acrylic_lar );
    new G4LogicalBorderSurface("Argon-->Acrylic4", physicalWorld , physicalAcrylical4, surface_acrylic_lar ); 
    new G4LogicalBorderSurface("Acrylic1-->FC1v", physicalAcrylical1, physicalFC1_v , surface_cryo_lar );
    new G4LogicalBorderSurface("Acrylic1-->FC2v", physicalAcrylical1, physicalFC2_v , surface_cryo_lar );
    new G4LogicalBorderSurface("Acrylic1-->FC3v", physicalAcrylical1, physicalFC3_v , surface_cryo_lar );
    new G4LogicalBorderSurface("Acrylic1-->FC4v", physicalAcrylical1, physicalFC4_v , surface_cryo_lar );
    new G4LogicalBorderSurface("Acrylic1-->FC1h", physicalAcrylical1, physicalFC1_h , surface_cryo_lar );
    new G4LogicalBorderSurface("Acrylic1-->FC2h", physicalAcrylical1, physicalFC2_h , surface_cryo_lar );
    new G4LogicalBorderSurface("Acrylic1-->FC3h", physicalAcrylical1, physicalFC3_h , surface_cryo_lar );
    new G4LogicalBorderSurface("Acrylic1-->FC4h", physicalAcrylical1, physicalFC4_h , surface_cryo_lar );
    new G4LogicalBorderSurface("Acrylic2-->FC1v", physicalAcrylical2, physicalFC1_v , surface_cryo_lar );
    new G4LogicalBorderSurface("Acrylic2-->FC2v", physicalAcrylical2, physicalFC2_v , surface_cryo_lar );
    new G4LogicalBorderSurface("Acrylic1-->FC3v", physicalAcrylical2, physicalFC3_v , surface_cryo_lar );
    new G4LogicalBorderSurface("Acrylic2-->FC4v", physicalAcrylical2, physicalFC4_v , surface_cryo_lar );
    new G4LogicalBorderSurface("Acrylic2-->FC1h", physicalAcrylical2, physicalFC1_h , surface_cryo_lar );
    new G4LogicalBorderSurface("Acrylic2-->FC2h", physicalAcrylical2, physicalFC2_h , surface_cryo_lar );
    new G4LogicalBorderSurface("Acrylic2-->FC3h", physicalAcrylical2, physicalFC3_h , surface_cryo_lar );
    new G4LogicalBorderSurface("Acrylic2-->FC4h", physicalAcrylical2, physicalFC4_h , surface_cryo_lar );
    new G4LogicalBorderSurface("Acrylic3-->FC1v", physicalAcrylical3, physicalFC1_v , surface_cryo_lar );
    new G4LogicalBorderSurface("Acrylic3-->FC2v", physicalAcrylical3, physicalFC2_v , surface_cryo_lar );
    new G4LogicalBorderSurface("Acrylic3-->FC3v", physicalAcrylical3, physicalFC3_v , surface_cryo_lar );
    new G4LogicalBorderSurface("Acrylic3-->FC4v", physicalAcrylical3, physicalFC4_v , surface_cryo_lar );
    new G4LogicalBorderSurface("Acrylic3-->FC1h", physicalAcrylical3, physicalFC1_h , surface_cryo_lar );
    new G4LogicalBorderSurface("Acrylic3-->FC2h", physicalAcrylical3, physicalFC2_h , surface_cryo_lar );
    new G4LogicalBorderSurface("Acrylic3-->FC3h", physicalAcrylical3, physicalFC3_h , surface_cryo_lar );
    new G4LogicalBorderSurface("Acrylic3-->FC4h", physicalAcrylical3, physicalFC4_h , surface_cryo_lar );
    new G4LogicalBorderSurface("Acrylic4-->FC1v", physicalAcrylical4, physicalFC1_v , surface_cryo_lar );
    new G4LogicalBorderSurface("Acrylic4-->FC2v", physicalAcrylical4, physicalFC2_v , surface_cryo_lar );
    new G4LogicalBorderSurface("Acrylic4-->FC3v", physicalAcrylical4, physicalFC3_v , surface_cryo_lar );
    new G4LogicalBorderSurface("Acrylic4-->FC4v", physicalAcrylical4, physicalFC4_v , surface_cryo_lar );
    new G4LogicalBorderSurface("Acrylic4-->FC1h", physicalAcrylical4, physicalFC1_h , surface_cryo_lar );
    new G4LogicalBorderSurface("Acrylic4-->FC2h", physicalAcrylical4, physicalFC2_h , surface_cryo_lar );
    new G4LogicalBorderSurface("Acrylic4-->FC3h", physicalAcrylical4, physicalFC3_h , surface_cryo_lar );
    new G4LogicalBorderSurface("Acrylic4-->FC4h", physicalAcrylical4, physicalFC4_h , surface_cryo_lar );
   
    //-----------------   Inserting PEN -------------------------------

    auto config_PEN = config["PEN"];
    auto density_PEN = config_PEN["density"].get<double>()*g/cm3;
    G4cout << "PEN Density: " << density_PEN*cm3/g << " g/cm3" << std::endl;
    std::vector<G4int>  natoms_PEN = config_PEN["n_elements"].get<std::vector<G4int>>();
    std::vector<G4String> elements_PEN = config_PEN["elements"].get<std::vector<G4String>>();
    int n_pen = 3;
    G4Material* PEN_mat = new G4Material("My_PEN",density_PEN,n_pen);
    G4cout << "PEN Elements: " << std::endl;
    for (size_t i=0; i<n_pen; i++) 
    {
        G4Element* elem = fNistManager->FindOrBuildElement(elements_PEN[i]);
        PEN_mat->AddElement(elem,natoms_PEN[i]);
        G4cout <<  elements_PEN[i] << " -- " << natoms_PEN[i] << std::endl;   
    }
    G4cout << "PEN Absorption Length: " << std::endl;
    std::string file_PEN_abs = config_PEN["abs_length"].get<string>();
    std::ifstream f_pen_abs(std::string("../configuration/")+file_PEN_abs);
    json data_pen = json::parse(f_pen_abs);
    size_t n_abs_pen = data_pen.size();
    std::vector<G4double> energies_abs_pen(n_abs_pen);
    std::vector<G4double> abs_pen(n_abs_pen);
    for (size_t i = 0; i < n_abs_pen; i++) 
    {
        energies_abs_pen[i] = data_pen[i]["E"].get<double>()*eV;  
        abs_pen[i] = data_pen[i]["l"].get<double>()*mm;
        G4cout <<   energies_abs_pen[i]/eV << " eV -- " << abs_pen[i]/m << " m" <<std::endl;   
    }

    G4cout << "PEN Emissium Spectrum: " << std::endl;
    std::string file_PEN_em = config_PEN["emissium"].get<string>();
    std::ifstream f_pen_em(std::string("../configuration/")+file_PEN_em);
    json data_pen_em = json::parse(f_pen_em);
    size_t n_em_pen = data_pen_em.size();
    std::vector<G4double> energies_em_pen(n_em_pen);
    std::vector<G4double> em_pen(n_em_pen);
    for (size_t i = 0; i < n_em_pen; i++) 
    {
        energies_em_pen[i] = data_pen_em[i]["E"].get<double>()*eV;  
        em_pen[i] = data_pen_em[i]["l"].get<double>();
        G4cout <<   energies_em_pen[i]/eV << " eV -- " << em_pen[i] <<std::endl;   
    }

    auto mpt_PEN = new G4MaterialPropertiesTable();
    auto refraction_index_pen = config_PEN["refraction_index"].get<double>();
    auto time_constant_pen=config_PEN["time_constant"].get<double>()*ns;
    auto eff_pen=config_PEN["efficiency"].get<double>();
    auto ray_pen=config_PEN["rayleigh_scattering"].get<double>()*um;

    G4cout << "PEN Refraction Index: " << refraction_index_pen << std::endl;
    G4cout << "PEN Efficiency " << eff_pen << std::endl;
    G4cout << "PEN Time constant: " << time_constant_pen/ns << " ns"<< std::endl;
    G4cout << "PEN Birks Constant: " << config_PEN["birks"].get<double>() << " cm/MeV" << std::endl;
    G4cout << "PEN Mean Excitation Energy: " << config_PEN["mean_excitation"].get<double>() << " eV" << std::endl;
    G4cout << "PEN Rayleigh: " << ray_pen/m << " m" << std::endl;

    mpt_PEN->AddProperty("RINDEX",{0.1*eV,15*eV},{refraction_index_pen,refraction_index_pen} , 2);
    mpt_PEN->AddProperty("WLSABSLENGTH", energies_abs_pen, abs_pen,n_abs_pen);
    mpt_PEN->AddProperty("WLSCOMPONENT", energies_em_pen, em_pen, n_em_pen);
    mpt_PEN->AddProperty("RAYLEIGH", {0.1*eV,15*eV} , {ray_pen, ray_pen}, 2);
    mpt_PEN->AddConstProperty("WLSTIMECONSTANT",time_constant_pen);
    mpt_PEN->AddConstProperty("WLSMEANNUMBERPHOTONS",eff_pen);
    PEN_mat->SetMaterialPropertiesTable(mpt_PEN);
    PEN_mat->GetIonisation()->SetBirksConstant(config_PEN["birks"].get<double>()*cm/MeV);
    PEN_mat->GetIonisation()->SetMeanExcitationEnergy(config_PEN["mean_excitation"].get<double>()*eV);
    
    double pen_X = acrylic_X-2*acry_thickness;
    double pen_Y = acrylic_Y;
    double pen_Z = acrylic_Z-2*acry_thickness;
    double pen_thickness = config_PEN["thickness"].get<double>()*cm;
    auto pen_wall_long = new G4Box("PEN_long", 0.5*pen_X, 0.5*pen_Y, 0.5*pen_thickness);
    auto pen_wall_end = new G4Box("PEN_end",  0.5*pen_thickness, 0.5*pen_Y, 0.5*pen_Z);

    auto logical_pen_wall_long = new G4LogicalVolume(pen_wall_long, PEN_mat , "PEN_long");
    auto logical_pen_wall_end = new G4LogicalVolume(pen_wall_end, PEN_mat , "PEN_end");

    G4VPhysicalVolume* physicalPEN1 = new G4PVPlacement(
    0,
    G4ThreeVector(
        0,
        acrylic_Ypos,
        cryostat_sizeZ/2
        - d_cryoZ
        - 0.5*FCv_sizeZ
        - cut_value1
        - acry_thickness
        - 0.5*pen_thickness
    ),
    logical_pen_wall_long,
    "PEN_long",
    logicWorld,
    true,
    1,
    checkOverlaps
    );

    G4VPhysicalVolume* physicalPEN2 = new G4PVPlacement(
        0,
        G4ThreeVector(
            0,
            acrylic_Ypos,
            -(cryostat_sizeZ/2
            - d_cryoZ
            - 0.5*FCv_sizeZ
            - cut_value1
            - acry_thickness
            - 0.5*pen_thickness)
        ),
        logical_pen_wall_long,
        "PEN_long",
        logicWorld,
        true,
        2,
        checkOverlaps
    );

    G4VPhysicalVolume* physicalPEN3 = new G4PVPlacement(
        0,
        G4ThreeVector(
            cryostat_sizeX/2
            - d_cryoX
            - 0.5*FCv_sizeX
            - cut_value1
            - acry_thickness
            - 0.5*pen_thickness,
            acrylic_Ypos,
            0
        ),
        logical_pen_wall_end,
        "PEN_end",
        logicWorld,
        true,
        1,
        checkOverlaps
    );

    G4VPhysicalVolume* physicalPEN4 = new G4PVPlacement(
        0,
        G4ThreeVector(
            -(cryostat_sizeX/2
            - d_cryoX
            - 0.5*FCv_sizeX
            - cut_value1
            - acry_thickness
            - 0.5*pen_thickness),
            acrylic_Ypos,
            0
        ),
        logical_pen_wall_end,
        "PEN_end",
        logicWorld,
        true,
        2,
        checkOverlaps
    );
    
    G4VisAttributes* visPEN = new G4VisAttributes(G4Colour(1, 1, 0.1)); 
    visPEN->SetForceSolid(true); 
    logical_pen_wall_long->SetVisAttributes(visPEN);
    logical_pen_wall_end->SetVisAttributes(visPEN);

    // SURFACES between PEN, acrylic and argon
    G4OpticalSurface* surface_PEN_lar = new G4OpticalSurface("surface_PEN_lar");
    surface_PEN_lar->SetModel(unified);
    surface_PEN_lar->SetType(dielectric_dielectric);
    surface_PEN_lar->SetFinish(polished);
    G4MaterialPropertiesTable* mpt_PENLar_Surface = new G4MaterialPropertiesTable();
    mpt_PENLar_Surface->AddProperty("REFLECTIVITY",{0.1*eV,15*eV}, {1,1},2);
    surface_PEN_lar->SetMaterialPropertiesTable(mpt_PENLar_Surface);

    G4OpticalSurface* surface_PEN_acry = new G4OpticalSurface("surface_PEN_acry");
    surface_PEN_acry->SetModel(unified);
    surface_PEN_acry->SetType(dielectric_dielectric);
    surface_PEN_acry->SetFinish(polished);
    G4MaterialPropertiesTable* mpt_PENAcry_Surface = new G4MaterialPropertiesTable();
    mpt_PENAcry_Surface->AddProperty("REFLECTIVITY",{0.1*eV,15*eV}, {1,1},2);
    surface_PEN_acry->SetMaterialPropertiesTable(mpt_PENAcry_Surface);

    new G4LogicalBorderSurface("Acrylic1 --> PEN1", physicalAcrylical1, physicalPEN1 , surface_PEN_acry );
    new G4LogicalBorderSurface("PEN1 --> Acrylic1", physicalPEN1, physicalAcrylical1 , surface_PEN_acry );
    new G4LogicalBorderSurface("Acrylic1 --> PEN2", physicalAcrylical1, physicalPEN2 , surface_PEN_acry );
    new G4LogicalBorderSurface("PEN2 --> Acrylic1", physicalPEN2, physicalAcrylical1 , surface_PEN_acry );
    new G4LogicalBorderSurface("Acrylic1 --> PEN3", physicalAcrylical1, physicalPEN3 , surface_PEN_acry );
    new G4LogicalBorderSurface("PEN3 --> Acrylic1", physicalPEN3, physicalAcrylical1 , surface_PEN_acry );
    new G4LogicalBorderSurface("Acrylic1 --> PEN4", physicalAcrylical1, physicalPEN4 , surface_PEN_acry );
    new G4LogicalBorderSurface("PEN4 --> Acrylic1", physicalPEN4, physicalAcrylical1 , surface_PEN_acry );
    new G4LogicalBorderSurface("Acrylic2 --> PEN1", physicalAcrylical2, physicalPEN1 , surface_PEN_acry );
    new G4LogicalBorderSurface("PEN1 --> Acrylic2", physicalPEN1, physicalAcrylical2 , surface_PEN_acry );
    new G4LogicalBorderSurface("Acrylic2 --> PEN2", physicalAcrylical2, physicalPEN2 , surface_PEN_acry );
    new G4LogicalBorderSurface("PEN2 --> Acrylic2", physicalPEN2, physicalAcrylical2 , surface_PEN_acry );
    new G4LogicalBorderSurface("Acrylic2 --> PEN3", physicalAcrylical2, physicalPEN3 , surface_PEN_acry );
    new G4LogicalBorderSurface("PEN3 --> Acrylic2", physicalPEN3, physicalAcrylical2 , surface_PEN_acry );
    new G4LogicalBorderSurface("Acrylic2 --> PEN4", physicalAcrylical2, physicalPEN4 , surface_PEN_acry );
    new G4LogicalBorderSurface("PEN4 --> Acrylic2", physicalPEN4, physicalAcrylical2 , surface_PEN_acry ); 
    new G4LogicalBorderSurface("Acrylic3 --> PEN1", physicalAcrylical3, physicalPEN1 , surface_PEN_acry );
    new G4LogicalBorderSurface("PEN1 --> Acrylic3", physicalPEN1, physicalAcrylical3 , surface_PEN_acry );
    new G4LogicalBorderSurface("Acrylic3 --> PEN2", physicalAcrylical3, physicalPEN2 , surface_PEN_acry );
    new G4LogicalBorderSurface("PEN2 --> Acrylic3", physicalPEN2, physicalAcrylical3 , surface_PEN_acry );
    new G4LogicalBorderSurface("Acrylic3 --> PEN3", physicalAcrylical3, physicalPEN3 , surface_PEN_acry );
    new G4LogicalBorderSurface("PEN3 --> Acrylic3", physicalPEN3, physicalAcrylical3 , surface_PEN_acry );
    new G4LogicalBorderSurface("Acrylic3 --> PEN4", physicalAcrylical3, physicalPEN4 , surface_PEN_acry );
    new G4LogicalBorderSurface("PEN4 --> Acrylic3", physicalPEN4, physicalAcrylical3 , surface_PEN_acry );
    new G4LogicalBorderSurface("Acrylic4 --> PEN1", physicalAcrylical4, physicalPEN1 , surface_PEN_acry );
    new G4LogicalBorderSurface("PEN1 --> Acrylic4", physicalPEN1, physicalAcrylical4 , surface_PEN_acry );
    new G4LogicalBorderSurface("Acrylic4 --> PEN2", physicalAcrylical4, physicalPEN2 , surface_PEN_acry );
    new G4LogicalBorderSurface("PEN2 --> Acrylic4", physicalPEN2, physicalAcrylical4 , surface_PEN_acry );
    new G4LogicalBorderSurface("Acrylic4 --> PEN3", physicalAcrylical4, physicalPEN3 , surface_PEN_acry );
    new G4LogicalBorderSurface("PEN3 --> Acrylic4", physicalPEN3, physicalAcrylical4 , surface_PEN_acry );
    new G4LogicalBorderSurface("Acrylic4 --> PEN4", physicalAcrylical4, physicalPEN4 , surface_PEN_acry );
    new G4LogicalBorderSurface("PEN4 --> Acrylic4", physicalPEN4, physicalAcrylical4 , surface_PEN_acry );

    new G4LogicalBorderSurface("LAr --> PEN1", physicalWorld, physicalPEN1 , surface_PEN_lar );
    new G4LogicalBorderSurface("LAr --> PEN2", physicalWorld, physicalPEN2 , surface_PEN_lar );
    new G4LogicalBorderSurface("LAr --> PEN3", physicalWorld, physicalPEN3 , surface_PEN_lar );
    new G4LogicalBorderSurface("LAr --> PEN4", physicalWorld, physicalPEN4 , surface_PEN_lar );
    new G4LogicalBorderSurface("PEN1 -->LAr", physicalPEN1 , physicalWorld, surface_PEN_lar );
    new G4LogicalBorderSurface("PEN2 -->LAr", physicalPEN2 , physicalWorld, surface_PEN_lar );
    new G4LogicalBorderSurface("PEN3 -->LAr", physicalPEN3 , physicalWorld, surface_PEN_lar );
    new G4LogicalBorderSurface("PEN4 -->LAr", physicalPEN4 , physicalWorld, surface_PEN_lar );

    // ----- Creating liquid argon REGION inside FC with lower light yield due to the electric field -----

    G4Material* lar_mat_inside = fNistManager->FindOrBuildMaterial(config_cryostat["filling"].get<string>());

    auto mpt_LAr2 = new G4MaterialPropertiesTable();
    mpt_LAr2->AddProperty("RINDEX", energies_r, rindex, n_rindex);
    mpt_LAr2->AddProperty("ABSLENGTH", absE_LAr, absLen_LAr, n_abs);
    mpt_LAr2->AddProperty("SCINTILLATIONCOMPONENT1",scintE_LAr,scint,n_scint);
    mpt_LAr2->AddProperty("SCINTILLATIONCOMPONENT2",scintE_LAr,scint,n_scint);
    mpt_LAr2->AddConstProperty("SCINTILLATIONYIELD",scint_lar["LY"].get<double>()/MeV);
    mpt_LAr2->AddConstProperty("SCINTILLATIONTIMECONSTANT1",scint_lar["Slow_Comp"].get<double>()*ns);
    mpt_LAr2->AddConstProperty("SCINTILLATIONTIMECONSTANT2",scint_lar["Fast_Comp"].get<double>()*ns);
    mpt_LAr2->AddConstProperty("SCINTILLATIONYIELD1", scint_lar["Slow_Amp"].get<double>());
    mpt_LAr2->AddConstProperty("SCINTILLATIONYIELD2", scint_lar["Fast_Amp"].get<double>());
    mpt_LAr2->AddConstProperty("RESOLUTIONSCALE",scint_lar["Resolution"].get<double>());
    mpt_LAr2->AddProperty("RAYLEIGH", energies_ray, ray_lar, n_rayleigh);
    lar_mat_inside->SetMaterialPropertiesTable(mpt_LAr2);

    double Lar_inside_X = pen_X-2*pen_thickness;
    double Lar_inside_Y = pen_Y;
    double Lar_inside_Z = pen_Z-2*pen_thickness;

    auto inside_argon = new G4Box("inside_argon", 0.5*Lar_inside_X, 0.5*Lar_inside_Y, 0.5*Lar_inside_Z);
    G4LogicalVolume* logicInsideArgon = new G4LogicalVolume(inside_argon,lar_mat_inside,"inside_argon");
    G4VPhysicalVolume* physicalInsideArgon = new G4PVPlacement(0,G4ThreeVector(),logicInsideArgon,"inside_argon",logicWorld,false,0,checkOverlaps);

    // Re-making some surfaces but with the argon inside

    new G4LogicalBorderSurface("LiquidArgonInside-->FC1v", physicalInsideArgon, physicalFC1_v , surface_cryo_lar );
    new G4LogicalBorderSurface("LiquidArgonInside-->FC2v", physicalInsideArgon, physicalFC2_v , surface_cryo_lar );
    new G4LogicalBorderSurface("LiquidArgonInside-->FC3v", physicalInsideArgon, physicalFC3_v , surface_cryo_lar );
    new G4LogicalBorderSurface("LiquidArgonInside-->FC4v", physicalInsideArgon, physicalFC4_v , surface_cryo_lar );
    new G4LogicalBorderSurface("LiquidArgonInside-->FC1h", physicalInsideArgon, physicalFC1_h , surface_cryo_lar );
    new G4LogicalBorderSurface("LiquidArgonInside-->FC2h", physicalInsideArgon, physicalFC2_h , surface_cryo_lar );
    new G4LogicalBorderSurface("LiquidArgonInside-->FC3h", physicalInsideArgon, physicalFC3_h , surface_cryo_lar );
    new G4LogicalBorderSurface("LiquidArgonInside-->FC4h", physicalInsideArgon, physicalFC4_h , surface_cryo_lar );
    new G4LogicalBorderSurface("LAr_Inside --> PEN1", physicalInsideArgon, physicalPEN1 , surface_PEN_lar );
    new G4LogicalBorderSurface("LAr_Inside --> PEN2", physicalInsideArgon, physicalPEN2 , surface_PEN_lar );
    new G4LogicalBorderSurface("LAr_Inside --> PEN3", physicalInsideArgon, physicalPEN3 , surface_PEN_lar );
    new G4LogicalBorderSurface("LAr_Inside --> PEN4", physicalInsideArgon, physicalPEN4 , surface_PEN_lar );
    new G4LogicalBorderSurface("PEN1 -->LAr_Inside", physicalPEN1 , physicalInsideArgon, surface_PEN_lar );
    new G4LogicalBorderSurface("PEN2 -->LAr_Inside", physicalPEN2 , physicalInsideArgon, surface_PEN_lar );
    new G4LogicalBorderSurface("PEN3 -->LAr_Inside", physicalPEN3 , physicalInsideArgon, surface_PEN_lar );
    new G4LogicalBorderSurface("PEN4 -->LAr_Inside", physicalPEN4 , physicalInsideArgon, surface_PEN_lar );
    new G4LogicalBorderSurface("LiquidArgonInside-->AnodeTop", physicalInsideArgon, physicalTopAnode , surface_anode_lar );
    new G4LogicalBorderSurface("LiquidArgonInside-->AnodeBottom", physicalInsideArgon, physicalBottomAnode , surface_anode_lar );

    // ---- Creating Cathode grid -------- 
    auto config_Cathode = config["Cathode_Grid"];
    double cathode_X = config_Cathode["length_x"].get<double>()*cm; //pen_X - 2*pen_thickness;
    double cathode_Z = config_Cathode["length_z"].get<double>()*cm; //pen_Z - 2*pen_thickness;
    double cathode_hole_X = config_Cathode["hole_x"].get<double>()*cm;
    double cathode_hole_Y = config_Cathode["hole_y"].get<double>()*cm;
    double cathode_thickness = config_Cathode["thickness"].get<double>()*cm;
    G4cout << "Cathode Size: " << cathode_X/m << " m X " << cathode_Z/m << " m" << std::endl;
    G4cout << "Cathode Hole Size: " << cathode_hole_X/m << " m X " << cathode_hole_Y/m << " m" << std::endl;

    auto cathode = new G4Box("Cathode", 0.5*cathode_X, 0.5*cathode_thickness, 0.5*cathode_Z);
    auto logicalCathode = new G4LogicalVolume(cathode, FC_mat, "Cathode");
    G4VPhysicalVolume* physicalCathode = new G4PVPlacement(0, G4ThreeVector(0,0,0), logicalCathode, "Cathode", logicInsideArgon, false, 0,true);

    // ----------- Inserting VIKUITI over and under cathode -------------
    
    auto reflector = new G4Box("Reflector", 0.5*cathode_X, 0.5*Reflector_thickness, 0.5*cathode_Z);
    G4LogicalVolume* logicReflector = new G4LogicalVolume(reflector,Reflector_mat,"Reflector");
    G4VisAttributes* visReflector = new G4VisAttributes(G4Colour(0.,0.5,0.9));
    visReflector->SetForceSolid(true);
    logicReflector->SetVisAttributes(visReflector); 

    /*G4VPhysicalVolume* physicalReflectorTop = new G4PVPlacement(0,G4ThreeVector(0,(cryostatThickness+Reflector_thickness)/2,0)
        ,logicReflector,"My_Reflector",logicInsideArgon,true,1,checkOverlaps); */

    G4VPhysicalVolume* physicalReflectorBottom = new G4PVPlacement(0,G4ThreeVector(0,-(cathode_thickness+Reflector_thickness)/2,0)
        ,logicReflector,"My_Reflector",logicInsideArgon,true,2,checkOverlaps);

      //-------------- Inserting PEN over on top of both Reflector ----------

    auto Pen_cathode = new G4Box("PEN_cathode", 0.5*cathode_X, 0.5*pen_thickness, 0.5*cathode_Z);
    G4LogicalVolume* logicPENCathode = new G4LogicalVolume(Pen_cathode,PEN_mat,"PEN_Cathode");

    logicPENCathode->SetVisAttributes(visPEN); 


    /* G4VPhysicalVolume* physicalPENTop = new G4PVPlacement(0,G4ThreeVector(0,(cathode_thickness+2*Reflector_thickness+pen_thickness)/2,0)
        ,logicPENCathode,"PEN_CATHODE",logicInsideArgon,true,1,checkOverlaps); */

    G4VPhysicalVolume* physicalPENBottom = new G4PVPlacement(0,G4ThreeVector(0,-(cathode_thickness+2*Reflector_thickness+pen_thickness)/2,0)
        ,logicPENCathode,"PEN_CATHODE",logicInsideArgon,true,2,checkOverlaps);

    // ---- Creating Cathode grid -------- PART 2 --- ADDING HOLES

    auto cathodehole = new G4Box("Cathode_Hole_argon", 0.5*cathode_hole_X, 0.5*cathode_thickness, 0.5*cathode_hole_Y);
    auto logicalHole = new G4LogicalVolume(cathodehole, lar_mat_inside, "Cathode_Hole_argon");
    
    double cathode_internal_distance = config_Cathode["internal_distance"].get<double>()*cm;
    double cathode_external_distance = config_Cathode["external_distance"].get<double>()*cm;
    double cathode_block_X = config_Cathode["block_x"].get<double>()*cm;
    double cathode_block_Z = config_Cathode["block_z"].get<double>()*cm;

    int n_blocks_x = static_cast<int>(std::round(cathode_X / cathode_block_X));
    int n_blocks_z = static_cast<int>(std::round(cathode_Z / cathode_block_Z));

    double x_inner = 0.5*(cathode_hole_X + cathode_external_distance);
    double x_outer = x_inner + cathode_hole_X + cathode_internal_distance;

    double z_inner = 0.5*(cathode_hole_Y + cathode_external_distance);
    double z_outer = z_inner + cathode_hole_Y + cathode_internal_distance;

    std::vector<double> cathode_x_local = {-x_outer, -x_inner, x_inner, x_outer};
    std::vector<double> cathode_z_local = {-z_outer, -z_inner, z_inner, z_outer};

    std::vector<double> cathode_x_positions;
    for (int ibx = 0; ibx < n_blocks_x; ibx++)
    {
        double block_center_x = (ibx - (n_blocks_x - 1)/2.0) * cathode_block_X;

        for (double x_local : cathode_x_local)
            cathode_x_positions.push_back(block_center_x + x_local);
    }

    std::vector<double> cathode_z_positions;
    for (int ibz = 0; ibz < n_blocks_z; ibz++)
    {
        double block_center_z = (ibz - (n_blocks_z - 1)/2.0) * cathode_block_Z;

        for (double z_local : cathode_z_local)
            cathode_z_positions.push_back(block_center_z + z_local);
    }

    auto hole_percentage = config_PEN["hole_size_percentage"].get<double>();
    G4cout << "PEN/VIKUITI Hole in Cathode: " << hole_percentage << std::endl;
    double pen_hole_X = sqrt(hole_percentage)*cathode_hole_X;
    double pen_hole_Y = sqrt(hole_percentage)*cathode_hole_Y;


    G4VisAttributes* visArgon = new G4VisAttributes(G4Colour(0.7,0.85,0.9));
    visArgon->SetVisibility(true);
    visArgon->SetForceSolid(true);
    logicalHole->SetVisAttributes(visArgon); 
    


    for (size_t ix = 0; ix < cathode_x_positions.size(); ix++)
    {
        for (size_t iz = 0; iz < cathode_z_positions.size(); iz++)
        {
            double x_pos = cathode_x_positions[ix];
            double z_pos = cathode_z_positions[iz];

            G4ThreeVector pos(x_pos, 0, z_pos);

            int number = ix * cathode_z_positions.size() + iz;

            auto physicalHole = new G4PVPlacement(0,pos,logicalHole,"Cathode_Hole_argon",logicalCathode,true,number,false);

            new G4LogicalBorderSurface("LAr_in_Hole_" + std::to_string(number) + "-->Cathode", physicalHole, physicalCathode, surface_cryo_lar );
            new G4LogicalBorderSurface("LAr_in_Hole_" + std::to_string(number) + "-->PEN2", physicalHole, physicalPENBottom, surface_PEN_lar );
            new G4LogicalBorderSurface("PEN2-->LAr_in_Hole_" + std::to_string(number) , physicalPENBottom,  physicalHole, surface_PEN_lar );

        }
    }

    new G4LogicalBorderSurface("LiquidArgon-->Cathode", physicalWorld, physicalCathode , surface_cryo_lar );

    new G4LogicalBorderSurface("LiquidArgon-->ReflectorB", physicalWorld, physicalReflectorBottom , surface_vikuiti_lar );
    new G4LogicalBorderSurface("LiquidArgonInside-->Cathode", physicalInsideArgon, physicalCathode , surface_cryo_lar );
    new G4LogicalBorderSurface("LiquidArgonInside-->ReflectorB", physicalInsideArgon, physicalReflectorBottom , surface_vikuiti_lar );
    new G4LogicalBorderSurface("ReflectorB-->Cathode",  physicalReflectorBottom ,physicalCathode, surface_vikuiti_cryo );
    new G4LogicalBorderSurface("PEN1 --> Cathode Grid", physicalPEN1 , physicalCathode, surface_cryo_lar );
    new G4LogicalBorderSurface("PEN2 --> Cathode Grid", physicalPEN2 , physicalCathode, surface_cryo_lar );
    new G4LogicalBorderSurface("PEN3 --> Cathode Grid", physicalPEN3 , physicalCathode, surface_cryo_lar );
    new G4LogicalBorderSurface("PEN4 --> Cathode Grid", physicalPEN4 , physicalCathode, surface_cryo_lar );



    new G4LogicalBorderSurface("PEN1 --> ReflectorG", physicalPEN1 , physicalReflectorBottom, surface_vikuiti_lar );
    new G4LogicalBorderSurface("PEN2 --> ReflectorG", physicalPEN2 , physicalReflectorBottom, surface_vikuiti_lar );
    new G4LogicalBorderSurface("PEN3 --> ReflectorG", physicalPEN3 , physicalReflectorBottom, surface_vikuiti_lar );
    new G4LogicalBorderSurface("PEN4 --> ReflectorG", physicalPEN4 , physicalReflectorBottom, surface_vikuiti_lar );

    new G4LogicalBorderSurface("PENBottom --> ReflectorBottom", physicalPENBottom , physicalReflectorBottom, surface_vikuiti_lar );
    new G4LogicalBorderSurface("PENBottom --> LAr", physicalPENBottom , physicalWorld, surface_PEN_lar );
    new G4LogicalBorderSurface("LAr --> PENBottom",  physicalWorld, physicalPENBottom , surface_PEN_lar );
    new G4LogicalBorderSurface("PENBottom --> LArInside", physicalPENBottom , physicalInsideArgon, surface_PEN_lar );
    new G4LogicalBorderSurface("LArInside --> PENBottom",  physicalInsideArgon, physicalPENBottom , surface_PEN_lar );

    G4VisAttributes* visCathode = new G4VisAttributes(G4Colour(0.9,0.9,0.9));
    visCathode->SetForceSolid(true);
    //visCathode->SetVisibility(false);
    logicalCathode->SetVisAttributes(visCathode); 

    // ------------- Assembling the SiPMs  ---------------------

    auto config_SiPM = config["SiPM"];

    G4Material* SiPM_mat     = fNistManager->FindOrBuildMaterial("G4_Si");
    G4double rindex_SiPM    = config_SiPM["refraction_index"].get<double>();
    G4double abslength_SiPM = config_SiPM["abs_length"].get<double>()*cm;
    G4double size_SiPM = config_SiPM["size"].get<double>()*cm;
    G4double thickness_SiPM = config_SiPM["thickness"].get<double>()*cm;
    G4double dv_SiPM = config_SiPM["distance_vertical"].get<double>()*cm;
    G4double dh_SiPM = config_SiPM["distance_horizontal"].get<double>()*cm;
    G4double SiPM_reflective = config_SiPM["reflectance"].get<double>();

    std::vector<std::vector<std::string>> matrix_pos = config_SiPM["map_ldu"].get<std::vector<std::vector<std::string>>>(); //[i][j] = [linha][coluna]
    G4cout << "SiPM Refraction Index: " << rindex_SiPM << std::endl;
    G4cout << "SiPM Absorption Length: " << abslength_SiPM/m << " m" << std::endl;
    G4cout << "SiPM Size Single Unity Matriz: " << size_SiPM/m << " m" << std::endl;
    G4cout << "SiPM Thickness: " << thickness_SiPM/m << " m" << std::endl;
    G4cout << "SiPM Vertical Distance: " << dv_SiPM/m << " m" << std::endl;
    G4cout << "SiPM Horizontal Distance: " << dh_SiPM/m << " m" << std::endl;
    G4cout << "SiPM Reflective: " << SiPM_reflective << std::endl;

    int size=2;
    G4double* RIndex_array   = new G4double[size];
    G4double* abs_array   = new G4double[size];

    for (int i=0; i<size; i++) 
    {     
        RIndex_array[i]   = rindex_SiPM; 
        abs_array[i] = abslength_SiPM;
    }

    G4double photonEnergy[2] = {0.1*eV , 15*eV}; 
    G4MaterialPropertiesTable* mpt_SiPM = new G4MaterialPropertiesTable();
    mpt_SiPM->AddProperty("RINDEX",photonEnergy,RIndex_array,size);
    mpt_SiPM->AddProperty("ABSLENGTH",photonEnergy,abs_array,size);
    SiPM_mat->SetMaterialPropertiesTable(mpt_SiPM);

    G4Box* sipm_template = new G4Box("sipm" , size_SiPM/2 , size_SiPM/2 , thickness_SiPM/2);
    sipmvis_logical = new G4LogicalVolume(sipm_template, SiPM_mat , "sipm");
    sipmuv_logical = new G4LogicalVolume(sipm_template, SiPM_mat , "sipm");

    // Surface properties
    auto sipmSurface = new G4OpticalSurface("SiPM_Surface");
    sipmSurface->SetType(dielectric_metal);
    //sipmSurface->SetType(dielectric_dielectric);    
    sipmSurface->SetFinish(polished);
    sipmSurface->SetModel(unified);
    G4MaterialPropertiesTable* mpt_surface = new G4MaterialPropertiesTable();
    G4double reflectivity[2] = {config_SiPM["reflectance"].get<double>() , config_SiPM["reflectance"].get<double>()};       
    G4double transmitance[2] = {1.0 , 1.0};    
    mpt_surface->AddProperty("TRANSMITTANCE", photonEnergy, transmitance, 2);
    mpt_surface->AddProperty("REFLECTIVITY", photonEnergy, reflectivity, 2);
    sipmSurface->SetMaterialPropertiesTable(mpt_surface);

    auto sipmSurfaceVikuiti = new G4OpticalSurface("SiPM_Surface_Vikuiti");
    sipmSurfaceVikuiti->SetType(dielectric_metal);
    sipmSurfaceVikuiti->SetFinish(polished);
    sipmSurfaceVikuiti->SetModel(unified);
    G4MaterialPropertiesTable* mpt_surface2 = new G4MaterialPropertiesTable();
    G4double reflectivity2[2] = {0.0 , 0.0};   
    mpt_surface2->AddProperty("REFLECTIVITY", photonEnergy,reflectivity2, 2);
    sipmSurfaceVikuiti->SetMaterialPropertiesTable(mpt_surface2);

    std::vector<std::vector<G4ThreeVector>> base_shifts(4, std::vector<G4ThreeVector>(4, G4ThreeVector(0,0,0)));
    std::vector<std::vector<G4ThreeVector>> base_shifts_end(4, std::vector<G4ThreeVector>(4, G4ThreeVector(0,0,0)));
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            base_shifts[i][j]=G4ThreeVector(-1.5*size_SiPM+i*size_SiPM,-1.5*size_SiPM+j*size_SiPM,0);
        }
    }
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            base_shifts_end[i][j]=G4ThreeVector(0,-1.5*size_SiPM+j*size_SiPM,-1.5*size_SiPM+i*size_SiPM);
        }
    }

    std::vector<G4PVPlacement*> sipm_physicals;
    sipm_physicals.clear();
    int n_ldu_x = cryostat_sizeX/dh_SiPM;
    int n_ldu_y = cryostat_sizeY/dv_SiPM;

    int cont=0;
    for (int ix = 0; ix < n_ldu_x; ix++)
    {
        for (int iy = 0; iy < n_ldu_y; iy++) 
        {
            double x_pos = (ix - (n_ldu_x)/2.0 + 0.5)*(dh_SiPM);
            double y_pos = (iy - (n_ldu_y)/2.0 + 0.5)*(dv_SiPM);
        
            // Se only_half=true, coloca SiPM apenas no volume inferior (y < 0)
            if (only_half==1 && y_pos >= 0)
                continue;

            
            for (int side = -1; side <= 1; side += 2)
            {
                double z_pos = side * (cryostat_sizeZ/2 - Reflector_thickness- thickness_SiPM/2);               
                G4ThreeVector base(x_pos, y_pos, z_pos);
                for(int i=0;i<4;i++)
                {
                    for(int j=0;j<4;j++)
                    {
                        if(matrix_pos[i][j]=="V")
                        {
                            sipm_physicals.push_back(new G4PVPlacement(0,base+base_shifts[i][j],sipmvis_logical,"sipmVIS:" + std::to_string(i) + "_" + std::to_string(j)
                                ,logicWorld,true,cont,false));
                        }
                        if(matrix_pos[i][j]=="U")
                        {
                            sipm_physicals.push_back(new G4PVPlacement(0,base+base_shifts[i][j],sipmuv_logical,"sipmUV:" + std::to_string(i) + "_" + std::to_string(j)
                                ,logicWorld,true,cont,false));
                        }    
                        new G4LogicalBorderSurface("World_to_SiPM", physicalWorld, sipm_physicals.back(), sipmSurface);
                        new G4LogicalBorderSurface("Vikuti_to_SiPM1", physicalVikuitiLong1, sipm_physicals.back(), sipmSurfaceVikuiti);
                        new G4LogicalBorderSurface("Vikuti_to_SiPM2", physicalVikuitiLong2, sipm_physicals.back(), sipmSurfaceVikuiti);
                    }
                }
                cont++;
            }
        }
    }
    
    auto rotY90 = new G4RotationMatrix();
    rotY90->rotateY(90.*deg);
    int n_ldu_z = (cryostat_sizeZ)/dh_SiPM-1;
    for (int iz = 0; iz < n_ldu_z; iz++)
    {
        for (int iy = 0; iy < n_ldu_y; iy++) 
        {
            double z_pos = (iz - (n_ldu_z)/2.0 + 0.5)*(dh_SiPM);
            double y_pos = (iy - (n_ldu_y)/2.0 + 0.5)*(dv_SiPM);

            // Se only_half=true, coloca SiPM apenas no volume inferior (y < 0)
            if (only_half==1 && y_pos >= 0)
                continue;

            
            for (int side = -1; side <= 1; side += 2)
            {
                double x_pos = side * (cryostat_sizeX/2 - Reflector_thickness -thickness_SiPM/2);               
                G4ThreeVector base(x_pos, y_pos, z_pos);
                for(int i=0;i<4;i++)
                {
                    for(int j=0;j<4;j++)
                    {
                        if(matrix_pos[i][j]=="V")
                        {
                            sipm_physicals.push_back(new G4PVPlacement(rotY90,base+base_shifts_end[i][j],sipmvis_logical,"sipmVIS:" + std::to_string(i) + "_" + std::to_string(j)
                                ,logicWorld,true,cont,false));
                        }
                        if(matrix_pos[i][j]=="U")
                        {
                            sipm_physicals.push_back(new G4PVPlacement(rotY90,base+base_shifts_end[i][j],sipmuv_logical,"sipmUV:" + std::to_string(i) + "_" + std::to_string(j)
                                ,logicWorld,true,cont,false));
                        }    
                        new G4LogicalBorderSurface("World_to_SiPM", physicalWorld, sipm_physicals.back(), sipmSurface);
                        new G4LogicalBorderSurface("Vikuti_to_SiPM1", physicalVikuitiEnd1, sipm_physicals.back(), sipmSurfaceVikuiti);
                        new G4LogicalBorderSurface("Vikuti_to_SiPM2", physicalVikuitiEnd2, sipm_physicals.back(), sipmSurfaceVikuiti);

                    }
                }
                cont++;
            }
        }
    }

    // Visible properties
    G4VisAttributes* visSiPM = new G4VisAttributes(G4Colour(0.01, 1, 0.1)); 
    visSiPM->SetForceSolid(true); 
    sipmvis_logical->SetVisAttributes(visSiPM);

    G4VisAttributes* visSiPMuv = new G4VisAttributes(G4Colour(0.6, 0, 0.6)); 
    visSiPMuv->SetForceSolid(true); 
    sipmuv_logical->SetVisAttributes(visSiPMuv);

    return physicalWorld;
  
}

void DetectorConstruction::ConstructSDandField()
{
    G4SDManager* SDman = G4SDManager::GetSDMpointer();

    if (SDman->FindSensitiveDetector("SensitiveDetector", false) == nullptr) 
    {
        auto* sensitiveModule = new SensitiveDetector("SensitiveDetector", "HitCollection");
        SDman->AddNewDetector(sensitiveModule);

        sipmvis_logical->SetSensitiveDetector(sensitiveModule);
        sipmuv_logical->SetSensitiveDetector(sensitiveModule);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
