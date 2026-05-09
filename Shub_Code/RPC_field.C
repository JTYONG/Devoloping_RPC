#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <chrono>
#include <cstdlib>
#include <filesystem>           // For handling directories
namespace fs = std::filesystem; // Simplify namespace usage
#include <TApplication.h>

#include "Garfield/SolidBox.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/MediumConductor.hh"
#include "Garfield/MediumPlastic.hh"
#include "Garfield/ComponentNeBem3d.hh"
#include "Garfield/ViewGeometry.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/Random.hh"

using namespace std;
using namespace Garfield;

int main(int argc, char *argv[])
{

    TApplication app("app", &argc, argv);

    MediumMagboltz gas;
    gas.SetComposition("Ar", 90., "CO", 10.);
    gas.SetTemperature(293.15);
    gas.SetPressure(760.);
    gas.Initialise(true);
    gas.EnablePenningTransfer();
    // const std::string path = std::getenv("GARFIELD_INSTALL");
    // gas.LoadIonMobility(path + "/share/Garfield/Data/IonMobility_Ar+_Ar.txt");
    // gas.LoadIonMobility(path + "/share/Garfield/Data/IonMobility_CO2+_CO2.txt");

    MediumConductor graphite;
    MediumConductor copper;
    MediumPlastic mylar;
    MediumPlastic glass;

    glass.SetDielectricConstant(10.0);
    mylar.SetDielectricConstant(3.1);

    // Geometry.
    GeometrySimple geo;
    geo.SetMedium(&gas);

    static constexpr double GasGapThickness = 0.025;             // cm - total gas gap thickness            250 micron = 0.25 mm = 0.025 cm
    static constexpr double GasGapCenterZ = 0.0 / 10.0;          // cm - centered at Z=0                    y = 0 mm =0 cm
    static constexpr double AnodeCathodeThickness = 0.02 / 10.0; // cm - thickness of HV layers             20 micron = 0.02 mm = 0.002 cm
    static constexpr double StripThickness = 0.05 / 10.0;        // cm - thickness of copper strips         50 micron = 0.005 cm
    static constexpr double StripWidthX = 30.0 / 10;             // cm - readout strip width                30 mm = 3 cm
    static constexpr double DetectorSizeX = 30.0 / 10.0;         // cm - detector size in X                 30 mm = 3 cm
    static constexpr double DetectorSizeY = 30.0 / 10.0;         // cm - detector size in Y                 30 mm = 3 cm
    static constexpr double HoneyCombThickness = 0;              // cm - honey comb layer thickness         0 mm = 0 cm (Remove Honeycomb)
    static constexpr double MylarThickness = 0.035;              // cm - mylar layer thickness              350 micron = 0.35 mm = 0.035 cm
    static constexpr double ResistiveGlassThickness = 0.12;      // cm - resistive glass layer thickness  1.2 mm = 0.12 cm
    static constexpr double AnodeVoltage = 5500.0;               // V - ANODE at +5.5 kV
    static constexpr double CathodeVoltage = -5500.0;            // V - CATHODE at -5.5 kV

    // Positions

    static constexpr double topglassZ = GasGapCenterZ + GasGapThickness / 2.0 + ResistiveGlassThickness / 2.0;
    static constexpr double bottomglassZ = GasGapCenterZ - GasGapThickness / 2.0 - ResistiveGlassThickness / 2.0;
    static constexpr double topgraphiteZ = topglassZ + ResistiveGlassThickness / 2.0 + AnodeCathodeThickness / 2.0;
    static constexpr double bottomgraphiteZ = bottomglassZ - ResistiveGlassThickness / 2.0 - AnodeCathodeThickness / 2.0;
    static constexpr double topmylarZ = topgraphiteZ + AnodeCathodeThickness / 2.0 + MylarThickness / 2.0;
    static constexpr double bottommylarZ = bottomgraphiteZ - AnodeCathodeThickness / 2.0 - MylarThickness / 2.0;
    static constexpr double topreadoutZ = topmylarZ + MylarThickness / 2.0 + StripThickness / 2.0;
    static constexpr double bottomreadoutZ = bottommylarZ - MylarThickness / 2.0 - StripThickness / 2.0;

    SolidBox topglass(0, 0, topglassZ, DetectorSizeX / 2.0, DetectorSizeY / 2.0, ResistiveGlassThickness / 2.0); // top glass
    topglass.SetBoundaryDielectric();
    geo.AddSolid(&topglass, &glass);

    SolidBox bottomglass(0, 0, bottomglassZ, DetectorSizeX / 2.0, DetectorSizeY / 2.0, ResistiveGlassThickness / 2.0); // top glass
    bottomglass.SetBoundaryDielectric();
    geo.AddSolid(&bottomglass, &glass);

    SolidBox topgraphite(0, 0, topgraphiteZ, DetectorSizeX / 2.0, DetectorSizeY / 2.0, AnodeCathodeThickness / 2.0); // top glass
    topgraphite.SetBoundaryPotential(AnodeVoltage);
    geo.AddSolid(&topgraphite, &graphite);

    SolidBox bottomgraphite(0, 0, bottomgraphiteZ, DetectorSizeX / 2.0, DetectorSizeY / 2.0, AnodeCathodeThickness / 2.0); // top glass
    bottomgraphite.SetBoundaryPotential(CathodeVoltage);
    geo.AddSolid(&bottomgraphite, &graphite);

    SolidBox topmylar(0, 0, topmylarZ, DetectorSizeX / 2.0, DetectorSizeY / 2.0, MylarThickness / 2.0); // top glass
    topmylar.SetBoundaryDielectric();
    geo.AddSolid(&topmylar, &mylar);

    SolidBox bottommylar(0, 0, bottommylarZ, DetectorSizeX / 2.0, DetectorSizeY / 2.0, MylarThickness / 2.0); // top glass
    bottommylar.SetBoundaryDielectric();
    geo.AddSolid(&bottommylar, &mylar);

    SolidBox topreadout(0, 0, topreadoutZ, DetectorSizeX / 2.0, DetectorSizeY / 2.0, StripThickness / 2.0); // top glass
    topreadout.SetBoundaryPotential(0);
    geo.AddSolid(&topreadout, &copper);

    SolidBox bottomreadout(0, 0, bottomreadoutZ, DetectorSizeX / 2.0, DetectorSizeY / 2.0, StripThickness / 2.0); // top glass
    bottomreadout.SetBoundaryPotential(0);
    geo.AddSolid(&bottomreadout, &copper);

    ComponentNeBem3d nebem;
    nebem.SetGeometry(&geo);
    nebem.SetNumberOfThreads(16);
    nebem.SetTargetElementSize(0.0002); // 2 micron
    nebem.SetMinMaxNumberOfElements(5, 20);
    // nebem.SetFastVolOptions(1, 1, 1);
    nebem.UseSVDInversion();
    // nebem.UseLUInversion();
    nebem.Initialise();

    // Section to view Geometry and field profile

    // Plot device geometry in 2D
    // ViewGeometry geomView2d;
    // geomView2d.SetGeometry(&geo);
    // geomView2d.SetArea(-5, -5, -1, 5, 5, 1);
    // geomView2d.SetPlane(0, 1, 0, 0, 0, 0.0);
    // geomView2d.Plot2d();

    // Plot device geometry in 3D
    // ViewGeometry geomView3d;
    // geomView3d.SetGeometry(&geo);
    // geomView3d.Plot();

    // ViewField PotProfileView;
    // PotProfileView.SetComponent(&nebem);
    // PotProfileView.PlotProfile(0.0, 0.0, bottomreadoutZ - StripThickness / 2.0, 0.0, 0.0, topreadoutZ + StripThickness / 2.0);

    // ViewField EfieldProfileView;
    // EfieldProfileView.SetComponent(&nebem);
    // EfieldProfileView.PlotProfile(0.0, 0.0, bottomreadoutZ - StripThickness / 2.0, 0.0, 0.0, topreadoutZ + StripThickness / 2.0, "e");

    ViewField PotProfileView1;
    PotProfileView1.SetComponent(&nebem);
    PotProfileView1.PlotProfile(-DetectorSizeX / 2.0, 0.0, 0.0, DetectorSizeX / 2.0, 0.0, 0);

    ViewField EfieldProfileView1;
    EfieldProfileView1.SetComponent(&nebem);
    EfieldProfileView1.PlotProfile(-DetectorSizeX / 2.0, 0.0, 0.0, DetectorSizeX / 2.0, 0.0, 0, "e");

    // // // Plot potential contour
    // ViewField fieldView;
    // fieldView.SetComponent(&nebem);
    // fieldView.SetArea(-DetectorSizeX, -DetectorSizeY, bottomreadoutZ - StripThickness / 2.0, DetectorSizeX, DetectorSizeY, topreadoutZ + StripThickness / 2.0);
    // fieldView.SetPlane(0, 0, 1, 0, 0, 0.0);
    // fieldView.PlotContour("e");

    fs::create_directory("Data");
    std::ofstream efield1;
    std::string filename1 = "Data/efield_center_z.txt";
    efield1.open(filename1);

    std::ofstream efield2;
    std::string filename2 = "Data/efield_center_x_z0.txt";
    efield2.open(filename2);

    const double z1 = GasGapCenterZ + GasGapThickness / 2.0 - 0.001; // 10 micron below the top glass bottom surface
    std::ofstream efield3;
    std::string filename3 = "Data/efield_center_x_z1.txt";
    efield3.open(filename3);

    const double z2 = GasGapCenterZ - GasGapThickness / 2.0 + 0.001; // 10 micron above the bottom glass top surface
    std::ofstream efield4;
    std::string filename4 = "Data/efield_center_x_z2.txt";
    efield4.open(filename4);
    std::cout << "z1: " << z1 << "  z2: " << z2 << std::endl;

    Medium *medium = nullptr; 

    double ex = 0., ey = 0., ez = 0., e = 0., v = 0.;
    int status = 0;

    // through Z axis
    for (double i = bottomreadoutZ; i <= topreadoutZ; i += 0.0002)
    {
        nebem.ElectricField(0, 0, i, ex, ey, ez, v, medium, status);
        e = sqrt((ex * ex + ey * ey + ez * ez));
        efield1 << i << " " << e << " " << v << "\n";
    }
    efield1.close();
    std::printf("Calculated Field and Potential along z-axis.\n");

    // through x axis
    // at z=0
    for (double i = -DetectorSizeX / 2.0; i <= DetectorSizeX / 2.0; i += 0.003) // 1000 data points
    {
        nebem.ElectricField(i, 0, 0, ex, ey, ez, v, medium, status);
        e = sqrt((ex * ex + ey * ey + ez * ez));
        efield2 << i << " " << e << " " << v << "\n";
    }
    efield2.close();
    std::printf("Calculated Field and Potential along x-axis at z=0.\n");

    // at z=Z1
    for (double i = -DetectorSizeX / 2.0; i <= DetectorSizeX / 2.0; i += 0.003) // 1000 data points
    {
        nebem.ElectricField(i, 0, z1, ex, ey, ez, v, medium, status);
        e = sqrt((ex * ex + ey * ey + ez * ez));
        efield3 << i << " " << e << " " << v << "\n";
    }
    efield3.close();
    std::printf("Calculated Field and Potential along x-axis at z=z1.\n");

    // at z=Z2
    for (double i = -DetectorSizeX / 2.0; i <= DetectorSizeX / 2.0; i += 0.003) // 1000 data points
    {
        nebem.ElectricField(i, 0, z2, ex, ey, ez, v, medium, status);
        e = sqrt((ex * ex + ey * ey + ez * ez));
        efield4 << i << " " << e << " " << v << "\n";
    }
    efield4.close();
    std::printf("Calculated Field and Potential along x-axis at z=z2.\n");

    app.Run(true);
}
