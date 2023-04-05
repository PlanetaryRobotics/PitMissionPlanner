#include "config.h"
#include "argh.h"
#include "toml.hpp"
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <array>
#include <filesystem>

void parseConfigFile(PlannerConfiguration &config, const std::string& fileName) {
    const std::string fileDir = std::filesystem::path(fileName).parent_path();

    const toml::table tbl = toml::parse_file(fileName);
    config.numVantages   = tbl["Basic"]["NumVantages"].value_or(config.numVantages);
    config.mapPitch      = tbl["Basic"]["MapPitch"].value_or(config.mapPitch);
    config.outputDir     = tbl["Basic"]["OutputDir"].value_or(config.outputDir) + "/";
    config.drawAnimation = tbl["Basic"]["DrawAnimation"].value_or(config.drawAnimation);

    config.landingSiteX = tbl["Lander"]["LocationXY"][0].value_or(config.landingSiteX);
    config.landingSiteY = tbl["Lander"]["LocationXY"][1].value_or(config.landingSiteY);
    config.landerHeight = tbl["Lander"]["Height"].value_or(config.landerHeight);

    config.roverHeight                 = tbl["Rover"]["Height"].value_or(config.roverHeight);
    config.roverSpeed                  = tbl["Rover"]["Speed"].value_or(config.roverSpeed);
    config.roverFOV                    = tbl["Rover"]["FOV"].value_or(config.roverFOV);
    config.roverLateralSlopeLimit      = tbl["Rover"]["LateralSlopeLimit"].value_or(config.roverLateralSlopeLimit);
    config.roverLongitudinalSlopeLimit = tbl["Rover"]["LongitudinalSlopeLimit"].value_or(config.roverLongitudinalSlopeLimit);
    config.roverPointTurnSlopeLimit    = tbl["Rover"]["PointTurnSlopeLimit"].value_or(config.roverPointTurnSlopeLimit);

    config.meshFile             = fileDir + "/" + tbl["Advanced"]["MeshFile"].value_or(config.meshFile);
    config.numProbes            = tbl["Advanced"]["NumProbes"].value_or(config.numProbes);
    config.numCandidates        = tbl["Advanced"]["NumCandidates"].value_or(config.numCandidates);
    config.maxVisAngle          = tbl["Advanced"]["MaxVisAngle"].value_or(config.maxVisAngle);
    config.maxVisRange          = tbl["Advanced"]["MaxVisRange"].value_or(config.maxVisRange);
    config.minVantageSeparation = tbl["Advanced"]["MinVantageSeparation"].value_or(config.minVantageSeparation);

    config.distCostMultiplier  = tbl["Advanced"]["DistCostMultiplier"].value_or(config.distCostMultiplier);
    config.slopeCostMultiplier = tbl["Advanced"]["SlopeCostMultiplier"].value_or(config.slopeCostMultiplier);
    config.turnCostMultiplier  = tbl["Advanced"]["TurnCostMultiplier"].value_or(config.turnCostMultiplier);
    config.heuristicMultiplier = tbl["Advanced"]["HeuristicMultiplier"].value_or(config.heuristicMultiplier);
}
void parseCommandLine(PlannerConfiguration &config, int argc, char* argv[]) {
    argh::parser cmdl(argc, argv);

    if( cmdl[{ "-h", "--help", "--Help" }] ) {
        fmt::print("usage: planranger [options] [config.toml]\n");
        fmt::print("  options:\n");
        fmt::print("    [Basic]\n");
        fmt::print("      -h,   --Help                    Show this help message.\n");
        fmt::print("      -nv,  --NumVantages             The number of vantages to include in the route.\n");
        fmt::print("      -p,   --MapPitch                The resolution (meters) to use for generating maps.\n");
        fmt::print("      -o,   --OutputDir               A directory to store the planner output.\n");
        fmt::print("   -anim,   --DrawAnimation           Generate an animation of the rover's route.\n");
        fmt::print("\n");

        fmt::print("    [Lander]\n");
        fmt::print("      -x,   --LandingSiteX            The x-coordinate of the landing site.\n");
        fmt::print("      -y,   --LandingSiteY            The y-coordinate of the landing site.\n");
        fmt::print("      -lh,  --LanderHeight            The height of the comms antenna on the lander.\n");
        fmt::print("\n");

        fmt::print("    [Rover]\n");
        fmt::print("      -rh,   --RoverHeight            The height of the comms antenna on the rover.\n");
        fmt::print("      -rs,   --RoverSpeed             The rover's driving speed in m/s .\n");
        fmt::print("      -fov,  --RoverFOV               The rover camera's field-of-view in degrees.\n");
        fmt::print("\n");

        fmt::print("      -lat,  --LateralSlopeLimit      The rover's max lateral slope in degrees.\n");
        fmt::print("      -lon,  --LongitudinalSlopeLimit The rover's max longitudinal slope in degrees.\n");
        fmt::print("      -turn, --PointTurnSlopeLimit    The steepest slope (in degrees) on which the rover can turn.\n");
        fmt::print("\n");

        fmt::print("    [Advanced]\n");
        fmt::print("      -mesh,   --MeshFile             A *.ply file representing the terrain.\n");
        fmt::print("      -np,     --NumProbes            The number of visibility probes to use for vantage selection.\n");
        fmt::print("      -nc,     --NumCandidates        The number of vantage candidates to consider for vantage selection.\n");
        fmt::print("      -va,     --MaxVisAngle          The maximum incidence angle (in degrees) that is considered a valid view hit.\n");
        fmt::print("      -vr,     --MaxVisRange          The maximum distance (in meters) beyond which a view hit is not counted.\n");
        fmt::print("      -minsep, --MinVantageSeparation The minimum distance between any two vantages.\n");
        fmt::print("\n");

        fmt::print("      -dmult,  --DistCostMultiplier   The unitless multiplier used to weight the distance cost term in the planning cost function.\n");
        fmt::print("      -smult,  --SlopeCostMultiplier  The unitless multiplier used to weight the slope cost term in the planning cost function.\n");
        fmt::print("      -tmult,  --TurnCostMultiplier   The unitless multiplier used to weight the turning cost term in the planning cost function.\n");
        fmt::print("      -hmult,  --HeuristicMultiplier  The unitless multiplier used to inflate the A* heuristic for faster, less-optimal planning.\n");
        fmt::print("\n");

        exit(0);
    }

    // If the first positional argument is a toml file, use it to configure the planner.
    if( !cmdl[1].empty() && cmdl[1].find(".toml") != std::string::npos ) {
        parseConfigFile(config, cmdl[1]);
    }

    // Basic
    cmdl({"NumVantages",      "nv"}, config.numVantages)                 >> config.numVantages;
    cmdl({"MapPitch",          "p"}, config.mapPitch)                    >> config.mapPitch;
    cmdl({"OutputDir",         "o"}, config.outputDir)                   >> config.outputDir;
    config.outputDir += "/";
    config.drawAnimation |= cmdl[{"DrawAnimation", "anim"}];

    // Lander
    cmdl({"LandingSiteX",      "x"}, config.landingSiteX)                >> config.landingSiteX;
    cmdl({"LandingSiteY",      "y"}, config.landingSiteY)                >> config.landingSiteY;
    cmdl({"LanderHeight",     "lh"}, config.landerHeight)                >> config.landerHeight;

    // Rover
    cmdl({"RoverHeight",      "rh"}, config.roverHeight)                 >> config.roverHeight;
    cmdl({"RoverSpeed",       "rs"}, config.roverSpeed)                  >> config.roverSpeed;
    cmdl({"RoverFOV",        "fov"}, config.roverFOV)                    >> config.roverFOV;

    cmdl({"LateralSlopeLimit",       "lat"}, config.roverLateralSlopeLimit)      >> config.roverLateralSlopeLimit;
    cmdl({"LongitudinalSlopeLimit",  "lon"}, config.roverLongitudinalSlopeLimit) >> config.roverLongitudinalSlopeLimit;
    cmdl({"PointTurnSlopeLimit",    "turn"}, config.roverPointTurnSlopeLimit)    >> config.roverPointTurnSlopeLimit;

    // Advanced
    cmdl({"MeshFile",               "mesh"}, config.meshFile)             >> config.meshFile;
    cmdl({"NumProbes",                "np"}, config.numProbes)            >> config.numProbes;
    cmdl({"NumCandidates",            "nc"}, config.numCandidates)        >> config.numCandidates;
    cmdl({"MaxVisAngle",              "va"}, config.maxVisAngle)          >> config.maxVisAngle;
    cmdl({"MaxVisRange",              "vr"}, config.maxVisRange)          >> config.maxVisRange;
    cmdl({"MinVantageSeparation", "minsep"}, config.minVantageSeparation) >> config.minVantageSeparation;

    cmdl({"DistCostMultiplier",       "dmult"}, config.distCostMultiplier)  >> config.distCostMultiplier;
    cmdl({"SlopeCostMultiplier",      "smult"}, config.slopeCostMultiplier) >> config.slopeCostMultiplier;
    cmdl({"TurnCostMultiplier",       "tmult"}, config.turnCostMultiplier)  >> config.turnCostMultiplier;
    cmdl({"HeuristicMultiplier",      "hmult"}, config.heuristicMultiplier) >> config.heuristicMultiplier;
}