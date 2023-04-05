#pragma once

#ifndef CONFIG_H
#include <array>
#include <string>
#define CONFIG_H
struct PlannerConfiguration {
  // Basic
  int numVantages = 18;         //
  double mapPitch = 1.0;        // meters
  std::string outputDir = "./"; //
  bool drawAnimation = false;

  // Lander
  double landingSiteX = 100; // meters
  double landingSiteY = 100; // meters
  double landerHeight = 1.0; // meters

  // Rover
  double roverHeight = 1.0;                  // meters
  double roverSpeed = 0.01;                  // m/s
  double roverFOV = 180;                     // degrees
  double roverLateralSlopeLimit = 12.0;      // degrees
  double roverLongitudinalSlopeLimit = 20.0; // degrees
  double roverPointTurnSlopeLimit = 5.0;     // degrees

  // Advanced
  std::string meshFile = "../meshes/lmp-interior.ply"; //

  int numProbes = 10000;     //
  int numCandidates = 10000; //

  double maxVisAngle = 55;            // degrees
  double maxVisRange = 300;           // meters
  double minVantageSeparation = 20.0; // meters

  std::array<double, 10> visMultipliers = {1.0, 1.0, 0.8, 0.5, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1}; //

  double distCostMultiplier = 1.0;   //
  double slopeCostMultiplier = 16.0; //
  double turnCostMultiplier = 1.0;   //
  double heuristicMultiplier = 1.0;  //
};
void parseConfigFile(PlannerConfiguration &config, const std::string &fileName);
void parseCommandLine(PlannerConfiguration &config, int argc, char *argv[]);

#endif // CONFIG_H