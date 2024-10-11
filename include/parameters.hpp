#pragma once

#include <cmath>
#include "tools/toml.hpp"


// Definable parameters
// ********************

const toml::table config_sim = toml::parse_file("./simulation_params.toml");

const std::string models_path = config_sim.at_path("models_path").as_string()->get();

// Times
const double worldTimeStep = 0.01;
const double timeStep = 0.1;
const double singleSupportDuration = config_sim.at_path("singleSupportDuration").as_floating_point()->get();
const double doubleSupportDuration = config_sim.at_path("doubleSupportDuration").as_floating_point()->get();
const double predictionTime = 1.0;
const int prev = 20;


const double comTargetHeight = config_sim.at_path("comTargetHeight").as_floating_point()->get();
// const double kSwingFoot = 0.05;

// Fixed step parameters
//const double stepx = 0.06;
//const double stepy = 0.1;

// Constraints
const double thetaMax = 0.30;
const double footConstraintSquareWidth = 0.1;

// Used in the code
// ****************

const double omega = sqrt(9.81/comTargetHeight);
const int N = (int)round(predictionTime/timeStep);
const int S = (int)round(singleSupportDuration/timeStep);
const int D = (int)round(doubleSupportDuration/timeStep);
const int M = 4; //ceil(N/(S+D));
const int timesRatio = (int)(timeStep/worldTimeStep); //ceil(N/(S+D));
const int singleSupportSamples = (int)round(singleSupportDuration/worldTimeStep);
const int doubleSupportSamples = (int)round(doubleSupportDuration/worldTimeStep);

//// IS-MPC parameters
//// *****************
//const double deltaXMax = 0.25;
//const double deltaYIn = 0.15;
//const double deltaYOut = 0.28;
//// Cost function weights
//const double qZd = 1;
//const double qVx = 0;
//const double qVy = 0;
//const double qZ = 0;
//const double qF = 10000;//10000000000;


const double robot_mass = 38.0549;
