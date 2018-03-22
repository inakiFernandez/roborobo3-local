/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */



#include "Increment/include/IncrementSharedData.h"

double IncrementSharedData::gSigmaRef = 0.0; // reference value of sigma
int IncrementSharedData::gPopulationSize = 1; // size of population per robot
bool IncrementSharedData::gClearPopulation; // empty population after replacement?
bool IncrementSharedData::gStoreOwn; // store own genome in local population?
int IncrementSharedData::gEvaluationTime = 0; // how long a controller will be evaluated on a robot

int IncrementSharedData::gSelectionMethod = 0; // Random by default

// global variable local to file -- TODO: move specific properties loader in dedicated WorldObserver
bool IncrementSharedData::gPropertiesLoaded = false;

int IncrementSharedData::gControllerType = 0; //MLP, Perceptron, Elman
int IncrementSharedData::gNbHiddenLayers = 1;
int IncrementSharedData::gNbNeuronsPerHiddenLayer = 5;
int IncrementSharedData::gNeuronWeightRange = 800;
bool IncrementSharedData::gWithBias = false;

bool IncrementSharedData::gCommunicationBySensors = true;
int IncrementSharedData::gFitness = -1;

std::string IncrementSharedData::gOutGenomeFile = "";
bool IncrementSharedData::gSaveGenome = true;
