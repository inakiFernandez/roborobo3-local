/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */



#include "SocialBehavior/include/SocialBehaviorSharedData.h"

double SocialBehaviorSharedData::gSigmaMin = 0.0;
double SocialBehaviorSharedData::gProbaMutation = 0.0;
double SocialBehaviorSharedData::gUpdateSigmaStep = 0.0;
double SocialBehaviorSharedData::gSigmaRef = 0.0; // reference value of sigma
double SocialBehaviorSharedData::gSigmaMax = 0.0; // maximal value of sigma
int SocialBehaviorSharedData::gEvaluationTime = 0; // how long a controller will be evaluated on a robot

bool SocialBehaviorSharedData::gSynchronization = true;

bool SocialBehaviorSharedData::gEnergyRequestOutput = 1;

double SocialBehaviorSharedData::gMonitorPositions;

bool SocialBehaviorSharedData::gPropertiesLoaded = false; // global variable local to file -- TODO: move specific properties loader in dedicated WorldObserver

int SocialBehaviorSharedData::gNbHiddenLayers = 1;
int SocialBehaviorSharedData::gNbNeuronsPerHiddenLayer = 5;
int SocialBehaviorSharedData::gNeuronWeightRange = 800;

bool SocialBehaviorSharedData::gSnapshots = true; // take snapshots
int SocialBehaviorSharedData::gSnapshotsFrequency = 50; // every N generations

int SocialBehaviorSharedData::gControllerType = -1; // cf. header for description

bool SocialBehaviorSharedData::gLimitGenomeTransmission = false; // default: do not limit.
int SocialBehaviorSharedData::gMaxNbGenomeTransmission = 65535; // default: arbitrarily set to 65535.

int SocialBehaviorSharedData::gSelectionMethod = 0; // default: random selection

int SocialBehaviorSharedData::gNotListeningStateDelay = 0;    // -1: infinite ; 0: no delay ; >0: delay
int SocialBehaviorSharedData::gListeningStateDelay = -1;      // -1: infinite ; 0: no delay ; >0: delay (ignored if gNotListeningStateDelay=-1)

bool SocialBehaviorSharedData::gLogGenome = false;

double SocialBehaviorSharedData::gIndividualMutationRate = 1.0;

int SocialBehaviorSharedData::gMutationOperator = 1; // 0: uniform, 1: gaussian

double SocialBehaviorSharedData::gSigma = 0.01; // 0.01 is just some random value.
