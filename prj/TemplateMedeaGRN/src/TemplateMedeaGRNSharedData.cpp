/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */



#include "TemplateMedeaGRN/include/TemplateMedeaGRNSharedData.h"

double TemplateMedeaGRNSharedData::gSigmaMin = 0.0;
double TemplateMedeaGRNSharedData::gProbaMutation = 0.0;
double TemplateMedeaGRNSharedData::gUpdateSigmaStep = 0.0;
double TemplateMedeaGRNSharedData::gSigmaRef = 0.0; // reference value of sigma
double TemplateMedeaGRNSharedData::gSigmaMax = 0.0; // maximal value of sigma
int TemplateMedeaGRNSharedData::gEvaluationTime = 0; // how long a controller will be evaluated on a robot

bool TemplateMedeaGRNSharedData::gSynchronization = true;

bool TemplateMedeaGRNSharedData::gEnergyRequestOutput = 1;

double TemplateMedeaGRNSharedData::gMonitorPositions;

bool TemplateMedeaGRNSharedData::gPropertiesLoaded = false; // global variable local to file -- TODO: move specific properties loader in dedicated WorldObserver



bool TemplateMedeaGRNSharedData::gSnapshots = true; // take snapshots
int TemplateMedeaGRNSharedData::gSnapshotsFrequency = 50; // every N generations

int TemplateMedeaGRNSharedData::gControllerType = -1; // cf. header for description

bool TemplateMedeaGRNSharedData::gLimitGenomeTransmission = false; // default: do not limit.
int TemplateMedeaGRNSharedData::gMaxNbGenomeTransmission = 65535; // default: arbitrarily set to 65535.

int TemplateMedeaGRNSharedData::gSelectionMethod = 0; // default: random selection


bool TemplateMedeaGRNSharedData::gLogGenome = false;

double TemplateMedeaGRNSharedData::gIndividualMutationRate = 1.0;

int TemplateMedeaGRNSharedData::gMutationOperator = 1; // 0: uniform, 1: gaussian

int TemplateMedeaGRNSharedData::gNbRegulatory = 0; //1



double TemplateMedeaGRNSharedData::gCrossoverProb = 0.0;
double TemplateMedeaGRNSharedData::gMutateProb = 0.8;

double TemplateMedeaGRNSharedData::gSelPressure = 1.0;
int TemplateMedeaGRNSharedData::gFitness = 0;

double TemplateMedeaGRNSharedData::gCommunicationRange = 30;
bool TemplateMedeaGRNSharedData::gCommunicationOnRadius = false;

int TemplateMedeaGRNSharedData::gMatingOperator = 0;
int TemplateMedeaGRNSharedData::gBroadcastTime = 1;
