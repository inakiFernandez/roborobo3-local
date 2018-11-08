/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */



#include "SocialBehaviorANN/include/SocialBehaviorANNSharedData.h"

double SocialBehaviorANNSharedData::gSigmaMin = 0.0;
double SocialBehaviorANNSharedData::gProbaMutation = 0.0;
double SocialBehaviorANNSharedData::gUpdateSigmaStep = 0.0;
double SocialBehaviorANNSharedData::gSigmaRef = 0.0; // reference value of sigma
double SocialBehaviorANNSharedData::gSigmaMax = 0.0; // maximal value of sigma
int SocialBehaviorANNSharedData::gEvaluationTime = 0; // how long a controller will be evaluated on a robot

bool SocialBehaviorANNSharedData::gSynchronization = true;

bool SocialBehaviorANNSharedData::gEnergyRequestOutput = 1;

double SocialBehaviorANNSharedData::gMonitorPositions;

bool SocialBehaviorANNSharedData::gPropertiesLoaded = false; // global variable local to file -- TODO: move specific properties loader in dedicated WorldObserver



bool SocialBehaviorANNSharedData::gSnapshots = true; // take snapshots
int SocialBehaviorANNSharedData::gSnapshotsFrequency = 50; // every N generations

int SocialBehaviorANNSharedData::gControllerType = -1; // cf. header for description

bool SocialBehaviorANNSharedData::gLimitGenomeTransmission = false; // default: do not limit.
int SocialBehaviorANNSharedData::gMaxNbGenomeTransmission = 65535; // default: arbitrarily set to 65535.

int SocialBehaviorANNSharedData::gSelectionMethod = 0; // default: random selection


bool SocialBehaviorANNSharedData::gLogGenome = false;

double SocialBehaviorANNSharedData::gIndividualMutationRate = 1.0;

int SocialBehaviorANNSharedData::gMutationOperator = 1; // 0: uniform, 1: gaussian

int SocialBehaviorANNSharedData::gNbRegulatory = 0; //1

bool SocialBehaviorANNSharedData::gWithBias = false;

double SocialBehaviorANNSharedData::gCrossoverProb = 0.0;
double SocialBehaviorANNSharedData::gMutateProb = 0.8;

double SocialBehaviorANNSharedData::gSelPressure = 1.0;
int SocialBehaviorANNSharedData::gFitness = 0;
int SocialBehaviorANNSharedData::gPopSize = 10;

double SocialBehaviorANNSharedData::gCommunicationRange = 30;
bool SocialBehaviorANNSharedData::gCommunicationOnRadius = false;
bool   SocialBehaviorANNSharedData::gDoMeasureDiv = true; //false;

int SocialBehaviorANNSharedData::gMatingOperator = 0;
int SocialBehaviorANNSharedData::gBroadcastTime = 1;
int SocialBehaviorANNSharedData::gMaturationTime = 5;

bool SocialBehaviorANNSharedData::gIsCentralized;


int SocialBehaviorANNSharedData::freqMeasureBehav = 25;
int SocialBehaviorANNSharedData::gNbInputsBehavior = 50;
std::vector<std::vector<double> > SocialBehaviorANNSharedData::gInputsBehavior; // = std::vector<std::vector<double> >() ;

void SocialBehaviorANNSharedData::initInputsBehavior(int n,int in)
{
  std::vector<double> inputs;

    for(int i=0; i < n; i++)
    {
        inputs.clear();

        for(int j = 0; j < in; j++)
        {
            float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
            inputs.push_back(r);
        }
        inputs.push_back(1.0);
        SocialBehaviorANNSharedData::gInputsBehavior.push_back(inputs);
    }
}
