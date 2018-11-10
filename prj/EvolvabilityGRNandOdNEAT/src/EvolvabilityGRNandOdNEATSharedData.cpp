/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */



#include "EvolvabilityGRNandOdNEAT/include/EvolvabilityGRNandOdNEATSharedData.h"

double EvolvabilityGRNandOdNEATSharedData::gSigmaMin = 0.0;
double EvolvabilityGRNandOdNEATSharedData::gProbaMutation = 0.0;
double EvolvabilityGRNandOdNEATSharedData::gUpdateSigmaStep = 0.0;
double EvolvabilityGRNandOdNEATSharedData::gSigmaRef = 0.0; // reference value of sigma
double EvolvabilityGRNandOdNEATSharedData::gSigmaMax = 0.0; // maximal value of sigma
int EvolvabilityGRNandOdNEATSharedData::gEvaluationTime = 0; // how long a controller will be evaluated on a robot

bool EvolvabilityGRNandOdNEATSharedData::gSynchronization = true;

bool EvolvabilityGRNandOdNEATSharedData::gEnergyRequestOutput = 1;

double EvolvabilityGRNandOdNEATSharedData::gMonitorPositions;

bool EvolvabilityGRNandOdNEATSharedData::gPropertiesLoaded = false; // global variable local to file -- TODO: move specific properties loader in dedicated WorldObserver



bool EvolvabilityGRNandOdNEATSharedData::gSnapshots = true; // take snapshots
int EvolvabilityGRNandOdNEATSharedData::gSnapshotsFrequency = 50; // every N generations

int EvolvabilityGRNandOdNEATSharedData::gControllerType = -1; // cf. header for description

bool EvolvabilityGRNandOdNEATSharedData::gLimitGenomeTransmission = false; // default: do not limit.
int EvolvabilityGRNandOdNEATSharedData::gMaxNbGenomeTransmission = 65535; // default: arbitrarily set to 65535.

int EvolvabilityGRNandOdNEATSharedData::gSelectionMethod = 0; // default: random selection


bool EvolvabilityGRNandOdNEATSharedData::gLogGenome = false;

double EvolvabilityGRNandOdNEATSharedData::gIndividualMutationRate = 1.0;

int EvolvabilityGRNandOdNEATSharedData::gMutationOperator = 1; // 0: uniform, 1: gaussian

int EvolvabilityGRNandOdNEATSharedData::gNbRegulatory = 0; //1

bool EvolvabilityGRNandOdNEATSharedData::gWithBias = false;

double EvolvabilityGRNandOdNEATSharedData::gCrossoverProb = 0.0;
double EvolvabilityGRNandOdNEATSharedData::gMutateProb = 0.8;

double EvolvabilityGRNandOdNEATSharedData::gSelPressure = 1.0;
int EvolvabilityGRNandOdNEATSharedData::gFitness = 0;
int EvolvabilityGRNandOdNEATSharedData::gPopSize = 10;

double EvolvabilityGRNandOdNEATSharedData::gCommunicationRange = 30;
bool EvolvabilityGRNandOdNEATSharedData::gCommunicationOnRadius = false;
bool   EvolvabilityGRNandOdNEATSharedData::gDoMeasureDiv = true; //false;

int EvolvabilityGRNandOdNEATSharedData::gMatingOperator = 0;
int EvolvabilityGRNandOdNEATSharedData::gBroadcastTime = 1;
int EvolvabilityGRNandOdNEATSharedData::gMaturationTime = 5;

bool EvolvabilityGRNandOdNEATSharedData::gIsCentralized;

std::string EvolvabilityGRNandOdNEATSharedData::gExpName;
int EvolvabilityGRNandOdNEATSharedData::freqMeasureBehav = 25;
int EvolvabilityGRNandOdNEATSharedData::gNbInputsBehavior = 50;
std::vector<std::vector<double> > EvolvabilityGRNandOdNEATSharedData::gInputsBehavior; // = std::vector<std::vector<double> >() ;

void EvolvabilityGRNandOdNEATSharedData::initInputsBehavior(int n,int in)
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
        if(EvolvabilityGRNandOdNEATSharedData::gControllerType == 0)
        {
            //inputs.push_back(1.0);
        }
        if(EvolvabilityGRNandOdNEATSharedData::gControllerType == 1)
        {
            inputs.push_back(1.0);
        }
        EvolvabilityGRNandOdNEATSharedData::gInputsBehavior.push_back(inputs);
    }
}
