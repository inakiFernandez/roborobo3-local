/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */



#ifndef SOCIALBEHAVIORANNSHAREDDATA_H
#define SOCIALBEHAVIORANNSHAREDDATA_H

#include "Utilities/Misc.h"
#include <math.h>
class SocialBehaviorANNSharedData {
    
public:
    
    // -----
    
    static double gSigmaMin; //! used with gDynamicSigma defined to true
    static double gUpdateSigmaStep; //!step used in the drecrease or increas of the value of sigma
    static double gSigmaRef; //! reference value of sigma
    static double gSigmaMax; //! maximal value of sigma
    static double gProbaMutation; //! probability of transmitting the current genome mutated with sigma ref
    static int gEvaluationTime; //! theoretical duration of a generation (ie. maximum time a controller will be evaluated on a robot)
    static int gIteration; //! used by every class to know what is the current iteration step of roborobo
    static bool gSynchronization; //!If set to false, a robot will restart its controller as soon as it has no more energy. If set to true, the robot without energy will wait and reload its controller at the same time as every other robots.
    
    static bool gEnergyRequestOutput; // does the robot can modulate its energy request (when being given some) ?
    
    static double gMonitorPositions; //! used in WorldObserver. Compute and log all necessary information for monitoring position and orientation wrt. center.
    
    static bool gPropertiesLoaded;
    
    
    static bool gSnapshots; // take snapshots
    static int gSnapshotsFrequency; // every N generations
    // controller type. Depending on the experiment it can be ANN( ex. 0: MLP, 1: Perceptron, 2: Elman)
    //or GRN, etc
    static int gControllerType;
    static bool gWithBias;

    static int gMaxNbGenomeTransmission; // max. nb of genome transmission per individual per lifetime.
    static bool gLimitGenomeTransmission; // enable/disable limiting transmission
    
    static int gSelectionMethod; // check PrjMateSelectionController::loadNewGenome() for authorized values.
    
    static int gNotListeningStateDelay; // -1: infinite ; 0: no delay ; >0: delay
    static int gListeningStateDelay;    // -1: infinite ; 0: no delay ; >0: delay (ignored if gNotListeningStateDelay=-1)
    
    static bool gLogGenome;
    
    static double gIndividualMutationRate;
    
    static int gMutationOperator;
    
    static double gSigma;

    static int gNbRegulatory;

    static double gCrossoverProb;
    static double gMutateProb;

    static double gSelPressure;
    static int gFitness;
    static int gPopSize;
    static double gCommunicationRange;
    static bool gCommunicationOnRadius;
    static bool gDoMeasureDiv;

    static int gMatingOperator;
    static int gBroadcastTime;
    static int gMaturationTime;
    // -----
    static bool gIsCentralized;

    static int freqMeasureBehav;

    static int gNbInputsBehavior;
    static std::vector<std::vector<double> > gInputsBehavior;
    static void initInputsBehavior(int n,int in);
    
    
};


#endif
