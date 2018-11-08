/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */

#include "Observers/AgentObserver.h"
#include "Observers/WorldObserver.h"
#include "EvolvabilityGRNandOdNEAT/include/EvolvabilityGRNandOdNEATWorldObserver.h"
#include "EvolvabilityGRNandOdNEAT/include/EvolvabilityGRNandOdNEATController.h"
#include "World/World.h"
#include <fstream>
EvolvabilityGRNandOdNEATWorldObserver::EvolvabilityGRNandOdNEATWorldObserver( World* world ) : WorldObserver( world )
{
    _world = world;
    
    // ==== loading project-specific properties
    
    gProperties.checkAndGetPropertyValue("gSigmaRef",&EvolvabilityGRNandOdNEATSharedData::gSigmaRef,true);
    gProperties.checkAndGetPropertyValue("gSigmaMin",&EvolvabilityGRNandOdNEATSharedData::gSigmaMin,true);
    gProperties.checkAndGetPropertyValue("gSigmaMax",&EvolvabilityGRNandOdNEATSharedData::gSigmaMax,true);
    
    gProperties.checkAndGetPropertyValue("gUpdateSigmaStep",&EvolvabilityGRNandOdNEATSharedData::gUpdateSigmaStep,true);
    gProperties.checkAndGetPropertyValue("gEvaluationTime",&EvolvabilityGRNandOdNEATSharedData::gEvaluationTime,true);
    gProperties.checkAndGetPropertyValue("gSynchronization",&EvolvabilityGRNandOdNEATSharedData::gSynchronization,true);
    
    gProperties.checkAndGetPropertyValue("gEnergyRequestOutput",&EvolvabilityGRNandOdNEATSharedData::gEnergyRequestOutput,false);
    
    gProperties.checkAndGetPropertyValue("gMonitorPositions",&EvolvabilityGRNandOdNEATSharedData::gMonitorPositions,true);
    
    
    gProperties.checkAndGetPropertyValue("gSnapshots",&EvolvabilityGRNandOdNEATSharedData::gSnapshots,false);
    gProperties.checkAndGetPropertyValue("gSnapshotsFrequency",&EvolvabilityGRNandOdNEATSharedData::gSnapshotsFrequency,false);
    
    gProperties.checkAndGetPropertyValue("gControllerType",&EvolvabilityGRNandOdNEATSharedData::gControllerType,true);
    
    gProperties.checkAndGetPropertyValue("gMaxNbGenomeTransmission",&EvolvabilityGRNandOdNEATSharedData::gMaxNbGenomeTransmission,true);
    gProperties.checkAndGetPropertyValue("gLimitGenomeTransmission",&EvolvabilityGRNandOdNEATSharedData::gLimitGenomeTransmission,true);
    gProperties.checkAndGetPropertyValue("gSelectionMethod",&EvolvabilityGRNandOdNEATSharedData::gSelectionMethod,true);
    
    
    gProperties.checkAndGetPropertyValue("gLogGenome",&EvolvabilityGRNandOdNEATSharedData::gLogGenome,false);
    
    gProperties.checkAndGetPropertyValue("gIndividualMutationRate",&EvolvabilityGRNandOdNEATSharedData::gIndividualMutationRate,false);

    gProperties.checkAndGetPropertyValue("gMutationOperator",&EvolvabilityGRNandOdNEATSharedData::gMutationOperator,false);
    gProperties.checkAndGetPropertyValue("gFitness",&EvolvabilityGRNandOdNEATSharedData::gFitness,true);

    gProperties.checkAndGetPropertyValue("gSelPressure",&EvolvabilityGRNandOdNEATSharedData::gSelPressure,true);
    gProperties.checkAndGetPropertyValue("gPopSize",&EvolvabilityGRNandOdNEATSharedData::gPopSize,true);


    gProperties.checkAndGetPropertyValue("gMutateProb",&EvolvabilityGRNandOdNEATSharedData::gMutateProb,true);
    gProperties.checkAndGetPropertyValue("gCrossoverProb",&EvolvabilityGRNandOdNEATSharedData::gCrossoverProb,true);
    gProperties.checkAndGetPropertyValue("gNbRegulatory",&EvolvabilityGRNandOdNEATSharedData::gNbRegulatory,true);


    gProperties.checkAndGetPropertyValue("gCommunicationRange",&EvolvabilityGRNandOdNEATSharedData::gCommunicationRange,true);
    gProperties.checkAndGetPropertyValue("gCommunicationOnRadius",&EvolvabilityGRNandOdNEATSharedData::gCommunicationOnRadius,true);

    gProperties.checkAndGetPropertyValue("gMatingOperator",&EvolvabilityGRNandOdNEATSharedData::gMatingOperator,true);
    gProperties.checkAndGetPropertyValue("gBroadcastTime",&EvolvabilityGRNandOdNEATSharedData::gBroadcastTime,true);    

    gProperties.checkAndGetPropertyValue("mutate_link_weights_prob",&Helper::mutateLinkWeightsProb,true);
    gProperties.checkAndGetPropertyValue("mutate_add_node_prob",&Helper::mutateAddNodeProb,true);
    gProperties.checkAndGetPropertyValue("mutate_add_link_prob",&Helper::mutateAddLinkProb,true);
    gProperties.checkAndGetPropertyValue("recur_only_prob",&Helper::recurOnlyProb,true);
    gProperties.checkAndGetPropertyValue("newstructure_tries",&Helper::newStructureTries,true);


    //TODO read other properties (cf shareddata)
    if ( !gRadioNetwork)
    {
        std::cout << "Error : gRadioNetwork must be true." << std::endl;
        exit(-1);
    }
    
    // * iteration and generation counters
    _generationItCount = -1;
    _generationCount = -1;
    int nbIn = -1;
    switch (EvolvabilityGRNandOdNEATSharedData::gFitness) {
    case 0:
        nbIn = 8;
        break;
    case 1:
        nbIn = 16;
        break;
    case 2:
        nbIn = 16;
        break;
    default:
        std::cerr << "[ERROR] Wrong or non implemented task" << std::endl;
        exit(-1);
        break;
    }
EvolvabilityGRNandOdNEATSharedData::initInputsBehavior(EvolvabilityGRNandOdNEATSharedData::gNbInputsBehavior, nbIn);
}

EvolvabilityGRNandOdNEATWorldObserver::~EvolvabilityGRNandOdNEATWorldObserver()
{
    // nothing to do.
}

void EvolvabilityGRNandOdNEATWorldObserver::reset()
{
}

void EvolvabilityGRNandOdNEATWorldObserver::step()
{
    _generationItCount++;
    
    if( _generationItCount == EvolvabilityGRNandOdNEATSharedData::gEvaluationTime+1 ) // switch to next generation.
    {
        // update iterations and generations counters
        _generationItCount = 0;
        _generationCount++;
    }
    //Local broadcast with distance: compute inter-robot distance matrix
    if(EvolvabilityGRNandOdNEATSharedData::gCommunicationOnRadius)
    {
    _robotDistances = std::vector<std::vector<double> >(gNumberOfRobots, std::vector<double>(gNumberOfRobots,-1));
        for ( int i = 0 ; i != gNumberOfRobots ; i++ )
        {
            EvolvabilityGRNandOdNEATController* c = (dynamic_cast<EvolvabilityGRNandOdNEATController*>(gWorld->getRobot(i)->getController()));
            double x = c->getWorldModel()->getXReal();
            double y = c->getWorldModel()->getYReal();

            for ( int j = i ; j != gNumberOfRobots ; j++ )
            {
                EvolvabilityGRNandOdNEATController* c2 = (dynamic_cast<EvolvabilityGRNandOdNEATController*>(gWorld->getRobot(j)->getController()));
                double x2 = c2->getWorldModel()->getXReal();
                double y2 = c2->getWorldModel()->getYReal();
                _robotDistances[i][j] = sqrt((x-x2)*(x-x2) + (y-y2)*(y-y2)); //.push_back(sqrt((x-x2)*(x-x2) + (y-y2)*(y-y2)));
                _robotDistances[j][i] = _robotDistances[i][j];
                //std::cout << _robotDistances[i][j] << "      ";
            }
            //std::cout << std::endl;

        }
        //std::cout << EvolvabilityGRNandOdNEATSharedData::gCommunicationRange << std::endl << std::endl << std::endl;
    }
    updateMonitoring();
    if (gWorld->getIterations() == (gMaxIt - 1))
    {
        std::ofstream log_file( //todo log in experiment folder
                "logOffspring.txt", std::ios_base::out | std::ios_base::app );
        int countGeneration = 0;
        for(auto it = _offspringStats.begin(); it != _offspringStats.end();it++)
        {
            log_file << "Generation Starts\n" << countGeneration++ << std::endl;
            for(auto it2 = (*it).begin(); it2 != (*it).end(); it2++ )
            {
                log_file <<  (*it2).second.toString()<< std::endl;
            }
            log_file <<  "Generation Ends\n\n\n\n\n\n\n\n" << std::endl;
        }
    }

    updateEnvironment();
    
}


void EvolvabilityGRNandOdNEATWorldObserver::updateEnvironment()
{
    // example: moving landmarks
    /*
    if ( gWorld->getIterations() % 2000 == 0 )
        for ( int i = 0 ; i != gLandmarks.size() ; i++ )
        {
            Point2d* position = new Point2d( 200+rand()%(gAreaWidth-400) , 200+rand()%(gAreaHeight-400) );
            gLandmarks[i].setPosition(*position);
        }
    */
}

void EvolvabilityGRNandOdNEATWorldObserver::updateMonitoring()
{
    // * Log at end of each generation

    if( //gWorld->getIterations() % EvolvabilityGRNandOdNEATSharedData::gEvaluationTime == 1 ||
            gWorld->getIterations() % EvolvabilityGRNandOdNEATSharedData::gEvaluationTime == EvolvabilityGRNandOdNEATSharedData::gEvaluationTime-1 ) // beginning(+1) *and* end of generation. ("==1" is required to monitor the outcome of the first iteration)
    {
        monitorPopulation();
    }
    //Simulate centralized algorithm by maintaining panmictic communication
    if(EvolvabilityGRNandOdNEATSharedData::gIsCentralized)
    {
        //Send around the end of the evaluation
        if((gWorld->getIterations() + 10) % EvolvabilityGRNandOdNEATSharedData::gEvaluationTime == 0)
        //if( _lifeIterationCount == EvolvabilityGRNandOdNEATSharedData::gEvaluationTime-10)
        {
            for (int i = 0 ; i != gNumberOfRobots ; i++ )
            {
                EvolvabilityGRNandOdNEATController* c1 = (dynamic_cast<EvolvabilityGRNandOdNEATController*>
                                          (gWorld->getRobot(i)->getController()));

                for (int j = 0 ; j != gNumberOfRobots ; j++ )
                {
                    EvolvabilityGRNandOdNEATController* c2 = (dynamic_cast<EvolvabilityGRNandOdNEATController*>
                                              (gWorld->getRobot(j)->getController()));
                    c2->storeGenomeHelper(c1->getCurrentGenome());//,c1->getCurrentId(), c1->getFitness(), c1->getBehavior());
                }
            }
        }
    }
    // * Every N generations, take a video (duration: one generation time)

    if ( EvolvabilityGRNandOdNEATSharedData::gSnapshots )
    {
        if ( ( gWorld->getIterations() ) % ( EvolvabilityGRNandOdNEATSharedData::gEvaluationTime * EvolvabilityGRNandOdNEATSharedData::gSnapshotsFrequency ) == 0 )
        {
            if ( gVerbose )
                std::cout << "[START] Video recording: generation #" << (gWorld->getIterations() / EvolvabilityGRNandOdNEATSharedData::gEvaluationTime ) << ".\n";
            gTrajectoryMonitorMode = 0;
            initTrajectoriesMonitor();
        }
        else
            if ( ( gWorld->getIterations() ) % ( EvolvabilityGRNandOdNEATSharedData::gEvaluationTime * EvolvabilityGRNandOdNEATSharedData::gSnapshotsFrequency ) == EvolvabilityGRNandOdNEATSharedData::gEvaluationTime - 1 )
            {
                if ( gVerbose )
                    std::cout << "[STOP]  Video recording: generation #" << (gWorld->getIterations() / EvolvabilityGRNandOdNEATSharedData::gEvaluationTime ) << ".\n";
                saveTrajectoryImage();
            }
    }    
}

void EvolvabilityGRNandOdNEATWorldObserver::monitorPopulation( bool localVerbose )
{
    // * monitoring: count number of active agents.
    
    int activeCount = 0;
    for ( int i = 0 ; i != gNumberOfRobots ; i++ )
    {
        if ( (dynamic_cast<EvolvabilityGRNandOdNEATController*>(gWorld->getRobot(i)->getController()))->getWorldModel()->isAlive() == true )
            activeCount++;
    }
    
    if ( gVerbose && localVerbose )
    {
        std::cout << "[gen:" << (gWorld->getIterations()/EvolvabilityGRNandOdNEATSharedData::gEvaluationTime) << ";it:" << gWorld->getIterations() << ";pop:" << activeCount << "]\n";
    }
    // Logging here
    double sumFitness = 0.0, sumAvgLocalPopFitness = 0.0;
    double sumItems = 0.0;
    double sumCollisions = 0.0;
    double gatheredGenomes = 0.0;
    double distance = 0.0;
    double nbUnits = 0.0;

    for ( int i = 0 ; i != gNumberOfRobots ; i++ )
    {

        EvolvabilityGRNandOdNEATController* c =(dynamic_cast<EvolvabilityGRNandOdNEATController*>(gWorld->getRobot(i)->getController()));
         sumFitness += c->getFitness();
         sumAvgLocalPopFitness += c-> getAvgPopFitness();
         sumItems += c->getCollectedItems();
         sumCollisions += c->getNbCollisions() - c->getCollectedItems();
         gatheredGenomes += c->getGenomesList().size();
         distance += c->getDSumTravelled();
         nbUnits += c->getNbUnits();
    }
    switch (EvolvabilityGRNandOdNEATSharedData::gFitness)
    {
        case 1:
        {
            sumFitness *= EvolvabilityGRNandOdNEATSharedData::gEvaluationTime;
            break;
        }
    }
    //Space-separated values on stdout:
    // Timestep AvgFitness LocalPopAvgFitnessAveragedOverRobots
    // GlobalDiversity AvgCollectedItems AvgCollisions AvgGatheredGenomes
    // AverageCoveredDistances AvgNumberOfComputingUnits
    std::cout << gWorld->getIterations() << " ";
    //<< (sumFitness  / gNumberOfRobots) / EvolvabilityGRNandOdNEATSharedData::gEvaluationTime

    // divided by two because each item gives 1 fitness point to both agents
    //std::cout << sumFitness / 2 << " " << sumAvgLocalPopFitness / 2 << std::endl;

    std::cout << sumFitness/(double) gNumberOfRobots <<
                 " " << sumAvgLocalPopFitness/(double) gNumberOfRobots;
    if(EvolvabilityGRNandOdNEATSharedData::gDoMeasureDiv)
    {
        std::cout << " " << computeGlobalDiversity();
    }
                        //" " << computeInterRobotDiversity();
    std::cout << " "
    << sumItems / gNumberOfRobots
    << " "
    << sumCollisions / gNumberOfRobots
    << " "
    << gatheredGenomes / gNumberOfRobots

                 //  " " << fullLaps <<
              //   " " << towardDoor0 <<
              //   " " << towardDoor1
    << "       "
    << distance / gNumberOfRobots
    << " "
    << nbUnits / gNumberOfRobots //nb neurones ou proteines (répartition des tailles de réseau par quartiles de fitness)
    << " ";
    // nbOffspring + fitness pour mesurer le kendall-tau-b pression de sélection

   /* if(EvolvabilityGRNandOdNEATSharedData::gDoMeasureDiv)
   {

        for ( int i = 0 ; i != gNumberOfRobots ; i++ )
        {

            EvolvabilityGRNandOdNEATController* c =(dynamic_cast<EvolvabilityGRNandOdNEATController*>(gWorld->getRobot(i)->getController()));
             //std::cout << c-> getDoorPassages() << " ";
             std::cout << " " << c->computeIntraRobotDiversity() << " ";
        }
    }*/
    std::cout << std::endl;
    _offspringStats.push_back(this->computeNbOffspringAndFitness());
    //std::cout //<< "[gen:"
    //        << (gWorld->getIterations()/EvolvabilityGRNandOdNEATSharedData::gEvaluationTime)
            //<< ", fitness: "
           // << " "
   //         << sumFitness / gNumberOfRobots
            //<< ", gatheredGenomes: "


    // Logging, population-level: alive
    //std::string sLog = std::string("") + std::to_string(gWorld->getIterations()) + ",pop,alive," + std::to_string(activeCount) + "\n";
    //gLogManager->write(sLog);
    //gLogManager->flush();

}
double EvolvabilityGRNandOdNEATWorldObserver::computeGlobalDiversity()
{
    double result = 0.0;
    std::vector<double> pairwiseDistances;
    if(gNumberOfRobots == 1)
        return 0.0;
    for (int i = 0 ; i < gNumberOfRobots ; i++ )
    {
        EvolvabilityGRNandOdNEATController* c1 = (dynamic_cast<EvolvabilityGRNandOdNEATController*>
                                  (gWorld->getRobot(i)->getController()));

        for (int j = i + 1 ; j < gNumberOfRobots ; j++ )
        {
            EvolvabilityGRNandOdNEATController* c2 = (dynamic_cast<EvolvabilityGRNandOdNEATController*>
                                      (gWorld->getRobot(j)->getController()));
            pairwiseDistances.push_back(c1->computeBehavDistance(c1->getBehavior(),c2->getBehavior()));
        }
    }
    for(auto d : pairwiseDistances)
    {
        result += d;
    }
    result = result/pairwiseDistances.size();
    //Average of pairwise distances between active behaviors
    return result;

}
std::map<GCIndividual,Stats> EvolvabilityGRNandOdNEATWorldObserver::computeNbOffspringAndFitness()
{
    //For computing selection pressure kendall-tau correlation coefficient
    //between number of offspring and other traits (fitness, gathered genomes,distance, collisions,items)
    std::map<GCIndividual,Stats> result;
    for (int i = 0 ; i < gNumberOfRobots ; i++ )
    {
        EvolvabilityGRNandOdNEATController* c1 = (dynamic_cast<EvolvabilityGRNandOdNEATController*>
                                  (gWorld->getRobot(i)->getController()));

        //if not the first generation
        if(gWorld->getIterations()> EvolvabilityGRNandOdNEATSharedData::gEvaluationTime)
        {
            this->_offspringStats.back()[c1->getPreviousMother()].incrementOffspring(); ;
           // std::cout << gWorld->getIterations() << ":" << this->getOffspringStats().back()[c1->getPreviousMother()].numberOffspring << std::endl;
        }
        result[c1->getGenome().id] = c1->getStats();
    }
    return result;

}
