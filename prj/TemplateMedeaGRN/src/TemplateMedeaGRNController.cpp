/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */

#include "TemplateMedeaGRN/include/TemplateMedeaGRNController.h"
#include "TemplateMedeaGRN/include/TemplateMedeaGRNWorldObserver.h"

#include "World/World.h"
#include "Utilities/Misc.h"
#include <math.h>
#include <string>

#include <neuralnetworks/MLP.h>
#include <neuralnetworks/Perceptron.h>
#include <neuralnetworks/Elman.h>



#include "grgen/classic.hpp"
#include "grgen/real.hpp"
#include "grgen/common.hpp"
#include "grgen/json.hpp"
#include "grgen/protein.hpp"
#include "grgen/grn.hpp"
#include <gaga/gaga.hpp>


//using namespace Neural;

TemplateMedeaGRNController::TemplateMedeaGRNController( RobotWorldModel *wm )
{
    _wm = wm;
    
    //_grn = NULL;
    // evolutionary engine
    _g.id.robot_id = _wm->_id;
    _g.id.gene_id = 0;
    _genomesList.clear();
    resetRobot();


    // behaviour
    _iteration = 0;
    _g.birthdate = 0;
    _g.collectedItems = 0;
    _g.nbCollisions = 0;
    _g.generations = 0;
    _g.nbFitnessUpdates = 0;
    _lifetime = -1;
    _nbGenomeTransmission = 0;

    _wm->updateLandmarkSensor(); // wrt closest landmark
    
    _wm->setRobotLED_colorValues(255, 0, 0);
    
}

TemplateMedeaGRNController::~TemplateMedeaGRNController()
{
    //delete &_grn;
    //_grn = NULL;
}

void TemplateMedeaGRNController::reset()
{
    //TODO reset GRN
}

void TemplateMedeaGRNController::step()
{
    _iteration++;
    // * step evolution
    stepEvolution();

    // * step controller
    stepBehaviour();

    double distance_sensor;
    double coef_obstacle = 1.0;

    double trans = _wm->_desiredTranslationalValue / gMaxTranslationalSpeed;
    double rot = _wm->_desiredRotationalVelocity  / gMaxRotationalSpeed;
    double deltaFitness;
    //Navigation fitness instant increment
    for(int i  = 0; i < _wm->_cameraSensorsNb; i++)
    {
        distance_sensor = _wm->getDistanceValueFromCameraSensor(i) /
                _wm->getCameraSensorMaximumDistanceValue(i);
        if (distance_sensor < coef_obstacle)
            coef_obstacle = distance_sensor;
    }

    //Fitness measurement and update
    switch (TemplateMedeaGRNSharedData::gFitness)
    {
    case 0:


        /*std::cout << trans //<< ", " << rot
                  << ", "
                  << coef_obstacle << std::endl;*/
        //deltaFitness: fitness contribution at current time-step
        //abs(trans) in [0,1], (1 - abs(rot)) in [0,1], coeffObstacle is in [0,1]
        deltaFitness = fabs(trans) * (1 - fabs(rot))
                * coef_obstacle;
        //if(gWorld->getRobot(_wm->getId())->isCollision())
        //    deltaFitness -= 0.001;
        //Incrementally averaging deltas
        _g.fitness += deltaFitness; // ((_g.fitness) * _lifetime + deltaFitness) / (_lifetime +1);
        break;
    case 1:
        //TODO Counting items: already done in agent observer
        break;
    case 2:
        //TODO Counting items cooperatively:  already done in agent observer
        break;
    case 3:
        //TODO Moving from zone to zone in a deceptive maze, done in agent observer, except last step
        break;
    default:
        break;
    }
    //if(gWorld->getRobot(_wm->getId())->isCollision())
    if(coef_obstacle < 0.1) //Collision
    {
        _g.nbCollisions +=1; //-= 0.001;
    }
    _g.nbFitnessUpdates++;
    _wm->setRobotLED_colorValues(0, 0, 0);

}

// ################ ######################## ################
// ################ BEHAVIOUR METHOD(S)      ################
// ################ ######################## ################

void TemplateMedeaGRNController::stepBehaviour()
{
    _wm->updateLandmarkSensor(); // update with closest landmark

    // ---- Build inputs ----
    
    std::vector<double>* inputs = new std::vector<double>(_nbInputs);
    int inputToUse = 0;
    int objectId = -1;
    // distance sensors
    for(int i  = 0; i < _wm->_cameraSensorsNb; i++)
    {
        objectId = _wm->getObjectIdFromCameraSensor(i);
        if ( ! PhysicalObject::isInstanceOf(objectId) )
        {
            (*inputs)[inputToUse] = 1.0 - _wm->getDistanceValueFromCameraSensor(i) / _wm->getCameraSensorMaximumDistanceValue(i);
            inputToUse++;
        }
        else
        {
            (*inputs)[inputToUse] = 0.0;
            inputToUse++;
        }

    }
   if((TemplateMedeaGRNSharedData::gFitness == 1) || (TemplateMedeaGRNSharedData::gFitness == 2))
   {
        int type = 1; //item to forage
        for(int i  = 0; i < _wm->_cameraSensorsNb; i++)
        {
            objectId = _wm->getObjectIdFromCameraSensor(i);
            // input: physical object?
            //sensing distance to energy item, 0.0 [not 1.0] if not energy item
            if ( PhysicalObject::isInstanceOf(objectId) )
            {
                if ( type == gPhysicalObjects[objectId - gPhysicalObjectIndexStartOffset]->getType() )
                  (*inputs)[inputToUse] = 1.0 - _wm->getDistanceValueFromCameraSensor(i) /
                      _wm->getCameraSensorMaximumDistanceValue(i);
                else
                    (*inputs)[inputToUse] = 0.0;
                inputToUse++;
            }
            else
            {
                //Not a physical object. But: should still fill in the inputs 0.0
                (*inputs)[inputToUse] = 0.0;
                inputToUse++;
            }

        }
   }


    // floor sensor
    //(*inputs)[inputToUse++] = (double)_wm->getGroundSensor_redValue()/255.0;
    //(*inputs)[inputToUse++] = (double)_wm->getGroundSensor_greenValue()/255.0;
    //(*inputs)[inputToUse++] = (double)_wm->getGroundSensor_blueValue()/255.0;
    
    // landmark (targeted landmark depends on g_skill)
    if ( gLandmarks.size() > 0 )
    {
         (*inputs)[inputToUse++] = _wm->getLandmarkDirectionAngleValue();
         (*inputs)[inputToUse++] = _wm->getLandmarkDistanceValue();
    }
    
    // energy level
   /* if ( gEnergyLevel )
    {
        (*inputs)[inputToUse++] = _wm->getEnergyLevel() / gEnergyMax;
    }*/

    // ---- set inputs, step, and read out ----
    setInputs(_g.controller,*inputs);


    /*for(unsigned int i=0; i< (*inputs).size();i++)
    {
        std::cout << _inputNames[i] << " : " << _g.controller.getProteinConcentration(_inputNames[i],GRN<RealC>::ProteinType_t::input) << std::endl;
    }*/


    _g.controller.step(); //_grn.step();_grn.step();_grn.step();_grn.step();

    std::vector<double> outputs = getOutputs(_g.controller);
    //std::cout << outputs[0] << ", " << outputs[1] << std::endl;
    double lw = outputs[0]; // * 2 - 1;
    double rw = outputs[1]; // * 2 - 1;


    //lw = 0.0; rw = 0.0;

    //TODO out in differential to (translational and rotational)
    _wm->_desiredTranslationalValue = (lw + rw) / 2;
    _wm->_desiredRotationalVelocity = (lw - rw) / 2;

    //TODO add noise to outputs

    // normalize to motor interval values
    _wm->_desiredTranslationalValue = _wm->_desiredTranslationalValue * gMaxTranslationalSpeed;
    _wm->_desiredRotationalVelocity = _wm->_desiredRotationalVelocity * gMaxRotationalSpeed;
    
    delete (inputs);
}

void TemplateMedeaGRNController::setInputs(GRN<RealC> &g, std::vector<double> in)
{
    //g.setInputConcentration("Bias",0.5);
    //TODO check if normalized in [0,1]

    for(unsigned int i=0; i< in.size();i++)
    {
        g.setInputConcentration(_inputNames[i], in[i]);
        //std::cout << _inputNames[i] << " : " << g.getProteinConcentration(_inputNames[i],GRN<RealC>::ProteinType_t::input) << std::endl;
    }

}
std::vector<double> TemplateMedeaGRNController::getOutputs(GRN<RealC> g)
{
    //std::cout << "LW+ :" << std::to_string(g.getOutputConcentration("LW+")) << ", LW- :" << std::to_string(g.getOutputConcentration("LW-"))
    //         << ", RW+ :" << std::to_string(g.getOutputConcentration("RW+")) << ", RW- :" << std::to_string(g.getOutputConcentration("LW-")) << std::endl;
    std::vector<double> result;

    //******************************************************************************************************************************
    //result.push_back(g.getOutputConcentration("LW"));
    //result.push_back(g.getOutputConcentration("RW"));
    //return result;
    //******************************************************************************************************************************

    double sumL = (g.getOutputConcentration("LW+") + g.getOutputConcentration("LW-"));
    double diffL = (g.getOutputConcentration("LW+") - g.getOutputConcentration("LW-"));

    if(sumL != 0)
        result.push_back( diffL / sumL  );
    else
        if(diffL >= 0)
            result.push_back(1.0);
        else
            result.push_back(-1.0);

    double sumR = (g.getOutputConcentration("RW+") + g.getOutputConcentration("RW-"));
    double diffR = (g.getOutputConcentration("RW+") - g.getOutputConcentration("RW-"));

    if(sumR != 0)
        result.push_back( diffR / sumR  );
    else
        if(diffR > 0)
            result.push_back(1.0);
        else
            result.push_back(-1.0);

    return result;
}

void TemplateMedeaGRNController::createGRN()
{
    //TODO    //if (_grn != NULL )   //    delete _grn;

    switch ( TemplateMedeaGRNSharedData::gControllerType )
    {
        case 0:
        {
            // Real-coded GRN
            _g.controller = GRN<RealC>(); //TODO create and init GRN
            //_nbInputs, _nbOutputs


            //_grn.addRandomProtein(GRN<RealC>::ProteinType_t::input, "Bias");
            for (int i = 0; i < _wm->_cameraSensorsNb; ++i)
            {
                _g.controller.addRandomProtein(GRN<RealC>::ProteinType_t::input, "S" + std::to_string(i));
                //_g.controller.addProtein(GRN<RealC>::ProteinType_t::input,"S" + std::to_string(i),Protein<3>({0.1,0.0,0.0},0.0));
                if((TemplateMedeaGRNSharedData::gFitness == 1)
                        || (TemplateMedeaGRNSharedData::gFitness == 2))
                    _g.controller.addRandomProtein(GRN<RealC>::ProteinType_t::input, "I" + std::to_string(i));
            }
            // left wheel differential coding
            _g.controller.addRandomProtein(GRN<RealC>::ProteinType_t::output, "LW+");
            _g.controller.addRandomProtein(GRN<RealC>::ProteinType_t::output, "LW-");
            // right wheel differential coding
            _g.controller.addRandomProtein(GRN<RealC>::ProteinType_t::output, "RW+");
            _g.controller.addRandomProtein(GRN<RealC>::ProteinType_t::output, "RW-");


            //_g.controller.addRandomProtein(GRN<RealC>::ProteinType_t::output, "LW");
            //_g.controller.addRandomProtein(GRN<RealC>::ProteinType_t::output, "RW");
            //_g.controller.addProtein(GRN<RealC>::ProteinType_t::output,"LW",Protein<3>({0.5,0.1,0.7},0.0));
            //_g.controller.addProtein(GRN<RealC>::ProteinType_t::output,"RW",Protein<3>({0.6,0.7,0.1},0.0));

            _g.controller.randomReguls(TemplateMedeaGRNSharedData::gNbRegulatory);
            _g.controller.randomParams();

            _g.controller.updateSignatures();
            for(int i = 0; i < 25; i++) _g.controller.step();
            //std::cout << _g.controller.serialize() << std::endl;
            break;
        }
        default: // default: no controller
            std::cerr << "[ERROR] gController type unknown (value: " << TemplateMedeaGRNSharedData::gControllerType << ").\n";
            exit(-1);
    };
    gProperties.checkAndGetPropertyValue("gModifRate",&_g.controller.config.MODIF_RATE,true);
    gProperties.checkAndGetPropertyValue("gAddRate",&_g.controller.config.ADD_RATE,true);
    gProperties.checkAndGetPropertyValue("gDeleteRate",&_g.controller.config.DEL_RATE,true);
}


// ################ ######################## ################
// ################ EVOLUTION ENGINE METHODS ################
// ################ ######################## ################

void TemplateMedeaGRNController::stepEvolution()
{
    _lifetime++;
    // * broadcasting genome : robot broadcasts its genome to all neighbors (contact-based wrt proximity sensors)
    //Maybe TODO radio networks
    if(doBroadcast())
    {
        broadcastGenome();
    }

    // * lifetime ended: replace genome (if possible)
    if( gWorld->getIterations() > 1 && gWorld->getIterations() % TemplateMedeaGRNSharedData::gEvaluationTime == 0 )
    {
        loadNewGenome();
        _nbGenomeTransmission = 0;
        reset();
    }
    else
    {
        _dSumTravelled = _dSumTravelled + getEuclidianDistance( _wm->getXReal(), _wm->getYReal(), _Xinit, _Yinit ); //remark: incl. squareroot.
    }
    
    //TODO log genome
    //std::string sLog = std::string("");   //gLogManager->write(sLog);

}
bool TemplateMedeaGRNController::doBroadcast()
{
    if(((double)rand() / RAND_MAX) <= getBroadcastRate()
            && gWorld->getIterations() % TemplateMedeaGRNSharedData::gEvaluationTime > TemplateMedeaGRNSharedData::gMaturationTime )
        return true;
    else
        return false;
}
double TemplateMedeaGRNController::getBroadcastRate()
{
    double result = 0.0;
    bool doRank = (TemplateMedeaGRNSharedData::gMatingOperator == 2);
    if(TemplateMedeaGRNSharedData::gMatingOperator > 0) //Rate proportional to fitness or rank
    {
        int rank = 0;
        double totalFitness = _g.getFitness();
        for(auto i: _genomesList)
        {
            if(_g.getFitness() > i.second.getFitness())
                rank++;
            totalFitness += i.second.getFitness();
        }
        if(doRank)
            result = (double)rank/(_genomesList.size() + 1);
        else
            result = _g.getFitness() / totalFitness;
    }
    else //Always broadcast
    {
        if(TemplateMedeaGRNSharedData::gMatingOperator == -2)
        {
          if((gWorld->getIterations() % TemplateMedeaGRNSharedData::gBroadcastTime) == 0 )
          {
              result = 1.0;
          }
          else
          {
              result = 0.0;
          }
        }
        else
        {
            result = 1.0;
        }
    }
    return result;

}

bool TemplateMedeaGRNController::storeGenome(Genome g) //(GRN<RealC> gReceived, GC senderId, int senderBirthdate, double fitness)
{
    //std::map<int,int>::const_iterator it = _birthdateList.find(senderBirthdate);
    //TODO filter on receive
    //TODO limited-size local pop?
    std::map< GC , Genome >::iterator it = _genomesList.find(g.id);
    if(_genomesList.end() == it)
    {
        _genomesList[g.id] = g;
    }
    else
    {
        _genomesList[g.id].fitness = g.fitness;
        _genomesList[g.id].nbCollisions = g.nbCollisions;
        _genomesList[g.id].collectedItems = g.collectedItems;
        _genomesList[g.id].nbFitnessUpdates = g.nbFitnessUpdates;
    }
    //_fitnessList[senderId] = fitness; //TODO multiObjective
    //_birthdateList[senderId] = senderBirthdate;
    return true;
}


void TemplateMedeaGRNController::resetRobot()
{
    _nbInputs = 0;
    //_nbInputs = 1;

    //_nbInputs = ( PhysicalObjectFactory::getNbOfTypes()+3+1 ) * _wm->_cameraSensorsNb; // nbOfTypes + ( isItAnAgent? + isItSameGroupId? + agentAngleDifference?) + isItAWall?
    
    _nbInputs += _wm->_cameraSensorsNb; // + 3; // proximity sensors + ground sensor (3 values)
    for(int i=0; i<_wm->_cameraSensorsNb; i++)
    {
        _inputNames.push_back("S" + std::to_string(i));
    }

    //If task=collect, add object sensors
    if ((TemplateMedeaGRNSharedData::gFitness == 1)
            || (TemplateMedeaGRNSharedData::gFitness == 2))
    {
        // gathering object distance
        _nbInputs +=  _wm->_cameraSensorsNb;
        for(int i=0; i<_wm->_cameraSensorsNb; i++)
        {
            _inputNames.push_back("I" + std::to_string(i));
        }
        //std::cout << "Created item inputs" << std::endl;

        //TOACTIVATE if color effector _nbOutputs += 1;
    }
    if ( gLandmarks.size() > 0 )
    {
        _nbInputs += 2; // incl. landmark (angle,dist)
        _inputNames.push_back("AngleL");
        _inputNames.push_back("DistL");
    }

    _nbOutputs = 2;
    _outputNames.push_back("LW");
    _outputNames.push_back("RW");

    _nbRegulatory = TemplateMedeaGRNSharedData::gNbRegulatory;


    createGRN(); //updateSignature and warmup called inside


    // initialize robot
    _g.fitness = 0.0;
    _g.collectedItems = 0;
    _g.nbCollisions= 0;
    _g.nbFitnessUpdates = 0;
    _g.generations = 1;


    _genomesList.clear();
    //_fitnessList.clear();
    //_birthdateList.clear();
    
}


void TemplateMedeaGRNController::broadcastGenome()
{
    //TODO communication through distance matrix (in WorldObserver)
    //TODO communication mating

    if(TemplateMedeaGRNSharedData::gCommunicationOnRadius)
    {
        std::vector<std::vector<double> > distanceMatrix = (dynamic_cast<TemplateMedeaGRNWorldObserver*>(gWorld->getWorldObserver()))->getRobotDistances();
        for(int i = 0; i < gNumberOfRobots ; i++)
        {

            if(distanceMatrix[_wm->getId()][i] < TemplateMedeaGRNSharedData::gCommunicationRange
                    && i != _wm->getId())
            {
                TemplateMedeaGRNController* targetRobotController = dynamic_cast<TemplateMedeaGRNController*>(gWorld->getRobot(i)->getController());
                GRN<RealC> copy = GRN<RealC>(_g.controller);
                Genome copyG;
                copyG.controller = copy;
                copyG.fitness = _g.fitness;
                copyG.nbCollisions = _g.nbCollisions;
                copyG.collectedItems = _g.collectedItems;
                copyG.nbFitnessUpdates = _g.nbFitnessUpdates;
                copyG.birthdate = _g.birthdate;
                copyG.id = _g.id;
                copyG.generations = _g.generations;
                bool success = targetRobotController->storeGenome(copyG);

                if ( success == true )
                    _nbGenomeTransmission++;

            }
        }
    }
    else
    {
        // remarque \todo: limiting genome transmission is sensitive to sensor order
        for( int i = 0 ; i < _wm->_cameraSensorsNb
             && ( TemplateMedeaGRNSharedData::gLimitGenomeTransmission == false
                  || ( TemplateMedeaGRNSharedData::gLimitGenomeTransmission == true
                       && _nbGenomeTransmission < TemplateMedeaGRNSharedData::gMaxNbGenomeTransmission ) ); i++)
        {
            int targetIndex = _wm->getObjectIdFromCameraSensor(i);

            if ( targetIndex >= gRobotIndexStartOffset && targetIndex < gRobotIndexStartOffset+gNumberOfRobots
                 && targetIndex != _wm->getId())   // sensor ray bumped into a robot : communication is possible
            {
                targetIndex = targetIndex - gRobotIndexStartOffset; // convert image registering index into robot id.

                TemplateMedeaGRNController* targetRobotController = dynamic_cast<TemplateMedeaGRNController*>(gWorld->getRobot(targetIndex)->getController());

                if ( !targetRobotController )
                {
                    std::cerr << "Error from robot " << _wm->getId() << " : the observer of robot " << targetIndex << " is not compatible." << std::endl;
                    exit(-1);
                }

                GRN<RealC> copy = GRN<RealC>(_g.controller);
                Genome copyG;
                copyG.controller = copy;
                copyG.fitness = _g.fitness;
                copyG.nbCollisions = _g.nbCollisions;
                copyG.collectedItems = _g.collectedItems;
                copyG.nbFitnessUpdates = _g.nbFitnessUpdates;
                copyG.birthdate = _g.birthdate;
                copyG.id = _g.id;
                copyG.generations = _g.generations;
                bool success = targetRobotController->storeGenome(copyG);
                // other agent stores my genome. Contaminant strategy.

                if ( success == true )
                    _nbGenomeTransmission++;
            }
        }
    }
}

void TemplateMedeaGRNController::loadNewGenome()
{
        storeGenome(_g); //,_genomeId, _birthdate, _fitness);
        //With selfinsemination always list>0
        //GRN<RealC>
        Genome parent1, parent2, offspring;
        bool isNew = false;
        switch ( TemplateMedeaGRNSharedData::gSelectionMethod )
        {
            case 0:
               parent1 = selectTournament(TemplateMedeaGRNSharedData::gSelPressure);
                break;
            default:
                std::cerr << "[ERROR] unknown selection method (gSelectionMethod = " << TemplateMedeaGRNSharedData::gSelectionMethod << ")\n";
                exit(-1);
        }


        if ( (double)rand()/RAND_MAX < TemplateMedeaGRNSharedData::gCrossoverProb)
        {
            switch ( TemplateMedeaGRNSharedData::gSelectionMethod )
            {
                case 0:
                    parent2 = selectTournament(TemplateMedeaGRNSharedData::gSelPressure);
                    break;
                default:
                    std::cerr << "[ERROR] unknown selection method (gSelectionMethod = " << TemplateMedeaGRNSharedData::gSelectionMethod << ")\n";
                    exit(-1);
            }
            offspring.controller = parent1.controller.crossover(parent2.controller);
            offspring.fitness = 0.0;
            offspring.collectedItems = 0;
            offspring.nbCollisions= 0;
            offspring.nbFitnessUpdates = 0;
            offspring.generations = 1;
            offspring.birthdate = gWorld->getIterations();
            //offspring.id.gene_id++;
            isNew = true;
        }
        else
        {
            offspring = parent1;
        }
        if ( (double)rand()/RAND_MAX < TemplateMedeaGRNSharedData::gMutateProb)
        {
           offspring.controller.mutate();
           offspring.fitness = 0.0;
           offspring.collectedItems = 0;
           offspring.nbCollisions= 0;
           offspring.nbFitnessUpdates = 0;
           offspring.generations = 1;
           offspring.birthdate = gWorld->getIterations();
           //offspring.id.gene_id++;
           isNew = true;
        }
        if(!isNew)
            offspring.generations++;
        /*cout << "************************************************************************************************************************************************************************************************************" << std::endl;
        std::cout << _grn.serialize() << std::endl;
        cout << "************************************************************************************************************************************************************************************************************" << std::endl;
        std::cout << offspring.serialize() << std::endl;
        cout << "************************************************************************************************************************************************************************************************************" << std::endl;
        */
        offspring.id.gene_id++;

        //TODO destroy previous grn
        _g =  offspring;
        _g.id.robot_id = _wm->getId();
        //_g.controller.setParam(0, 1.0);
        //_g.controller.setParam(1, 1.0);

        //_g.controller.updateSignatures();
        for(int i = 0; i < 25; i++) _g.controller.step();

        _genomesList.clear(); //TODO destroy GRN objects inside

        //_fitnessList.clear();
        //_birthdateList.clear();
        _Xinit = _wm->getXReal();
        _Yinit = _wm->getYReal();
        _dSumTravelled = 0;
        //_fitness = 0.0;
    }

//GRN<RealC>
Genome TemplateMedeaGRNController::selectTournament(double sp)
{
    /* size of the tournament */
    int size = _genomesList.size();
    int inspected = sp * (double) size;

    /* shuffle indexes */
    std::vector<GC> v;
    for(auto i: _genomesList)
        v.push_back(i.first);
    std::random_shuffle(v.begin(), v.end());

    /* get the best from the inspected */
    double max_fit =  _genomesList[v[0]].getFitness();
    GC    best_g  =  v[0];

    for (int i=1 ; i<inspected; i++)
    {
        double f  = _genomesList[v[i]].getFitness();
        if(f > max_fit)
        {
            max_fit = f;
            best_g = v[i] ;
        }
    }

    return _genomesList[best_g];
}
