/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */

#include "SocialBehaviorANN/include/SocialBehaviorANNController.h"
#include "SocialBehaviorANN/include/SocialBehaviorANNWorldObserver.h"

#include "odneatgc/helper.h"
#include "odneatgc/gene.h"

#include "World/World.h"
#include "Utilities/Misc.h"
#include <math.h>
#include <string>
#include <set>
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

SocialBehaviorANNController::SocialBehaviorANNController( RobotWorldModel *wm )
{
    _wm = wm;
     _g.nn = NULL;
     // weight boundaries NN
     _minValue = -5.0;
     _maxValue = 5.0;

     //_g.controller = NULL;
    // evolutionary engine
    _g.id.robot_id = _wm->_id;
    _genomeId= 0;
    _g.id.gene_id = _genomeId;

    _genomesList.clear();
    resetRobot();


    // behaviour
    _iteration = 0;
    _g.birthdate = 0;
    _g.collectedItems = 0;
    _g.nbCollisions = 0;
    _g.generations = 0;
    _g.nbFitnessUpdates = 1;

    /*
    std::cout << "Begin" << std::endl;

    std::vector<double> inputs; // = new std::vector<double>(_nbInputs);
    inputs.clear();
    std::vector<double> outputs;
    for(int k = 0; k < 3; k++)
    {
        for(int i = 0; i < 2; i++)
        {

            for(int j  = 0; j < 17; j++)
            {
                inputs.push_back(((double)rand()/RAND_MAX)); //1.0); //
            }
            std::cout << "I:   " << std::flush;
            for(int j  = 0; j < 17; j++)
            {
                std::cout << inputs[j] << ", ";
            }
            std::cout << std::endl;
            switch ( SocialBehaviorANNSharedData::gControllerType )
            {
                case 0:
                {
                    setInputs(_g.grn,inputs);
                    _g.grn.step();
                    outputs = getOutputs(_g.grn);
                    break;
                }
                case 1:
                {
                    _g.nn->load_sensors (&((inputs)[0]));
                    if (!(_g.nn->activate ()))
                    { std::cerr << "[ERROR]"<< std::endl; exit (-1);}
                    for (auto out_iter  = _g.nn->outputs.begin();
                         out_iter != _g.nn->outputs.end();
                         out_iter++)
                    {
                        outputs.push_back((*out_iter)->activation);
                    }
                    break;
                }
                default:
                {
                    std::cerr << "Wrong controller type" << std::endl;
                    exit(-1);
                }
            }
            std::cout << "O:   ";
            for(auto it = outputs.begin(); it != outputs.end(); it++)
            {
                std::cout << (*it) << ", "; //std::endl;
            }
            std::cout << std::endl;
            std::cout << std::endl;
            inputs.clear();
            outputs.clear();
        }
        switch ( SocialBehaviorANNSharedData::gControllerType )
        {
            case 0:
            {
                _g.grn.mutate();
                _g.grn.updateSignatures();
                for(int i = 0; i < 25; i++) _g.grn.step();
                break;
            }
            case 1:
            {
                _g.nnGenome->mutate_link_weights(SocialBehaviorANNSharedData::gSigmaRef);
                _g.nn = _g.nnGenome->genesis();
                 break;
            }
            default:
            {
                std::cerr << "Wrong controller type" << std::endl;
                exit(-1);
            }
        }

    }
    exit(-1);
    */


    _lifetime = -1;
    _nbGenomeTransmission = 0;

    _wm->updateLandmarkSensor(); // wrt closest landmark
    
    _wm->setRobotLED_colorValues(255, 0, 0);
    
}

SocialBehaviorANNController::~SocialBehaviorANNController()
{
    if(SocialBehaviorANNSharedData::gControllerType == 0)
    {

        //delete &_grn;
        //_grn = NULL;
    }
    else if (SocialBehaviorANNSharedData::gControllerType == 1)
    {

        delete _g.nn;
        _g.nn = NULL;

    }
    else
    {
        //default
        std::cerr << "[ERROR] Wrong type of controller (not 0 or 1), but: " <<SocialBehaviorANNSharedData::gControllerType << std::endl;
        exit(-1);

    }
    //delete _g;

    _genomesList.clear();
    //_tabu.clear();
}

void SocialBehaviorANNController::reset()
{
    if(SocialBehaviorANNSharedData::gControllerType == 0)
    {
    }
    else if (SocialBehaviorANNSharedData::gControllerType == 1)
    {


    }

    /*_g.nbFitnessUpdates = 0;
    //_energy = SocialBehaviorANNSharedData::gDefaultInitialEnergy;
    _g.fitness =  0 ; //_energy;
    if (SocialBehaviorANNSharedData::gFitness == 1)
    {
        _g.collectedItems = 0;
    }
    //_genome->species = -1;
    _g.nbFitnessUpdates ++;
    _g.birthdate = gWorld->getIterations();
    _lifetime = -1;
    createController();*/

    //recompute_all_species();
    //NOT Adding new genome (even if fitness is not yet measured)
    //TODO fix add to pop
    //GenomeData copy = _g;
    //copy.nnGenome = _g.nnGenome->duplicate();
    //copy.nn = copy.nnGenome->genesis();
    //copy.grn = GRN(_g.grn);
    //storeGenome(copy);//,_energy

}

void SocialBehaviorANNController::step()
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
    switch (SocialBehaviorANNSharedData::gFitness)
    {
    case 0:
       //deltaFitness: fitness contribution at current time-step
        //abs(trans) in [0,1], (1 - abs(rot)) in [0,1], coeffObstacle is in [0,1]
        deltaFitness = fabs(trans) * (1 - fabs(rot))
                * coef_obstacle;
        //if(gWorld->getRobot(_wm->getId())->isCollision()) deltaFitness -= 0.001;
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

void SocialBehaviorANNController::stepBehaviour()
{
    _wm->updateLandmarkSensor(); // update with closest landmark

    // ---- Build inputs ----
    
    //std::vector<double>* inputs = new std::vector<double>(_nbInputs);
    //inputs->clear();
    std::vector<double> inputs;
    int inputToUse = 0;
    int objectId = -1;
    // distance sensors
    for(int i  = 0; i < _wm->_cameraSensorsNb; i++)
    {
        objectId = _wm->getObjectIdFromCameraSensor(i);
        if ( ! PhysicalObject::isInstanceOf(objectId) )
        {
            //(*inputs)[inputToUse] = 1.0 - _wm->getDistanceValueFromCameraSensor(i) / _wm->getCameraSensorMaximumDistanceValue(i);
            inputs.push_back(1.0 - _wm->getDistanceValueFromCameraSensor(i) / _wm->getCameraSensorMaximumDistanceValue(i));
            inputToUse++;
        }
        else
        {
            //(*inputs)[inputToUse] = 0.0;
            inputs.push_back(0.0);
            inputToUse++;
        }

    }
   if((SocialBehaviorANNSharedData::gFitness == 1)
           || (SocialBehaviorANNSharedData::gFitness == 2))
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
                {
                    //(*inputs)[inputToUse] = 1.0 - _wm->getDistanceValueFromCameraSensor(i) /
                    //  _wm->getCameraSensorMaximumDistanceValue(i);
                    inputs.push_back(1.0 - _wm->getDistanceValueFromCameraSensor(i) /
                                       _wm->getCameraSensorMaximumDistanceValue(i));
                }
                else
                {
                    //(*inputs)[inputToUse] = 0.0;
                    inputs.push_back(0.0);
                }
                inputToUse++;
            }
            else
            {
                //Not a physical object. But: should still fill in the inputs 0.0
                inputs.push_back(0.0);
                //(*inputs)[inputToUse] = 0.0;
                inputToUse++;
            }

        }
   }
   //Bias
   //(*inputs)[inputToUse++] = 1.0;
    inputs.push_back(1.0);

    // floor sensor
    //(*inputs)[inputToUse++] = (double)_wm->getGroundSensor_redValue()/255.0;
    //(*inputs)[inputToUse++] = (double)_wm->getGroundSensor_greenValue()/255.0;
    //(*inputs)[inputToUse++] = (double)_wm->getGroundSensor_blueValue()/255.0;
    
    // landmark (targeted landmark depends on g_skill)
    //if ( gLandmarks.size() > 0 )
    //{
     //    (*inputs)[inputToUse++] = _wm->getLandmarkDirectionAngleValue();
      //   (*inputs)[inputToUse++] = _wm->getLandmarkDistanceValue();
    //}
    
    // energy level
   /* if ( gEnergyLevel )
    {
        (*inputs)[inputToUse++] = _wm->getEnergyLevel() / gEnergyMax;
    }*/
   /*for(int i = 0; i < _nbInputs; i++)
   {
       std::cout << inputs[i] << ", ";
   }
    std::cout << std::endl;*/

    // ---- set inputs, step, and read out ----
    std::vector<double> outputs;
    if(SocialBehaviorANNSharedData::gControllerType == 0)
    {
        setInputs(_g.grn,inputs);

        _g.grn.step();

        outputs = getOutputs(_g.grn);

    }
    else if (SocialBehaviorANNSharedData::gControllerType == 1)
    {
        //std::cout << "step" << std::endl;

        _g.nn->load_sensors (&((inputs)[0]));
        //std::cout << "step2" << std::endl;

        if (!(_g.nn->activate ()))
        {
            std::cerr << "[ERROR] Activation of ANN not correct: genome R"
               << _g.id.robot_id << ", G" << _g.id.gene_id << std::endl;
            //save_genome();
            exit (-1);
        }
        //std::cout << "step3" << std::endl;

        for (auto out_iter  = _g.nn->outputs.begin();
             out_iter != _g.nn->outputs.end();
             out_iter++)
        {
            outputs.push_back((*out_iter)->activation);
        }
    }
    else
    {
        //default
        std::cerr << "[ERROR] Wrong type of controller (not 0 or 1), but: " <<SocialBehaviorANNSharedData::gControllerType << std::endl;
        exit(-1);
    }

    double lw, rw;
    //std::cout << outputs[0] << ", " << outputs[1] << std::endl;
    if (SocialBehaviorANNSharedData::gControllerType == 0)
    {
        lw = outputs[0] * 2 - 1;
        rw = outputs[1] * 2 - 1;
    }
    else if (SocialBehaviorANNSharedData::gControllerType == 1)
    {
        lw = outputs[0];
        rw = outputs[1];
    }

    //std::cout << lw << ", " << rw << std::endl;

    //lw = 0.0; rw = 0.0;

    //TODO out in differential to (translational and rotational)
    _wm->_desiredTranslationalValue = (lw + rw) / 2;
    _wm->_desiredRotationalVelocity = (lw - rw) / 2;

    //TODO add noise to outputs

    // normalize to motor interval values
    _wm->_desiredTranslationalValue = _wm->_desiredTranslationalValue * gMaxTranslationalSpeed;
    _wm->_desiredRotationalVelocity = _wm->_desiredRotationalVelocity * gMaxRotationalSpeed;
    
    //delete (inputs);
}

void SocialBehaviorANNController::setInputs(GRN<RealC> &g, std::vector<double> in)
{
    //TODO check if normalized in [0,1]
    for(unsigned int i=0; i< _nbInputs;i++)//_inputNames.size();i++)
    {
        //std::cout << i << "," << in[i] << "," << _inputNames[i] << std::endl;
        g.setInputConcentration(_inputNames[i], in[i]);
    }

}
std::vector<double> SocialBehaviorANNController::getOutputs(GRN<RealC> g)
{
    //std::cout << "LW+ :" << std::to_string(g.getOutputConcentration("LW+")) << ", LW- :" << std::to_string(g.getOutputConcentration("LW-"))
    //         << ", RW+ :" << std::to_string(g.getOutputConcentration("RW+")) << ", RW- :" << std::to_string(g.getOutputConcentration("LW-")) << std::endl;
    std::vector<double> result;

    //******************************************************************************************************************************
    //result.push_back(g.getOutputConcentration("LW"));
    //result.push_back(g.getOutputConcentration("RW"));
    //return result;
    //******************************************************************************************************************************

    /*
    //differential wheel coding
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
    */
    //direct wheel coding
    result.push_back(g.getOutputConcentration("LW"));
    result.push_back(g.getOutputConcentration("RW"));
    return result;
}

void SocialBehaviorANNController::createController() //GRN()
{
    //TODO    //if (_grn != NULL )   //    delete _grn;
    GC genomeId;
    genomeId.robot_id = _wm->_id;
    genomeId.gene_id = 0;
    switch ( SocialBehaviorANNSharedData::gControllerType )
    {
        case 0:
        {
            // Real-coded GRN
            _g.grn = GRN<RealC>();
            //_nbInputs, _nbOutputs


            _g.grn.addRandomProtein(GRN<RealC>::ProteinType_t::input, "Bias");
            for (int i = 0; i < _wm->_cameraSensorsNb; ++i)
            {
                _g.grn.addRandomProtein(GRN<RealC>::ProteinType_t::input, "S" + std::to_string(i));
                if((SocialBehaviorANNSharedData::gFitness == 1)
                        || (SocialBehaviorANNSharedData::gFitness == 2))
                    _g.grn.addRandomProtein(GRN<RealC>::ProteinType_t::input, "I" + std::to_string(i));
            }
            // left wheel direct coding
            _g.grn.addRandomProtein(GRN<RealC>::ProteinType_t::output, "LW");
            // left wheel differential coding
            //_g.grn.addRandomProtein(GRN<RealC>::ProteinType_t::output, "LW+");
            //_g.grn.addRandomProtein(GRN<RealC>::ProteinType_t::output, "LW-");
            //right wheel direct coding
            _g.grn.addRandomProtein(GRN<RealC>::ProteinType_t::output, "RW");
            // right wheel differential coding
            //_g.grn.addRandomProtein(GRN<RealC>::ProteinType_t::output, "RW+");
            //_g.grn.addRandomProtein(GRN<RealC>::ProteinType_t::output, "RW-");


            _g.grn.randomReguls(SocialBehaviorANNSharedData::gNbRegulatory);
            _g.grn.randomParams();

            _g.grn.updateSignatures();
            for(int i = 0; i < 25; i++) _g.grn.step();

            gProperties.checkAndGetPropertyValue("gModifRate",&_g.grn.config.MODIF_RATE,true);
            gProperties.checkAndGetPropertyValue("gAddRate",&_g.grn.config.ADD_RATE,true);
            gProperties.checkAndGetPropertyValue("gDeleteRate",&_g.grn.config.DEL_RATE,true);

            break;
        }
        case 1:
        {
        // Inputs, outputs, 0 hidden neurons, fully connected.
        //Initial Genes=>common to all agents, thus identified by a common historical marker
            //TODO for the struct Genome

            _g.nnGenome = new Genome (genomeId,_nbInputs, _nbOutputs); //,SocialBehaviorANNSharedData::gNbRegulatory);
            //Initialize with braitenberg
            std::vector<double> w;
            w = {1.0, 0.5, -1.5, -0.5, -1.0, 0.0, 0.0, 0.0,     //obst. sensor to left wh.
                 -1.0, -0.5, 0.5, 0.5, 1.0, 0.0, 0.0, 0.0,      //item sensor to left wh.
                  0.25,
                 -1.0, -0.5, -1.5, 0.5, 1.0, 0.0, 0.0, 0.0,       //obst. sensor to right wh.
                 1.0, 0.5, 0.5, -0.5, -1.0, 0.0, 0.0, 0.0,     //item sensor to right wh.
                 0.25
                };
            int i =0;
            for(auto &curgene : _g.nnGenome->genes)
            {
                curgene->lnk->weight =w[i]*15;
                i++;
            }
            //_g.nnGenome->mutate_link_weights(1.0);

            //Fully connected
            _g_count = _nbInputs * _nbOutputs + 1;
            _n_count = _nbInputs + _nbOutputs + 1;

            if (_g.nn != NULL)
                delete _g.nn;
            _g.nn= _g.nnGenome->genesis();
            break;
        }
        default: // default: no controller
            std::cerr << "[ERROR] gController type unknown (value: " << SocialBehaviorANNSharedData::gControllerType << ").\n";
            exit(-1);
    };
    if(SocialBehaviorANNSharedData::gDoMeasureDiv)
    {
        _g.behavior = computeFunctionalControllerBehavior(_g);
    }

}


// ################ ######################## ################
// ################ EVOLUTION ENGINE METHODS ################
// ################ ######################## ################

void SocialBehaviorANNController::stepEvolution()
{
    _lifetime++;
    // * broadcasting genome : robot broadcasts its genome to all neighbors (contact-based wrt proximity sensors)
    //Maybe TODO radio networks
    if(doBroadcast())
    {
        broadcastGenome();
    }

    // * lifetime ended: replace genome (if possible)
    if( gWorld->getIterations() > 1 && gWorld->getIterations() % SocialBehaviorANNSharedData::gEvaluationTime == 0 )
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
bool SocialBehaviorANNController::doBroadcast()
{
    if(!SocialBehaviorANNSharedData::gIsCentralized)
    {
        if(((double)rand() / RAND_MAX) <= getBroadcastRate()
                && (gWorld->getIterations() % SocialBehaviorANNSharedData::gEvaluationTime) > SocialBehaviorANNSharedData::gMaturationTime )
            return true;
        else
            return false;
    }
}
double SocialBehaviorANNController::getBroadcastRate()
{
    double result = 0.0;
    bool doRank = (SocialBehaviorANNSharedData::gMatingOperator == 2);
    if(SocialBehaviorANNSharedData::gMatingOperator > 0) //Rate proportional to fitness or rank
    {
        int rank = 1;
        double totalFitness = _g.getFitness();
        for(auto i: _genomesList)
        {
            if(_g.getFitness() >= i.second.getFitness())
                rank++;
            totalFitness += i.second.getFitness();
        }
        if(doRank)
        {
            result = (double)rank/((double)_genomesList.size() + 1.0);
        }
        else
        {
            result = (_g.getFitness() + 1) / (totalFitness + 1);
        }
    }
    else
    {
        if(SocialBehaviorANNSharedData::gMatingOperator == -2)
        {
          if((gWorld->getIterations() % SocialBehaviorANNSharedData::gBroadcastTime) == 0 )
          {
              result = 1.0;
          }
          else
          {
              result = 0.0;
          }
        }
        else  //Always broadcast
        {
            result = 1.0;
        }
    }
    return result;

}

bool SocialBehaviorANNController::storeGenome(GenomeData g) //(GRN<RealC> gReceived, GCIndividual senderId, int senderBirthdate, double fitness)
{
    //std::map<int,int>::const_iterator it = _birthdateList.find(senderBirthdate);
    //TODO filter on receive. maybe based on behavioral distance
    std::map< GCIndividual , GenomeData >::iterator it = _genomesList.find(g.id);
    bool stored = false;
    if(_genomesList.end() != it)
    {
        _genomesList[g.id].fitness = g.fitness;
        //_genomesList[g.id].behavior = g.behavior;
        _genomesList[g.id].nbCollisions = g.nbCollisions;
        _genomesList[g.id].collectedItems= g.collectedItems;
        _genomesList[g.id].nbFitnessUpdates = g.nbFitnessUpdates;
        //TODO delete g
        if(SocialBehaviorANNSharedData::gControllerType == 0)
        {
            //delete g.grn;
        }
        else if(SocialBehaviorANNSharedData::gControllerType == 1)
        {
            delete g.nnGenome;
            delete g.nn;
        }
        stored = true;
    }
    else
    {
        if(_genomesList.size() < SocialBehaviorANNSharedData::gPopSize)
        {
            _genomesList[g.id] = g;
            stored = true;
        }
        else
        {
            double minFitness = std::numeric_limits<double>::max();
            struct GCIndividual indexWorse;
            indexWorse.robot_id = -1; indexWorse.gene_id = -1;
            std::map<GCIndividual, GenomeData>::iterator it = _genomesList.begin();
            for(; it != _genomesList.end(); it++)
            {
                if ((*it).second.fitness < minFitness)
                {
                    minFitness = (*it).second.fitness;
                    indexWorse = (*it).first;
                }
            }
            if (minFitness <= g.fitness)
            {
                if(SocialBehaviorANNSharedData::gControllerType == 0)
                {
                    //delete _genomesList[indexWorse].grn;
                }
                else if(SocialBehaviorANNSharedData::gControllerType == 1)
                {
                    delete _genomesList[indexWorse].nnGenome;
                    delete _genomesList[indexWorse].nn;
                }
                _genomesList.erase(indexWorse);
                _genomesList[g.id] = g;

            }

        }
    }    
    //TODO multiObjective
    return stored;
}


void SocialBehaviorANNController::resetRobot()
{
    _nbInputs = 0;

    //_nbInputs = ( PhysicalObjectFactory::getNbOfTypes()+3+1 ) * _wm->_cameraSensorsNb; // nbOfTypes + ( isItAnAgent? + isItSameGroupId? + agentAngleDifference?) + isItAWall?
    
    _nbInputs += _wm->_cameraSensorsNb; // + 3; // proximity sensors + ground sensor (3 values)
    for(int i=0; i<_wm->_cameraSensorsNb; i++)
    {
        _inputNames.push_back("S" + std::to_string(i));
    }

    //If task=collect, add object sensors
    if ((SocialBehaviorANNSharedData::gFitness == 1)
            || (SocialBehaviorANNSharedData::gFitness == 2))
    {
        // gathering object distance
        _nbInputs +=  _wm->_cameraSensorsNb;
        for(int i=0; i<_wm->_cameraSensorsNb; i++)
        {
            _inputNames.push_back("I" + std::to_string(i));
        }

        //TOACTIVATE if color effector _nbOutputs += 1;
    }
    /*if ( gLandmarks.size() > 0 )
    {
        _nbInputs += 2; // incl. landmark (angle,dist)
        _inputNames.push_back("AngleL");
        _inputNames.push_back("DistL");
    }*/
    //Energy
    //_nbInputs += 1;
    //Bias
    _inputNames.push_back("Bias");
    _nbInputs += 1;

    _nbOutputs = 2;

    /*
    _outputNames.push_back("LW+");
    _outputNames.push_back("RW-");
    _outputNames.push_back("LW-");
    _outputNames.push_back("RW+");
    */
    _outputNames.push_back("LW");
    _outputNames.push_back("RW");
    if(SocialBehaviorANNSharedData::gControllerType == 0)
    {
        _nbRegulatory = SocialBehaviorANNSharedData::gNbRegulatory;
        //createGRN(); //updateSignature and warmup called inside
    }
    else if(SocialBehaviorANNSharedData::gControllerType == 1)
    {


    }
    createController();
    // initialize robot
    _g.fitness = 0.0;
    _g.collectedItems = 0;
    _g.nbCollisions= 0;
    _g.nbFitnessUpdates = 1;
    _g.generations = 1;


    _genomesList.clear(); //TODO use same data structure for pop
    //_fitnessList.clear();
    //_birthdateList.clear();
    
}


void SocialBehaviorANNController::broadcastGenome()
{
    //TODO communication through distance matrix (in WorldObserver)
    //TODO communication mating

    if(SocialBehaviorANNSharedData::gCommunicationOnRadius)
    {
        std::vector<std::vector<double> > distanceMatrix = (dynamic_cast<SocialBehaviorANNWorldObserver*>(gWorld->getWorldObserver()))->getRobotDistances();
        for(int i = 0; i < gNumberOfRobots ; i++)
        {

            if(distanceMatrix[_wm->getId()][i] < SocialBehaviorANNSharedData::gCommunicationRange
                    && i != _wm->getId())
            {
                SocialBehaviorANNController* targetRobotController = dynamic_cast<SocialBehaviorANNController*>(gWorld->getRobot(i)->getController());

                GRN<RealC> copy;
                Genome* copyNNG = NULL;
                Network* copyNN = NULL;
                if(SocialBehaviorANNSharedData::gControllerType == 0)
                {
                    copy = GRN<RealC>(_g.grn);
                }
                else if(SocialBehaviorANNSharedData::gControllerType == 1)
                {
                    copyNNG = _g.nnGenome->duplicate();
                    copyNN = copyNNG->genesis();
                }
                GenomeData copyG;
                copyG.grn = copy;
                copyG.nnGenome = copyNNG;
                copyG.nn = copyNN;
                copyG.fitness = _g.fitness;
                copyG.behavior = _g.behavior;
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
             && ( SocialBehaviorANNSharedData::gLimitGenomeTransmission == false
                  || ( SocialBehaviorANNSharedData::gLimitGenomeTransmission == true
                       && _nbGenomeTransmission < SocialBehaviorANNSharedData::gMaxNbGenomeTransmission ) ); i++)
        {
            int targetIndex = _wm->getObjectIdFromCameraSensor(i);

            if ( targetIndex >= gRobotIndexStartOffset && targetIndex < gRobotIndexStartOffset+gNumberOfRobots
                 && targetIndex != _wm->getId())   // sensor ray bumped into a robot : communication is possible
            {
                targetIndex = targetIndex - gRobotIndexStartOffset; // convert image registering index into robot id.

                SocialBehaviorANNController* targetRobotController = dynamic_cast<SocialBehaviorANNController*>(gWorld->getRobot(targetIndex)->getController());

                if ( !targetRobotController )
                {
                    std::cerr << "Error from robot " << _wm->getId() << " : the observer of robot " << targetIndex << " is not compatible." << std::endl;
                    exit(-1);
                }
                GRN<RealC> copy;
                Genome* copyNNG = NULL;
                Network* copyNN = NULL;
                if(SocialBehaviorANNSharedData::gControllerType == 0)
                {
                    copy = GRN<RealC>(_g.grn);
                }
                else if(SocialBehaviorANNSharedData::gControllerType == 1)
                {
                    copyNNG = _g.nnGenome->duplicate();
                    copyNN = copyNNG->genesis();
                }

                GenomeData copyG;
                copyG.grn = copy;
                copyG.nnGenome = copyNNG;
                copyG.nn = copyNN;
                copyG.fitness = _g.fitness;
                copyG.behavior = _g.behavior;
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
//TODO add behavior descriptors and diversity measuring here and in world observer.
//TODO: evolvability how to measure, define measure etc
void SocialBehaviorANNController::loadNewGenome()
{
        storeGenome(_g);
        //With selfinsemination always list>0
        //TODO GRN and NN cases with similar operators
        //(pay attention to comparison on fair terms)
        GenomeData parent1, parent2, offspring;
        bool isNew = false;
        _genomeId++;
        switch ( SocialBehaviorANNSharedData::gSelectionMethod )
        {
            case 0:
               parent1 = selectTournament(SocialBehaviorANNSharedData::gSelPressure);
                break;
            default:
                std::cerr << "[ERROR] unknown selection method (gSelectionMethod = " << SocialBehaviorANNSharedData::gSelectionMethod << ")\n";
                exit(-1);
        }


        if ( ((double)rand()/RAND_MAX < SocialBehaviorANNSharedData::gCrossoverProb) && (_genomesList.size() > 1))
        {
            switch ( SocialBehaviorANNSharedData::gSelectionMethod )
            {
                case 0:
                    parent2 = selectTournament(SocialBehaviorANNSharedData::gSelPressure);
                    break;
                default:
                    std::cerr << "[ERROR] unknown selection method (gSelectionMethod = " << SocialBehaviorANNSharedData::gSelectionMethod << ")\n";
                    exit(-1);
            }

            if(SocialBehaviorANNSharedData::gControllerType == 0)
            {
                offspring.grn = parent1.grn.crossover(parent2.grn);
            }
            else if (SocialBehaviorANNSharedData::gControllerType == 1)
            {
                GC newId = {_wm->getId(),_genomeId};
                offspring.nnGenome = parent1.nnGenome->mate(parent2.nnGenome,newId,parent1.fitness,parent2.fitness);
            }
            offspring.fitness = 0.0;
            offspring.collectedItems = 0;
            offspring.nbCollisions= 0;
            offspring.nbFitnessUpdates = 1;
            offspring.generations = 1;
            offspring.birthdate = gWorld->getIterations();
            //offspring.id.gene_id++;
            isNew = true;
        }
        else
        {
            offspring = parent1;
        }
        if ( (double)rand()/RAND_MAX < SocialBehaviorANNSharedData::gMutateProb)
        {

            if(SocialBehaviorANNSharedData::gControllerType == 0)
            {
                offspring.grn.mutate();
            }
            else if (SocialBehaviorANNSharedData::gControllerType == 1)
            {
                GC newId = {_wm->getId(),_genomeId};
                offspring.nnGenome->mutate(SocialBehaviorANNSharedData::gSigmaRef,_wm->getId(),newId,_n_count,_g_count);
                offspring.nn = offspring.nnGenome->genesis();
            }

           offspring.fitness = 0.0;
           offspring.collectedItems = 0;
           offspring.nbCollisions= 0;
           offspring.nbFitnessUpdates = 1;
           offspring.generations = 1;
           offspring.birthdate = gWorld->getIterations();
           //offspring.id.gene_id++;
           isNew = true;
        }

        if(isNew)
        {

            if (SocialBehaviorANNSharedData::gControllerType == 0)
            {
                //_g.controller.setParam(0, 1.0);
                //_g.controller.setParam(1, 1.0);
                offspring.grn.updateSignatures();
                for(int i = 0; i < 25; i++) offspring.grn.step();
                //delete _g.grn;
            }
            else if (SocialBehaviorANNSharedData::gControllerType == 1)
            {
                offspring.nn = offspring.nnGenome->genesis();
                //delete _g.nnGenome;
                //delete _g.nn;
            }
            offspring.id = {_wm->getId(),_genomeId};
            if(SocialBehaviorANNSharedData::gDoMeasureDiv)
            {
                offspring.behavior = computeFunctionalControllerBehavior(offspring);
            }
        offspring.id.gene_id++;
        }
        else
        {
            //offspring = _g;
            _genomesList.erase(_genomesList.find(_g.id));
            offspring.fitness = 0.0;
            offspring.collectedItems = 0;
            offspring.nbCollisions = 0;
            offspring.nbFitnessUpdates = 1;
            offspring.generations++;
        }

        //delete _g.grn;
        //delete _g.nnGenome; //TODO do delete nn?
        //TODO destroy previous grn
        _g =  offspring;
        _g.id.robot_id = _wm->getId();
        if(SocialBehaviorANNSharedData::gControllerType == 0)
        {
            gProperties.checkAndGetPropertyValue("gModifRate",&_g.grn.config.MODIF_RATE,true);
            gProperties.checkAndGetPropertyValue("gAddRate",&_g.grn.config.ADD_RATE,true);
            gProperties.checkAndGetPropertyValue("gDeleteRate",&_g.grn.config.DEL_RATE,true);
        }
        std::vector<GCIndividual> v;
        for(auto i: _genomesList)
            v.push_back(i.first);
        for(int i=0; i< _genomesList.size();i++)
        {
            if(SocialBehaviorANNSharedData::gControllerType == 0)
            {
                //delete _genomesList[v[i]].grn;

            }
            else if (SocialBehaviorANNSharedData::gControllerType ==1)
            {
                //delete _genomesList[v[i]].nnGenome;
                //delete _genomesList[v[i]].nn;
            }
        }
        //std::cout << "afterdelete" << std::endl;
        _genomesList.clear(); //TODO leak? destroy GRN and NN objects inside

        //_fitnessList.clear();
        //_birthdateList.clear();
        _Xinit = _wm->getXReal();
        _Yinit = _wm->getYReal();
        _dSumTravelled = 0;
        //_fitness = 0.0;
    }

GenomeData SocialBehaviorANNController::selectTournament(double sp)
{
    /* size of the tournament */
    int size = _genomesList.size();
    int inspected = sp * (double) size;

    /* shuffle indexes */
    std::vector<GCIndividual> v;
    for(auto i: _genomesList)
        v.push_back(i.first);
    std::random_shuffle(v.begin(), v.end());

    /* get the best from the inspected */
    double max_fit =  _genomesList[v[0]].getFitness();    
    GCIndividual    best_g  =  v[0];
    /*
    std::cout << best_g.robot_id << "-"
              << best_g.gene_id << ": "
              << _genomesList[v[0]].getFitness()* SocialBehaviorANNSharedData::gEvaluationTime
              << std::endl;
    */
    for (int i=1 ; i<inspected; i++)
    {
        double f  = _genomesList[v[i]].getFitness();        
        //std::cout << v[i].robot_id << "-" << v[i].gene_id << ": " << f* SocialBehaviorANNSharedData::gEvaluationTime << std::endl;

        if(f > max_fit)
        {
            max_fit = f;
            best_g = v[i] ;
        }
    }
    return _genomesList[best_g];
}
double SocialBehaviorANNController::computeBehavDistance(std::vector< std::vector<double> > b1,std::vector< std::vector<double> > b2)
{
    double result = 0.0;

    for(unsigned int i = 0 ;  i < b1.size() ; i++)
    {
        double sampleDistance = 0.0;
        //Compute Euclidian distance
        for(unsigned int j = 0; j < b1[i].size(); j++)
        {
            sampleDistance += (b1[i][j] - b2[i][j]) * (b1[i][j] - b2[i][j]);
        }
        result+=sqrt(sampleDistance);
    }
    result =  result / b1.size();
    return result;
}
double SocialBehaviorANNController::computeIntraRobotDiversity()
{
    double result = 0.0;
    //TODO tocheck
    for(auto it = _genomesList.begin();it!=_genomesList.end();it++)
    {
        for(auto it2 = _genomesList.begin();it2!=_genomesList.end();it2++)
        {
            if(it!=it2)
            {
                result += computeBehavDistance(getBehavior(std::get<0>(*it)),getBehavior(std::get<0>(*it2)));
            }
        }
    }
    return result/(_genomesList.size() * _genomesList.size() - 2);
}

double SocialBehaviorANNController::computeCurrentVSLocalPopDiversity()
{
    double result = 0.0;
    //TODO tocheck. todo getbehavior
    for(auto it = _genomesList.begin();it!=_genomesList.end();it++)
    {
        result += computeBehavDistance(getBehavior(),computeFunctionalControllerBehavior(it->second));//&std::get<1>(*it)));
    }
    return result/_genomesList.size();
}
std::vector<std::vector<double> > SocialBehaviorANNController::computeFunctionalControllerBehavior(GenomeData g) //std::vector<double> g
{
    std::vector<double> outputs;
    std::vector<double> inputs;
    std::vector<std::vector<double> > result;
    //TODO maybe copy and recreate????
    if(SocialBehaviorANNSharedData::gControllerType == 0)
    {

    }
    else if (SocialBehaviorANNSharedData::gControllerType == 1)
    {

    }
    else
    {
        //default
        std::cerr << "[ERROR] Wrong type of controller (not 0 or 1), but: " <<SocialBehaviorANNSharedData::gControllerType << std::endl;
        exit(-1);

    }
    for(int i=0; i < SocialBehaviorANNSharedData::gNbInputsBehavior; i++)
    {
        inputs = SocialBehaviorANNSharedData::gInputsBehavior[i];
        if(SocialBehaviorANNSharedData::gControllerType == 0)
        {
            setInputs(_g.grn, inputs);
            _g.grn.step();
            outputs = getOutputs(_g.grn);
        }
        else if (SocialBehaviorANNSharedData::gControllerType == 1)
        {
            _g.nn->load_sensors (&(inputs[0]));
            if (!(_g.nn->activate ()))
            {
                std::cerr << "[ERROR] Activation of ANN not correct: genome R"
                   << _g.id.robot_id << ", G" << _g.id.gene_id << std::endl;
                //save_genome();
                exit (-1);
            }
            for (auto out_iter  = _g.nn->outputs.begin();
                 out_iter != _g.nn->outputs.end();
                 out_iter++)
                outputs.push_back((*out_iter)->activation);
        }

        result.push_back(outputs);
    }
    //delete neuralNet;
    return result;
}
std::vector<std::vector<double> > SocialBehaviorANNController::getBehavior(GCIndividual id)
{
    std::vector<std::vector<double> > result;
    result = _genomesList[id].behavior; //_behaviorsList[id];
    return result;
}
