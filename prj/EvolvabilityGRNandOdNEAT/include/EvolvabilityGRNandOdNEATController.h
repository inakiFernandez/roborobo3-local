/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */



#ifndef EVOLVABILITYGRNANDODNEATCONTROLLER_H
#define EVOLVABILITYGRNANDODNEATCONTROLLER_H

#include "RoboroboMain/common.h"
#include "RoboroboMain/roborobo.h"
#include "Utilities/Graphics.h"
#include "Controllers/Controller.h"
#include "WorldModels/RobotWorldModel.h"
#include "EvolvabilityGRNandOdNEAT/include/EvolvabilityGRNandOdNEATAgentObserver.h"
#include <neuralnetworks/NeuralNetwork.h>
#include "odNeat/include/odNeatAgentObserver.h"
#include "odNeat/include/odNeatSpecies.h"
#include <odneatgc/network.h>
#include <odneatgc/genome.h>

#include "grgen/classic.hpp"
#include "grgen/real.hpp"
#include "grgen/common.hpp"
#include "grgen/json.hpp"
#include "grgen/protein.hpp"
#include "grgen/grn.hpp"

#include <gaga/gaga.hpp>

#include <iomanip>

using namespace Neural;
struct GCIndividual
{
    int robot_id;
    int gene_id;
    bool operator<(const GCIndividual &o)  const {
        if(gene_id == o.gene_id)
            return robot_id < o.robot_id;

        return gene_id < o.gene_id;
    }
    bool operator==(const GCIndividual &o)  const {
        return gene_id == o.gene_id && robot_id == o.robot_id;
    }
    friend std::ostream& operator<<(std::ostream& os, const GCIndividual& gene_clock);
};

struct GenomeData
{
    GCIndividual id;
    GRN<RealC> grn;
    Genome *nnGenome;
    Network *nn;
    std::vector<std::vector<double>> behavior;
    int birthdate;
    double fitness;
    int collectedItems;
    int nbCollisions;
    double cumulatedFitness;
    int nbFitnessUpdates;
    int generations;
    double getFitness()
    {
        double result = -1;
        if(EvolvabilityGRNandOdNEATSharedData::gFitness == 0)
            {
                result = fitness / nbFitnessUpdates;

            }
        else if(EvolvabilityGRNandOdNEATSharedData::gFitness == 1)
            {
                result = (double)collectedItems / nbFitnessUpdates; // - 0.001 * nbCollisions;

            }
        else
        {
            //default

                std::cerr << "Wrong fitness or not implemented: " << EvolvabilityGRNandOdNEATSharedData::gFitness << std::endl;
                exit(-1);
        }
        return result;
    }
};

class EvolvabilityGRNandOdNEATController : public Controller
{
private:
    int _iteration;
    //int _birthdate; // evaluation when this controller was initialized.
    int _lifetime;

    bool _isListening;
    int _notListeningDelay;
    int _listeningDelay;
    
    int _nbGenomeTransmission;
    
    //std::vector<double> _parameters;
    //std::string _nnType;
    //double _fitness;

    //GRN<RealC> _grn;
    //GCIndividual _genomeId;
    //TODO update genome ID
    GenomeData _g;

    void createController();
    
    //bool _isAlive; // agent stand still if not.
    //bool _isNewGenome;

    //GRN<RealC>
    GenomeData selectTournament(double sp);
    
    void mutate();

    void stepBehaviour();
    void stepEvolution();
    
    void broadcastGenome();
    void loadNewGenome();
    
    // evolutionary engine
    std::map<GCIndividual, GenomeData > _genomesList;
    //std::map<GCIndividual, double > _fitnessList;
    //std::map<GCIndividual,int> _birthdateList; // store the birthdate of the received controllers (useful for monitoring).
    // GRN
    unsigned int _nbInputs;
    std::vector<std::string> _inputNames;
    std::vector<std::string> _outputNames;
    unsigned int _nbOutputs;
    unsigned int _nbRegulatory;
    unsigned int _nbHiddenLayers;
    std::vector<unsigned int>* _nbNeuronsPerHiddenLayer;
    //gene (link) and neuron local counters
    int _g_count; int _n_count;
    int _genomeId;

    // logging purpose
    double _Xinit;
    double _Yinit;
    double _dSumTravelled;
    
    void setInputs(GRN<RealC> &g, std::vector<double> in);
    std::vector<double> getOutputs(GRN<RealC> g);
    bool storeGenome(GenomeData g);
    //GRN<RealC> genome, GCIndividual senderId, int senderBirthdate, double fitness);
    void resetRobot();

    double getBroadcastRate();

    bool doBroadcast();

    void logCurrentState();
    std::vector<std::vector<double> > _behavior;
    void updateBehavior(std::vector<double> descriptor)
    {
        _behavior.push_back(descriptor);
    }

    void cleanBehavior()
    {
        for(auto v : _behavior)
        {
            v.clear();
        }
        _behavior.clear();
    }
    bool doUpdateBehavior()
    {
        return ((_lifetime  % EvolvabilityGRNandOdNEATSharedData::freqMeasureBehav) == 0 );
    }
public:
    void updateFitness(double v){_g.fitness+= v;}
    double getFitness(){return _g.getFitness();}
    int getCollectedItems(){return _g.collectedItems;}
    int getNbCollisions(){return _g.nbCollisions;}
    double getDSumTravelled(){return _dSumTravelled;}
    double getNbUnits()
    {
        if(EvolvabilityGRNandOdNEATSharedData::gControllerType == 0)
        {
            return _g.grn.getProteinSize(ProteinType::regul);
        }
        else if(EvolvabilityGRNandOdNEATSharedData::gControllerType == 1)
        {
            double result = 0.0;
            for(auto it = _g.nnGenome->nodes.begin(); it != _g.nnGenome->nodes.end(); it++)
            {
                if((*it)->gen_node_label == HIDDEN)
                    result +=1.0;
            }
            return result;
        }
    }
    std::map<GCIndividual, GenomeData > getGenomesList(){return _genomesList;}
    EvolvabilityGRNandOdNEATController(RobotWorldModel *wm);
    ~EvolvabilityGRNandOdNEATController();
    double _minValue, _maxValue;
    void reset();
    void step();
    
    GenomeData getGenome()
    {
        return _g;
    }
    GenomeData getCurrentGenome()
    {
        return _g;
    }
    void collectItem()
    {
        _g.collectedItems++;
    }

    int getBirthdate() { return _g.birthdate; }
    double getAvgPopFitness()
    {
        double result = 0.0;
        for (auto it = _genomesList.begin(); it != _genomesList.end();++it)
            result += (*it).second.fitness;
        if (_genomesList.size() != 0)
            return result / _genomesList.size();
        else
            return 0.0;
    }
    std::vector<std::vector <double> > getBehavior()
    {
        return _g.behavior;
    }
    double computeIntraRobotDiversity();
    double computeCurrentVSLocalPopDiversity();
    double computeBehavDistance(std::vector< std::vector<double> > b1,std::vector< std::vector<double> > b2);
    void storeGenomeHelper(GenomeData genome)
    {
        storeGenome(genome);
    }
    std::map<GCIndividual, std::vector<std::vector<double> > > _behaviorsList;

    std::vector<std::vector <double> > computeFunctionalControllerBehavior(GenomeData g); //Genome* std::vector<double> g);
    std::vector<std::vector<double> > getBehavior(GCIndividual id);
};


#endif

