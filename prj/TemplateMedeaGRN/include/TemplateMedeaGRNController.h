/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */



#ifndef TEMPLATEMEDEAGRNCONTROLLER_H
#define TEMPLATEMEDEAGRNCONTROLLER_H

#include "RoboroboMain/common.h"
#include "RoboroboMain/roborobo.h"
#include "Utilities/Graphics.h"
#include "Controllers/Controller.h"
#include "WorldModels/RobotWorldModel.h"
#include "TemplateMedeaGRN/include/TemplateMedeaGRNAgentObserver.h"
#include <neuralnetworks/NeuralNetwork.h>


#include "grgen/classic.hpp"
#include "grgen/real.hpp"
#include "grgen/common.hpp"
#include "grgen/json.hpp"
#include "grgen/protein.hpp"
#include "grgen/grn.hpp"

#include <gaga/gaga.hpp>

#include <iomanip>

using namespace Neural;
struct GC
{
    int robot_id;
    int gene_id;
    bool operator<(const GC &o)  const {
        if(gene_id == o.gene_id)
            return robot_id < o.robot_id;

        return gene_id < o.gene_id;
    }
    bool operator==(const GC &o)  const {
        return gene_id == o.gene_id && robot_id == o.robot_id;
    }
    friend std::ostream& operator<<(std::ostream& os, const GC& gene_clock);
};

struct Genome
{
    GC id;
    GRN<RealC> controller;
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
        if(TemplateMedeaGRNSharedData::gFitness == 0)
            {
                result = fitness / nbFitnessUpdates;

            }
        else if(TemplateMedeaGRNSharedData::gFitness == 1)
            {
                result = (double)collectedItems / nbFitnessUpdates; // - 0.001 * nbCollisions;

            }
        else
        {
            //default

                std::cerr << "Wrong fitness or not implemented: " << TemplateMedeaGRNSharedData::gFitness << std::endl;
                exit(-1);
        }
        return result;
    }
};

class TemplateMedeaGRNController : public Controller
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
    //GC _genomeId;
    //TODO update genome ID
    Genome _g;

    void createGRN();
    
    //bool _isAlive; // agent stand still if not.
    //bool _isNewGenome;

    //GRN<RealC>
    Genome selectTournament(double sp);
    
    void mutate();

    void stepBehaviour();
    void stepEvolution();
    
    void broadcastGenome();
    void loadNewGenome();
    
    // evolutionary engine
    std::map<GC, Genome > _genomesList;
    //std::map<GC, double > _fitnessList;
    //std::map<GC,int> _birthdateList; // store the birthdate of the received controllers (useful for monitoring).
    // GRN
    unsigned int _nbInputs;
    std::vector<std::string> _inputNames;
    std::vector<std::string> _outputNames;
    unsigned int _nbOutputs;
    unsigned int _nbRegulatory;
    
    // logging purpose
    double _Xinit;
    double _Yinit;
    double _dSumTravelled;
    
    void setInputs(GRN<RealC> &g, std::vector<double> in);
    std::vector<double> getOutputs(GRN<RealC> g);
    bool storeGenome(Genome g);
    //GRN<RealC> genome, GC senderId, int senderBirthdate, double fitness);
    void resetRobot();

    double getBroadcastRate();

    bool doBroadcast()
    {
        if(((double)rand() / RAND_MAX) <= getBroadcastRate())
            return true;
        else
            return false;
    }

    void logCurrentState();
    
public:
    void updateFitness(double v){_g.fitness+= v;}
    double getFitness(){return _g.getFitness();}
    int getCollectedItems(){return _g.collectedItems;}
    int getNbCollisions(){return _g.nbCollisions;}
    std::map<GC, Genome > getGenomesList(){return _genomesList;}
    TemplateMedeaGRNController(RobotWorldModel *wm);
    ~TemplateMedeaGRNController();
    
    void reset();
    void step();
    
    Genome getGenome()
    {
        return _g;
    }
    void collectItem()
    {
        _g.collectedItems++;
    }

    int getBirthdate() { return _g.birthdate; }
    
};


#endif

