/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */



#ifndef COLLECT2CONTROLLER_H
#define COLLECT2CONTROLLER_H

#include "RoboroboMain/common.h"
#include "RoboroboMain/roborobo.h"
#include "Utilities/Graphics.h"
#include "Controllers/Controller.h"
#include "WorldModels/RobotWorldModel.h"
#include "Collect2/include/Collect2AgentObserver.h"
#include <neuralnetworks/NeuralNetwork.h>

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

class Collect2Controller : public Controller
{
private:
    int _iteration;
    int _birthdate; // evaluation when this controller was initialized.
    double _currentFitness;

    std::vector<double> _parameters;
    std::string _nnType;
    std::vector<int> _nbHiddenNeuronsPerLayer;
    std::vector<int> _nbBiaisNeuronsPerLayer;
    NeuralNetwork* nn;

    void createNN();
    
    //bool _isAlive; // agent stand still if not.
    bool _isNewGenome;
    std::vector<std::vector<double> > behavior;
    
    void selectRandomGenome();
    void selectBestGenome();
    void selectRankBasedGenome();
    void selectTournament(double sp);

    void mutate(float sigma);

    void stepBehaviour();
    void stepEvolution();
    
    void broadcastGenome();
    void loadNewGenome();
    
    unsigned int computeRequiredNumberOfWeights();

    bool getNewGenomeStatus() { return _isNewGenome; }
    void setNewGenomeStatus( bool __status ) { _isNewGenome = __status; }
    
    // evolutionary engine
    std::vector<double> _genome; // current genome in evaluation
    GC _genomeId;
    std::vector<double> _previousGenome; // 1+1-online-ES surviving genome

    std::vector<double> _currentGenome;
    float _currentSigma;
    int _lifetime;
    
    // ANN
    double _minValue;
    double _maxValue;
    unsigned int _nbInputs;
    unsigned int _nbOutputs;
    unsigned int _nbHiddenLayers;
    std::vector<unsigned int>* _nbNeuronsPerHiddenLayer;
    
    void storeGenome(std::vector<double> genome, GC senderId, double fitness);
    void storeOwnGenome();
    void resetRobot();
    int _countDoorPassages;
    double _minDistNextDoor;

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
        return ((_lifetime  % Collect2SharedData::freqMeasureBehav) == 0 );
    }
    
public:

    Collect2Controller(RobotWorldModel *wm);
    ~Collect2Controller();

    double computeIntraRobotDiversity();
    double computeCurrentVSLocalPopDiversity();
    double computeBehavDistance(std::vector< std::vector<double> > b1,std::vector< std::vector<double> > b2);
    void reset();
    void step();
    void updateFitness(double delta);
    int getBirthdate() { return _birthdate; }
    int getNbInputs(){return _nbInputs;}
    double getFitness()
    {
        if(Collect2SharedData::gFitness==3)
        {
            return getDoorPassages() * 1000.0 + (1000.0 - getMinDistanceToNextDoor());
        }
        return _currentFitness;
    }
    std::vector<double> getCurrentGenome()
    {
        return _currentGenome;
    }
    GC getCurrentId()
    {
        return _genomeId;
    }

    void storeGenomeHelper(std::vector<double> genome, GC senderId, double fitness)
    {
        storeGenome(genome,senderId,fitness);
    }

    int getDoorPassages()
    {
        return _countDoorPassages;
    }
    int passDoor()
    {
        return _countDoorPassages++;
    }

    double getAvgPopFitness()
    {
        double result = 0.0;
        for (auto it = _fitnessList.begin(); it != _fitnessList.end();++it)
            result += (*it).second;
        if (_fitnessList.size() != 0)
            return result / _fitnessList.size();
        else
            return 0.0;
    }
    void logGenome(std::string s);
    std::map<GC, std::vector<double> > _genomesList;
    std::map<GC, double > _fitnessList;

    int _lastZone = -1;
    double getDistanceToNextDoor();
    void resetMinDistNextDoor()
    {
        _minDistNextDoor =  std::numeric_limits<double>::max();
    }
    double getMinDistanceToNextDoor()
    {
        double result = getDistanceToNextDoor();
        if(result < _minDistNextDoor)
        {
            _minDistNextDoor = result;
        }
        return result;
    }
    std::vector<std::vector <double> > getBehavior()
    {
        return _behavior;
    }
    std::vector<std::vector <double> > computeFunctionalControllerBehavior(std::vector<double> g);

};


#endif

