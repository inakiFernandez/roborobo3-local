// Gaga: lightweight simple genetic algorithm library
// Copyright (c) Jean Disset 2016, All rights reserved.

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library.

#ifndef GAMULTI_HPP
#define GAMULTI_HPP

/****************************************
 *       TO ENABLE PARALLELISATION
 * *************************************/
// before including this file,
// #define OMP if you want OpenMP parallelisation
// #define CLUSTER if you want MPI parralelisation
// #define CLUSTER if you want MPI parralelisation
#ifdef CLUSTER
#include <mpi.h>
#include <cstring>
#endif
#ifdef OMP
#include <omp.h>
#endif

#include <assert.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <algorithm>
#include <chrono>
#include <cstring>
#include <cstring>
#include <deque>
#include <fstream>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include "include/json.hpp"

#define PURPLE "\033[35m"
#define PURPLEBOLD "\033[1;35m"
#define BLUE "\033[34m"
#define BLUEBOLD "\033[1;34m"
#define GREY "\033[30m"
#define GREYBOLD "\033[1;30m"
#define YELLOW "\033[33m"
#define YELLOWBOLD "\033[1;33m"
#define RED "\033[31m"
#define REDBOLD "\033[1;31m"
#define CYAN "\033[36m"
#define CYANBOLD "\033[1;36m"
#define GREEN "\033[32m"
#define GREENBOLD "\033[1;32m"
#define NORMAL "\033[0m"

namespace GAGA {

using std::vector;
using std::string;
using std::unordered_set;
using std::map;
using std::unordered_map;
using std::cout;
using std::cerr;
using std::endl;
using fpType = std::vector<std::vector<double>>;  // footprints for novelty
using json = nlohmann::json;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;
using std::chrono::system_clock;

/******************************************************************************************
 *                                 GAGA LIBRARY
 *****************************************************************************************/
// This file contains :
// 1 - the Individual class template : an individual's generic representation, with its
// dna, fitnesses and behavior footprints (for novelty)
// 2 - the main GA class template

// About parallelisation :
// before including this file,
// #define OMP if you want OpenMP parallelisation
// #define CLUSTER if you want MPI parallelisation

/*****************************************************************************
 *                         INDIVIDUAL CLASS
 * **************************************************************************/
// A valid DNA class must have (see examples folder):
// DNA mutate()
// DNA crossover(DNA& other)
// static DNA random(int argc, char** argv)
// json& constructor
// void reset()
// json toJson()

template <typename DNA> struct Individual {
	DNA dna;
	map<string, double> fitnesses;  // map {"fitnessCriterName" -> "fitnessValue"}
	fpType footprint;               // individual's footprint for novelty computation
	string infos;                   // custom infos, description, whatever...
	bool evaluated = false;
	bool wasAlreadyEvaluated = false;
	double evalTime = 0.0;
	map<string, double> stats;  // custom stats

	Individual() {}
	explicit Individual(const DNA &d) : dna(d) {}

	explicit Individual(const json &o) {
		assert(o.count("dna"));
		dna = DNA(o.at("dna").dump());
		if (o.count("footprint")) footprint = o.at("footprint").get<fpType>();
		if (o.count("fitnesses")) fitnesses = o.at("fitnesses").get<decltype(fitnesses)>();
		if (o.count("infos")) infos = o.at("infos");
		if (o.count("evaluated")) evaluated = o.at("evaluated");
		if (o.count("alreadyEval")) wasAlreadyEvaluated = o.at("alreadyEval");
		if (o.count("evalTime")) evalTime = o.at("evalTime");
	}

	// Exports individual to json
	json toJSON() const {
		json o;
		o["dna"] = dna.serialize();
		o["fitnesses"] = fitnesses;
		o["footprint"] = footprint;
		o["infos"] = infos;
		o["evaluated"] = evaluated;
		o["alreadyEval"] = wasAlreadyEvaluated;
		o["evalTime"] = evalTime;
		return o;
	}

	// Exports a vector of individual to json
	static json popToJSON(const vector<Individual<DNA>> &p) {
		json o;
		json popArray;
		for (auto &i : p) popArray.push_back(i.toJSON());
		o["population"] = popArray;
		return o;
	}

	// Loads a vector of individual from json
	static vector<Individual<DNA>> loadPopFromJSON(const json &o) {
		assert(o.count("population"));
		vector<Individual<DNA>> res;
		json popArray = o.at("population");
		for (auto &ind : popArray) res.push_back(Individual<DNA>(ind));
		return res;
	}
};

/*********************************************************************************
 *                                 GA CLASS
 ********************************************************************************/
// DNA requirements : see Individual class;
//
// Evaluaor class requirements (see examples folder):
// constructor(int argc, char** argv)
// void operator()(const Individual<DNA>& ind)
// const string name
//
// TYPICAL USAGE :
//
// GA<DNAType, EvalType> ga;
// ga.setPopSize(400);
// return ga.start();

enum class SelectionMethod { paretoTournament, randomObjTournament };
template <typename DNA> class GA {
 protected:

	size_t popSize = 500;                 // nb of individuals in the population
	size_t nbElites = 1;                  // nb of elites to keep accross generations
	size_t nbSavedElites = 1;             // nb of elites to save
	string evaluatorName;                 // name of the given evaluator func
	double crossoverProba = 0.2;          // crossover probability
	double mutationProba = 0.5;           // mutation probablility
    //SelectionMethod selecMethod = SelectionMethod::paretoTournament;

	// for novelty:
	bool novelty = false;             // enable novelty
	double minNoveltyForArchive = 1;  // min novelty for being added to the general archive
	size_t KNN = 15;                  // size of the neighbourhood for novelty

	// for speciation:
	bool speciation = false;           // enable speciation
	double speciationThreshold = 0.2;  // min distance between two dna of same specie
	size_t minSpecieSize = 15;         // minimum specie size
	double minSpeciationThreshold = 0.03;
	double maxSpeciationThreshold = 0.5;
	double speciationThresholdIncrement = 0.01;
	std::function<double(const Individual<DNA> &, const Individual<DNA> &)>
	    indDistanceFunction = [](const auto &, const auto &) { return 0.0; };
	const unsigned int MAX_SPECIATION_TRIES = 100;
	vector<double> speciationThresholds;  // spec thresholds per specie

	/********************************************************************************
	 *                                 SETTERS
	 ********************************************************************************/
 public:
	using Iptr = Individual<DNA> *;
	using DNA_t = DNA;
	void setPopSize(size_t s) { popSize = s; }
	size_t getPopSize() { return popSize; }
	void setNbElites(size_t n) { nbElites = n; }
	size_t getNbElites() { return nbElites; }
	void setNbSavedElites(size_t n) { nbSavedElites = n; }
	void setCrossoverProba(double p) {
		crossoverProba = p <= 1.0 ? (p >= 0.0 ? p : 0.0) : 1.0;
	}
	double getCrossoverProba() { return crossoverProba; }
	void setMutationProba(double p) {
		mutationProba = p <= 1.0 ? (p >= 0.0 ? p : 0.0) : 1.0;
	}
	double getMutationProba() { return mutationProba; }

    // for novelty:
	void enableNovelty() { novelty = true; }
	void disableNovelty() { novelty = false; }
	bool noveltyEnabled() { return novelty; }
	void setKNN(size_t n) { KNN = n; }
	size_t getKNN() { return KNN; }
	void setMinNoveltyForArchive(double m) { minNoveltyForArchive = m; }
	double getMinNoveltyForArchive() { return minNoveltyForArchive; }

	// for speciation:
    /*void enableSpeciation() {
		nextGeneration = [this]() { speciationNextGen(); };
		speciation = true;
	}
	void disableSpeciation() {
		nextGeneration = [this]() { classicNextGen(); };
		speciation = false;
	}

	bool speciationEnabled() { return speciation; }
	void setMinSpeciationThreshold(double s) { minSpeciationThreshold = s; }
	double getMinSpeciationThreshold() { return minSpeciationThreshold; }
	void setMaxSpeciationThreshold(double s) { maxSpeciationThreshold = s; }
	double getMaxSpeciationThreshold() { return maxSpeciationThreshold; }
	void setSpeciationThreshold(double s) { speciationThreshold = s; }
	double getSpeciationThreshold() { return speciationThreshold; }
	void setSpeciationThresholdIncrement(double s) { speciationThresholdIncrement = s; }
	double getSpeciationThresholdIncrement() { return speciationThresholdIncrement; }
	void setMinSpecieSize(double s) { minSpecieSize = s; }
    double getMinSpecieSize() { return minSpecieSize; }*/
	void setIndDistanceFunction(
	    std::function<double(const Individual<DNA> &, const Individual<DNA> &)> f) {
		indDistanceFunction = f;
	}
	vector<vector<Iptr>> species;  // pointers to the individuals of the species

	////////////////////////////////////////////////////////////////////////////////////

	std::random_device rd;
	std::default_random_engine globalRand = std::default_random_engine(rd());

 protected:
	vector<Individual<DNA>>
	    archive;  // when novelty is enabled, we store the novel individuals there
	size_t currentGeneration = 0;
	bool customInit = false;

	// returns a reference (transforms pointer into reference)
	template <typename T> inline T &ref(T &obj) { return obj; }
	template <typename T> inline T &ref(T *obj) { return *obj; }
	template <typename T> inline const T &ref(const T &obj) { return obj; }
	template <typename T> inline const T &ref(const T *obj) { return *obj; }

 public:
	/*********************************************************************************
	 *                              CONSTRUCTOR
	 ********************************************************************************/
    GA(int ac, char **av)
    {
	}

    ~GA()
    {
	}

    void step()
    {
        /*stepBehavior();
        updateFitness();
        broadcast();

        if(endGeneration)
        {
           loadNewGenome();
        }*/



        //savePop();
        //novelty saveArchive();
        //saveParetoFront();
        //updateStats(totalTime); saveGenStats(); saveIndStats();

        ++currentGeneration;
    }


	/*********************************************************************************
	 *                            NEXT POP GETTING READY
	 ********************************************************************************/
    void classicNextGen()
    {
        /*evaluate();
        updateNovelty();
        parent=select();

        if(crossover)
        {
            parent2 = select();
            offspring = parent.crossover(parent2);
        }
        else
        {
            offspring = parent;
        }

        if(mutate)
        {
            offspring.mutate();
        }*/
    }

	// - evaluation de toute la pop, sans se soucier des espèces.
	// - choix des nouveaux représentants parmis les espèces précédentes (clonage)
	// - création d'une nouvelle population via selection/mutation/crossover intra-espece
	// - regroupement en nouvelles espèces en utilisant la distance aux représentants
	// (création d'une nouvelle espèce si distance < speciationThreshold)
	// - on supprime les espèces de taille < minSpecieSize
	// - on rajoute des individus en les faisant muter depuis une espèce aléatoire (nouvelle
	// espèce à chaque tirage) et en les rajoutants à cette espèce
	// - c'est reparti :D

    /*template <typename I>  // I is ither Individual<DNA> or Individual<DNA>*
	vector<Individual<DNA>> produceNOffsprings(size_t n, vector<I> &popu,
                                               size_t nElites = 0)
    {
        std::uniform_real_distribution<double> d(0.0, 1.0);
		vector<Individual<DNA>> nextGen;
		nextGen.reserve(n);

		auto selection = getSelectionMethod<vector<I>>();

		size_t nCross = crossoverProba * (n - s);
            auto *p0 = selection(popu);
			auto *p1 = selection(popu);
			Individual<DNA> offspring(p0->dna.crossover(p1->dna));
			nextGen[i] = offspring;

        size_t nMut = mutationProba * (n - s);

            nextGen[i] = *selection(popu);
			nextGen[i].dna.mutate();

            return nextGen;
	}

    bool paretoDominates(const Individual<DNA> &a, const Individual<DNA> &b) const
    {
		for (auto &o : a.fitnesses) {
			if (!isBetter(o.second, b.fitnesses.at(o.first))) return false;
		}
		return true;
	}

	vector<Individual<DNA> *> getParetoFront(
        const std::vector<Individual<DNA> *> &ind) const
    {
		// naive algorithm. Should be ok for small ind.size()
		assert(ind.size() > 0);
		vector<Individual<DNA> *> pareto;
        for (size_t i = 0; i < ind.size(); ++i)
        {
			bool dominated = false;
			for (auto &j : pareto) {
				if (paretoDominates(*j, *ind[i])) {
					dominated = true;
					break;
				}
        }
			if (!dominated) {
				for (size_t j = i + 1; j < ind.size(); ++j) {
					if (paretoDominates(*ind[j], *ind[i])) {
						dominated = true;
						break;
					}
				}
				if (!dominated) {
					pareto.push_back(ind[i]);
				}
			}
		}
		return pareto;
    }*/

 protected:
	/*********************************************************************************
	 *                          NOVELTY RELATED METHODS
	 ********************************************************************************/
	// Novelty works with footprints. A footprint is just a vector of vector of doubles.
	// It is recommended that those doubles are within a same order of magnitude.
	// Each vector<double> is a "snapshot": it represents the state of the evaluation of
	// one individual at a certain time. Thus, a complete footprint is a combination
	// of one or more snapshot taken at different points in the
	// simulation (a vector<vector<double>>).
	// Snapshot must be of same size accross individuals.
	// Footprint must be set in the evaluator (see examples)

    /*static double getFootprintDistance(const fpType &f0, const fpType &f1)
    {
		assert(f0.size() == f1.size());
		double d = 0;
		for (size_t i = 0; i < f0.size(); ++i) {
			for (size_t j = 0; j < f0[i].size(); ++j) {
				d += std::pow(f0[i][j] - f1[i][j], 2);
			}
		}
		return sqrt(d);
    }*/

	// computeAvgDist (novelty related)
	// returns the average distance of a footprint fp to its k nearest neighbours
	// in an archive of footprints
    /*static double computeAvgDist(size_t K, const vector<Individual<DNA>> &arch,
                                 const fpType &fp)
    {
		double avgDist = 0;
		if (arch.size() > 1) {
			size_t k = arch.size() < K ? static_cast<size_t>(arch.size()) : K;
			vector<Individual<DNA>> knn;
			knn.reserve(k);
			vector<double> knnDist;
			knnDist.reserve(k);
			std::pair<double, size_t> worstKnn = {getFootprintDistance(fp, arch[0].footprint),
			                                      0};  // maxKnn is the worst among the knn
			for (size_t i = 0; i < k; ++i) {
				knn.push_back(arch[i]);
				double d = getFootprintDistance(fp, arch[i].footprint);
				knnDist.push_back(d);
				if (d > worstKnn.first) {
					worstKnn = {d, i};
				}
			}
			for (size_t i = k; i < arch.size(); ++i) {
				double d = getFootprintDistance(fp, arch[i].footprint);
				if (d < worstKnn.first) {  // this one is closer than our worst knn
					knn[worstKnn.second] = arch[i];
					knnDist[worstKnn.second] = d;
					worstKnn.first = d;
					// we update maxKnn
					for (size_t j = 0; j < knn.size(); ++j) {
						if (knnDist[j] > worstKnn.first) {
							worstKnn = {knnDist[j], j};
						}
					}
				}
			}
			assert(knn.size() == k);
			for (size_t i = 0; i < knn.size(); ++i) {
				assert(getFootprintDistance(fp, knn[i].footprint) == knnDist[i]);
				avgDist += knnDist[i];
			}
			avgDist /= static_cast<double>(knn.size());
		}
		return avgDist;
    }*/
    /*void updateNovelty()
    {
        auto savedArchiveSize = archive.size();
        //UpdateArchive	archive.push_back(ind);

		std::pair<Individual<DNA> *, double> best = {&population[0], 0};
		vector<Individual<DNA>> toBeAdded;
        for (auto &ind : population)
        {
			double avgD = computeAvgDist(KNN, archive, ind.footprint);
			bool added = false;
			if (avgD > minNoveltyForArchive) {
				toBeAdded.push_back(ind);
				added = true;
			}
			if (avgD > best.second) best = {&ind, avgD};

			ind.fitnesses["novelty"] = avgD;
		}
		archive.resize(savedArchiveSize);
		archive.insert(std::end(archive), std::begin(toBeAdded), std::end(toBeAdded));

    }*/

/* public:

	void saveArchive() {
		json o = Individual<DNA>::popToJSON(archive);
		o["evaluator"] = evaluatorName;
		std::stringstream baseName;
		baseName << folder << "/gen" << currentGeneration;
		mkdir(baseName.str().c_str(), 0777);
		std::stringstream fileName;
		fileName << baseName.str() << "/archive" << currentGeneration << ".pop";
		std::ofstream file;
		file.open(fileName.str());
		file << o.dump();
		file.close();
    }*/
};
}  // namespace GAGA
#endif
