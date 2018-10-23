#ifndef _ODNEATGCGENOME_H_
#define _ODNEATGCGENOME_H_

#include <vector>
#include "odneatgc/innov.h"
#include <algorithm>
#include "odneatgc/helper.h"

namespace ODNEATGC
{
  class Gene; class NNode; class Network;

  enum mutator { GAUSSIAN = 0,};
  //----------------------------------------------------------------------- 
  //A Genome is the primary source of genotype information used to create
  //a phenotype.  It contains 2 major constituents:
  //  1) A list of NNodes => NNode specifications
  //  2) A list of Genes with Links that point to Traits from (1)
  //(3) Is the primary source of innovation in the evolutionary Genome.     
  //    Each Gene in (3) has a marker telling when it arose historically.   
  //    Thus, these Genes can be used to speciate the population, and the   
  //    list of Genes provide an evolutionary history of innovation and     
  //    link-building.
  
  class Genome 
  {

  public:
    GC genome_id; GC mom_id; GC dad_id;
    int species;

    int nbFitnessUpdates;//average on how many energy measures on current robot
    std::vector<NNode*> nodes; //List of NNodes for the Network
    std::vector<Gene*> genes; //List of innovation-tracking genes

    Network *phenotype; //Allows Genome to be matched with its Network
    
    //Special constructor which spawns off an input file
    //This constructor assumes that some routine has already read in GENOMESTART
    Genome(GC id, std::ifstream &iFile);
    
    //Special constructor that creates a Genome:
    //Fully linked, no hidden nodes with default innovation numbers
    Genome(GC idGenome,int num_in,int num_out);
    //Fully linked, 1 hidden layer with num_hid hidden nodes with default innovation numbers
    Genome(GC id,int num_in,int num_out, int num_hid);
    Genome(GC id, std::vector<NNode*> n, std::vector<Gene*> g);

    //Destructor kills off all lists
    ~Genome();
    
    //Generate a network phenotype from this Genome with specified id
    Network *genesis();
    
    // Dump this genome to specified file
    void print_to_file(std::ostream &outFile);    
    // Wrapper for print_to_file above
    void print_to_filename(const char *filename);
    
    // Duplicate this Genome to create a new one 
    Genome *duplicate();
    
    // ******* MUTATORS *******    
    //Launch all mutations and returns the mutated genome
    //Check out properties file, section NEAT - parameters for the probabilities
    Genome *mutate(float sigma, int idR, GC idNewGenome, int &nodeId, int &genecounter);
    
    // Add Gaussian noise to all linkweights with variance power ^ 2
    void mutate_link_weights(double power);
    int getNumberSynapses(Gene* gene);
    void initialize_link_weights();
    
    double capWeights(double w);
    
    // These last kinds of mutations return false if they fail
    //   They can fail under certain conditions,  being unable
    //   to find a suitable place to make the mutation.
    //   Generally, if they fail, they can be called again if desired. 
    
    // Mutate genome by adding a node
    bool mutate_add_node(int tries,int idR,int &nodeId, int &genecounter);
    
    // Mutate the genome by adding a new link between 2 random NNodes 
    bool mutate_add_link(int tries,int idR, int &genecounter);
    //Toggle the activation of "times" random links (enable to disabled or viceversa)
    void mutate_toggle_enable(unsigned int times);

    // ****** MATING METHODS ***** 
    
    // This method mates this Genome with another Genome g.  
    //   For every point in each Genome, where each Genome shares
    //   the innovation number, the Gene is chosen randomly from 
    //   either parent.  If one parent has an innovation absent in 
    //   the other, the baby will inherit the innovation 
    //   Interspecies mating leads to most genes being inherited.
    //   Otherwise, excess genes come from most fit parent.
    Genome *mate(Genome *g,GC genomeid,double fitness1, double fitness2);
    
    double dissimilarity(Genome *g);
    
  protected:
    //Inserts a NNode into a given ordered list of NNodes in order
    void node_insert(std::vector<NNode*> &nlist, NNode *n);
    
    //Adds a new gene that has been created through a mutation in the
    //*correct order* into the list of genes in the genome
    void add_gene(std::vector<Gene*> &glist,Gene *g);
  };
}// namespace ODNEATGC

#endif
