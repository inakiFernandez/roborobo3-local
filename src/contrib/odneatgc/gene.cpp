#include "odneatgc/gene.h"
#include "odneatgc/genome.h"
#include "odneatgc/innov.h"

#include <iostream>
#include <string.h>
#include <sstream>

using namespace ODNEATGC;

Gene::Gene(double w, NNode *inode, NNode *onode, bool recur, innov innovNum)
{
  lnk = new Link(w, inode, onode, recur);
  innovation_num = innovNum;
  
  enable = true;
}

Gene::Gene(Gene *g,NNode *inode,NNode *onode) 
{
  lnk=new Link((g->lnk)->weight,inode,onode,(g->lnk)->is_recurrent);
  innovation_num=g->innovation_num;
  enable=g->enable;
}

Gene::Gene(const char *argline, std::vector<NNode*> &nodes) 
{
  NNode *inode;
  NNode *onode;
  double weight;
  int recur;
  
  std::vector<NNode*>::iterator curnode;
  
  char inIdR[30], inGc[30];
  char outIdR[30], outGc[30];
  char geneIdR[30], geneGc[30];

  innov innovNodeIn, innovNodeOut;

  //Get the gene parameters
  
  std::stringstream ss(argline);
  ss >> inIdR >> inGc>> outIdR>> outGc>> weight >> recur >> geneIdR >> geneGc >> enable;

  innovNodeIn.idR= atoi(inIdR);
  innovNodeIn.gc= atoi(inGc);

  innovNodeOut.idR= atoi(outIdR);
  innovNodeOut.gc= atoi(outGc);

  innovation_num.idR = atoi(geneIdR);
  innovation_num.gc= atoi(geneGc);

  //Get a pointer to the input node
  curnode=nodes.begin();
  while(!(((*curnode)->node_id)==innovNodeIn))
    ++curnode;
  inode=(*curnode);
  
  //Get a pointer to the output node
  curnode=nodes.begin();
  while(!(((*curnode)->node_id)==innovNodeOut))
    ++curnode;
  onode=(*curnode);
  
  lnk=new Link(weight,inode,onode,recur);

}

Gene::Gene(const Gene& gene)
{
	innovation_num = gene.innovation_num;
	enable = gene.enable;

	lnk = new Link(*gene.lnk);
}

Gene::~Gene() {
	delete lnk;
}

void Gene::print_to_file(std::ostream &outFile) 
{
    outFile<<"gene ";
    outFile<<(lnk->in_node)->node_id.idR<<" ";
    outFile<<(lnk->in_node)->node_id.gc<<" ";
    outFile<<(lnk->out_node)->node_id.idR<<" ";
    outFile<<(lnk->out_node)->node_id.gc<<" ";
    outFile<<(lnk->weight)<<" ";
    outFile<<(lnk->is_recurrent)<<" ";
    outFile<<innovation_num.idR<<" ";
    outFile<<innovation_num.gc<<" ";
    outFile<<enable<<std::endl;

}
