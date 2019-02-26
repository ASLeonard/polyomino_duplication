#include "duplication_model.hpp"
#include <functional>

namespace simulation_params
{
  uint16_t population_size=1;

  uint32_t generation_limit=100,independent_trials=1,run_offset=0;

  double fitness_jump=2;
  
}

struct PopulationGenotype {
  static constexpr size_t PID_depth=25;
  Genotype genotype, subunits;
  std::vector<Phenotype_ID> pids;
  std::map<Phenotype_ID, std::array<bool,PID_depth> > PID_tracker;
  std::map<Phenotype_ID, std::pair<std::pair<size_t,size_t>, size_t> > PID_info;
  std::map<Phenotype_ID, std::vector<std::tuple<InteractionPair,size_t,size_t,double> > > PID_deet;

  PopulationGenotype(void) : subunits(simulation_params::n_tiles*4) {genotype.resize(simulation_params::n_tiles*4); RandomiseGenotype(genotype);};
  
};


void DimerModelTable(FitnessPhenotypeTable* pt);

/* Main evolution runners */
uint32_t DiscoverInteraction(bool self_interaction); 
uint32_t DecayInteraction(bool self_interaction, uint8_t gap);
void InteractionMetrics();

void EvolutionRunner();

void EvolveHomology(std::string run_details,bool self);

void EvolvingHomology();


void SetRuntimeConfigurations(int argc, char* argv[]);
//void UpdatePhylogenyTrackers(PopulationGenotype& PG, std::map<Phenotype_ID, std::map<std::pair<size_t,size_t>, size_t> >& Homology_tracker,size_t a, size_t b);


template<typename T, typename A>
void BinaryWriter(std::ofstream& bfile,std::vector<T,A>& vec) {
  bfile.write(reinterpret_cast<const char*>(vec.data()), vec.size()*sizeof(T));
  std::fill(vec.begin(),vec.end(),0);
}
