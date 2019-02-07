#include "duplication_model.hpp"
#include <functional>

namespace simulation_params
{
  uint16_t population_size=1;

  uint32_t generation_limit=100,independent_trials=1,run_offset=0;

  double fitness_jump=2;
  
}

struct PopulationGenotype {
  Genotype genotype;
  Phenotype_ID pid;
  PopulationGenotype(void) : genotype(simulation_params::n_tiles*4), pid{1,0} {RandomiseGenotype(genotype);};
  
};


void DimerModelTable(FitnessPhenotypeTable* pt);

/* Main evolution runners */
uint32_t DiscoverInteraction(bool self_interaction); 
uint32_t DecayInteraction(bool self_interaction, uint8_t gap);
void InteractionMetrics();

void EvolutionRunner();

void EvolvingHomology();


void SetRuntimeConfigurations(int argc, char* argv[]);


template<typename T, typename A>
void BinaryWriter(std::ofstream& bfile,const std::vector<T,A>& vec) {
  bfile.write(reinterpret_cast<const char*>(vec.data()), vec.size()*sizeof(T));
}
