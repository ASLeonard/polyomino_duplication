#include "duplication_model.hpp"
#include <functional>

namespace simulation_params
{
  extern uint16_t population_size;

  uint32_t generation_limit=100,independent_trials=1,run_offset=0;

  double fitness_jump=2;
  double mutation_rate=.1;
}

struct PopulationGenotype {
  Genotype genotype;
  Phenotype_ID pid;
  //PopulationGenotype(void) : genotype(simulation_params::n_tiles*4), pid{1,0} {std::generate(genotype.begin(),genotype.end(),[] () {return InterfaceAssembly::GenRandomSite();});};
  PopulationGenotype(void) : genotype(simulation_params::n_tiles*4), pid{1,0} {RandomiseGenotype(genotype);};
  
};


void DimerModelTable(FitnessPhenotypeTable* pt);

/* Main evolution runners */
uint32_t Evo1();
uint32_t Evo2();
uint32_t DEvo1(uint8_t gap);
uint32_t DEvo2(uint8_t gap);
uint8_t DEvo3(uint8_t gap);
void EvRu();
uint32_t Evo(uint8_t ttype); 
void EvolutionRunner();
void EvolutionRunner2();

void SetRuntimeConfigurations(int argc, char* argv[]);


template<typename T, typename A>
void BinaryWriter(std::ofstream& bfile,const std::vector<T,A>& vec) {
  bfile.write(reinterpret_cast<const char*>(vec.data()), vec.size()*sizeof(T));
}
