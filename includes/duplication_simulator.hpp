#include "duplication_model.hpp"
#include <functional>

namespace simulation_params
{
  extern uint16_t population_size;
  uint16_t fitness_period=100;
  uint32_t generation_limit=100,independent_trials=1,run_offset=0;
  bool random_initilisation=true;
  double fitness_jump=2,fitness_rise=10;
}

struct PopulationGenotype {
  BGenotype genotype;
  Phenotype_ID pid;
  PopulationGenotype(void) : genotype(simulation_params::n_tiles*4), pid{1,0} {};
};


void DimerModelTable(FitnessPhenotypeTable* pt);

/* Main evolution runners */
uint32_t Evo1();
uint32_t Evo2();
uint32_t DEvo1(uint8_t gap);
uint32_t DEvo2(uint8_t gap);
uint32_t Evo(uint8_t ttype); 
void EvolutionRunner();
void EvolutionRunner2();

//std::vector<uint16_t> RouletteWheelSelection(std::vector<double>& fitnesses);

void SetRuntimeConfigurations(int argc, char* argv[]);

