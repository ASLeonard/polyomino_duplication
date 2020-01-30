#include "duplication_model.hpp"
#include <functional>

//global parameters
namespace simulation_params
{
  uint16_t population_size = 1;
  uint32_t generation_limit = 100, independent_trials = 1, run_offset = 0;
  double fitness_jump = 2;
}

using PhenotypeEdgeInformation = std::vector<std::tuple<InteractionPair, uint8_t, uint8_t, double> >; 

//special structure for an individual to hold their evolution history and assembly details
struct PopulationGenotype {
  static constexpr size_t PID_depth=20;
  Genotype active_space, neutral_space;

  std::map<Phenotype_ID, std::map<InteractionPair,uint16_t> > pid_interactions;

  std::map<Phenotype_ID, std::array<bool,PID_depth> > PID_tracker;
  std::map<Phenotype_ID, std::tuple<uint32_t, uint16_t, PhenotypeEdgeInformation> > PID_details;
  std::set<Phenotype_ID> PID_lineage;
  std::vector<size_t> PID_hierarchy;

  PopulationGenotype() : neutral_space(simulation_params::n_tiles*4+4) {InterfaceAssembly::RandomiseGenotype(neutral_space); active_space.insert(active_space.end(),neutral_space.begin(),neutral_space.begin()+4); neutral_space.erase(neutral_space.begin(),neutral_space.begin()+4);};
};

//Associated functions to track evolutionary histories
void UpdatePhylogenyTrackers(PopulationGenotype& PG, std::vector<std::tuple<Phenotype_ID,uint32_t, uint16_t, PhenotypeEdgeInformation > >& Homology_tracker,uint32_t generation, uint16_t pop_index);


/* Main evolution runners */
void EvolutionRunner();
void EvolvePopulation(const std::string& run_details);

//?
uint32_t DiscoverInteraction(bool self_interaction,bool duplication=false); 
uint32_t DecayInteraction(bool self_interaction, uint8_t gap);
uint32_t DecayDup(uint8_t gap);
void InteractionMetrics();

void EvolveHomology(std::string& run_details,bool self);
void EvolvingHomology();

void SetRuntimeConfigurations(int argc, char* argv[]);

//easy write method to write a binary file and clear the vector
template<typename T, typename A>
void BinaryWriter(std::ofstream& bfile,std::vector<T,A>& vec) {
  bfile.write(reinterpret_cast<const char*>(vec.data()), vec.size()*sizeof(T));
  std::fill(vec.begin(),vec.end(),0);
}
