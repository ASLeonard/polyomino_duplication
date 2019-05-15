#pragma once
#include <climits>
#include <bitset>

#include "core_genotype.hpp"
#include "core_phenotype.hpp"
#include "core_evolution.hpp"

constexpr uint8_t interface_size=128;
using interface_type = std::bitset<interface_size>;

using Genotype = std::vector<interface_type>;



class InterfaceAssembly : public PolyominoAssembly<InterfaceAssembly> {

  protected:
  inline static std::array<double,interface_size+1> binding_probabilities{};
  
public:

  inline static thread_local auto GenRandomSite = []() {
    std::bitset<interface_size> bits;
    std::bernoulli_distribution d(.5);
    for(size_t n = 0; n < interface_size; ++n)
      bits[n] = d(RNG_Engine);
    return bits;
  };
  
  static double InteractionMatrix(const interface_type, const interface_type);
  static size_t Mutation(Genotype& genotype, bool duplication, bool insertion, bool deletion);
  static void SetBindingStrengths();
  static void PrintBindingStrengths();

  inline static double temperature=0,binding_threshold=1;
  inline static double mutation_rate=0,duplication_rate=0,insertion_rate=0,deletion_rate=0;
  inline static uint8_t samming_threshold=0;
  
};


Genotype StripMonomers(const Genotype genotype);
namespace simulation_params
{
  inline uint8_t model_type=1,n_tiles=1;
  inline double homologous_threshold=.1;
}

std::vector<uint8_t> CalculateHomology(const Genotype& genotype);
bool IsHomologous(const Genotype& genotype, const uint8_t T1, const uint8_t T2);

namespace interface_model
{   
  interface_type ReverseBits(interface_type v);
  uint8_t SammingDistance(interface_type face1,interface_type face2);

  /* ASSEMBLY */
  std::map<Phenotype_ID,uint16_t> PolyominoAssemblyOutcome(Genotype& binary_genome, FitnessPhenotypeTable* pt, std::map<Phenotype_ID, std::map<InteractionPair,uint16_t> >& pid_interactions);//std::map<Phenotype_ID, std::set<InteractionPair>>& pid_interactions);
  
}

void EvolutionRunner();
void RandomiseGenotype(Genotype& genotype);
void SplitActiveNeutralSpaces(Genotype& active, Genotype& neutral);
void EvolvePopulation(std::string run_details);


Genotype GenerateTargetGraph(std::map<uint8_t,std::vector<uint8_t>> edge_map,uint8_t graph_size);
void EnsureNeutralDisconnections(Genotype& genotype);

size_t GenotypeDuplication(Genotype& genotype);
void GenotypeInsertion(Genotype& genotype);
size_t GenotypeDeletion(Genotype& genotype);
