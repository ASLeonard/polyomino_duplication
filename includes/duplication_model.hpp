#pragma once
#include <climits>
#include <bitset>

#include "core_genotype.hpp"
#include "core_phenotype.hpp"
#include "core_evolution.hpp"

constexpr uint8_t interface_size=128;
using interface_type = std::bitset<interface_size>;

using BGenotype = std::vector<interface_type>;

template< size_t size>
typename std::bitset<size> random_bitset() {
    typename std::bitset<size> bits;
    std::bernoulli_distribution d(.5);
    for(size_t n = 0; n < size; ++n) 
        bits[n] = d(RNG_Engine);
    return bits;
}


class InterfaceAssembly : public PolyominoAssembly<InterfaceAssembly> {

  protected:
  inline static std::array<double,interface_size+1> binding_probabilities{};
  
public:

  inline static thread_local auto GenRandomSite = []() {return random_bitset<interface_size>();};
  
  static double InteractionMatrix(const interface_type, const interface_type);
  static void Mutation(BGenotype& genotype);

  static void SetBindingStrengths();
  static void PrintBindingStrengths();
  
};



namespace simulation_params
{
  extern uint8_t model_type,n_tiles,samming_threshold;
  extern uint16_t dissociation_time;
  extern double temperature,binding_threshold,mu_prob;
}

namespace interface_model
{   
  interface_type ReverseBits(interface_type v);
  uint8_t SammingDistance(interface_type face1,interface_type face2);

  /* ASSEMBLY */
  double PolyominoAssemblyOutcome(BGenotype& binary_genome, FitnessPhenotypeTable* pt,Phenotype_ID& pid,std::set<InteractionPair>& pid_interactions);
  
}


BGenotype GenerateTargetGraph(std::map<uint8_t,std::vector<uint8_t>> edge_map,uint8_t graph_size);
void EnsureNeutralDisconnections(BGenotype& genotype);

