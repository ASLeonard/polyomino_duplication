#pragma once
#include <climits>
#include <bitset>

#include "core_genotype.hpp"
#include "core_phenotype.hpp"
#include "core_evolution.hpp"

//compiler optional argument to compile with a fixed length, otherwise defaults to 64
#ifndef GCC_INTERFACE_LENGTH
#define GCC_INTERFACE_LENGTH 64
#endif

constexpr uint8_t interface_size = GCC_INTERFACE_LENGTH;
using interface_type = std::bitset<interface_size>;
using Genotype = std::vector<interface_type>;

//global parameters
namespace simulation_params
{
  inline uint8_t n_tiles = 1;
}

//specialised implementation of polyomino assembly model (from polyomino_core repository)
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
  static void Mutation(Genotype& genotype);
  static void SetBindingStrengths();
  static void PrintBindingStrengths();

  inline static double temperature = 0, binding_threshold = 1;
  inline static double mutation_rate = 0, duplication_rate = 0, deletion_rate = 0;
  inline static uint8_t samming_threshold = 0;
};

namespace interface_model
{   
  interface_type ReverseBits(interface_type v);
  uint8_t SammingDistance(interface_type face1,interface_type face2);

  /*MAIN ASSEMBLY IMPLEMENTATION*/
  std::map<Phenotype_ID,uint16_t> PolyominoAssemblyOutcome(Genotype& binary_genome, FitnessPhenotypeTable* pt, std::map<Phenotype_ID, std::map<InteractionPair,uint16_t> >& pid_interactions);
}

std::vector<uint8_t> CalculateHomology(const Genotype& genotype);

//Associated functions to maintain genotypes during evolution
Genotype StripMonomers(const Genotype& genotype);
void SplitActiveNeutralSpaces(Genotype& active, Genotype& neutral);
void EnsureNeutralDisconnections(Genotype& genotype);
