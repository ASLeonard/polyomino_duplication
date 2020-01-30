#include "duplication_model.hpp"
#include <functional>

//defines how two binding sites may interact
double InterfaceAssembly::InteractionMatrix(const interface_type face_1,const interface_type face_2) {
  return binding_probabilities[interface_model::SammingDistance(face_1,face_2)];
}

//set binding strengths according to function
void InterfaceAssembly::SetBindingStrengths() {
  samming_threshold = static_cast<uint8_t>(interface_size*(1-binding_threshold));
  //use magic value .999 to avoid floating errors on cutoff boundary
  if(interface_size*(1-binding_threshold) - samming_threshold > .999)
    ++samming_threshold;
  for(size_t i = 0; i <= samming_threshold; ++i)
      binding_probabilities[i] = std::pow(1-double(i)/interface_size,temperature);
}

void InterfaceAssembly::PrintBindingStrengths() {
  for(auto b : binding_probabilities)
    std::cout<<b<<std::endl;
}    

void InterfaceAssembly::Mutation(Genotype& genotype) {
  //Genotype undergoes deletion, return after if empty genotype
  if(std::bernoulli_distribution(deletion_rate)(RNG_Engine)) {
    const size_t del_gene= std::uniform_int_distribution<size_t>(0,genotype.size()/4-1)(RNG_Engine)/4;
    genotype.erase(genotype.begin()+4*del_gene,genotype.begin()+4*del_gene+4);
    if(genotype.empty())
      return;
  }

  //Genotype undergoes duplication, but duplicate monomers has no effect
  if(std::bernoulli_distribution(duplication_rate)(RNG_Engine)) {
    if(!InterfaceAssembly::GetActiveInterfaces(genotype).empty()) {
      const size_t dup_gene= std::uniform_int_distribution<size_t>(0,genotype.size()/4-1)(RNG_Engine);
      genotype.insert(genotype.end(),genotype.begin()+4*dup_gene,genotype.begin()+4*dup_gene+4);
    }
  }  

  //Genotype undergoes mutation, each element on each binary string with a chance to flip 
  for(interface_type& base : genotype)
    for(uint8_t nth_bit=0; nth_bit<interface_size; ++nth_bit)
      if(std::bernoulli_distribution(mutation_rate)(RNG_Engine))
        base.flip(nth_bit);
}

//calculates how sequence identical every pairwise interface in genotype is
std::vector<uint8_t> CalculateHomology(const Genotype& genotype) {
  if(genotype.size()<8)
    return {};
  std::vector<uint8_t> homology;
  homology.reserve(genotype.size()/8 * (genotype.size()-4));
  for(uint8_t tile_1 = 0; tile_1<genotype.size()-4; tile_1+=4)
    for(uint8_t tile_2 = tile_1+4; tile_2<genotype.size(); tile_2+=4)
      for(uint8_t face=0; face < 4; ++face)
        homology.emplace_back((genotype[tile_1+face]^genotype[tile_2+face]).count());
  return homology;
}

//removes all subunits in genotype that do not have at least one interacting edge
Genotype StripMonomers(const Genotype& genotype) {
  std::vector<bool> oligomers(genotype.size()/4);
  for(size_t f1=0;f1<genotype.size();++f1) {
    for(size_t f2=f1;f2<genotype.size();++f2) {
      if(interface_model::SammingDistance(genotype[f1],genotype[f2])<=InterfaceAssembly::samming_threshold) {
        oligomers[f1/4] = true;
        oligomers[f2/4] = true;
      }
    }
  }

  Genotype fresh;
  for(size_t nth = 0; nth < oligomers.size(); ++nth)
    if(!oligomers[nth])
      fresh.insert(fresh.end(),genotype.begin()+nth*4,genotype.begin()+nth*4+4);
  return fresh;
}


void SplitActiveNeutralSpaces(Genotype& active, Genotype& neutral) {
  //Reassess total genotype to find total interactions
  Genotype full = active;
  full.insert(full.end(), neutral.begin(), neutral.end());
  //strip non-interacting parts, store in neutral
  neutral = StripMonomers(InterfaceAssembly::StripNoncodingGenotype(full));
  active = full;
  full.insert(full.end(),neutral.begin(),neutral.end());

  //ensure the neutral space is the same size as requested
  if(neutral.size() < simulation_params::n_tiles*4) {
    const size_t n_edges = InterfaceAssembly::GetActiveInterfaces(active).size();
    Genotype spare(simulation_params::n_tiles*4-neutral.size());
    if(n_edges != InterfaceAssembly::GetActiveInterfaces(full).size())
      std::exit(1);
    //insert raw genotype subunit by subunit, ensuring no new edges are formed
    while(true) {
      Genotype temp = full;
      InterfaceAssembly::RandomiseGenotype(spare);
      temp.insert(temp.end(),spare.begin(),spare.end());
      if(InterfaceAssembly::GetActiveInterfaces(temp).size() == n_edges)
        break;
    }
    neutral.insert(neutral.end(),spare.begin(),spare.end());
  }
  //trim neutral space if too large
  while(neutral.size() > simulation_params::n_tiles*4)
    neutral.erase(neutral.begin(),neutral.begin()+4);
}

namespace interface_model
{
  //reverse binary string
  interface_type ReverseBits(interface_type v) {
    interface_type s;
    for(size_t i = 0; i < interface_size/2; ++i) {
        bool t = v[i];
        s[i] = v[interface_size-i-1];
        s[interface_size-i-1] = t;
    }
    return s;
  }
  //calculate Hamming distance of interface and a reverse of another
  uint8_t SammingDistance(interface_type face1,interface_type face2) {
    return (face1 ^ ReverseBits(~face2)).count();
  }

  std::map<Phenotype_ID,uint16_t> PolyominoAssemblyOutcome(Genotype& binary_genome,FitnessPhenotypeTable* pt,std::map<Phenotype_ID, std::map<InteractionPair,uint16_t> >& pid_interactions) {
    if(binary_genome.empty())
      return {{UNBOUND_pid,1}};
    
    Genotype genotype = binary_genome;
    const std::vector<std::pair<InteractionPair,double> > edges = InterfaceAssembly::GetActiveInterfaces(genotype);

    if(edges.empty())
      return {{Phenotype_ID{1,0},1}};

    Phenotype phen;
    std::vector<Phenotype_ID> Phenotype_IDs;
    Phenotype_IDs.reserve(pt->phenotype_builds);

    for(uint16_t nth=0;nth<pt->phenotype_builds;++nth) {
      auto [assembly_information,interacting_indices]=InterfaceAssembly::AssemblePolyomino(edges);
      switch(assembly_information.size()) {
      case 0: //unbound
        pt->ClearIncomplete();
        pid_interactions.clear();
        return {{UNBOUND_pid,1}};
      default:
        phen=GetPhenotypeFromGrid(assembly_information);
        Phenotype_IDs.emplace_back(pt->GetPhenotypeID(phen));
        for(const auto& inter : interacting_indices)
          ++pid_interactions[Phenotype_IDs.back()][inter];
      }        
      interacting_indices.clear();
    }
    
    pt->RelabelPIDs(Phenotype_IDs);
    pt->RelabelMaps(pid_interactions,true);
    pt->UpdateFitnesses();
    
    std::map<Phenotype_ID,uint16_t> ID_counter=pt->PhenotypeFrequencies(Phenotype_IDs);

    for(auto iter= pid_interactions.begin(); iter!= pid_interactions.end();) {
      if (ID_counter.find(iter->first)==ID_counter.end()) {
        iter = pid_interactions.erase(iter);
      } else {
        ++iter;
      }
    }
    return ID_counter;
  }
}
