#include "duplication_model.hpp"
#include <functional>


void InterfaceAssembly::SetBindingStrengths() {
  samming_threshold=static_cast<uint8_t>(interface_size*(1-binding_threshold));
  for(size_t i=0;i<=samming_threshold;++i)
      binding_probabilities[i]=std::pow(1-double(i)/interface_size,temperature);
}

void InterfaceAssembly::PrintBindingStrengths() {
  for(auto b : binding_probabilities)
    std::cout<<b<<std::endl;
}    

void InterfaceAssembly::Mutation(Genotype& genotype, bool duplication=false, bool insertion=false, bool deletion=false) {
  if(duplication && std::bernoulli_distribution(duplication_rate)(RNG_Engine))
    GenotypeDuplication(genotype);
  if(insertion && std::bernoulli_distribution(insertion_rate)(RNG_Engine))
    GenotypeInsertion(genotype);
  if(deletion && std::bernoulli_distribution(deletion_rate)(RNG_Engine))
    GenotypeDeletion(genotype);

  for(interface_type& base : genotype)
    for(uint8_t nth_bit=0; nth_bit<interface_size; ++nth_bit)
      if(std::bernoulli_distribution(mutation_rate)(RNG_Engine))
        base.flip(nth_bit);  
}

double InterfaceAssembly::InteractionMatrix(const interface_type face_1,const interface_type face_2) {
  return binding_probabilities[interface_model::SammingDistance(face_1,face_2)];
}

void GenotypeInsertion(Genotype& genotype) {
  const size_t N_edges = InterfaceAssembly::GetActiveInterfaces(genotype).size();
  Genotype new_gene(4);
  genotype.insert(genotype.end(),new_gene.begin(),new_gene.end());
  do {
    RandomiseGenotype(new_gene);
    std::move(new_gene.begin(),new_gene.end(),genotype.rbegin());
  } while (N_edges != InterfaceAssembly::GetActiveInterfaces(genotype).size());

}

void GenotypeDuplication(Genotype& genotype) {
  const size_t dup_gene= std::uniform_int_distribution<size_t>(0,genotype.size())(RNG_Engine)/4;
  genotype.insert(genotype.end(),genotype.begin()+4*dup_gene,genotype.begin()+4*dup_gene+4);
}

void GenotypeDeletion(Genotype& genotype) {
  if(genotype.size()<=4)
    return;
  const size_t del_gene= std::uniform_int_distribution<size_t>(0,genotype.size())(RNG_Engine)/4;
  genotype.erase(genotype.begin()+4*del_gene,genotype.begin()+4*del_gene+4);
}

namespace interface_model
{
  interface_type ReverseBits(interface_type v) {
    interface_type s;
    for(size_t i = 0; i < interface_size/2; ++i) {
        bool t = v[i];
        s[i] = v[interface_size-i-1];
        s[interface_size-i-1] = t;
    }
    return s;
  }
  
  uint8_t SammingDistance(interface_type face1,interface_type face2) {
    return (face1 ^ ReverseBits(~face2)).count();
  }

  double PolyominoAssemblyOutcome(Genotype& binary_genome,FitnessPhenotypeTable* pt,Phenotype_ID& pid,std::set<InteractionPair>& pid_interactions) {
    InterfaceAssembly::StripNoncodingGenotype(binary_genome);
    const std::vector<std::pair<InteractionPair,double> > edges = InterfaceAssembly::GetActiveInterfaces(binary_genome);


    std::vector<int8_t> assembly_information;
    Phenotype phen;
    std::vector<Phenotype_ID> Phenotype_IDs;
    Phenotype_IDs.reserve(pt->phenotype_builds);
    std::set<InteractionPair > interacting_indices;
    std::map<Phenotype_ID, std::set<InteractionPair> > phenotype_interactions;

    for(uint16_t nth=0;nth<pt->phenotype_builds;++nth) {
      assembly_information=InterfaceAssembly::AssemblePolyomino(edges,interacting_indices);
      if(assembly_information.size()>0) {
        phen=GetPhenotypeFromGrid(assembly_information);
        Phenotype_IDs.emplace_back(pt->GetPhenotypeID(phen));
        phenotype_interactions[Phenotype_IDs.back()].insert(interacting_indices.begin(),interacting_indices.end());
      }
      else
        Phenotype_IDs.emplace_back(0,0);
      interacting_indices.clear();
    }

    pt->RelabelPhenotypes(Phenotype_IDs,phenotype_interactions);
    
    std::map<Phenotype_ID,uint16_t> ID_counter=pt->PhenotypeFrequencies(Phenotype_IDs);

    if(!ID_counter.empty())
      pid=std::max_element(ID_counter.begin(),ID_counter.end(),[] (const auto & p1, const auto & p2) {return p1.second < p2.second;})->first;
    else
      pid=NULL_pid;
    

    //pid_interactions=phenotype_interactions[pid];

    /*
    if(simulation_params::model_type==1)
      return pt->SingleFitness(pid,ID_counter[pid]);
    if(simulation_params::model_type==2) {
      for(auto kv : ID_counter) {
        if(kv.first!=NULL_pid && kv.first!=pid)
          pid_interactions.merge(phenotype_interactions[kv.first]);
      }

    }
    */
    return pt->GenotypeFitness(ID_counter);
  }


}//end interface_model namespace

void RandomiseGenotype(Genotype& genotype) {
  do {
    std::generate(genotype.begin(),genotype.end(),InterfaceAssembly::GenRandomSite);
  }while(!InterfaceAssembly::GetActiveInterfaces(genotype).empty());
}

Genotype GenerateTargetGraph(std::map<uint8_t,std::vector<uint8_t>> edge_map,uint8_t graph_size) {
  const uint8_t total_edges=std::accumulate(edge_map.begin(),edge_map.end(),0,[](uint8_t size,const auto & p1) {return size+p1.second.size();});
  Genotype graph(graph_size);
  
  std::uniform_int_distribution<uint8_t> delta_ser(0,InterfaceAssembly::samming_threshold);
  std::vector<uint8_t> bits(interface_size);
  std::iota(bits.begin(),bits.end(),0);
  constexpr uint8_t shift_r=interface_size/2;
  
  do {
    RandomiseGenotype(graph); 
    for(auto edge : edge_map) {
      graph[edge.first]=InterfaceAssembly::GenRandomSite();
      for(uint8_t connector : edge.second) {
        /*!make perfectly self-interacting*/
        if(edge.first!=connector)
          graph[connector]=interface_model::ReverseBits(~graph[edge.first]);
        else 
          graph[connector]=(interface_model::ReverseBits(~(graph[edge.first]>>shift_r))>>shift_r) | ((graph[edge.first]>>shift_r)<<shift_r);

        std::shuffle(bits.begin(),bits.end(),RNG_Engine);
        const uint8_t delta_s = delta_ser(RNG_Engine)/((edge.first==connector) ? 2 : 1);
        for(uint8_t b=0; b<delta_s;++b)
          graph[connector] ^=(interface_type(1)<<bits[b]);
      }
    }
  }while(InterfaceAssembly::GetActiveInterfaces(graph).size()!=total_edges);
  return graph;
}
