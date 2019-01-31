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

size_t InterfaceAssembly::Mutation(Genotype& genotype, bool duplication=false, bool insertion=false, bool deletion=false) {
  size_t holder=0;
  if(duplication && std::bernoulli_distribution(duplication_rate)(RNG_Engine)) {
    holder+=1;
    GenotypeDuplication(genotype);
  }
  if(insertion && std::bernoulli_distribution(insertion_rate)(RNG_Engine))
    GenotypeInsertion(genotype);
  if(deletion && std::bernoulli_distribution(deletion_rate)(RNG_Engine)) {
    holder+=2;
    GenotypeDeletion(genotype);
  }

  for(interface_type& base : genotype)
    for(uint8_t nth_bit=0; nth_bit<interface_size; ++nth_bit)
      if(std::bernoulli_distribution(mutation_rate)(RNG_Engine))
        base.flip(nth_bit);
  return holder;
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

size_t GenotypeDuplication(Genotype& genotype) {
  const size_t dup_gene= std::uniform_int_distribution<size_t>(0,genotype.size())(RNG_Engine)/4;
  genotype.insert(genotype.end(),genotype.begin()+4*dup_gene,genotype.begin()+4*dup_gene+4);
  return dup_gene;
}

size_t GenotypeDeletion(Genotype& genotype) {
  const size_t del_gene= std::uniform_int_distribution<size_t>(0,genotype.size())(RNG_Engine)/4;
  genotype.erase(genotype.begin()+4*del_gene,genotype.begin()+4*del_gene+4);
  return del_gene;
}

bool IsHomologous(const Genotype& genotype, const uint8_t T1, const uint8_t T2) {
  uint16_t homology_count=0;
  for(uint8_t face = 0; face < 4; ++face)
    homology_count+=(genotype[T1*4+face]^genotype[T2*4+face]).count();
  return homology_count/(interface_size*4)>simulation_params::homologous_threshold;
}

Genotype StripMonomers(const Genotype genotype) {
  std::vector<bool> oligomers(genotype.size()/4);
  const std::vector<std::pair<InteractionPair,double> > edges = InterfaceAssembly::GetActiveInterfaces(genotype);
  for(auto& edge : edges) {
    oligomers[edge.first.first/4]=true;
    oligomers[edge.first.second/4]=true;
  }
  Genotype fresh;
  for(size_t nth = 0; nth < oligomers.size(); ++nth)
    if(oligomers[nth])
      fresh.insert(fresh.end(),genotype.begin()+nth*4,genotype.begin()+nth*4+4);
  return fresh;

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

  std::map<Phenotype_ID,uint16_t> PolyominoAssemblyOutcome(Genotype& binary_genome,FitnessPhenotypeTable* pt,std::vector<uint8_t>& homologies) { //std::set<InteractionPair>& pid_interactions __attribute__((unused))) {
    if(binary_genome.empty())
      return {{UNBOUND_pid,1}};
      
    //InterfaceAssembly::StripNoncodingGenotype(binary_genome);
    Genotype genotype = StripMonomers(binary_genome);
    if(genotype.empty())
      return {{Phenotype_ID{1,0},1}};

    const std::vector<std::pair<InteractionPair,double> > edges = InterfaceAssembly::GetActiveInterfaces(genotype);


    std::vector<int8_t> assembly_information;
    Phenotype phen;
    std::vector<Phenotype_ID> Phenotype_IDs;
    Phenotype_IDs.reserve(pt->phenotype_builds);
    std::set<InteractionPair > interacting_indices;
    std::map<Phenotype_ID, std::set<InteractionPair> > phenotype_interactions;

    for(uint16_t nth=0;nth<pt->phenotype_builds;++nth) {
      assembly_information=InterfaceAssembly::AssemblePolyomino(edges,interacting_indices);
      switch(assembly_information.size()) {
      case 0: //unbound
        return {{UNBOUND_pid,1}};
      case 3: //monomer
        //Phenotype_IDs.emplace_back(Phenotype_ID{1,0});
        //break;
        [[fallthrough]]
      default:
        phen=GetPhenotypeFromGrid(assembly_information);
        Phenotype_IDs.emplace_back(pt->GetPhenotypeID(phen));
        phenotype_interactions[Phenotype_IDs.back()].insert(interacting_indices.begin(),interacting_indices.end());
      }
      /*
      if(assembly_information.size()>0) {
        phen=GetPhenotypeFromGrid(assembly_information);
        Phenotype_IDs.emplace_back(pt->GetPhenotypeID(phen));
        phenotype_interactions[Phenotype_IDs.back()].insert(interacting_indices.begin(),interacting_indices.end());
      }
      else {
        pid=UNBOUND_pid;
  	return 0;
	//Phenotype_IDs.emplace_back(0,0);
      }
      */
        
      interacting_indices.clear();
    }

    pt->RelabelPhenotypes(Phenotype_IDs,phenotype_interactions);
    
    std::map<Phenotype_ID,uint16_t> ID_counter=pt->PhenotypeFrequencies(Phenotype_IDs,NULL_pid,true);
    if(ID_counter.find(Phenotype_ID{2,1})!=ID_counter.end()) {
      for(auto external_interaction : phenotype_interactions[Phenotype_ID{2,1}]) {
        //auto external_interaction=*phenotype_interactions[Phenotype_ID{2,1}].begin();
        uint8_t T1=external_interaction.first/4, R1=external_interaction.first%4, T2=external_interaction.second/4, R2=external_interaction.second%4;
        for(uint8_t faces=0; faces<4; ++faces) {
          homologies.emplace_back((genotype[(T1*4)+(R1+faces)%4]^genotype[(T2*4)+(R2+faces)%4]).count());
        }
      }
    }
    return ID_counter;
    /*
    
    if(simulation_params::model_type==1) {
      if(std::find(Phenotype_IDs.begin(),Phenotype_IDs.end(),Phenotype_ID{2,0})!=Phenotype_IDs.end() && std::find(Phenotype_IDs.begin(),Phenotype_IDs.end(),Phenotype_ID{2,1})!=Phenotype_IDs.end())
        pid=Phenotype_ID{2,2};
      else
        pid=std::max_element(ID_counter.begin(),ID_counter.end(),[] (const auto & p1, const auto & p2) {return p1.second < p2.second;})->first;
      
      if(ID_counter.find(Phenotype_ID{2,1})!=ID_counter.end()) {
        for(auto external_interaction : phenotype_interactions[Phenotype_ID{2,1}]) {
          //auto external_interaction=*phenotype_interactions[Phenotype_ID{2,1}].begin();
          uint8_t T1=external_interaction.first/4, R1=external_interaction.first%4, T2=external_interaction.second/4, R2=external_interaction.second%4;
          for(uint8_t faces=0; faces<4; ++faces) {
            homologies.emplace_back((genotype[(T1*4)+(R1+faces)%4]^genotype[(T2*4)+(R2+faces)%4]).count());
          }
        }
        /*
        if(phenotype_interactions[Phenotype_ID{2,1}].size()==1) {
          int x=1;
        }
        else {
          std::cout<<"hmm "<<phenotype_interactions[Phenotype_ID{2,1}].size()<<"\n";
          for(auto x : phenotype_interactions[Phenotype_ID{2,1}])
            std::cout<<+x.first<<" "<<+x.second<<" ";
          std::cout<<"\n";
        }
        */
    /*
      }
      
    }
    else {
      if(!ID_counter.empty())
        pid=std::max_element(ID_counter.begin(),ID_counter.end(),[] (const auto & p1, const auto & p2) {return p1.second < p2.second;})->first;
      else
        pid=NULL_pid;
    }
    */

    
    

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
    //return pt->GenotypeFitness(ID_counter);
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
