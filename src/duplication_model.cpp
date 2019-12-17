#include "duplication_model.hpp"
#include <functional>


void InterfaceAssembly::SetBindingStrengths() {
  samming_threshold=static_cast<uint8_t>(interface_size*(1-binding_threshold));
  if(interface_size*(1-binding_threshold) - samming_threshold > .999)
    ++samming_threshold;
  for(size_t i=0;i<=samming_threshold;++i)
      binding_probabilities[i]=std::pow(1-double(i)/interface_size,temperature);
}

void InterfaceAssembly::PrintBindingStrengths() {
  for(auto b : binding_probabilities)
    std::cout<<b<<std::endl;
}    

size_t InterfaceAssembly::Mutation(Genotype& genotype, bool duplication=false, bool insertion=false, bool deletion=false) {
  size_t holder=0;
  if(deletion && std::bernoulli_distribution(deletion_rate)(RNG_Engine)) {
    holder+=2;
    GenotypeDeletion(genotype);
  }
  if(genotype.empty())
    return 0;
  if(duplication && std::bernoulli_distribution(duplication_rate)(RNG_Engine)) {
    holder+=1;
    GenotypeDuplication(genotype);
  }
  if(insertion && std::bernoulli_distribution(insertion_rate)(RNG_Engine)) {
    holder+=1;
    GenotypeInsertion(genotype);
  }

  for(interface_type& base : genotype)
    for(uint8_t nth_bit=0; nth_bit<interface_size; ++nth_bit)
      if(std::bernoulli_distribution(mutation_rate)(RNG_Engine))
        base.flip(nth_bit);
  return holder;
}

double InterfaceAssembly::InteractionMatrix(const interface_type face_1,const interface_type face_2) {
  /*
  #warning "Temporary disable homomeric"
  if(face_1 == face_2)
    return 0;
  */
  return binding_probabilities[interface_model::SammingDistance(face_1,face_2)];
}

void GenotypeInsertion(Genotype& genotype) {
  genotype.reserve(genotype.size()+4);
  std::generate_n(std::back_inserter(genotype), 4, InterfaceAssembly::GenRandomSite);
}

size_t GenotypeDuplication(Genotype& genotype) {
  if(InterfaceAssembly::GetActiveInterfaces(genotype).empty())
    return 0;
  const size_t dup_gene= std::uniform_int_distribution<size_t>(0,genotype.size()/4-1)(RNG_Engine);
  genotype.insert(genotype.end(),genotype.begin()+4*dup_gene,genotype.begin()+4*dup_gene+4);
  return dup_gene;
}

size_t GenotypeDeletion(Genotype& genotype) {
  const size_t del_gene= std::uniform_int_distribution<size_t>(0,genotype.size()/4-1)(RNG_Engine)/4;
  genotype.erase(genotype.begin()+4*del_gene,genotype.begin()+4*del_gene+4);
  return del_gene;
}

bool IsHomologous(const Genotype& genotype, const uint8_t T1, const uint8_t T2) {
  uint16_t homology_count=0;
  for(uint8_t face = 0; face < 4; ++face)
    homology_count+=(genotype[T1*4+face]^genotype[T2*4+face]).count();
  return homology_count/(interface_size*4)>simulation_params::homologous_threshold;
}

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

Genotype StripMonomers2(const Genotype& genotype) {
  std::vector<bool> oligomers(genotype.size()/4);
  for(size_t f1=0;f1<genotype.size();++f1) {
    for(size_t f2=f1;f2<genotype.size();++f2) {
      if(interface_model::SammingDistance(genotype[f1],genotype[f2])<=InterfaceAssembly::samming_threshold) {
        oligomers[f1/4]=true;
        oligomers[f2/4]=true;
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
  Genotype full = active;
  full.insert(full.end(), neutral.begin(), neutral.end());
  //Genotype stripped=full;
  neutral = StripMonomers2(InterfaceAssembly::StripNoncodingGenotype(full));
  //neutral = StripMonomers2(neutral);
  active=full;
  full.insert(full.end(),neutral.begin(),neutral.end());
  //add spare space
  //while(neutral.size() < simulation_params::n_tiles*4)
  //  GenotypeInsertion(neutral);
  if(neutral.size() < simulation_params::n_tiles*4) {
    const size_t n_edges = InterfaceAssembly::GetActiveInterfaces(active).size();
    Genotype spare(simulation_params::n_tiles*4-neutral.size());
    if(n_edges != InterfaceAssembly::GetActiveInterfaces(full).size())
      std::exit(1);
    while(true) {
      Genotype temp = full;
      RandomiseGenotype(spare);
      temp.insert(temp.end(),spare.begin(),spare.end());
    if(InterfaceAssembly::GetActiveInterfaces(temp).size() == n_edges)
      break;
}
    neutral.insert(neutral.end(),spare.begin(),spare.end());
  }
  while(neutral.size() > simulation_params::n_tiles*4) {
    neutral.erase(neutral.begin(),neutral.begin()+4);
  }

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

  std::map<Phenotype_ID,uint16_t> PolyominoAssemblyOutcome(Genotype& binary_genome,FitnessPhenotypeTable* pt,std::map<Phenotype_ID, std::map<InteractionPair,uint16_t> >& pid_interactions) {
    if(binary_genome.empty())
      return {{UNBOUND_pid,1}};
      
    //InterfaceAssembly::StripNoncodingGenotype(binary_genome);
    //Genotype genotype = binary_genome;
    //Genotype genotype = StripMonomers(binary_genome);
    //Genotype genotype = StripMonomers(binary_genome);
    
    Genotype genotype = binary_genome;
    const std::vector<std::pair<InteractionPair,double> > edges = InterfaceAssembly::GetActiveInterfaces(genotype);

    if(edges.empty())
      return {{Phenotype_ID{1,0},1}};

    

    Phenotype phen;
    std::vector<Phenotype_ID> Phenotype_IDs;
    Phenotype_IDs.reserve(pt->phenotype_builds);

    //std::map<Phenotype_ID, std::set<InteractionPair> > phenotype_interactions;

    for(uint16_t nth=0;nth<pt->phenotype_builds;++nth) {
      auto [assembly_information,interacting_indices]=InterfaceAssembly::AssemblePolyomino(edges);
      switch(assembly_information.size()) {
      case 0: //unbound
        pt->ClearIncomplete();
        pid_interactions.clear();
        return {{UNBOUND_pid,1}};
        //case 3: //monomer
        //Phenotype_IDs.emplace_back(Phenotype_ID{1,0});
        //break;
        //[[fallthrough]]
      default:
        phen=GetPhenotypeFromGrid(assembly_information);
        Phenotype_IDs.emplace_back(pt->GetPhenotypeID(phen));
        for(const auto& inter : interacting_indices)
          ++pid_interactions[Phenotype_IDs.back()][inter];
        //pid_interactions[Phenotype_IDs.back()].insert(interacting_indices.begin(),interacting_indices.end());
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
