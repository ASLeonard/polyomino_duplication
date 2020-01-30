#include "duplication_simulator.hpp"
#include <iostream>

#ifndef FULL_WRITE
#define FULL_WRITE 0
#endif

#define STRINGIZE(x) #x
#define STRINGIZE_VALUE_OF(x) STRINGIZE(x)

#ifndef ROOT_FILE_PATH
#define ROOT_FILE_PATH
#endif

const std::string file_base_path = STRINGIZE_VALUE_OF(ROOT_FILE_PATH);



void InteractionMetrics() {
  const uint32_t N_runs=simulation_params::independent_trials;
  
  std::ofstream f_out(file_base_path+"Discovery_"+std::to_string(InterfaceAssembly::binding_threshold)+".BIN", std::ios::binary);
  std::vector<uint32_t> res_S(N_runs),res_A(N_runs);
/*
#pragma omp parallel for schedule(dynamic) 
  for(uint32_t r=0;r < N_runs;++r) {
    res_S[r]= DiscoverInteraction(true,false);
    res_A[r]= DiscoverInteraction(false,false);
  }    
  BinaryWriter(f_out,res_S);
  BinaryWriter(f_out,res_A); 
  return;
    */
  std::ofstream f_out2(file_base_path+"Decay_"+std::to_string(InterfaceAssembly::binding_threshold)+".BIN", std::ios::binary);
  for(uint8_t gap=0;gap<=InterfaceAssembly::samming_threshold/2;++gap) {
#pragma omp parallel for schedule(dynamic) 
    for(uint32_t r=0;r < N_runs;++r) {
      //if(gap%2==0 && InterfaceAssembly::samming_threshold%2==0)
      //  res_S[r] = DecayDup(true,gap/2);
      res_A[r]= DecayDup(gap);
    }
    //BinaryWriter(f_out2,res_S);
    BinaryWriter(f_out2,res_A);
  }  
}

void EvolvingHomology() {
  const uint32_t N_runs=simulation_params::independent_trials;
  std::ofstream f_out(file_base_path+"EHom_"+std::to_string(InterfaceAssembly::binding_threshold)+".BIN", std::ios::binary);
  for(uint8_t self=0;self<2;++self) {
#pragma omp parallel for schedule(dynamic) 
    for(uint32_t r=0;r < N_runs;++r) {
      std::vector<uint8_t> res;
      res.reserve(simulation_params::generation_limit*4);

      std::uniform_int_distribution<uint8_t> dis(0, interface_size-1);
      std::uniform_int_distribution<uint8_t> dis2(0, 7);

    
      Genotype g(8);
      InterfaceAssembly::RandomiseGenotype(g);
      if(self) {
        for(size_t n=0;n<interface_size/2;++n)
          g[0][interface_size-n]=~g[0][n];
        std::move(g.begin(),g.begin()+4,g.begin()+4);
    
      }
      else
        g[4]=interface_model::ReverseBits(~g[0]);
    
      for(uint32_t generation=1;generation<=simulation_params::generation_limit;++generation) {
        auto hom = CalculateHomology(g);
        res.insert(res.end(),hom.begin(),hom.end());
        g[dis2(RNG_Engine)].flip(dis(RNG_Engine));
      }
#pragma omp critical(wrout)
      {
        BinaryWriter(f_out,res);
      }
    }
  }
    
}

uint32_t DiscoverInteraction(bool self_interaction,bool duplicated) {
  Genotype genotype(2);
  InterfaceAssembly::RandomiseGenotype(genotype);
  genotype[duplicated]=genotype[0];

  std::uniform_int_distribution<size_t> index_dist(0, interface_size-1);
  for(uint32_t generation=1;generation<=simulation_params::generation_limit;++generation) {
    genotype[0].flip(index_dist(RNG_Engine));
    if(interface_model::SammingDistance(genotype[0],genotype[!self_interaction])<=InterfaceAssembly::samming_threshold)
      return generation;
  }

  return 0;
}

uint32_t DecayInteraction(bool self_interaction, uint8_t gap) {
  interface_type geno1=InterfaceAssembly::GenRandomSite(), geno2=interface_model::ReverseBits(~geno1);
  if(self_interaction)
    for(size_t n=0;n<interface_size/2;++n)
      geno1[interface_size-1-n]=~geno1[n];

  std::uniform_int_distribution<uint8_t> dis(0, interface_size-1);
  std::vector<uint8_t> bits(interface_size/(self_interaction?2:1));
  std::iota(bits.begin(),bits.end(),0);
  std::shuffle(bits.begin(),bits.end(),RNG_Engine);

  for(uint8_t b=0; b<gap;++b)
    geno1.flip(bits[b]);

  for(uint32_t generation=1;generation<=simulation_params::generation_limit;++generation) {
    geno1.flip(dis(RNG_Engine));
    if(interface_model::SammingDistance(geno1,self_interaction?geno1:geno2)>InterfaceAssembly::samming_threshold)
      return generation;
  }
  
  return 0;

}

uint32_t DecayDup(uint8_t gap) {
  interface_type geno1=InterfaceAssembly::GenRandomSite(), geno2=interface_model::ReverseBits(~geno1);
  for(size_t n=0;n<interface_size/2;++n)
    geno1[interface_size-1-n]=~geno1[n];

  std::uniform_int_distribution<uint8_t> dis(0, interface_size-1);
  std::vector<uint8_t> bits(interface_size/2);
  std::iota(bits.begin(),bits.end(),0);
  std::shuffle(bits.begin(),bits.end(),RNG_Engine);
  std::bernoulli_distribution B_d(0.5);

  for(uint8_t b=0; b<gap;++b)
    geno1.flip(bits[b]);
  geno2 = geno1;

  for(uint32_t generation=1;generation<=simulation_params::generation_limit;++generation) {
    if(B_d(RNG_Engine))
      geno1.flip(dis(RNG_Engine));
    else
      geno2.flip(dis(RNG_Engine));

    if(interface_model::SammingDistance(geno1,geno2)<=InterfaceAssembly::samming_threshold) {
      if(interface_model::SammingDistance(geno1,geno1)>InterfaceAssembly::samming_threshold && interface_model::SammingDistance(geno2,geno2)>InterfaceAssembly::samming_threshold)
        return 1;
  }
    else
      return 0;
  }
  
  return -1;

}

void EvolveHomology(std::string& run_details,bool self) {
  std::string file_simulation_details=run_details+".BIN";
  std::ofstream fout_homology(file_base_path+"Bomology"+file_simulation_details,std::ios::binary); 

  std::vector<double> population_fitnesses(simulation_params::population_size);
  std::vector<PopulationGenotype> evolving_population(simulation_params::population_size),reproduced_population;
  reproduced_population.resize(simulation_params::population_size);

  Genotype g(8);
  InterfaceAssembly::RandomiseGenotype(g);
  if(!self) {
    for(size_t n=0;n<interface_size/2;++n)
      g[0][interface_size-n-1]=~g[0][n];
    std::move(g.begin(),g.begin()+4,g.begin()+4);
  }
  else
    g[4]=interface_model::ReverseBits(~g[0]);

  for(auto& ep : evolving_population)
    ep.active_space=Genotype(g);

  std::uniform_int_distribution<uint8_t> dis(0, interface_size-1);
  std::uniform_int_distribution<uint8_t> dis2(0, 7);

  std::vector<uint8_t> res;
  res.reserve(simulation_params::generation_limit*4*simulation_params::population_size);
  
  for(uint32_t generation=0;generation<simulation_params::generation_limit;++generation) { /*! MAIN EVOLUTION LOOP */
    
    uint16_t nth_genotype=0;
    for(PopulationGenotype& evolving_genotype : evolving_population) { /*! GENOTYPE LOOP */
      auto hom = CalculateHomology(evolving_genotype.active_space);
      res.insert(res.end(),hom.begin(),hom.end());
           
      evolving_genotype.active_space[dis2(RNG_Engine)].flip(dis(RNG_Engine));
      population_fitnesses[nth_genotype]=0;
      const std::vector<std::pair<InteractionPair,double> > edges = InterfaceAssembly::GetActiveInterfaces(evolving_genotype.active_space);
      for(auto edge : edges) {
        if(edge.first==InteractionPair{0,4})
          population_fitnesses[nth_genotype]=1;
        if(edge.first!=InteractionPair{0,0} && edge.first!=InteractionPair{0,4} && edge.first!=InteractionPair{4,4}) {
          population_fitnesses[nth_genotype]=0;
          break;
        }
      }
      ++nth_genotype;
    } /*! END GENOTYPE LOOP */

    /*! SELECTION */
    uint16_t nth_repro=0;
    auto reproducing_selection=RouletteWheelSelection(population_fitnesses);
    for(uint16_t selected : reproducing_selection) {
      reproduced_population[nth_repro++]=evolving_population[selected];
    }
    evolving_population.swap(reproduced_population);    
  }
  BinaryWriter(fout_homology,res);  
}

void EvolutionRunner() {
  /*!PYTHON INFORMATION*/
  Phenotype::DETERMINISM_LEVEL=1;

  const uint16_t N_runs=simulation_params::independent_trials;
#pragma omp parallel for schedule(dynamic) 
  for(uint16_t r=0;r < N_runs;++r)
    EvolvePopulation("_Run"+std::to_string(r+simulation_params::run_offset));
}

void UpdatePhylogenyTrackers(PopulationGenotype& PG, std::vector<std::tuple<Phenotype_ID,uint32_t, uint16_t, PhenotypeEdgeInformation > >& Homology_tracker,uint32_t generation, uint16_t pop_index) {
  //shift all observations 1 generation down
  for(auto& kv : PG.PID_tracker) {
    kv.second.front()=false;
    std::rotate(kv.second.begin(),kv.second.begin()+1,kv.second.end());
  }

  //increment tracking for all pids, store original discovery for new ones
  for(const auto& pid_kv : PG.pid_interactions) {
    const auto pid = pid_kv.first;
    if(pid==UNBOUND_pid || pid==NULL_pid || pid==Phenotype_ID{1,0})
      continue;
    if(PG.PID_details.find(pid)==PG.PID_details.end()) {
      PhenotypeEdgeInformation edge_info;
      for(auto edge : PG.pid_interactions[pid]) {
        const auto hom = (PG.active_space[edge.first.first] ^ PG.active_space[edge.first.second]).count();
        const auto str = interface_model::SammingDistance(PG.active_space[edge.first.first], PG.active_space[edge.first.second]);
        edge_info.emplace_back(edge.first,hom,str,1);//edge.second/weighting);
      }

      PG.PID_details[pid]= std::make_tuple(generation,pop_index,edge_info);
    }

    PG.PID_tracker[pid].back()=true;
  }

  //remove outdated observations from tracking
  for(auto iter = PG.PID_tracker.begin(); iter!=PG.PID_tracker.end();) {
    if(std::find(iter->second.begin(),iter->second.end(),true)!=iter->second.end()) 
      ++iter;
    else {
      PG.PID_details.erase(iter->first);
      iter=PG.PID_tracker.erase(iter);
    }
  }
    
  //iterate through all currently tracked phenotypes
  for(const auto& [pid, pid_record] : PG.PID_tracker) {
   
    //this pid has survived long enough to be meaningful
    if(std::accumulate(pid_record.begin(),pid_record.end(),0)==PG.PID_depth) {

      //if this pid is not in the individuals lineage, add it and update global record
      if(PG.PID_lineage.find(pid)==PG.PID_lineage.end() && (PG.PID_lineage.empty() || pid.first>PG.PID_lineage.crbegin()->first)) {
        PG.PID_lineage.emplace(pid);

        auto combined_details=std::tuple_cat(std::tie(pid),PG.PID_details[pid]);
        auto iter=std::find(Homology_tracker.begin(),Homology_tracker.end(),combined_details);
        size_t hierarchy_index;
        if(iter==Homology_tracker.end()) {
          hierarchy_index=Homology_tracker.size();
          Homology_tracker.emplace_back(combined_details); 
        }
        else
          hierarchy_index = std::distance(Homology_tracker.begin(), iter);
        PG.PID_hierarchy.emplace_back(hierarchy_index);

      }
    }
  }
}


void EvolvePopulation(const std::string& run_details) {
  std::string file_simulation_details=run_details+".txt";

  std::ofstream fout_evo(file_base_path+"EvoRecord_Mu"+std::to_string(InterfaceAssembly::mutation_rate)+"_S"+std::to_string(InterfaceAssembly::binding_threshold)+"_D" + std::to_string(InterfaceAssembly::duplication_rate)+".txt",std::ios::app );
  
  std::string fname_phenotype(file_base_path+"PhenotypeTable"+file_simulation_details);

  std::ofstream fout_selection_history, fout_phenotype_IDs, fout_size, fout_interactions, fout_strength, fout_homology;

  if(FULL_WRITE) {
    fout_selection_history.open(file_base_path+"Selections"+file_simulation_details,std::ios::binary);
    fout_phenotype_IDs.open(file_base_path+"PIDs"+file_simulation_details,std::ios::out );
    fout_size.open(file_base_path+"Size"+file_simulation_details,std::ios::out);
    fout_interactions.open(file_base_path+"Interactions"+file_simulation_details,std::ios::out);
    fout_strength.open(file_base_path+"Strengths"+file_simulation_details,std::ios::binary);
    fout_homology.open(file_base_path+"Homology"+file_simulation_details,std::ios::binary);
  }  

  std::vector<std::tuple<Phenotype_ID, uint32_t,uint16_t, PhenotypeEdgeInformation > > evolution_record;
  std::set< std::vector<size_t> > global_sets;
    
  std::vector<double> population_fitnesses(simulation_params::population_size);
  std::vector<PopulationGenotype> evolving_population(simulation_params::population_size),reproduced_population;
  reproduced_population.resize(simulation_params::population_size);
  std::vector<uint16_t> binary_homologies(interface_size+1), binary_strengths(interface_size+1);
  std::map<std::pair<int,int>, std::vector<uint16_t>> binary_homology;

  for(int i =0; i<4;++i)
    for(int j=0; j<4; ++j)
      binary_homology[std::make_pair(i,j)] = binary_homologies;
  
  
  FitnessPhenotypeTable pt = FitnessPhenotypeTable();
  pt.fit_func=[](double s) {return s;};


  for(uint32_t generation=0;generation<simulation_params::generation_limit;++generation) { /*! MAIN EVOLUTION LOOP */

    uint16_t nth_genotype=0;
    for(PopulationGenotype& evolving_genotype : evolving_population) { /*! GENOTYPE LOOP */
      
      if(evolving_genotype.active_space.empty()) {
        std::cout<<"NULL GENOTYPE AT START OF GENERATION\ngen "<<generation<<" nth "<<nth_genotype<<"\n";
        return;
      }
      if(evolving_genotype.active_space.size()>100) {
        std::cout<<"TOO LONG GENOTYPE\ngen "<<generation<<" nth "<<nth_genotype<<"\n";
        return;
      }
      
      InterfaceAssembly::Mutation(evolving_genotype.active_space);

      evolving_genotype.pid_interactions.clear();

      SplitActiveNeutralSpaces(evolving_genotype.active_space,evolving_genotype.neutral_space);
      auto pid_map=interface_model::PolyominoAssemblyOutcome(evolving_genotype.active_space,&pt,evolving_genotype.pid_interactions);
      
      switch(pid_map.begin()->first.first) {
      case 255:
        population_fitnesses[nth_genotype]=0;
        break;
      case 1:
        population_fitnesses[nth_genotype]=1;
        break;
      default:
        population_fitnesses[nth_genotype]=pt.GenotypeFitness(pid_map);
      }      

      UpdatePhylogenyTrackers(evolving_genotype,evolution_record,generation,nth_genotype);
      if(!evolving_genotype.PID_hierarchy.empty())
        global_sets.emplace(evolving_genotype.PID_hierarchy);
      
      for(const auto& pid_kv : evolving_genotype.pid_interactions) {
        for(const auto& ints_kv : pid_kv.second) {
          auto IP = ints_kv.first;
          binary_homologies[(evolving_genotype.active_space[IP.first]^evolving_genotype.active_space[IP.second]).count()]+=ints_kv.second;
          binary_homology[std::minmax(IP.first%4,IP.second%4)][(evolving_genotype.active_space[IP.first]^evolving_genotype.active_space[IP.second]).count()]+=ints_kv.second;

          binary_strengths[interface_size-interface_model::SammingDistance(evolving_genotype.active_space[IP.first],evolving_genotype.active_space[IP.second])]+=ints_kv.second;
        }
      }
           
      ++nth_genotype;
      
      if(FULL_WRITE) {
        for(auto& kv : pid_map) 
          fout_phenotype_IDs << +kv.first.first << " " << +kv.first.second << " ";
        for(auto edge : InterfaceAssembly::GetActiveInterfaces(evolving_genotype.active_space))
          fout_interactions << +edge.first.first << " " << +edge.first.second << " ";
        fout_interactions << ",";
      
        fout_phenotype_IDs << ",";
        fout_size << evolving_genotype.active_space.size()/4 << " ";
      }
      
    } /*! END GENOTYPE LOOP */
    
    constexpr uint8_t BREAK_SIZE_LIMIT = 10;
    if(!FULL_WRITE && std::find_if(evolution_record.begin(), evolution_record.end(),[=](auto record){return std::get<0>(record).first >= BREAK_SIZE_LIMIT;})!=evolution_record.end())
      break;

    /*! SELECTION */
    uint16_t nth_repro=0;
    auto reproducing_selection=RouletteWheelSelection(population_fitnesses);
    for(uint16_t selected : reproducing_selection) {
      reproduced_population[nth_repro++]=evolving_population[selected];
    }
    evolving_population.swap(reproduced_population);
    

    if(FULL_WRITE) {
      BinaryWriter(fout_selection_history,reproducing_selection);
      fout_phenotype_IDs << "\n";
      fout_size << "\n";
      fout_interactions << "\n";
      for(int i =0; i<4;++i)
        for(int j=0; j<4; ++j)
          BinaryWriter(fout_homology,binary_homology[std::make_pair(i,j)]);
      BinaryWriter(fout_strength,binary_strengths);
    }

    
    
  } /* END EVOLUTION LOOP */

#pragma omp critical(evo_record_write)
  {
  for(const auto& kv : evolution_record) {
    const auto pid = std::get<0>(kv);
    fout_evo<<+pid.first<<" "<<+pid.second<<" "<<std::get<1>(kv)<<" "<<std::get<2>(kv)<<" ";
    for(const auto& pv : std::get<3>(kv))
      fout_evo<<+std::get<0>(pv).first<<" "<<+std::get<0>(pv).second<<" "<<+std::get<1>(pv)<<" "<<+std::get<2>(pv)<<" "<<std::get<3>(pv)<<" ";
    fout_evo<<"\n";
    }
  for(const auto& evo_route : global_sets) {
    for(const auto& element : evo_route)
      fout_evo << element << " ";
    fout_evo <<", ";
  }
  fout_evo << "\n\n";

  }
  
}

/********************/
/*******!MAIN!*******/
/********************/


int main(int argc, char* argv[]) {
  char run_option;
  if(argc<=1) {
    std::cout<<"Too few arguments"<<std::endl;
    return 0;
  }

  run_option=argv[1][1];
  SetRuntimeConfigurations(argc,argv);
  Genotype g(4);
  
  switch(run_option) {
  case 'B':
    EvolvingHomology();
    break;
  case 'E':
    EvolutionRunner();
    break;
  case 'M':
    InteractionMetrics();
    break;
  case '?':
    InterfaceAssembly::PrintBindingStrengths();
    break;
  case 'H':
  default:
    std::cout<<"Polyomino interface model\n**Simulation Parameters**\nN: number of tiles\nP: population size\nK: generation limit\nB: number of phenotype builds\n";
    std::cout<<"\n**Model Parameters**\nU: mutation probability (per interface)\nT: temperature\nI: unbound size factor\nA: misbinding rate\nM: Fitness factor\n";
    std::cout<<"\n**Run options**\nR: evolution without fitness\nE: evolution with fitness\n";
    break;
  }
  return 0;
}

void SetRuntimeConfigurations(int argc, char* argv[]) {
  if(argc<3 && argv[1][1]!='H')
    std::cout<<"Invalid Parameters"<<std::endl;
  else {
    for(uint8_t arg=2;arg<argc;arg+=2) {
      switch(argv[arg][1]) {
        /*! model basics */
      case 'N': simulation_params::n_tiles=std::stoi(argv[arg+1]);break;
      case 'P': simulation_params::population_size=std::stoi(argv[arg+1]);break;
      case 'G': simulation_params::generation_limit=std::stoi(argv[arg+1]);break;
      case 'D': simulation_params::independent_trials=std::stoi(argv[arg+1]);break;
      case 'V': simulation_params::run_offset=std::stoi(argv[arg+1]);break;
        
      case 'B': FitnessPhenotypeTable::phenotype_builds=std::stoi(argv[arg+1]);break;
      case 'X': FitnessPhenotypeTable::UND_threshold=std::stod(argv[arg+1]);break;
      case 'F': FitnessPhenotypeTable::fitness_factor=std::stod(argv[arg+1]);break;

      case 'M': InterfaceAssembly::mutation_rate=std::stod(argv[arg+1]);break;
      case 'J': InterfaceAssembly::duplication_rate=std::stod(argv[arg+1]);break;
      case 'K': InterfaceAssembly::insertion_rate=std::stod(argv[arg+1]);break;
      case 'L': InterfaceAssembly::deletion_rate=std::stod(argv[arg+1]);break;
      case 'Y': InterfaceAssembly::binding_threshold=std::stod(argv[arg+1]);break;
      case 'T': InterfaceAssembly::temperature=std::stod(argv[arg+1]);break;
        
      default: std::cout<<"Unknown Parameter Flag: "<<argv[arg][1]<<std::endl;
      }
    }
    InterfaceAssembly::SetBindingStrengths();

  }
}



