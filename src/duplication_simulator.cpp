#include "duplication_simulator.hpp"
#include <iostream>

constexpr bool BINARY_WRITE_FILES=false;
bool KILL_BACK_MUTATIONS=false;
const std::string file_base_path="//scratch//asl47//Data_Runs//Bulk_Data//";
const std::string shared_base_path="//rscratch//asl47//Duplication//EvoRecords//";
const std::string shared_base_path2="//rscratch//asl47//Duplication//Metrics//";
const std::map<Phenotype_ID,uint8_t> phen_stages{{{0,0},0},{{10,0},4},{{1,0},1},{{2,0},2},{{4,0},2},{{4,1},3},{{8,0},3},{{12,0},4},{{16,0},4}};


void InteractionMetrics() {
  const uint32_t N_runs=simulation_params::independent_trials;
  std::ofstream f_out(shared_base_path2+"Discovery_"+std::to_string(InterfaceAssembly::binding_threshold)+".BIN", std::ios::binary);
  std::vector<uint32_t> res_S(N_runs),res_A(N_runs);
#pragma omp parallel for schedule(dynamic) 
  for(uint32_t r=0;r < N_runs;++r) {
    res_S[r]= DiscoverInteraction(false,false);
    res_A[r]= DiscoverInteraction(false,true);
  }    
  BinaryWriter(f_out,res_S);
  BinaryWriter(f_out,res_A); 
  return;

  std::ofstream f_out2(file_base_path+"Decay_"+std::to_string(InterfaceAssembly::binding_threshold)+".BIN", std::ios::binary);
  for(uint8_t gap=0;gap<=InterfaceAssembly::samming_threshold;++gap) {
#pragma omp parallel for schedule(dynamic) 
    for(uint32_t r=0;r < N_runs;++r) {
      if(gap%2==0 && InterfaceAssembly::samming_threshold%2==0)
        res_S[r] = DecayInteraction(true,gap/2);
      res_A[r]= DecayInteraction(false,gap);
    }
    BinaryWriter(f_out2,res_S);
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
      RandomiseGenotype(g);
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
  RandomiseGenotype(genotype);
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
void EvolveHomology(std::string run_details,bool self) {
  std::string file_simulation_details=run_details+".BIN";
  std::ofstream fout_homology(file_base_path+"Bomology"+file_simulation_details,std::ios::binary); 

  std::vector<double> population_fitnesses(simulation_params::population_size);
  std::vector<PopulationGenotype> evolving_population(simulation_params::population_size),reproduced_population;
  reproduced_population.resize(simulation_params::population_size);

  Genotype g(8);
  RandomiseGenotype(g);
  if(!self) {
    for(size_t n=0;n<interface_size/2;++n)
      g[0][interface_size-n-1]=~g[0][n];
    std::move(g.begin(),g.begin()+4,g.begin()+4);
  }
  else
    g[4]=interface_model::ReverseBits(~g[0]);

  for(auto& ep : evolving_population)
    ep.genotype=Genotype(g);

  std::uniform_int_distribution<uint8_t> dis(0, interface_size-1);
  std::uniform_int_distribution<uint8_t> dis2(0, 7);

  std::vector<uint8_t> res;
  res.reserve(simulation_params::generation_limit*4*simulation_params::population_size);
  
  for(uint32_t generation=0;generation<simulation_params::generation_limit;++generation) { /*! MAIN EVOLUTION LOOP */
    
    uint16_t nth_genotype=0;
    for(PopulationGenotype& evolving_genotype : evolving_population) { /*! GENOTYPE LOOP */
      auto hom = CalculateHomology(evolving_genotype.genotype);
      res.insert(res.end(),hom.begin(),hom.end());
           
      evolving_genotype.genotype[dis2(RNG_Engine)].flip(dis(RNG_Engine));
      population_fitnesses[nth_genotype]=0;
      const std::vector<std::pair<InteractionPair,double> > edges = InterfaceAssembly::GetActiveInterfaces(evolving_genotype.genotype);
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

  Phenotype::DETERMINISM_LEVEL=3;
  KILL_BACK_MUTATIONS=true;
  const std::string py_exec = "python3 ";
  const std::string py_loc = "~/Documents/PolyDev/duplication/scripts/interface_analysis.py ";
  const std::string py_mode="internal "+std::to_string(simulation_params::model_type);
  
  const std::string py_CALL=py_exec + py_loc + py_mode + " "+std::to_string(BINARY_WRITE_FILES)+" ";
  const std::string python_params=" "+std::to_string(InterfaceAssembly::binding_threshold)+" "+std::to_string(InterfaceAssembly::temperature)+" "+std::to_string(InterfaceAssembly::mutation_rate)+" "+std::to_string(InterfaceAssembly::duplication_rate)+" "+std::to_string(InterfaceAssembly::insertion_rate)+" "+std::to_string(InterfaceAssembly::deletion_rate)+" "+std::to_string(simulation_params::population_size);

  const uint16_t N_runs=simulation_params::independent_trials;
#pragma omp parallel for schedule(dynamic) 
  for(uint16_t r=0;r < N_runs;++r) {
    //EvolveHomology("_Run"+std::to_string(r+simulation_params::run_offset),r%2);
    EvolvePopulation("_Run"+std::to_string(r+simulation_params::run_offset));
    /*!PYTHON CALL*/
    //std::system((py_CALL+std::to_string(r)+python_params).c_str());
    /*!PYTHON CALL*/

    
  }
}

void UpdatePhylogenyTrackers(PopulationGenotype& PG, std::map<Phenotype_ID, std::tuple<uint32_t, uint16_t, PhenotypeEdgeInformation > >& Homology_tracker,uint32_t generation, uint16_t pop_index) {
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
        const auto hom = (PG.subunits[edge.first.first] ^ PG.subunits[edge.first.second]).count();
        const auto str = interface_model::SammingDistance(PG.subunits[edge.first.first], PG.subunits[edge.first.second]);
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
    
  //Record homology for "successful" observations
  if(PG.PID_tracker.size()==1 && std::accumulate(PG.PID_tracker.begin()->second.begin(),PG.PID_tracker.begin()->second.end(),0)==PG.PID_depth) {
    const auto pid= PG.PID_tracker.begin()->first;
    if(PG.PID_lineage.find(pid)==PG.PID_lineage.end()) {

    }
    if(Homology_tracker.find(pid)==Homology_tracker.end())
      //return;
      Homology_tracker[pid] = PG.PID_details[pid];
  }
  
}



void EvolvePopulation(std::string run_details) {
  std::string file_simulation_details=run_details+".txt";
    
  std::ofstream fout_strength(file_base_path+"Strengths"+file_simulation_details,std::ios::binary);
  std::ofstream fout_zomology(file_base_path+"Zomology"+file_simulation_details,std::ios::binary);
  std::ofstream fout_evo(shared_base_path+"EvoRecord_Mu"+std::to_string(InterfaceAssembly::mutation_rate)+"_S"+std::to_string(InterfaceAssembly::binding_threshold)+(InterfaceAssembly::duplication_rate > InterfaceAssembly::insertion_rate ? "_D" + std::to_string(InterfaceAssembly::duplication_rate) : "_I" + std::to_string(InterfaceAssembly::insertion_rate))+".txt",std::ios::app );
  
  std::string fname_phenotype(file_base_path+"PhenotypeTable"+file_simulation_details);
  
  std::ofstream fout_selection_history, fout_phenotype_IDs, fout_size, fout_interactions2;

  const bool FULL_WRITE=false;
  if(FULL_WRITE) {
    fout_selection_history.open(file_base_path+"Selections"+file_simulation_details,std::ios::binary);
    fout_phenotype_IDs.open(file_base_path+"PIDs"+file_simulation_details,std::ios::out );
    fout_size.open(file_base_path+"Size"+file_simulation_details,std::ios::out);
    fout_interactions2.open(file_base_path+"Binteractions"+file_simulation_details,std::ios::out);
    
  }

  
  

  std::map<Phenotype_ID, std::tuple< uint32_t,uint16_t, PhenotypeEdgeInformation > > evolution_record;
  //std::set
    
  std::vector<double> population_fitnesses(simulation_params::population_size);
  std::vector<PopulationGenotype> evolving_population(simulation_params::population_size),reproduced_population;
  reproduced_population.resize(simulation_params::population_size);

  std::vector<uint16_t> binary_homologies(interface_size+1), binary_strengths(interface_size+1);
  
  Phenotype::DETERMINISM_LEVEL=1;
  FitnessPhenotypeTable pt = FitnessPhenotypeTable();
  pt.fit_func=[](double s) {return s;};

   
  if(simulation_params::model_type==17) {
    pt.known_phenotypes[1].emplace_back(Phenotype(1,1,{1}));
    pt.known_phenotypes[2].emplace_back(Phenotype(2,1,{1,3}));
    pt.known_phenotypes[2].emplace_back(Phenotype(2,1,{1,5}));
    pt.phenotype_fitnesses[1].emplace_back(1);
    pt.phenotype_fitnesses[2].emplace_back(4);
    pt.phenotype_fitnesses[2].emplace_back(4);
    pt.FIXED_TABLE=true;
  }
  
  
  //Genotype assembly_genotype;
  //std::vector<uint8_t> pIDs(simulation_params::population_size*2);
  //Phenotype_ID prev_ev;
 
  
  for(uint32_t generation=0;generation<simulation_params::generation_limit;++generation) { /*! MAIN EVOLUTION LOOP */


    uint16_t nth_genotype=0;
    for(PopulationGenotype& evolving_genotype : evolving_population) { /*! GENOTYPE LOOP */
      
      if(evolving_genotype.genotype.empty()) {
        std::cout<<"NULL GENOTYPE AT START OF GENERATION\ngen "<<generation<<" nth "<<nth_genotype<<"\n";
        return;
      }
      InterfaceAssembly::Mutation(evolving_genotype.genotype,1,1,1);
            
      evolving_genotype.subunits=evolving_genotype.genotype;

      evolving_genotype.pid_interactions.clear();
      auto pid_map=interface_model::PolyominoAssemblyOutcome(evolving_genotype.subunits,&pt,evolving_genotype.pid_interactions);
      
      const std::vector<std::pair<InteractionPair,double> > edges = InterfaceAssembly::GetActiveInterfaces(evolving_genotype.subunits);
      
      
      //evolving_genotype.subunits=assembly_genotype;
      switch(pid_map.begin()->first.first) {
      case 255:
        population_fitnesses[nth_genotype]=0;
        break;
      case 1:
        population_fitnesses[nth_genotype]=1;
        break;
      default:
        if(evolving_genotype.genotype.size()>(8*pid_map.rbegin()->first.first))
          population_fitnesses[nth_genotype]=0; //oversized genotype
        else
          population_fitnesses[nth_genotype]=pt.GenotypeFitness(pid_map);
      }

      //Add pids to genotype
      //evolving_genotype.pids.clear();
      //for(auto& kv : pid_map)
      //  evolving_genotype.pids.emplace_back(kv.first);
      
      

      UpdatePhylogenyTrackers(evolving_genotype,evolution_record,generation,nth_genotype);
      
      for(const auto& pid_kv : evolving_genotype.pid_interactions) {
        for(const auto& ints_kv : pid_kv.second) {
          auto IP = ints_kv.first;
          binary_homologies[(evolving_genotype.subunits[IP.first]^evolving_genotype.subunits[IP.second]).count()]+=ints_kv.second;
          binary_strengths[interface_size-interface_model::SammingDistance(evolving_genotype.subunits[IP.first],evolving_genotype.subunits[IP.second])]+=ints_kv.second;
        }
      }
           
      ++nth_genotype;
      
      if(FULL_WRITE) {
        for(auto& kv : pid_map) 
          fout_phenotype_IDs<<+kv.first.first<<" "<<+kv.first.second<<" ";
        for(auto edge : edges)
          fout_interactions2<<+edge.first.first<<" "<<+edge.first.second<<" ";
        fout_interactions2<<",";
      
        fout_phenotype_IDs<<",";
        fout_size<<evolving_genotype.subunits.size()/4<<" ";
      }
      
    } /*! END GENOTYPE LOOP */


    /*! SELECTION */
    uint16_t nth_repro=0;
    auto reproducing_selection=RouletteWheelSelection(population_fitnesses);
    for(uint16_t selected : reproducing_selection) {
      reproduced_population[nth_repro++]=evolving_population[selected];
    }
    evolving_population.swap(reproduced_population);
    
    
    //BinaryWriter(fout_zomology,binary_homologies);
    //BinaryWriter(fout_strength,binary_strengths);

    if(FULL_WRITE) {
      BinaryWriter(fout_selection_history,reproducing_selection);
      fout_phenotype_IDs<<"\n";
      fout_size<<"\n";
      fout_interactions2<<"\n";
    }



    
    
  } /* END EVOLUTION LOOP */

#pragma omp critical(evo_record_write)
  {
  for(auto kv : evolution_record) {
    fout_evo<<+kv.first.first<<" "<<+kv.first.second<<" "<<std::get<0>(kv.second)<<" "<<std::get<1>(kv.second)<<" ";
      for(auto pv : std::get<2>(kv.second))
        fout_evo<<+std::get<0>(pv).first<<" "<<+std::get<0>(pv).second<<" "<<+std::get<1>(pv)<<" "<<+std::get<2>(pv)<<" "<<std::get<3>(pv)<<" ";
      fout_evo<<"\n";
    }
  }
  
  //pt.PrintTable(fname_phenotype);
  
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
      case 'A': simulation_params::model_type=std::stoi(argv[arg+1]);break; 
      case 'H': simulation_params::homologous_threshold=std::stod(argv[arg+1]);break;
        
      case 'B': FitnessPhenotypeTable::phenotype_builds=std::stoi(argv[arg+1]);break;
      case 'X': FitnessPhenotypeTable::UND_threshold=std::stod(argv[arg+1]);break;
      case 'F': FitnessPhenotypeTable::fitness_factor=std::stod(argv[arg+1]);break;

        /*! run configurations */
      

        /*! simulation specific */
        
        //DONE IN INIT FILE
      case 'M': InterfaceAssembly::mutation_rate=std::stod(argv[arg+1]);break;
      case 'J': InterfaceAssembly::duplication_rate=std::stod(argv[arg+1]);break;
      case 'K': InterfaceAssembly::insertion_rate=std::stod(argv[arg+1]);break;
      case 'L': InterfaceAssembly::deletion_rate=std::stod(argv[arg+1]);break;
      case 'Y': InterfaceAssembly::binding_threshold=std::stod(argv[arg+1]);break;
      case 'T': InterfaceAssembly::temperature=std::stod(argv[arg+1]);break;
        
        //case 'Q': FitnessPhenotypeTable::fitness_jump=std::stod(argv[arg+1]);break;

      default: std::cout<<"Unknown Parameter Flag: "<<argv[arg][1]<<std::endl;
      }
    }
    InterfaceAssembly::SetBindingStrengths();

  }
}



