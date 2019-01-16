#include "duplication_simulator.hpp"
#include <iostream>

constexpr bool BINARY_WRITE_FILES=false;
bool KILL_BACK_MUTATIONS=false;
const std::string file_base_path="//scratch//asl47//Data_Runs//Bulk_Data//";
const std::map<Phenotype_ID,uint8_t> phen_stages{{{0,0},0},{{10,0},4},{{1,0},1},{{2,0},2},{{4,0},2},{{4,1},3},{{8,0},3},{{12,0},4},{{16,0},4}};

namespace simulation_params {
  uint16_t population_size=100;
  double fitness_factor=1;
}



void EvRu() {
  const uint32_t N_runs=simulation_params::independent_trials;
  std::ofstream f_out(file_base_path+"Decor"+std::to_string(simulation_params::binding_threshold)+".BIN", std::ios::binary);
  std::vector<uint8_t> res(N_runs);
  for(uint8_t gap=0; gap<= simulation_params::samming_threshold/2;gap++) {
#pragma omp parallel for schedule(dynamic) 
    for(uint32_t r=0;r < N_runs;++r) {
      do {
        res[r]= DEvo3(gap);
      }while(res[r]==2);
    }    
    BinaryWriter(f_out,res);
  }
}
  


void EvolutionRunnerz() {
  const uint32_t N_runs=simulation_params::independent_trials;
  std::ofstream f_out(file_base_path+"Discovs"+std::to_string(simulation_params::binding_threshold)+".BIN");

  
#pragma omp parallel for schedule(dynamic) 
  for(uint32_t r=0;r < N_runs;++r) {
    uint32_t gen=0;
    if(simulation_params::samming_threshold%2==0)
      gen = Evo1();
#pragma omp critical(file_write)
    f_out<<+gen<<" ";
    
  }
  
  f_out<<"\n";
#pragma omp parallel for schedule(dynamic) 
  for(uint32_t r=0;r < N_runs;++r) {
    uint32_t gen = Evo2();
#pragma omp critical(file_write)
    f_out<<+gen<<" ";
  }
  f_out<<"\n";
}

void EvolutionRunner2() {
  const uint32_t N_runs=simulation_params::independent_trials;
  std::ofstream f_out(file_base_path+"Decays"+std::to_string(simulation_params::binding_threshold)+".BIN");
  for(uint8_t gap=0;gap<=simulation_params::samming_threshold;++gap) {
#pragma omp parallel for schedule(dynamic) 
    for(uint32_t r=0;r < N_runs;++r) {
      uint32_t gen=0;
      if(gap%2==0 && simulation_params::samming_threshold%2==0)
        gen = DEvo1(gap/2);
#pragma omp critical(file_write)
      f_out<<+gen<<" ";
    
    }
  
    f_out<<"\n";
#pragma omp parallel for schedule(dynamic) 
    for(uint32_t r=0;r < N_runs;++r) {
      uint32_t gen = DEvo2(gap);
#pragma omp critical(file_write)
      f_out<<+gen<<" ";
    }
    f_out<<"\n";
  }
}

uint32_t Evo1() {
  interface_type geno=InterfaceAssembly::GenRandomSite();
  std::uniform_int_distribution<uint8_t> dis(0, interface_size-1);
  for(uint32_t generation=1;generation<=simulation_params::generation_limit;++generation) {
    geno.flip(dis(RNG_Engine));
    if(interface_model::SammingDistance(geno,geno)<=simulation_params::samming_threshold)
      return generation;
  }
  return 0;
}
uint32_t Evo2() {
  interface_type  geno1,geno2;
  do {
    geno1=InterfaceAssembly::GenRandomSite();
    geno2=InterfaceAssembly::GenRandomSite();
  }while(interface_model::SammingDistance(geno1,geno2)<=simulation_params::samming_threshold);

  std::uniform_int_distribution<uint8_t> dis(0, interface_size-1);
  for(uint32_t generation=1;generation<=simulation_params::generation_limit;++generation) {
    geno1.flip(dis(RNG_Engine));
    if(interface_model::SammingDistance(geno1,geno2)<=simulation_params::samming_threshold)
      return generation;
  }
  return 0;
}

uint32_t DEvo1(uint8_t gap) {
  interface_type geno=InterfaceAssembly::GenRandomSite();
  for(size_t n=0;n<interface_size/2;++n)
    geno[interface_size-n]=~geno[n];

  std::uniform_int_distribution<uint8_t> dis(0, interface_size-1);
  std::vector<uint8_t> bits(interface_size/2);
  std::iota(bits.begin(),bits.end(),0);
  std::shuffle(bits.begin(),bits.end(),RNG_Engine);

  for(uint8_t b=0; b<gap;++b)
    geno.flip(bits[b]);

  
  for(uint32_t generation=1;generation<=simulation_params::generation_limit;++generation) {
     geno.flip(dis(RNG_Engine));
    if(interface_model::SammingDistance(geno,geno)>simulation_params::samming_threshold)
      return generation;
  }
  return 0;
}

uint32_t DEvo2(uint8_t gap) {

  interface_type geno1=InterfaceAssembly::GenRandomSite();
  interface_type geno2=interface_model::ReverseBits(~geno1);

  std::uniform_int_distribution<uint8_t> dis(0, interface_size-1);
  std::vector<uint8_t> bits(interface_size);
  std::iota(bits.begin(),bits.end(),0);
  std::shuffle(bits.begin(),bits.end(),RNG_Engine);

  for(uint8_t b=0; b<gap;++b)
    geno1.flip(bits[b]);

  
  for(uint32_t generation=1;generation<=simulation_params::generation_limit;++generation) {
    geno2.flip(dis(RNG_Engine));
    if(interface_model::SammingDistance(geno1,geno2)>simulation_params::samming_threshold)
      return generation;
  }
  return 0;
}

uint8_t DEvo3(uint8_t gap) {

  interface_type geno1=InterfaceAssembly::GenRandomSite();
  for(size_t n=0;n<interface_size/2;++n)
    geno1[interface_size-n]=~geno1[n];
 
  
  std::uniform_int_distribution<uint8_t> dis(0, interface_size-1);
  std::vector<uint8_t> bits(interface_size);
  std::iota(bits.begin(),bits.end(),0);
  std::shuffle(bits.begin(),bits.end(),RNG_Engine);
  std::bernoulli_distribution bern(0.5);
  for(uint8_t b=0; b<gap;++b)
    geno1.flip(bits[b]);
  interface_type geno2=geno1;
  

  for(uint32_t generation=1;generation<=100000;++generation) {
    if(InterfaceAssembly::GetActiveInterfaces(std::vector<interface_type>{geno1,geno2}).size()==1) {
      return interface_model::SammingDistance(geno1,geno2)<=simulation_params::samming_threshold;
    }

    if(bern(RNG_Engine))
      geno1.flip(dis(RNG_Engine));
    else
      geno2.flip(dis(RNG_Engine));


  }
  return 2;

}



void EvolutionRunner() {
  /*!PYTHON INFORMATION*/
  const std::string py_exec = "python3 ";
  const std::string py_loc = "~/Documents/PolyDev/polyomino_interfaces/scripts/interface_analysis.py ";
  const std::string py_mode="internal "+std::to_string(simulation_params::model_type);
  
  const std::string py_CALL=py_exec + py_loc + py_mode + " "+std::to_string(BINARY_WRITE_FILES)+" ";
  const std::string python_params=" "+std::to_string(simulation_params::binding_threshold)+" "+std::to_string(simulation_params::temperature)+" "+std::to_string(simulation_params::mutation_rate)+" "+std::to_string(simulation_params::fitness_factor)+" "+std::to_string(simulation_params::population_size);

  const uint16_t N_runs=simulation_params::independent_trials;
#pragma omp parallel for schedule(dynamic) 
  for(uint16_t r=0;r < N_runs;++r) {
    EvolvePopulation("_Run"+std::to_string(r+simulation_params::run_offset));
    /*!PYTHON CALL*/
    //std::system((py_CALL+std::to_string(r)+python_params).c_str());
    /*!PYTHON CALL*/

    
  }
}

void EvolvePopulation(std::string run_details) {
  std::string file_simulation_details=run_details+".txt";
    
  std::ofstream fout_strength(file_base_path+"Strengths"+file_simulation_details,std::ios::out);
  std::ofstream fout_phenotype(file_base_path+"PhenotypeTable"+file_simulation_details,std::ios::out);  
  std::ofstream fout_selection_history(file_base_path+"Selections"+file_simulation_details,std::ios::out);    
  std::ofstream fout_phenotype_IDs(file_base_path+"PIDs"+file_simulation_details,std::ios::out );
  
  
  std::vector<double> population_fitnesses(simulation_params::population_size);
  std::vector<PopulationGenotype> evolving_population(simulation_params::population_size),reproduced_population;
  reproduced_population.resize(simulation_params::population_size);

  
  FitnessPhenotypeTable pt = FitnessPhenotypeTable();
  pt.fit_func=[](double s) {return s*s;};


  std::set<InteractionPair> pid_interactions;
  Genotype assembly_genotype;
  Phenotype_ID prev_ev;

  
  
  for(uint32_t generation=0;generation<simulation_params::generation_limit;++generation) { /*! MAIN EVOLUTION LOOP */


    uint16_t nth_genotype=0;
    for(PopulationGenotype& evolving_genotype : evolving_population) { /*! GENOTYPE LOOP */
      
      InterfaceAssembly::Mutation(evolving_genotype.genotype);
      
            
      //const std::vector<std::pair<InteractionPair,double> > edges = InterfaceAssembly::GetActiveInterfaces(evolving_genotype.genotype);

      assembly_genotype=evolving_genotype.genotype;
      prev_ev=evolving_genotype.pid;


      population_fitnesses[nth_genotype]=interface_model::PolyominoAssemblyOutcome(assembly_genotype,&pt,evolving_genotype.pid,pid_interactions);
      ++nth_genotype;



      for(auto x : pid_interactions)
        fout_strength<<+x.first<<" "<<+x.second<<" "<<+interface_model::SammingDistance(assembly_genotype[x.first],assembly_genotype[x.second])<<".";
      fout_strength<<",";
      fout_phenotype_IDs << +evolving_genotype.pid.first <<" "<<+evolving_genotype.pid.second<<" ";
      
    } /*! END GENOTYPE LOOP */

    /*! SELECTION */
    uint16_t nth_repro=0;
    for(uint16_t selected : RouletteWheelSelection(population_fitnesses)) {
      reproduced_population[nth_repro++]=evolving_population[selected];
      fout_selection_history<<+selected<<" ";
    }
    evolving_population.swap(reproduced_population);



    fout_selection_history<<"\n";
    fout_phenotype_IDs<<"\n";
    fout_strength<<"\n";
    
    
  } /* END EVOLUTION LOOP */
  pt.PrintTable(fout_phenotype);
  
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
    EvRu();
    break;
  case 'E':
    EvolutionRunner();
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
      case 'K': simulation_params::generation_limit=std::stoi(argv[arg+1]);break;
      case 'B': FitnessPhenotypeTable::phenotype_builds=std::stoi(argv[arg+1]);break;
      case 'X': FitnessPhenotypeTable::UND_threshold=std::stod(argv[arg+1]);break;

        /*! run configurations */
      case 'D': simulation_params::independent_trials=std::stoi(argv[arg+1]);break;
      case 'V': simulation_params::run_offset=std::stoi(argv[arg+1]);break;

        /*! simulation specific */
        
        //DONE IN INIT FILE
      case 'M': simulation_params::mutation_rate=std::stod(argv[arg+1]);break;
      case 'Y': simulation_params::binding_threshold=std::stod(argv[arg+1]);break;
      case 'T': simulation_params::temperature=std::stod(argv[arg+1]);break;
        

      case 'A': simulation_params::model_type=std::stoi(argv[arg+1]);break;   
        
      case 'F': simulation_params::fitness_factor=std::stod(argv[arg+1]);break;
      case 'J': simulation_params::fitness_jump=std::stod(argv[arg+1]);break;

      default: std::cout<<"Unknown Parameter Flag: "<<argv[arg][1]<<std::endl;
      }
    }
    InterfaceAssembly::SetBindingStrengths();

  }
}



