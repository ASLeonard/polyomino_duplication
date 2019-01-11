#include "duplication_simulator.hpp"
#include <iostream>

constexpr bool BINARY_WRITE_FILES=false;
bool KILL_BACK_MUTATIONS=false;
const std::string file_base_path="//rscratch//asl47//Discs//";
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
  


void EvolutionRunner() {
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



uint32_t Evo(uint8_t ttype) {


  FitnessPhenotypeTable pt = FitnessPhenotypeTable();
  DimerModelTable(&pt);
  pt.known_phenotypes[2][0].tiling={1,ttype};

  std::vector<double> population_fitnesses(simulation_params::population_size);
  std::vector<PopulationGenotype> evolving_population(simulation_params::population_size),reproduced_population(simulation_params::population_size);
  
  for(auto& species : evolving_population) {
    species.genotype.resize(2);
    InterfaceAssembly::RandomiseGenotype(species.genotype);
  }
  
  std::set<InteractionPair> pid_interactions;  
  for(uint32_t generation=0;generation<simulation_params::generation_limit;++generation) { /*! MAIN EVOLUTION LOOP */
    uint16_t nth_genotype=0;
    for(PopulationGenotype& evolving_genotype : evolving_population) { /*! GENOTYPE LOOP */
      InterfaceAssembly::Mutation(evolving_genotype.genotype);                       
      BGenotype assembly_genotype=evolving_genotype.genotype;      
      population_fitnesses[nth_genotype]=interface_model::PolyominoAssemblyOutcome(assembly_genotype,&pt,evolving_genotype.pid,pid_interactions);
      ++nth_genotype;
      if(evolving_genotype.pid.first==2)
        return generation;

    } /*! END GENOTYPE LOOP */

    /*! SELECTION */
    uint16_t nth_repro=0;
    for(uint16_t selected : RouletteWheelSelection(population_fitnesses)) {
      reproduced_population[nth_repro++]=evolving_population[selected];
    }
    evolving_population.swap(reproduced_population);
  }
  return simulation_params::generation_limit;
}



void DimerModelTable(FitnessPhenotypeTable* pt) {
  pt->FIXED_TABLE=true;
  pt->known_phenotypes[1].emplace_back(Phenotype{1,1, {1}});
  pt->known_phenotypes[2].emplace_back(Phenotype{2,1, {1,3}});
  pt->phenotype_fitnesses[1].emplace_back(1);
  pt->phenotype_fitnesses[2].emplace_back(2);

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
      case 'R': simulation_params::random_initilisation=std::stoi(argv[arg+1])>0;break;

        /*! simulation specific */
        
        //DONE IN INIT FILE
      case 'M': simulation_params::mu_prob=std::stod(argv[arg+1]);break;
      case 'Y': simulation_params::binding_threshold=std::stod(argv[arg+1]);break;
      case 'T': simulation_params::temperature=std::stod(argv[arg+1]);break;
        

      case 'A': simulation_params::model_type=std::stoi(argv[arg+1]);break;   
      case 'H': simulation_params::dissociation_time=std::stoi(argv[arg+1]);break;
        
      case 'F': simulation_params::fitness_factor=std::stod(argv[arg+1]);break;
      case 'J': simulation_params::fitness_jump=std::stod(argv[arg+1]);break;
        
      case 'O': simulation_params::fitness_period=std::stod(argv[arg+1]);break;
      case 'G': simulation_params::fitness_rise=std::stod(argv[arg+1]);break;
        
      default: std::cout<<"Unknown Parameter Flag: "<<argv[arg][1]<<std::endl;
      }
    }
    InterfaceAssembly::SetBindingStrengths();

  }
}



