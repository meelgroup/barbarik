/*****************************************************************************************[Main.cc]
Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#include <errno.h>

#include <signal.h>
#include <zlib.h>

#include "utils/System.h"
#include "utils/ParseUtils.h"
#include "utils/Options.h"
#include "core/Dimacs.h"
#include "core/Solver.h"

#include <vector>
#include <map>
#include <string>
#include <sstream>
#include "mtl/Sort.h"
#include <algorithm>

//#include <boost/random/linear_congruential.hpp>
//#include <boost/random/uniform_int.hpp>
//#include <boost/random/uniform_real.hpp>
//#include <boost/random/variate_generator.hpp>
//#include <boost/random.hpp>

//using namespace boost;
using namespace Minisat;

//=================================================================================================


void printStats(Solver& solver)
{
    double cpu_time = cpuTime();
    //double mem_used = memUsedPeak();
	    double mem_used = 0;
    printf("restarts              : %"PRIu64"\n", solver.starts);
    printf("conflicts             : %-12"PRIu64"   (%.0f /sec)\n", solver.conflicts   , solver.conflicts   /cpu_time);
    printf("decisions             : %-12"PRIu64"   (%4.2f %% random) (%.0f /sec)\n", solver.decisions, (float)solver.rnd_decisions*100 / (float)solver.decisions, solver.decisions   /cpu_time);
    printf("propagations          : %-12"PRIu64"   (%.0f /sec)\n", solver.propagations, solver.propagations/cpu_time);
    printf("conflict literals     : %-12"PRIu64"   (%4.2f %% deleted)\n", solver.tot_literals, (solver.max_literals - solver.tot_literals)*100 / (double)solver.max_literals);
    if (mem_used != 0) printf("Memory used           : %.2f MB\n", mem_used);
    printf("CPU time              : %g s\n", cpu_time);
}


static Solver* solver;
// Terminate by notifying the solver and back out gracefully. This is mainly to have a test-case
// for this feature of the Solver as it may take longer than an immediate call to '_exit()'.
static void SIGINT_interrupt(int signum) { solver->interrupt(); }

// Note that '_exit()' rather than 'exit()' has to be used. The reason is that 'exit()' calls
// destructors and may cause deadlocks if a malloc/free function happens to be running (these
// functions are guarded by locks for multithreaded use).
static void SIGINT_exit(int signum) {
    printf("\n"); printf("*** INTERRUPTED ***\n");
    if (solver->verbosity > 0){
        printStats(*solver);
        printf("\n"); printf("*** INTERRUPTED ***\n"); }
    _exit(1); }

	/*
struct bar : std::unary_function<unsigned, unsigned> {
    boost::mt19937 &_state;
    unsigned operator()(unsigned i) {
        boost::uniform_int<> rng(0, i - 1);
        return rng(_state);
    }
    bar(boost::mt19937 &state) : _state(state) {}
};
	
void shuffle(std::vector<unsigned> &vec, boost::mt19937 &state)
{
    bar rand(state);
    std::random_shuffle(vec.begin(), vec.end(), rand);
}
*/
double logsumexp(std::vector <double> nums) {
  double max_exp = nums[0], sum = 0.0;
  size_t ct = nums.size();
  size_t i;

  for (i = 1 ; i < ct ; i++)
    if (nums[i] > max_exp)
      max_exp = nums[i];

  for (i = 0; i < ct ; i++)
    sum += exp(nums[i] - max_exp);

  return log(sum) + max_exp;
}


double compute_solution_weight(std::vector <bool> sol, double * weights, std::vector <unsigned > var_order, int l)
{
double log_prod = 0.0;
for (int i = 0; i < l; i++)
	{
	unsigned v = var_order[i];
	if (weights[v]!=-1)			// it is a weighted variable
		if (sol[v]==true)
			log_prod += log(weights[v]);
		else
			log_prod += log(1.0-weights[v]);
	}
return log_prod;
}

double compute_normalization_weight(std::vector < std::vector <bool> > sol , double * weights, std::vector <unsigned > var_order, int l)
{
std::vector <double> val;
val.resize(sol.size());
for (int i = 0; i < (int) sol.size(); i++)
	val[i]=compute_solution_weight(sol[i],weights,var_order,l);
return logsumexp(val);
}

//=================================================================================================
// Main:


int main(int argc, char** argv)
{
    try {
        setUsageHelp("USAGE: %s [options] <input-file> <result-output-file>\n\n  where input may be either in plain or gzipped DIMACS.\n");
        // printf("This is MiniSat 2.0 beta\n");
        
#if defined(__linux__)
        fpu_control_t oldcw, newcw;
        _FPU_GETCW(oldcw); newcw = (oldcw & ~_FPU_EXTENDED) | _FPU_DOUBLE; _FPU_SETCW(newcw);
        printf("WARNING: for repeatability, setting FPU to use double precision\n");
#endif
        // Extra options:
        //
        IntOption    verb   ("MAIN", "verb",   "Verbosity level (0=silent, 1=some, 2=more).", 1, IntRange(0, 2));
        IntOption    cpu_lim("MAIN", "cpu-lim","Limit on CPU time allowed in seconds.\n", INT32_MAX, IntRange(0, INT32_MAX));
        IntOption    mem_lim("MAIN", "mem-lim","Limit on memory usage in megabytes.\n", INT32_MAX, IntRange(0, INT32_MAX));
        
		
				// Sampling options:
		
		BoolOption   random_variable_ordering    ("SAMPLING", "rvar",    "Random variable ordering or lex", false);
		IntOption    nsamples ("SAMPLING", "nsamples","Number of sampling iterations.\n", 10, IntRange(0, 300000000));
		IntOption    ks ("SAMPLING", "k","Number of samples per level.\n", 50, IntRange(0, 100000000));
		
        parseOptions(argc, argv, true);

        Solver S;
        double initial_time = cpuTime();

        S.verbosity = verb;
        
        solver = &S;
        // Use signal handlers that forcibly quit until the solver will be able to respond to
        // interrupts:
        signal(SIGINT, SIGINT_exit);
        signal(SIGXCPU,SIGINT_exit);

        // Set limit on CPU-time:
        if (cpu_lim != INT32_MAX){
            rlimit rl;
            getrlimit(RLIMIT_CPU, &rl);
            if (rl.rlim_max == RLIM_INFINITY || (rlim_t)cpu_lim < rl.rlim_max){
                rl.rlim_cur = cpu_lim;
                if (setrlimit(RLIMIT_CPU, &rl) == -1)
                    printf("WARNING! Could not set resource limit: CPU-time.\n");
            } }

        // Set limit on virtual memory:
        if (mem_lim != INT32_MAX){
            rlim_t new_mem_lim = (rlim_t)mem_lim * 1024*1024;
            rlimit rl;
            getrlimit(RLIMIT_AS, &rl);
            if (rl.rlim_max == RLIM_INFINITY || new_mem_lim < rl.rlim_max){
                rl.rlim_cur = new_mem_lim;
                if (setrlimit(RLIMIT_AS, &rl) == -1)
                    printf("WARNING! Could not set resource limit: Virtual memory.\n");
            } }
        
        if (argc == 1)
            printf("Reading from standard input... Use '--help' for help.\n");
        
        gzFile in = (argc == 1) ? gzdopen(0, "rb") : gzopen(argv[1], "rb");
        if (in == NULL)
            printf("ERROR! Could not open file: %s\n", argc == 1 ? "<stdin>" : argv[1]), exit(1);
        
        if (S.verbosity > 0){
            printf("============================[ Problem Statistics ]=============================\n");
            printf("|                                                                             |\n"); }
        
		double weighted_problem = false;
		
		double * VarWeights = parse_weighted_DIMACS(in,S);
        if (VarWeights!=0)
			weighted_problem = true;
			
		//parse_DIMACS(in, S);
        gzclose(in);
        FILE* res = (argc >= 3) ? fopen(argv[2], "wb") : NULL;
        
        if (S.verbosity > 0){
            printf("|  Number of variables:  %12d                                         |\n", S.nVars());
            printf("|  Number of clauses:    %12d                                         |\n", S.nClauses()); }
        
        double parsed_time = cpuTime();
        if (S.verbosity > 0){
            printf("|  Parse time:           %12.2f s                                       |\n", parsed_time - initial_time);
            printf("|                                                                             |\n"); }
 
        // Change to signal-handlers that will only notify the solver and allow it to terminate
        // voluntarily:
        signal(SIGINT, SIGINT_interrupt);
        signal(SIGXCPU,SIGINT_interrupt);
       
	   /*
        if (!S.simplify()){
            if (res != NULL) fprintf(res, "UNSAT\n"), fclose(res);
            if (S.verbosity > 0){
                printf("===============================================================================\n");
                printf("Solved by unit propagation\n");
                printStats(S);
                printf("\n"); }
            printf("UNSATISFIABLE\n");
            exit(20);
        }        
		*/
		
        vec<Lit> dummy;
        lbool ret = S.solveLimited(dummy);

        if (S.verbosity > 0){
            printStats(S);
            printf("\n"); }
        printf(ret == l_True ? "SATISFIABLE\n" : ret == l_False ? "UNSATISFIABLE\n" : "INDETERMINATE\n");
        if (res != NULL){
            if (ret == l_True){
                fprintf(res, "SAT\n");
                for (int i = 0; i < S.nVars(); i++)
                    if (S.model[i] != l_Undef)
                        fprintf(res, "%s%s%d", (i==0)?"":" ", (S.model[i]==l_True)?"":"-", i+1);
                fprintf(res, " 0\n");
								
            }else if (ret == l_False)
                fprintf(res, "UNSAT\n");
            else
                fprintf(res, "INDET\n");
            fclose(res);
        }
		
		int var_num = S.nVars();
		std::vector <std::vector <std::vector <bool> > > AllSamples;
		std::vector <std::vector <bool> > InputSamples;
		std::vector <std::vector <bool> > OutputSamples;
		std::vector <double > multipliers;
		std::vector <double > avg_multipliers;
		avg_multipliers.resize(var_num+1,0.0);
				
		std::vector <double > sampled_log_z;
		
		std::vector <unsigned > var_ordering;		// variable ordering
		// generate a random variable ordering
		var_ordering.resize(var_num);
		for (int hh =0;hh<var_num;hh++)
			var_ordering[hh]=hh;
		
		int actual_number_of_samples = 0;
		//bool random_variable_ordering = false;
		
		//boost::mt19937 rng;         // produces randomness out of thin air                                 
		srand(S.random_seed);
		
		int k = ks;
		//std::random_shuffle(var_ordering.begin(),var_ordering.end());
		if (weighted_problem && verb>1)
		{
		for (int hh =0;hh<var_num;hh++)
			printf("%f,",VarWeights[hh]);
		}	
		//shuffle(var_ordering, rng);			//!!!!!!!!!!!!
		//var_ordering[0]=100;
		//var_ordering[100]=0;
		
		
		//int k = 200;

		S.verbosity = 0;
		std::map <std::string, int > counts;
		//int nsamples = 10;
		
		for (int ss =1;ss<=nsamples;ss++)
		{
		multipliers.clear();
		//AllSamples.clear();
		
		if (random_variable_ordering)
			{
			std::random_shuffle(var_ordering.begin(),var_ordering.end());
			//shuffle(var_ordering, rng);
			/*
			for (int hh =0;hh<var_num-1;hh++)
				{
				int l = (rand() % var_num-hh);
				l = l +hh;
				unsigned b = var_ordering[hh];
				var_ordering[hh]=var_ordering[l];
				var_ordering[l]=b;
				}
				*/
			}
		if (ret == l_True){						// sampling code
			if (verb>1)
				printf("Proceeding with sampling\n");			
			InputSamples.resize(1);
			for (int level =0;level<S.nVars()+1;level++)
				{				
				OutputSamples.clear();
				int prev_samp_num = (int) InputSamples.size();
				
				if (weighted_problem)
				{
					// get array of probabilities
					double WNormalization = compute_normalization_weight(InputSamples, VarWeights, var_ordering, level);
					printf("%f\n",WNormalization);
				}
				
				if (verb>0)
					printf("----------------- Level %d (%d)-----------------------------\n",level,prev_samp_num);
				for (int j=0;j<std::min(prev_samp_num,k);j++)
				{
					// sample s_j from Samples							
					vec<Lit>    cur_model;
					if (level>0)
						{
						//int l = irand(seed,InputSamples.size());
						int l = rand() % InputSamples.size();
						
						//boost::uniform_int<> roll(0,InputSamples.size()-1);
						//int l = roll(rng);  
						
						cur_model.growTo(level);
						for (int i = 0; i < level; i++) {
								
								if (i==level-1)							// change last								
								{
								cur_model[i] = mkLit(var_ordering[i],InputSamples[l][var_ordering[i]]);
								if (verb>1)
									printf("x_%d = %d\n",var_ordering[i]+1,!toInt(InputSamples[l][var_ordering[i]]));
								}
								else
								{
								cur_model[i] = mkLit(var_ordering[i],!InputSamples[l][var_ordering[i]]);
								if (verb>1)
									printf("x_%d = %d\n",var_ordering[i]+1,toInt(InputSamples[l][var_ordering[i]]));
								}								
								}
						if (verb>1)
							printf("\n");
						OutputSamples.push_back(InputSamples[l]);
						InputSamples.erase(InputSamples.begin()+l);		// sampling without replacement
						}
					
					// check if there is another one
					lbool ret2 = S.solveLimited(cur_model);
					
					if (ret2 == l_True){
					
					std::vector <bool> solution;
					solution.resize(S.nVars());
					for (int i = 0; i < S.nVars(); i++)
						if (S.model[i] != l_Undef)
							{
							solution[i]=(S.model[i]==l_True);
							//printf("%d,",toInt(solution[i]));
							}
					OutputSamples.push_back(solution);

					}				
				}
				if (weighted_problem)
					{
					double AvgDelta = 0.0;
					if (level==0)
						AvgDelta = 1.0;
					else
					{
						for (int i = 0; i < (int) OutputSamples.size(); i++)
							{
							unsigned v = var_ordering[level-1];
							if (VarWeights[v]!=-1)			// it is a weighted variable
								{
								if (OutputSamples[i][v]==true)
									AvgDelta += VarWeights[v];
								else
									AvgDelta += 1.0-VarWeights[v];
								}
							else
								AvgDelta += 1.0;
							}
					}
					//printf("qverqage delta %e:\n",AvgDelta);
					multipliers.push_back( AvgDelta / std::min(prev_samp_num,k) );
					}
				else
					multipliers.push_back( (double) (OutputSamples.size()+0.0) / std::min(prev_samp_num,k) );
				//AllSamples.push_back(InputSamples);
				InputSamples = OutputSamples;
				}
			
			actual_number_of_samples = actual_number_of_samples +OutputSamples.size();
			// output the samples
			if (verb>0)
			{
			printf("Outputting samples:\n");
			for (int l = 0; l < (int) OutputSamples.size(); l++)
				{
				//int z = 0;
				std::string bit_representation;
				for (int i = 0; i < var_num; i++)
					{
						std::stringstream sss;//create a stringstream
						sss << toInt(OutputSamples[l][i]);//add number to the stream
						bit_representation = bit_representation +sss.str();
						if (i!=var_num-1)
							{
							if (verb>0)						
							printf("%d,",toInt(OutputSamples[l][i]));							
							}
						else
							{
							if (verb>0)						
							printf("%d\n",toInt(OutputSamples[l][i]));							
							}
					}						
				counts[bit_representation] = counts[bit_representation]+1;
				}
			}	
			// output the estimated multipliers
			if (verb>=2)
				{
				printf("Outputting mutlipliers:\n");
				for (int l = 0; l < (int) multipliers.size(); l++)
					printf("%f,",multipliers[l]);
				printf("\n");	
				}
			/*	
			double logZ = 0.0;			// log-partition function
				for (int l = 0; l < (int) multipliers.size(); l++)
					logZ = logZ + log(multipliers[l]);
			if (verb>0)
				printf("Log partition function: %f\n",logZ);
			sampled_log_z.push_back(logZ);
			*/
			
			double avg_log_z_hat = 0.0;
			double logZ = 0.0;			// log-partition function
				for (int l = 0; l < (int) multipliers.size(); l++)
					{
					if (!random_variable_ordering)		// we can average the multipliers						
						{
							avg_multipliers[l] = (avg_multipliers[l] * (ss-1)+multipliers[l])/ss;
							avg_log_z_hat = avg_log_z_hat + log(avg_multipliers[l]);
						}
					logZ = logZ + log(multipliers[l]);					
					}
			if (verb>0)
				printf("Log partition function: %f\n",logZ);
			sampled_log_z.push_back(logZ);
			
			double ss_time = cpuTime();
			double log_z_hat =0.0;
			if (random_variable_ordering)
				log_z_hat = logsumexp(sampled_log_z)-log(sampled_log_z.size());		
			else
				log_z_hat = avg_log_z_hat;
			printf("z <%e, %e > , time= %12.2f, samp=%d\n",exp(log_z_hat),exp(log_z_hat),ss_time-initial_time,ss);
			fflush (stdout);
        }
		}
		
		
		
		double chi = 0.0;
		double e = (1.0/counts.size())*actual_number_of_samples;
		std::map<std::string,int>::iterator it;
		for ( it=counts.begin() ; it != counts.end(); it++ )
			{
			double o = (*it).second+0.0;
			printf("%d\n",(*it).second);
			chi = chi + pow(o-e,2)/e;
			}
		printf("Different : %d\n",counts.size());
		printf("Chi-square : %f\n",chi);
		
		double log_z_hat = logsumexp(sampled_log_z)-log(nsamples);
		printf("Estimated log-z: %f\n",log_z_hat);
		printf("Estimated Z: %e\n",exp(log_z_hat));
		
		
		
#ifdef NDEBUG
        exit(ret == l_True ? 10 : ret == l_False ? 20 : 0);     // (faster than "return", which will invoke the destructor for 'Solver')
#else
        return (ret == l_True ? 10 : ret == l_False ? 20 : 0);
#endif
    } catch (OutOfMemoryException&){
        printf("===============================================================================\n");
        printf("INDETERMINATE\n");
        exit(0);
    }
}
