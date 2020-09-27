#include <vector>
#include <iostream>
#include <ctime>
#include <time.h>
#include <math.h>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <stdexcept>

bool init_done = false;

const int SCREEN_SIZE = 600;
const int LATTICE_SIZE = 200;
const float SCALE_RATE = (float)SCREEN_SIZE / (float)LATTICE_SIZE;
std::vector< std::vector<double> > antibiotic;
std::vector< std::vector<double> > vulnerable;
std::vector< std::vector<double> > inhibit;
std::vector< std::vector<double> > kill;

//ART: 3,0,3,0.2,0.02,

std::string RUN_FILE_NAME = "TEST";

// PARAMS.
int KILL_RADIUS = 2;
int INHIBIT_RADIUS = 2;
int GROW_RADIUS = 2;
//float MARGIN = 0.1;
float KILL_MARGIN = 0.2;
float INHIBIT_MARGIN = 0.01;
float MUTATE_SIZE = 0.01;
int VERSION = 0;
//const int KILL_RADIUS = 3;
//const int INHIBIT_RADIUS = 1;
//const int GROW_RADIUS = 3;
//const float MARGIN = 0.003;
//const float MUTATE_SIZE = 0.05;
float SPECIES_BIN = 0.05;
int SPECIES_PER_DIM = (int)(1.0/SPECIES_BIN) + 1;
int NUM_SPECIES = pow(SPECIES_PER_DIM, 3);

int EPOCHS = 1000;
int unchanged_counter = 0;
int NUM_RUNS_COUNTER = 0;
std::string BASE_LOCATION = "";

int SEED = 13;

// recording.
std::vector< std::vector<int> > species_distributions;

// methods.
void init();
void reset();
int g_mod(int a, int b);
void mutate(int i, int j);
bool in_margin(float value, float target, float margin);
void mark_kill_neighbours();
void mark_inhibit_neighbours();
void initiate_kill_list();
void spread();
void update();
void record_species_distributions();
void save_species_distributions();
float run_step(
	int kill_radius, int inhibit_radius, int grow_radius, 
	float kill_margin, float inhibit_margin, float mutate_size,
	int version);

std::vector<int> int_range_inclusive(int min, int max);
std::vector<float> float_range_inclusive(float min, float max, int count);
template<typename T>
void print_vector(std::vector<T> v, std::string name);
void run_param_sweep(int num_args, char** args);
void run_params_from_doc_sweep(int num_args, char** args);

std::vector<float> get_all_run_stats();
void compute_shannon_diversity(std::vector<int> histogram, 
        float& shannon_diversity, float& shannon_equit_index);
float mean(std::vector<float> v);
float variance(std::vector<float> v, float mean);
void compute_mean_var_shannon_diversity_over_time(
        std::vector< std::vector<int> > data,
        float& mean_shannon_diversity, float& var_shannon_diversity,
        float& mean_shannon_equit_index, float& var_shannon_equit_index);
std::vector<float> compute_shannon_diversity_over_time(
        std::vector< std::vector<int> > data);
void save_species_distributions();
void save_all_run_stats();

void init()
{
    std::cout << "INIT\n";
    if(!init_done)
    {
	    //srand(time(NULL));
	    srand(SEED);

	    antibiotic = std::vector< std::vector<double> >(LATTICE_SIZE);
	    vulnerable = std::vector< std::vector<double> >(LATTICE_SIZE);
	    inhibit = std::vector< std::vector<double> >(LATTICE_SIZE);
	    kill = std::vector< std::vector<double> >(LATTICE_SIZE);
	    for(int i = 0; i < LATTICE_SIZE; ++i)
	    {
		antibiotic.at(i) = std::vector<double>(LATTICE_SIZE);
		vulnerable.at(i) = std::vector<double>(LATTICE_SIZE);
		inhibit.at(i) = std::vector<double>(LATTICE_SIZE);
		kill.at(i) = std::vector<double>(LATTICE_SIZE);

	    }

	    // recording data setup.
	    species_distributions = std::vector< std::vector<int> >();
	    init_done = true;
    }

    reset();

}

void reset()
{
    // recording data setup.
    species_distributions.clear();

    // reset everything.
    for(int i = 0; i < LATTICE_SIZE; ++i)
    {
        for(int j = 0; j < LATTICE_SIZE; ++j)
        {
            antibiotic.at(i).at(j) = -1;
            vulnerable.at(i).at(j) = -1;
            inhibit.at(i).at(j) = -1;
            kill.at(i).at(j) = -1;
        }
    }

    // randomise organisms.
    // for(int i = 0; i < LATTICE_SIZE; ++i)
    // {
    //     for(int j = 0; j < LATTICE_SIZE; ++j)
    //     {
    //         antibiotic.at(i).at(j) = (double)rand()/(double)RAND_MAX;
    //         vulnerable.at(i).at(j) = (double)rand()/(double)RAND_MAX;
    //         inhibit.at(i).at(j) = (double)rand()/(double)RAND_MAX;
    //         float a = antibiotic.at(i).at(j);
    //         float v = vulnerable.at(i).at(j);
    //         float h = inhibit.at(i).at(j);
    //         //std::cout << a << "," << v << "," << h << "\n";
    //     }
    // }
    
    // centre organism.
    int i = LATTICE_SIZE/2;
    int j = LATTICE_SIZE/2;
    antibiotic.at(i).at(j) = (double)rand()/(double)RAND_MAX;
    vulnerable.at(i).at(j) = (double)rand()/(double)RAND_MAX;
    inhibit.at(i).at(j) = (double)rand()/(double)RAND_MAX;
    
    //// RPS.
    //int i = LATTICE_SIZE/4;
    //int j = LATTICE_SIZE/4;
    //antibiotic.at(i).at(j) = 0.1;
    //vulnerable.at(i).at(j) = 0.4;
    //inhibit.at(i).at(j) = 0.7;
    //int i2 = LATTICE_SIZE/4*3;
    //int j2 = LATTICE_SIZE/4;
    //antibiotic.at(i2).at(j2) = 0.7;
    //vulnerable.at(i2).at(j2) = 0.1;
    //inhibit.at(i2).at(j2) = 0.4;
    //int i3 = LATTICE_SIZE/2;
    //int j3 = LATTICE_SIZE/4*3;
    //antibiotic.at(i3).at(j3) = 0.4;
    //vulnerable.at(i3).at(j3) = 0.7;
    //inhibit.at(i3).at(j3) = 0.1;
}

int g_mod(int a, int b)
{
    return ((a % b) + b) % b;
}

void mutate(int i, int j)
{
    float a = antibiotic.at(i).at(j) + ((double)rand()/(double)RAND_MAX) * 2.0*MUTATE_SIZE - MUTATE_SIZE;
    float v = vulnerable.at(i).at(j) + ((double)rand()/(double)RAND_MAX) * 2.0*MUTATE_SIZE - MUTATE_SIZE;
    float h = inhibit.at(i).at(j) + ((double)rand()/(double)RAND_MAX) * 2.0*MUTATE_SIZE - MUTATE_SIZE;
//    float a = antibiotic.at(i).at(j);
//    float v = vulnerable.at(i).at(j);
//    float h = inhibit.at(i).at(j);
//    float p = 0.1;
//    if((double)rand()/(double)RAND_MAX < p)
//	    a = (double)rand()/(double)RAND_MAX;
//    if((double)rand()/(double)RAND_MAX < p)
//	    v = (double)rand()/(double)RAND_MAX;
//    if((double)rand()/(double)RAND_MAX < p)
//	    h = (double)rand()/(double)RAND_MAX;
    if(a < 0.0)
        a += 1.0;
    else if(a > 1.0)
        a -= 1.0;
    if(v < 0.0)
        v += 1.0;
    else if(v > 1.0)
        v -= 1.0;
    if(h < 0.0)
        h += 1.0;
    else if(h > 1.0)
        h -= 1.0;
    antibiotic.at(i).at(j) = a;
    vulnerable.at(i).at(j) = v;
    inhibit.at(i).at(j) = h;
    //std::cout << a << ", " << v << ", " << h << "\n";
}

bool in_margin(float value, float target, float margin)
{
    float lower = target - margin;
    float upper = target + margin;
    //std::cout << "LOW " << lower << "\n";
    //std::cout << "UP " << upper << "\n";
    //std::cout << "VALUE " << value << "\n";
    //if(lower < 0.0)
    //    lower += 1.0;
    //    upper += 1.0;
    //if(upper > 1.0)
    //    upper -= 1.0;
    //    lower -= 1.0;
    return (value > lower && value < upper) || (value + 1.0 > lower && value + 1.0 < upper) || (value - 1.0 > lower && value - 1.0 < upper);
    //return value > lower && value < upper;
}

void mark_kill_neighbours()
{
    for(int i = 0; i < LATTICE_SIZE; ++i)
    {
        for(int j = 0; j < LATTICE_SIZE; ++j)
        {
            if(antibiotic.at(i).at(j) != -1)
            {
                float a = antibiotic.at(i).at(j);
                //TODO: square lattice for now.
                for(int k = i - KILL_RADIUS; k < i + KILL_RADIUS + 1; ++k)
                {
                    for(int m = j - KILL_RADIUS; m < j + KILL_RADIUS + 1; ++m)
                    {
                        int k_wrap = g_mod(k, LATTICE_SIZE);
                        int m_wrap = g_mod(m, LATTICE_SIZE);
                        float v = vulnerable.at(k_wrap).at(m_wrap);
                        if((i != k or j != m) and in_margin(a, v, KILL_MARGIN))
                            kill.at(k_wrap).at(m_wrap) = a;
                    }
                }
            }
        }
    }
}

void mark_inhibit_neighbours()
{
    for(int i = 0; i < LATTICE_SIZE; ++i)
    {
        for(int j = 0; j < LATTICE_SIZE; ++j)
        {
            if(antibiotic.at(i).at(j) != -1)
            {
                float h = inhibit[i][j];
                //TODO: square lattice for now.
                for(int k = i - INHIBIT_RADIUS; k < i + INHIBIT_RADIUS + 1; ++k)
                {
                    for(int m = j - INHIBIT_RADIUS; m < j + INHIBIT_RADIUS + 1; ++m)
                    {
                        int k_wrap = g_mod(k, LATTICE_SIZE);
                        int m_wrap = g_mod(m, LATTICE_SIZE);
                        if((i != k or j != m) and kill.at(k_wrap).at(m_wrap) != -1)
                        {
                            float a = kill.at(k_wrap).at(m_wrap);
                            if(in_margin(a, h, INHIBIT_MARGIN))
                                kill.at(k_wrap).at(m_wrap) = -1; //turn-off kill.
                        }
                    }
                }
            }
        }
    }
}

void initiate_kill_list()
{
    for(int i = 0; i < LATTICE_SIZE; ++i)
    {
        for(int j = 0; j < LATTICE_SIZE; ++j)
        {
            if(kill.at(i).at(j) != -1)
            {
                antibiotic.at(i).at(j) = -1;
                vulnerable.at(i).at(j) = -1;
                inhibit.at(i).at(j) = -1;
                kill.at(i).at(j) = -1;
            }
        }
    }
}

void spread()
{
    bool unchanged = true;
    for(int i = 0; i < LATTICE_SIZE; ++i)
    {
        for(int j = 0; j < LATTICE_SIZE; ++j)
        {
            if(antibiotic.at(i).at(j) != -1)
            {
                // choose a random point within the range.
                int r_i = rand() % (2 * GROW_RADIUS + 1) + (i - GROW_RADIUS);
                int r_j = rand() % (2 * GROW_RADIUS + 1) + (j - GROW_RADIUS);
                r_i = g_mod(r_i, LATTICE_SIZE);
                r_j = g_mod(r_j, LATTICE_SIZE);
                if(antibiotic.at(r_i).at(r_j) == -1) // if empty...
                {
                    antibiotic.at(r_i).at(r_j) = antibiotic.at(i).at(j);
                    vulnerable.at(r_i).at(r_j) = vulnerable.at(i).at(j);
                    inhibit.at(r_i).at(r_j) = inhibit.at(i).at(j);
                    mutate(r_i, r_j);
		    unchanged = false;
                }
            }
        }
    }
    if(unchanged)
        unchanged_counter++;
    else
        unchanged_counter = 0;
}

void update()
{
    mark_kill_neighbours();

    mark_inhibit_neighbours();

    initiate_kill_list();

    spread();
}

//void hash_species(int a, int v, int i)
//{
//    //int h = (int)(a/SPECIES_BIN) + (SPECIES_PER_DIM+1)*(int)(v/SPECIES_PER_BIN) + (SPECIES_PER_DIM+(SPECIES_PER_DIM+1)*SPECIES_PER_DIM+1)*(int)(i/SPECIES_PER_BIN);
//    int h = (int)(a/SPECIES_BIN) + (SPECIES_PER_DIM+1)*(int)(v/SPECIES_PER_BIN) + ((SPECIES_PER_DIM+SPECIES_PER_DIM+1)*SPECIES_PER_DIM+1)*(int)(i/SPECIES_PER_BIN);
//    return h
//}

void record_species_distributions()
{
    int unstackedDistr[SPECIES_PER_DIM][SPECIES_PER_DIM][SPECIES_PER_DIM];
    for(int i = 0; i < SPECIES_PER_DIM; ++i) {
        for(int j = 0; j < SPECIES_PER_DIM; ++j) {
            for(int k = 0; k < SPECIES_PER_DIM; ++k) {
                unstackedDistr[i][j][k] = 0;
            }
        }
    }
	for(int k = 0; k < LATTICE_SIZE; ++k)
	{
		for(int j = 0; j < LATTICE_SIZE; ++j)
		{
			if(antibiotic.at(k).at(j) != -1)
			{
				//int s_id = 1 * (int)(antibiotic.at(k).at(j) / SPECIES_BIN) + 10 * (int)(antibiotic.at(k).at(j) / SPECIES_BIN) + 100 * (int)(antibiotic.at(k).at(j) / SPECIES_BIN);
				//int s_id = 1 * (int)(antibiotic.at(k).at(j) / SPECIES_BIN) + (SPECIES_PER_DIM+1) * (int)(antibiotic.at(k).at(j) / SPECIES_BIN) +  (SPECIES_PER_DIM+(SPECIES_PER_DIM+1)*SPECIES_PER_DIM+1)* (int)(antibiotic.at(k).at(j) / SPECIES_BIN);
                //int s_id = hash_species(antibiotic.at(k).at(j), vulnerable.at(k).at(j), inhibit.at(k).at(j));
				//std::cout << s_id << "\n";
				//distr.at(s_id) += 1;
                int x = (int)(antibiotic.at(k).at(j) / SPECIES_BIN);
                int y = (int)(vulnerable.at(k).at(j) / SPECIES_BIN);
                int z = (int)(inhibit.at(k).at(j) / SPECIES_BIN);
                //std::cout << x << ", " << y << ", " << z << "\n";
                //std::cout << SPECIES_PER_DIM << "\n";
                assert (x >= 0 && x < SPECIES_PER_DIM);
                assert (y >= 0 && y < SPECIES_PER_DIM);
                assert (z >= 0 && z < SPECIES_PER_DIM);
                unstackedDistr[x][y][z] += 1;
			}
		}
	}

	std::vector<int> distr;
    for(int i = 0; i < SPECIES_PER_DIM; ++i) {
        for(int j = 0; j < SPECIES_PER_DIM; ++j) {
            for(int k = 0; k < SPECIES_PER_DIM; ++k) {
                distr.push_back(unstackedDistr[i][j][k]);
            }
        }
    }
    assert (distr.size() == sizeof(unstackedDistr) / sizeof(int));

	species_distributions.push_back(distr);
}

std::vector<float> get_all_run_stats()
{
    std::vector<float> run = std::vector<float>{
        (float)KILL_RADIUS, (float)INHIBIT_RADIUS, (float)GROW_RADIUS,
        (float)KILL_MARGIN, (float)INHIBIT_MARGIN, (float)MUTATE_SIZE,
        (float)VERSION, (float)SPECIES_BIN, (float)NUM_SPECIES, (float)EPOCHS,
        (float)species_distributions.size()
    };
    float mean_sd, var_sd;
    float mean_eqi, var_eqi;
    compute_mean_var_shannon_diversity_over_time(
        species_distributions,
        mean_sd, var_sd,
        mean_eqi, var_eqi);
    run.push_back(mean_sd);
    run.push_back(var_sd);
    run.push_back(mean_eqi);
    run.push_back(var_eqi);

    return run;
}

void compute_shannon_diversity(std::vector<int> histogram, 
        float& shannon_diversity, float& shannon_equit_index)
{
    // compute the mass of the histogram.
    int total = 0;
    for(int h : histogram)
        total += h;
    shannon_diversity = 0.0f;
    for(int i = 0; i < histogram.size(); ++i)
    {
        float p = 0.0;
        if(total > 0.0f)
            p = (float)histogram.at(i) / (float)total;
        
        assert (p >= 0.0);
        if(p > 0.0)
            shannon_diversity -= p * log2(p);
    }

    shannon_equit_index = shannon_diversity / log2(NUM_SPECIES);
}

float mean(std::vector<float> v)
{
    float m = 0.0f;
    for(float i : v)
        m += i;
    m /= v.size();
    return m;
}

float variance(std::vector<float> v, float mean)
{
    float var = 0.0f;
    for(int i : v)
        var += pow((i - mean), 2);
    var /= v.size()-1;
    return var;
}

std::vector<float> compute_shannon_diversity_over_time(
        std::vector< std::vector<int> > data)
{
    std::vector<float> shannon_diversities(data.size());
    for(int i = 0; i < data.size(); ++i)
    {
        float sd, eqi;
        compute_shannon_diversity(data.at(i), sd, eqi);
        shannon_diversities.at(i) = sd;
    }
    return shannon_diversities;
}

void compute_mean_var_shannon_diversity_over_time(
        std::vector< std::vector<int> > data,
        float& mean_shannon_diversity, float& var_shannon_diversity,
        float& mean_shannon_equit_index, float& var_shannon_equit_index)
{
    std::vector<float> shannon_diversities(data.size());
    std::vector<float> shannon_equit_indices(data.size());
    for(int i = 0; i < data.size(); ++i)
    {
        float sd, eqi;
        compute_shannon_diversity(data.at(i), sd, eqi);
        shannon_diversities.at(i) = sd;
        shannon_equit_indices.at(i) = eqi;
    }

    mean_shannon_diversity = mean(shannon_diversities);
    var_shannon_diversity = variance(shannon_diversities, mean_shannon_diversity);
    mean_shannon_equit_index = mean(shannon_equit_indices);
    var_shannon_equit_index = variance(shannon_equit_indices, mean_shannon_equit_index);
}

void save_species_distributions()
{
	std::ostringstream os;
	os << BASE_LOCATION << "k" << KILL_RADIUS << "_i" << INHIBIT_RADIUS << "_d" << GROW_RADIUS << "_m" << MUTATE_SIZE << "_epk" << KILL_MARGIN << "_epi" << INHIBIT_MARGIN << "_v" << VERSION << ".csv";
	std::string name = os.str();
	std::cout << "SAVING AT: " << name << "\n";
	std::fstream fout;
	fout.open(name, std::ios::out);
	for(int k = 0; k < species_distributions.size(); ++k)
	{
		for(int j = 0; j < NUM_SPECIES; ++j)
		{
			fout << species_distributions.at(k).at(j) << ", ";
		}

		fout << "\n";
	}
	fout.close();
}

void save_shannon_diversity(std::vector<float> shannon_diversity)
{
	std::ostringstream os;
	os << BASE_LOCATION << "k" << KILL_RADIUS << "_i" << INHIBIT_RADIUS << "_d" << GROW_RADIUS << "_m" << MUTATE_SIZE << "_epk" << KILL_MARGIN << "_epi" << INHIBIT_MARGIN << "_v" << VERSION << ".csv";
	std::string name = os.str();
	std::cout << "SAVING AT: " << name << "\n";
	std::fstream fout;
	fout.open(name, std::ios::out);
	for(int k = 0; k < shannon_diversity.size(); ++k)
	{
        fout << shannon_diversity.at(k) << ", ";
	}
	fout.close();
}

void save_all_run_stats()
{
    std::ostringstream os;

    std::cout << "SAVING AT: " << RUN_FILE_NAME << "\n";
    
    std::fstream fout;
    fout.open(RUN_FILE_NAME, std::ios::out | 
            (NUM_RUNS_COUNTER == 0 ? std::ios::trunc : std::ios::app));
    std::vector<float> run_data = get_all_run_stats();
    for(int i = 0; i < run_data.size(); ++i)
    {
        fout << run_data.at(i) << (i < run_data.size() - 1 ? ", " : "");
    }

    fout << "\n";
}

float run_step(
	int kill_radius, int inhibit_radius, int grow_radius, 
	float kill_margin, float inhibit_margin, float mutate_size,
	int version)
{
    KILL_RADIUS = kill_radius;
    INHIBIT_RADIUS = inhibit_radius;
    GROW_RADIUS = grow_radius;
    KILL_MARGIN = kill_margin;
    INHIBIT_MARGIN = inhibit_margin;
    MUTATE_SIZE = mutate_size;
    VERSION = version;
    unchanged_counter = 0;

    init();

    bool headless = true;


    // comment & uncomment if headless (TOOD: make cleaner).

    double sec_delay = 0.1;

    int counter = 0;

    clock_t begin_time = clock();

    while(headless)
    {
        update();
	    record_species_distributions();

        if(counter >= EPOCHS || unchanged_counter >= 5)
        {
            //save_species_distributions();
            save_all_run_stats();
            //std::vector<float> sd = compute_shannon_diversity_over_time(species_distributions);
            //save_shannon_diversity(sd);
            break;
        }

        counter += 1;
    }

    float time_diff = (clock() - begin_time)/(double)CLOCKS_PER_SEC;
    std::cout << "TIME: " << time_diff << "\n";
    return time_diff;
}

std::vector<int> int_range_inclusive(int min, int max)
{
    std::vector<int> values(max-min+1);
    for(int i = min; i < max+1; ++i)
        values.at(i-min) = i;
    assert (values.at(0) == min);
    assert (values.at(max-min) == max);
    return values;
}

std::vector<float> float_range_inclusive(float min, float max, int count)
{
    std::vector<float> values(count);
    for(float m = min, i = 0; i < count; m += (max - min)/(float)(count-1), ++i)
        values.at(i) = m;
    for(float f : values)
        std::cout << f << ", ";
    assert(abs(values.at(0) - min) < 0.0001);
    assert(abs(values.at(count-1) - max) < 0.0001);
    return values;
}

template<typename T>
void print_vector(std::vector<T> v, std::string name)
{
    std::cout << name << ": ";
    for(T i : v)
        std::cout << i << ",";
    std::cout << "\n";
}

void run_param_sweep(int num_args, char** args)
{
    int i = 2;
    std::vector<int> k_sweep = int_range_inclusive(atoi(args[i]), atoi(args[i + 1]));
    i += 2;
    std::cout << "FIRST\n";
    std::vector<int> i_sweep = int_range_inclusive(atoi(args[i]), atoi(args[i + 1]));
    i += 2;
    std::vector<int> d_sweep = int_range_inclusive(atoi(args[i]), atoi(args[i + 1]));
    i += 2;
    std::vector<float> m_sweep = float_range_inclusive(atof(args[i]), atof(args[i + 1]), atoi(args[i + 2]));
    i += 3;
    std::vector<float> epk_sweep = float_range_inclusive(atof(args[i]), atof(args[i + 1]), atoi(args[i + 2]));
    i += 3;
    std::vector<float> epi_sweep = float_range_inclusive(atof(args[i]), atof(args[i + 1]), atoi(args[i + 2]));
    i += 3;
    //std::vector<int> version = int_range_inclusive(1, atoi(args[i]));
    std::vector<int> version = std::vector<int>{ atoi(args[i]) }; // only one version is supported right now.

    SEED = atoi(args[i + 1]);

    std::ostringstream os;
    os << BASE_LOCATION;
    os << "v2_";
    for(int i = 2; i < num_args; ++i)
    {
        os << args[i] << "_";
    }
    os << ".csv";
    RUN_FILE_NAME = os.str();

    std::cout << "RUNNING PARAMETERS: \n";
    print_vector<int>(k_sweep, "kill radius");
    print_vector<int>(i_sweep, "inhibit radius");
    print_vector<int>(d_sweep, "grow radius");
    print_vector<float>(m_sweep, "mutation");
    print_vector<float>(epk_sweep, "kill margin");
    print_vector<float>(epi_sweep, "inhibit margin");
    print_vector<int>(version, "version (repeats)");

    int total = k_sweep.size() * i_sweep.size() * d_sweep.size() * m_sweep.size() * epk_sweep.size() * epi_sweep.size() * version.size();
    NUM_RUNS_COUNTER = 0;
    float avg_time = 0.0;
    std::cout << "TOTAL RUNS: " << total << "\n";

    for(int v : version) {
        for(int k : k_sweep) {
            for(int i : i_sweep) {
                for(int d : d_sweep) {
                    for(float m : m_sweep) {
                        for(float ek : epk_sweep) {
                            for(float ei : epi_sweep) {
                                std::cout << "START: " << k << ", " << i << ", " << d << ", " << m << ", " << ek << ", " << ei << ", " << v << "\n";
                                float time_complete = run_step(k, i, d, ek, ei, m, v);
                                if(avg_time == 0.0)
                                    avg_time = time_complete;
                                else
                                    avg_time = 0.9*avg_time + 0.1*time_complete;
                                std::cout << "DONE: " << k << ", " << i << ", " << d << ", " << m << ", " << ek << ", " << ei << ", " << v << "\n";
                                NUM_RUNS_COUNTER++;
                                std::cout << "\t" << NUM_RUNS_COUNTER << "/" << total << "\n";
                                std::cout << "PRED. TIME LEFT: " << avg_time*(total-NUM_RUNS_COUNTER)/(60.0*60.0) << "h\n";
                            } } } } } } }
}

void run_params_from_doc(int numargs, char** args)
{
    // first arg is a file location. Next arg is the range of rows to pull from that file.
    if(numargs != 6)
    {
        std::cerr << "ARGS: " << numargs << ". Must have at least 5 arguments to run.\n";
        return;
    }
    std::string file_location = args[2];
    int rowLower = atoi(args[3]);
    int rowUpper = atoi(args[4]);

    std::fstream fin;
    fin.open(file_location, std::ios::in);
    std::string line;
    std::cout << "RUNNING: \n";
    std::vector<int> k_sweep;
    std::vector<int> i_sweep;
    std::vector<int> d_sweep;
    std::vector<float> m_sweep;
    std::vector<float> epk_sweep;
    std::vector<float> epi_sweep;
    int row_count = 0;
    getline(fin, line); // header.
    while(getline(fin, line))
    {
        if(row_count < rowLower || row_count >= rowUpper)
        {
            row_count++;
            continue;
        }

        std::stringstream line_stream(line);
        int col_count = 0;
        while(line_stream.good())
        {
            std::string value;
            getline(line_stream, value, ',');
            if(col_count == 0) {
                k_sweep.push_back(atoi(value.c_str()));
            }
            else if(col_count == 1) {
                i_sweep.push_back(atoi(value.c_str()));
            }
            else if(col_count == 2) {
                d_sweep.push_back(atoi(value.c_str()));
            }
            else if(col_count == 3) {
                m_sweep.push_back(atof(value.c_str()));
            }
            else if(col_count == 4) {
                epk_sweep.push_back(atof(value.c_str()));
            }
            else if(col_count == 5) {
                epi_sweep.push_back(atof(value.c_str()));
            }
            col_count++;
        }
        std::cout << line << "\n";
        row_count++;
    }
    fin.close();

    SEED = atoi(args[5]);

    std::ostringstream os;
    os << BASE_LOCATION;
    os << "v2_";
    for(int i = 3; i < numargs; ++i)
    {
        os << args[i] << "_";
    }
    os << ".csv";
    RUN_FILE_NAME = os.str();

    //int total = k_sweep.size() * i_sweep.size() * d_sweep.size() * m_sweep.size() * epk_sweep.size() * epi_sweep.size();
    int total = k_sweep.size();
    NUM_RUNS_COUNTER = 0;
    float avg_time = 0.0;
    std::cout << "TOTAL RUNS: " << total << "\n";

    // run the selected params.
    for(int e = 0; e < k_sweep.size(); ++e)
    {
        int k = k_sweep.at(e);
        int i = i_sweep.at(e);
        int d = d_sweep.at(e);
        float m = m_sweep.at(e);
        float ek = epk_sweep.at(e);
        float ei = epi_sweep.at(e);
        std::cout << "START: " << k << ", " << i << ", " << d << ", " << m << ", " << ek << ", " << ei << "\n";
        float time_complete = run_step(k, i, d, ek, ei, m, 0);
        //float time_complete = 0.0;
        if(avg_time == 0.0)
            avg_time = time_complete;
        else
            avg_time = 0.9*avg_time + 0.1*time_complete;
        std::cout << "DONE: " << k << ", " << i << ", " << d << ", " << m << ", " << ek << ", " << ei << "\n";
        NUM_RUNS_COUNTER++;
        std::cout << "\t" << NUM_RUNS_COUNTER << "/" << total << "\n";
        std::cout << "PRED. TIME LEFT: " << avg_time*(total-NUM_RUNS_COUNTER)/(60.0*60.0) << "h\n";
    }
}

int main(int num_args, char** args)
{
    std::cout << "NUM ARGS: " << num_args << "\n";
    for(int i = 0; i < num_args; ++i)
        std::cout << i << ": " << args[i] << "\n";

    BASE_LOCATION = args[1];
    run_params_from_doc(num_args, args);

    //if(num_args == 2)
    //{
    //    BASE_LOCATION = args[1];
    //    std::cout << "!- RUNNING NO PARAMETER SWEEP.\n";
    //    run_step(KILL_RADIUS, INHIBIT_RADIUS, GROW_RADIUS, KILL_MARGIN, INHIBIT_MARGIN, MUTATE_SIZE, VERSION);
    //}
    //else if(num_args == 19)
    //{
    //    BASE_LOCATION = args[1];
    //    std::cout << "!- RUNNING PARAMETER SWEEP.\n";
    //    run_param_sweep(num_args, args);
    //}
    //else
    //{
    //    throw std::invalid_argument("arguments: <filename> (opt: <k_l> <k_h> <i_l> <i_h> <g_l> <g_h> <m_l> <m_h> <m_c> <epk_l> <epk_h> <epk_c> <epi_l> <epi_h> <epi_c> <v> <seed>)");
    //}

    return 0;
}
