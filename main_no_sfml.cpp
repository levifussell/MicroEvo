#include <vector>
#include <iostream>
#include <ctime>
#include <time.h>
#include <math.h>
#include <fstream>
#include <sstream>

bool init_done = false;

const int SCREEN_SIZE = 600;
const int LATTICE_SIZE = 200;
const float SCALE_RATE = (float)SCREEN_SIZE / (float)LATTICE_SIZE;
std::vector< std::vector<double> > antibiotic;
std::vector< std::vector<double> > vulnerable;
std::vector< std::vector<double> > inhibit;
std::vector< std::vector<double> > kill;

//ART: 3,0,3,0.2,0.02,

// PARAMS.
int KILL_RADIUS = 2;
int INHIBIT_RADIUS = 1;
int GROW_RADIUS = 2;
//float MARGIN = 0.1;
float KILL_MARGIN = 0.1;
float INHIBIT_MARGIN = 0.01;
float MUTATE_SIZE = 0.05;
int VERSION = 0;
//const int KILL_RADIUS = 3;
//const int INHIBIT_RADIUS = 1;
//const int GROW_RADIUS = 3;
//const float MARGIN = 0.003;
//const float MUTATE_SIZE = 0.05;
float SPECIES_BIN = 0.05;
int NUM_SPECIES = pow((1.0/SPECIES_BIN), 3);

int EPOCHS = 500;
int unchanged_counter = 0;
std::string BASE_LOCATION = "";

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

void init()
{
    std::cout << "INIT\n";
    if(!init_done)
    {
	    //srand(time(NULL));
	    srand(13);

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

void record_species_distributions()
{
	std::vector<int> distr(NUM_SPECIES, 0);
	for(int k = 0; k < LATTICE_SIZE; ++k)
	{
		for(int j = 0; j < LATTICE_SIZE; ++j)
		{
			if(antibiotic.at(k).at(j) != -1)
			{
				int s_id = 1 * (int)(antibiotic.at(k).at(j) / SPECIES_BIN) + 10 * (int)(antibiotic.at(k).at(j) / SPECIES_BIN) + 100 * (int)(antibiotic.at(k).at(j) / SPECIES_BIN);
				//std::cout << s_id << "\n";
				distr.at(s_id) += 1;
			}
		}
	}

	species_distributions.push_back(distr);
}

void save_species_distributions()
{
	std::ostringstream os;
	os << BASE_LOCATION << "data/" << "k" << KILL_RADIUS << "_i" << INHIBIT_RADIUS << "_d" << GROW_RADIUS << "_m" << MUTATE_SIZE << "_epk" << KILL_MARGIN << "_epi" << INHIBIT_MARGIN << "_v" << VERSION << ".csv";
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

    clock_t start_time = clock();
    double sec_delay = 0.1;

    int counter = 0;

    clock_t begin_time = clock();

    while(headless)
    {
        //if((clock() - start_time) / (double)CLOCKS_PER_SEC > sec_delay)
        if(true)
        {
            //std::cout << "UPDATE\n";
            start_time = clock();

            update();
	    record_species_distributions();
        }

	if(counter >= EPOCHS || unchanged_counter >= 5)
	{
	    save_species_distributions();
	    break;
	}

	counter += 1;
    }

    float time_diff = (clock() - begin_time)/(double)CLOCKS_PER_SEC;
    std::cout << "TIME: " << time_diff << "\n";
    return time_diff;
}

int main(int num_args, char** args)
{

    BASE_LOCATION = args[1];

    bool param_sweep = true;
    if(num_args != 3)
        throw std::invalid_argument("TWO INPUTS REQUIRED.");
    //TODO: for now, 'k' is used to parallelise on clusters, but easy to make work with all params for later.
    std::vector<int> k_sweep{ std::stoi(args[2]) }; 
    //std::vector<int> k_sweep{1,2,3,5,10,20};
    std::vector<int> i_sweep{0,1,2,5,10,20};
    std::vector<int> d_sweep{1,2,3,5,10,20};
    std::vector<float> m_sweep{0.01,0.04,0.08,0.2};
    std::vector<float> epk_sweep{0.05,0.1,0.2,0.4};
    std::vector<float> epi_sweep{0.01,0.05,0.1,0.2};
    std::vector<int> version{0,1,2};

    int total = k_sweep.size() * i_sweep.size() * d_sweep.size() * m_sweep.size() * epk_sweep.size() * epi_sweep.size() * version.size();
    int count = 0;
    float avg_time = 0.0;

    if(param_sweep)
    {
	for(int v : version)
	{
		for(int k : k_sweep)
		{
			for(int i : i_sweep)
			{
				for(int d : d_sweep)
				{
					for(float m : m_sweep)
					{
						for(float ek : epk_sweep)
						{
							for(float ei : epi_sweep)
							{
								std::cout << "START: " << k << ", " << i << ", " << d << ", " << m << ", " << ek << ", " << ei << ", " << v << "\n";
								float time_complete = run_step(k, i, d, ek, ei, m, v);
								if(avg_time == 0.0)
									avg_time = time_complete;
								else
									avg_time = 0.9*avg_time + 0.1*time_complete;
								std::cout << "DONE: " << k << ", " << i << ", " << d << ", " << m << ", " << ek << ", " << ei << ", " << v << "\n";
								count += 1;
								std::cout << "\t" << count << "/" << total << "\n";
								std::cout << "PRED. TIME LEFT: " << avg_time*(total-count)/(60.0*60.0) << "h\n";
							}
						}
					}
				}
			}
		}
	}
    }
    else
    {
	std::cout << "RUNNING NO PARAMETER SWEEP.\n";
    	run_step(KILL_RADIUS, INHIBIT_RADIUS, GROW_RADIUS, KILL_MARGIN, INHIBIT_MARGIN, MUTATE_SIZE, VERSION);
    }

    return 0;
}
