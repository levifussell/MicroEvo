#include <SFML/Graphics.hpp>
#include <vector>
#include <iostream>
#include <ctime>
#include <time.h>
#include <math.h>
#include <fstream>
#include <sstream>

/*
 * COMPILE: g++ main.cpp -lsfml-audio -lsfml-graphics -lsfml-window -lsfml-system -o run.out
 * RUN: ./run.out data/<save-location>
 *
 * This version will have an interactive element prompting the user to choose the the type it is.
 * 
 */

bool init_done = false;

const int SCREEN_SIZE = 600;
const int LATTICE_SIZE = 200;
const float SCALE_RATE = (float)SCREEN_SIZE / (float)LATTICE_SIZE;
std::vector< std::vector<double> > antibiotic;
std::vector< std::vector<double> > vulnerable;
std::vector< std::vector<double> > inhibit;
std::vector< std::vector<double> > kill;

sf::Image drawTexture;

//ART: 3,0,3,0.2,0.02,

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
int NUM_SPECIES = pow((1.0/SPECIES_BIN), 3);

int RUN_CLASSIFICATION = -1;

int SEED = 13;

int EPOCHS = 5000;
int unchanged_counter = 0;
int NUM_RUNS_COUNTER = 0;
std::string BASE_LOCATION = "";
std::string RUN_FILE_NAME = "TEST";

// recording.
std::vector< std::vector<int> > species_distributions;

// methods.
void init();
void reset();
int g_mod(int a, int b);
sf::Color to_color(int i, int j);
void mutate(int i, int j);
bool in_margin(float value, float target, float margin);
void mark_kill_neighbours();
void mark_inhibit_neighbours();
void initiate_kill_list();
void spread();
void add_to_image();
void update();
void record_species_distributions();
void save_species_distributions();
float run_step(
	int kill_radius, int inhibit_radius, int grow_radius, 
	float kill_margin, float inhibit_margin, float mutate_size,
	int version);
std::vector<float> get_all_run_stats();
void save_all_run_stats();
bool check_classification(int keyCode);
bool terminate;

void init()
{
    std::cout << "INIT\n";
    if(!init_done)
    {
	    //srand(time(NULL));
	    srand(SEED);
	    drawTexture.create(LATTICE_SIZE, LATTICE_SIZE, sf::Color::Red);

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

sf::Color to_color(int i, int j)
{
    if(antibiotic.at(i).at(j) != -1)
    {
	sf::Uint8 a = (sf::Uint8)(antibiotic.at(i).at(j)*255.0);
	sf::Uint8 v = (sf::Uint8)(vulnerable.at(i).at(j)*255.0);
	sf::Uint8 h = (sf::Uint8)(inhibit.at(i).at(j)*255.0);
	//std::cout << (int)(antibiotic.at(i).at(j)*255) << "," << (int)(vulnerable.at(i).at(j)*255) << "," << (int)(inhibit.at(i).at(j)*255) << "\n";
        sf::Color c(a, v, h);
        return c;
    }
    else
    {
        sf::Color c(0.0, 0.0, 0.0);
        return c;
    }
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

void add_to_image()
{
    for(int i = 0; i < LATTICE_SIZE; ++i)
    {
        for(int j = 0; j < LATTICE_SIZE; ++j)
        {
            drawTexture.setPixel(i,j, to_color(i,j));
        }
    }
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
	//std::ostringstream os;
	//os << BASE_LOCATION << "data/" << "k" << KILL_RADIUS << "_i" << INHIBIT_RADIUS << "_d" << GROW_RADIUS << "_m" << MUTATE_SIZE << "_epk" << KILL_MARGIN << "_epi" << INHIBIT_MARGIN << "_v" << VERSION << ".csv";
	//os << BASE_LOCATION << "k" << KILL_RADIUS << "_i" << INHIBIT_RADIUS << "_d" << GROW_RADIUS << "_m" << MUTATE_SIZE << "_epk" << KILL_MARGIN << "_epi" << INHIBIT_MARGIN << "_v" << VERSION << ".csv";
	//std::string name = os.str();
	std::cout << "SAVING AT: " << RUN_FILE_NAME << "\n";
	std::fstream fout;
	fout.open(RUN_FILE_NAME, std::ios::out);
	for(int k = 0; k < species_distributions.size(); ++k)
	{
		for(int j = 0; j < NUM_SPECIES; ++j)
		{
			fout << species_distributions.at(k).at(j);
            if(j < NUM_SPECIES - 1)
                fout << ",";
		}

		fout << "\n";
	}
	fout.close();
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
        
        //assert (p >= 0.0);
        if(p > 0.0)
            shannon_diversity -= p * log2(p);
    }

    shannon_equit_index = shannon_diversity / log2(NUM_SPECIES);
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
void save_shannon_diversity(std::string experimentName) //std::vector<float> shannon_diversity)
{
    std::vector<float> shannon_diversity = compute_shannon_diversity_over_time(species_distributions);


	//std::ostringstream os;
	//os << RUN_FI << "k" << KILL_RADIUS << "_i" << INHIBIT_RADIUS << "_d" << GROW_RADIUS << "_m" << MUTATE_SIZE << "_epk" << KILL_MARGIN << "_epi" << INHIBIT_MARGIN << "_v" << VERSION << ".csv";
	//std::string name = os.str();
	std::cout << "SAVING AT: " << RUN_FILE_NAME << "\n";
	std::fstream fout;
    fout.open(RUN_FILE_NAME, std::ios::out | std::ios::trunc);
    fout << experimentName << "\n";
    //
	for(int k = 0; k < shannon_diversity.size(); ++k)
	{
        fout << shannon_diversity.at(k) << "\n";
	}
	fout.close();
}

std::vector<float> get_all_run_stats()
{
    std::vector<float> run = std::vector<float>{
        (float)KILL_RADIUS, (float)INHIBIT_RADIUS, (float)GROW_RADIUS,
        (float)KILL_MARGIN, (float)INHIBIT_MARGIN, (float)MUTATE_SIZE,
        (float)VERSION
    };

    run.push_back(RUN_CLASSIFICATION);

    return run;
}


bool check_classification(int keyCode)
{
    bool choiceMade = false;

    if(keyCode == sf::Keyboard::K)
    {
        RUN_CLASSIFICATION = 0;
        choiceMade = true;
    }
    else if(keyCode == sf::Keyboard::C)
    {
        RUN_CLASSIFICATION = 1;
        choiceMade = true;
    }
    else if(keyCode == sf::Keyboard::S)
    {
        RUN_CLASSIFICATION = 2;
        choiceMade = true;
    }
    else if(keyCode == sf::Keyboard::E)
    {
        RUN_CLASSIFICATION = 3;
        choiceMade = true;
    }
    else if(keyCode == sf::Keyboard::U)
    {
        RUN_CLASSIFICATION = 4;
        choiceMade = true;
    }
    else if(keyCode == sf::Keyboard::X)
    {
        RUN_CLASSIFICATION = 5;
        choiceMade = true;
    }

    return choiceMade;
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

    bool headless = false;

    std::ostringstream runNameStream;
    runNameStream  << "PARAMS: " << kill_radius << ", " << inhibit_radius << ", " << grow_radius << ", " << mutate_size << ", " << kill_margin << ", " << inhibit_margin << ", " << version;
    std::string runName = runNameStream.str();

    // comment & uncomment if headless (TOOD: make cleaner).
    sf::RenderWindow window(sf::VideoMode(SCREEN_SIZE, SCREEN_SIZE), "micro");
    //sf::RenderWindow window;

    clock_t start_time = clock();
    double sec_delay = 0.1;

    int counter = 0;

    clock_t begin_time = clock();

    start_time = clock();
    double start_delay = 0.0f;
    while(headless || window.isOpen())
    {
        //bool terminate = false;
        //if(!headless)
        //{
            sf::Event event;
            while(window.pollEvent(event))
            {
                if(event.type == sf::Event::Closed)
                {
                    window.close();
                }
            }
        //}

        //if((clock() - start_time) / (double)CLOCKS_PER_SEC > start_delay)
        if(true) // just to remember the placeholder of slowing down the sim.
        {
            //std::cout << "UPDATE\n";

            if((clock() - start_time) / (double)CLOCKS_PER_SEC > start_delay)
            {
                update();
            }

            //std::cout << "STEP: " << counter << "\n";

            if(!headless)
            {
                window.clear(sf::Color::Black);

                // draw here.
                add_to_image();
                sf::Texture tex;
                tex.loadFromImage(drawTexture);
                sf::Sprite sp;
                sp.setTexture(tex, true);
                sp.setScale(SCALE_RATE, SCALE_RATE);
                window.draw(sp);

                window.display();
            }
        }

        //if(counter >= EPOCHS || unchanged_counter >= 5)
        if(terminate)//(check_classification())
        {
            terminate = false;

            if(!headless)
                window.close();
            else
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
    if(num_args != 8)
    {
        std::cerr << "needs 7 arguments to run: <kill radius> <inhibit radius> <grow radius> <mutate size> <kill margin> <inhibit margin> <seed>\n";
        return 1;
    }

    int k = atoi(args[1]);
    int i = atoi(args[2]);
    int d = atoi(args[3]);
    float m = atof(args[4]);
    float ek = atof(args[5]);
    float ei = atof(args[6]);
    int v = atoi(args[7]);


    NUM_RUNS_COUNTER = 0;
    std::cout << "STARTING: " << k << ", " << i << ", " << d << ", " << m << ", " << ek << ", " << ei << ", " << v << "\n";

    SEED = v;
    float time_complete = run_step(k, i, d, ek, ei, m, v);

    return 0;
}
