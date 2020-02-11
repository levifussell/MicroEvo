#include <SFML/Graphics.hpp>
#include <vector>
#include <iostream>
#include <ctime>

const int SCREEN_SIZE = 300;
const int LATTICE_SIZE = 100;
const float SCALE_RATE = (float)SCREEN_SIZE / (float)LATTICE_SIZE;
std::vector< std::vector<double> > antibiotic;
std::vector< std::vector<double> > vulnerable;
std::vector< std::vector<double> > inhibit;
std::vector< std::vector<double> > kill;

sf::Image drawTexture;

// PARAMS.
const int KILL_RADIUS = 1;
const int INHIBIT_RADIUS = 1;
const int GROW_RADIUS = 3;
const float MARGIN = 0.3;
const float MUTATE_SIZE = 0.1;

const int EPOCHS = 10;

// methods.
int init();
int g_mod(int a, int b);
sf::Color to_color(int i, int j);
void mutate(int i, int j);
bool in_margin(float value, float target);
void mark_kill_neighbours();
void mark_inhibit_neighbours();
void initiate_kill_list();
void spread();
void add_to_image();
void update();

int init()
{
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
    int i = LATTICE_SIZE/2;
    int j = LATTICE_SIZE/2;
    antibiotic.at(i).at(j) = (double)rand()/(double)RAND_MAX;
    vulnerable.at(i).at(j) = (double)rand()/(double)RAND_MAX;
    inhibit.at(i).at(j) = (double)rand()/(double)RAND_MAX;

}

int g_mod(int a, int b)
{
    return ((a % b) + b) % b;
}

sf::Color to_color(int i, int j)
{
    if(antibiotic.at(i).at(j) != -1)
    {
        float a = antibiotic.at(i).at(j);
        float v = vulnerable.at(i).at(j);
        float h = inhibit.at(i).at(j);
        sf::Color c(a*255.0, v*255.0, h*255.0);
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
}

bool in_margin(float value, float target)
{
    float lower = target - MARGIN;
    float upper = target + MARGIN;
    if(lower < 0.0)
        lower += 1.0;
        upper += 1.0;
    if(upper > 1.0)
        upper -= 1.0;
        lower -= 1.0;
    //return (value > lower && value < upper) || (value + 1.0> lower && value + 1.0 < upper) || (value - 1.0 > lower && value - 1.0 < upper);
    return value > lower && value < upper;
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
                for(int k = i - KILL_RADIUS; k < i + KILL_RADIUS; ++k)
                {
                    for(int m = j - KILL_RADIUS; m < j + KILL_RADIUS; ++m)
                    {
                        int k_wrap = g_mod(k, LATTICE_SIZE);
                        int m_wrap = g_mod(m, LATTICE_SIZE);
                        float v = vulnerable.at(k_wrap).at(m_wrap);
                        if(in_margin(a, v))
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
                for(int k = i - INHIBIT_RADIUS; k < i + INHIBIT_RADIUS; ++k)
                {
                    for(int m = j - INHIBIT_RADIUS; m < j + INHIBIT_RADIUS; ++m)
                    {
                        int k_wrap = g_mod(k, LATTICE_SIZE);
                        int m_wrap = g_mod(m, LATTICE_SIZE);
                        if(kill.at(k_wrap).at(m_wrap) != -1)
                        {
                            float a = kill.at(k_wrap).at(m_wrap);
                            if(in_margin(a, h))
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
    for(int i = 0; i < LATTICE_SIZE; ++i)
    {
        for(int j = 0; j < LATTICE_SIZE; ++j)
        {
            if(antibiotic.at(i).at(j) != -1)
            {
                // choose a random point within the range.
                int r_i = rand() % (2 * GROW_RADIUS) + (i - GROW_RADIUS);
                int r_j = rand() % (2 * GROW_RADIUS) + (j - GROW_RADIUS);
                r_i = g_mod(r_i, LATTICE_SIZE);
                r_j = g_mod(r_j, LATTICE_SIZE);
                if(antibiotic.at(r_i).at(r_j) == -1) // if empty...
                {
                    antibiotic.at(r_i).at(r_j) = antibiotic.at(i).at(j);
                    vulnerable.at(r_i).at(r_j) = vulnerable.at(i).at(j);
                    inhibit.at(r_i).at(r_j) = inhibit.at(i).at(j);
                    mutate(r_i, r_j);
                }
            }
        }
    }
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

    //mark_inhibit_neighbours();

    initiate_kill_list();

    spread();
}

int main(int num_args, char** args)
{
    init();

    //sf::RenderWindow window(sf::VideoMode(LATTICE_SIZE, LATTICE_SIZE), "micro");
    sf::RenderWindow window(sf::VideoMode(SCREEN_SIZE, SCREEN_SIZE), "micro");

    clock_t start_time = clock();
    double sec_delay = 0.1;

    while(window.isOpen())
    {
        sf::Event event;
        while(window.pollEvent(event))
        {
            if(event.type == sf::Event::Closed)
                window.close();
        }

        if((clock() - start_time) / (double)CLOCKS_PER_SEC > sec_delay)
        //if(true)
        {
            std::cout << "UPDATE\n";
            start_time = clock();

            update();

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

    return 0;
}
