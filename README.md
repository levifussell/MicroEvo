# Spatial Microbial Community Dynamics Using a Continuous Species Interaction Model #

by Anshuman Swain<sup>\*</sup><sup>1</sup>, Levi Fussell<sup>\*</sup><sup>2</sup>, William F Fagan<sup>1</sup>

<sup>1</sup>Department of Biology, University of Maryland, College Park, MD, USA

<sup>2</sup>Institute of Perception, Action and Behaviour, University of Edinburgh, Edinburgh, UK

<sup>\*</sup>*contributed equally*

### Instructions for Running the Interactive Process. ###

For those interested in running some experiments for themselves, we have outlined the steps below for downloading the code base, compiling it, and running it locally on your machine. These steps focus on MacOS, but should also work for Linux.

(1) Clone the repository to a local path on your machine: `git clone https://github.com/levifussell/MicroEvo`.

(2) You will need to download SFML, which is the low-level graphics library I use for rendering. First, try it the way we describe here. If this doesn't work for you, then follow the instructions [here](https://www.sfml-dev.org/tutorials/2.5/start-osx.php).

(3) Getting 'g++' installed. Try the following different options:
    (a) Type `g++` into a terminal and it might prompt you on how to install it (or it happens to be installed already.
    (b) Install XCode (just google), and then you should have it on your terminal (test by typing `g++`).
    (c) Try one of the options described here: https://stackoverflow.com/questions/2122425/how-do-i-install-g-on-macos-x

(4) Install SFML according to here: https://www.sfml-dev.org/tutorials/2.5/start-osx.php (or navigate to your approrpiate OS).

(5) Compile the code using 'g++':
    (a) type and run in a terminal: `g++ main_param_select.cpp -lsfml-audio -lsfml-graphics -lsfml-window -lsfml-system -o run.out`,
    
or if that fails, `g++ --std=c++11 main_param_select.cpp -lsfml-audio -lsfml-graphics -lsfml-window -lsfml-system -o run.out`, 
    
or if that fails, `g++ --std=c++x0 main_param_select.cpp -lsfml-audio -lsfml-graphics -lsfml-window -lsfml-system -o run.out`,

(6) Run a test:
    (a) type and run in a terminal: `./run.out <kill radius> <inhibit radius> <grow radius> <mutate size> <kill margin> <inhibit margin> <seed>`
    (b) if the parameter values are reasonable, but the simulation is black, try changing the <seed>.

### Abstract ###

Comprehending the assembly and maintenance of microbial diversity in natural communities, despite the abundance of antagonistic interactions, is a major problem of interest in biology. A common framework to study the problem is through cyclic dominance games over pairwise interactions. Recent papers incorporating higher-order interactions in these models have successfully explained high diversity of microbes, especially in communities where antibiotic producing, sensitive, and resistant strains co-exist. But most of these models are based on a small number of discrete species, assume a notion of pure cyclic dominance, and focus on low mutation rate regimes, none of which best represents the highly interlinked, quickly evolving and continuous nature of microbial phenotypic space. Here, we present a model of species in a continuous space, with mutual higher order interactions, to examine the assembly and stability of microbial communities. Specifically, we focus on toxin production, vulnerability, and inhibition among the simulated species. We observe intricate interaction between certain parameters that generates highly divergent patterns of diversity and spatial community dynamics. We find that spatial properties are better predicted by species interaction constraints rather than mobility, and that community formation time, mobility, and mutation rate best explain the patterns of diversity.

