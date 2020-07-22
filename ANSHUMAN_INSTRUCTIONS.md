### INSTRUCTIONS FOR ANSHUMAN TO RUN THE INTERACTIVE PROCESS ###

Hi, Anshuman.

You will need to download SFML, which is the low-level graphics library I use for rendering. The instructions for MacOS can be found here: https://www.sfml-dev.org/tutorials/2.5/start-osx.php
But first, try it this way and see if it works. If it doesn't then follow the instructions on the cited page to install XCode, and then try my run command I describe below.

(1) Getting 'g++' installed. Try the following different options:
    (a) Type 'g++' into a terminal and hopefully it prompts you on how to install it (or it happens to be installed already.
    (b) Install XCode (just google), and then you should have it on your terminal (test by typing 'g++' (no parenthesis)).
    (c) Try one of the options described here: https://stackoverflow.com/questions/2122425/how-do-i-install-g-on-macos-x

(2) Install SFML according to here: https://www.sfml-dev.org/tutorials/2.5/start-osx.php

(2) Compile the code using 'g++':
    (a) type and run in a terminal: g++ main_interactive.cpp -lsfml-audio -lsfml-graphics -lsfml-window -lsfml-system -o run.out

(3) Run a test:
    (a) type and run in a terminal: ./run.out <save-location>
    (b) where <save-location> is the directory and name of the file you want to save.

(4) Run your sweep:
    (a) open 'main_interactive.cpp' and go to the bottom of the file. The inputs for the sweep parameters are highlighted.
    (b) then re-run step 4.

