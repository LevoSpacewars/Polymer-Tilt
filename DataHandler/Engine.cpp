


#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include "MainMenuHandler.h"
#include "GraphHandler.h"
#include "Compiler.h"
namespace fs = std::filesystem;
using namespace std;


int main(int argc, char*argv[]) {

    int choice = -1;
    if (argc >= 2)
    {
        for ( int i =1; i < argc; i++)
        {
            Compiler * compiler = new Compiler(string(argv[i]),true);

        }
        
        exit(0);
    }
    MainMenuHandler *menu = new MainMenuHandler();
    while (choice != 0)
    {
        

        choice = menu->restart();

        if(choice == 1);
            //GraphHandler * graph_h = new GraphHandler();
        if(choice == 2);
            Compiler * compiler = new Compiler(menu->getPath(), false);
    
    }

    

    //have the user what they want to graph | type h -N to give desc of the graph option
        //1. compare tilt graphs (select N graphs but cannot repeat)
        //2. make individual tilt graphs
        //3. make individual heatmap graphs
        //4...N other customizations with desc


    //have the user select or group what they want to compile
    
        //1. compile data for selection
        //2. try different comilation 
    
}


