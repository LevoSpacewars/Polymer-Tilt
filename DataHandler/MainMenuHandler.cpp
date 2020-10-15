#include "MainMenuHandler.h"

MainMenuHandler::MainMenuHandler()
{
    cout << "Please insert path to Simulation Directory:";
    cin >> c_path;

    
}

int MainMenuHandler::restart()
{
    string path = c_path;
    
    cout <<" what would you like to do? " << endl;
    cout <<"1. Gragh\n" << "2. Compile Data\n"<<"3. misc" <<endl;
    
    int number;
    cin >> number;
    return number;
    
   
}

string MainMenuHandler::getPath()
{
    return c_path;
}


