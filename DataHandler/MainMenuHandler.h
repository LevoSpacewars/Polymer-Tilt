

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
namespace fs = std::filesystem;
using namespace std;
enum Activity
{
    GRPAH,
    COMPILE 
};
class MainMenuHandler
{
private:

    string c_path = "";

    vector<string>* getSimulationDirectories(string path);
    bool contains(fs::path path, string file_type);
    /* data */
    
public:
    string getPath();

    int restart();
    MainMenuHandler(/* args */);
    ~MainMenuHandler();
};

