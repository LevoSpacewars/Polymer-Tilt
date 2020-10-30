extern "C"{
    #include "gsd.h"
}
#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <vector>
#include <cmath>
#include <filesystem>
namespace fs = std::filesystem;
using namespace std;


struct PolymerProfile
{
    float   sheerForceRange [2] = {0,0};
    int   df                  = 0;
    int     length              = 0;
    int     lines               = 0;
    float   tension             = 0;
    float   k_amplitude         = 0;
    float   p_amplitude         = 0;
    float   gamma               = 0;
    float   kbT                 = 0;
    float   dt                  = 0;
    int     sampleRate          = 0;
    int     runLength           = 0;
    float   boxdimx             = 0;
    float   boxdimy             = 0;
    string  read_direciton      = "";

};
struct HeatMapParameters {
    int x = 0;
    int y = 0;
    float width = 0;
    float height = 0;
    int rezx = 0;
    int rezy = 0;
};
class Compiler
{
private:

    float pi = 0.3;
    float pf = 0.6;
    struct PolymerProfile profile;
    struct gsd_handle handler;
    vector<float> output;
    vector<float> uoutput;
    vector<float> dx;
    vector<float> udx;
    vector<float> dy;
    vector<float> udy;
    vector<float> length;
    vector<float> ulength;

    vector<float> sheerTension;
    

    int compileData(string* fname,float interval);
    bool truncateFiles();
    int definePolymerProfile(string* parameter_file_Location,struct PolymerProfile * p);
    void unwrapData(float ** data, int p_n, int p_length,int step);
    float* calcAveragePosition(float ** data, int p_n, int p_length,int runLength);
    float* calcAverageDx(float **x, int p_n, int p_length);
    bool exportDensityFunction_avg(float** xa, float ** ya, int p_n, int p_length,float force_value, string path);
    bool exportDensityFunction_raw(float** x, float ** y, int p_n, int p_length,float force_value, string path);
    float* calcAverageLength(float ** x, float ** y, int p_n, int p_length);
    bool  writeProfileOutput(float** x, float ** y, int p_n, int p_length,float force_value, string path);
    bool  writeHeatMap(float** x, float ** y, int time_steps, int nParticles,float force_value,bool normalized, HeatMapParameters param, string path);
    float * calcSystemOutput(float ** dx, float ** length, float sheer_tension);
    float * calcSheerTension();
    int writeResults(string path);

    float interval = 0.5;
    int runLength;
    string current_path;

    vector<string>* getSimulationDirectories(string path);
    bool contains(fs::path path, string file_type);


public:

    Compiler(string path, bool singular);

    ~Compiler();

};