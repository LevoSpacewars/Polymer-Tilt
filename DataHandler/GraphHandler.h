// extern "C"{
//     #include "gsd.h"
// }
// #include <iostream>
// #include <vector>
// #include <string>
// #include <filesystem>
// namespace fs = std::filesystem;
// using namespace std;
// struct PolymerProfile
// {
//     float   sheerForceRange [2] = {0,0};
//     int   df                  = 0;
//     int     length              = 0;
//     int     lines               = 0;
//     float   tension             = 0;
//     float   k_amplitude         = 0;
//     float   p_amplitude         = 0;
//     float   gamma               = 0;
//     float   kbT                 = 0;
//     float   dt                  = 0;
//     int     sampleRate          = 0;
//     int     runLength           = 0;
//     float   boxdimx             = 0;
//     float   boxdimy             = 0;
//     string  read_direciton      = "";

// };
// struct Bounds {
//     int x = 0;
//     int y = 0;
//     float width = 0;
//     float height = 0;
// };

// class GraphHandler
// {
// private:
//     int mem_size =0;
//     PolymerProfile profile;
//     struct gsd_handle handler;
//     int runLength;
//     float interval =0;
//     float * datax;
//     float * datay;

//     int definePolymerProfile(string* parameter_file_Location,struct PolymerProfile * p);
//     vector<string> getSimulationDirectories(string path);
//     bool contains(fs::path path, string file_type);
//     int createDensityPlots(Bounds bound, int xrez, int yrez);
//     bool getRawData(string* name, float interval);
//     void unwrapData(float ** data, int p_n, int p_length,int step);

// public:
//     GraphHandler(string path);
//     ~GraphHandler();



// };


// //
// //TODO: create desnity plot? 
// //run different simulations?
// //try and get a 