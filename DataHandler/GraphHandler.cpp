// #include "GraphHandler.h"



// GraphHandler::GraphHandler(string path)
// {
//     vector<string>* paths = getSimulationDirectories(path);

//     for (int i = 0; i < paths->size(); i++)
//     {
//         cout << i<< ". " + paths->at(i)<<endl;
//     }
//     cout <<paths->size()<< ". ALL"<<endl;

//     cout << "choose which file to compile data for: ";
//     int index = -1;
//     cin >> index;
//     cin >> index;
    
        
//     string polyProfileName = paths->at(index) + "/_simulation_parameters.txt";
//     string polyDataName = paths->at(index) + "/trajectory.gsd";


//     this->definePolymerProfile(&polyProfileName, &profile);
//     this->getRawData(&polyDataName, interval);
// }


// bool GraphHandler::getRawData(string * filename, float interval)
// {
//      cout<<"compiling data"<<endl;
//     runLength = (int) (this->profile.runLength / this->profile.sampleRate);
//     int l_polymer = this->profile.length;
//     int n_polymers = this->profile.lines; 

//     int df = (int)(this->profile.df);
//     int t_step     = 0;
//     int line     = 0;
//     int n_runs   = (int)df;
//     int c_run    = 0;
//     int lstep    = 0;
//     int t_adj   =0;
//     float prev    = 0;
//     float prevt   = 0;
//     float L       = this->profile.boxdimx;
//     float currentt = 0;
//     float * raw_x = new float[l_polymer*n_polymers * ((int)((1-interval)*runLength)) * df];
//     float * raw_y = new float[l_polymer*n_polymers * ((int)((1-interval)*runLength)) * df];

//     int gsd_open_error = gsd_open(&this->handler,filename->c_str(), GSD_OPEN_READONLY);

//     auto e = gsd_find_chunk(&this->handler,0,"particles/position");
//     mem_size = l_polymer*n_polymers * ((int)((1-interval)*runLength));


//     for (int i = 0; i< df;i++)
//     {
//         //Begin by allocating the arrays needed for position data;
//         int current_run = i;

//         float * pos_x;
//         float * pos_y;


//         // pos_x = (float*)malloc(l_polymer*n_polymers * ((int)(1-interval))*runLength * sizeof(float));
//         // pos_y = (float*)malloc(l_polymer*n_polymers * ((int)(1-interval))*runLength * sizeof(float));
//         int mem_size = l_polymer*n_polymers * ((int)((1-interval)*runLength));
//         pos_x = new float[mem_size];
//         pos_y = new float[mem_size];

        
//         cout<<i<<endl;
//         //Next Sort through the data and "unwrap" the polymers from periodic boundary
//         // Begin by 
//         for(int j = int(runLength*interval); j < runLength; j++)
//         {
//             t_step = i*runLength + j;
//             t_adj = j - (int)(runLength*interval);

//             auto chunk_entry = gsd_find_chunk(&this->handler,t_step,"particles/position"); //retrives the chunk information from a time step from the gsd file
//             //cout<<"assigning raw data"<<endl;
//             //cout<< e->N * 3 * sizeof(float) <<endl;
//             float * raw_data = new float[e->N * 3]; // (number of particles) by (dimensions)
//             // cout<<"failed?"<<endl;
//             int errorch = gsd_read_chunk(&this->handler,raw_data, chunk_entry); // retrives data from chunk
            

//             int indext = 0;
//             int base_offset = t_adj * (l_polymer*n_polymers);
//             for (int v = e->N*((int)((interval)*e->M)); v < e->N*e->M; v=v+3) // 0,1,2 | ,3,4,5 |,6,7,8 ... // need to implement a proper starting index func
//             {
//                  pos_x[indext + base_offset] = (float) raw_data[v];
//                  pos_y[indext + base_offset] = (float) raw_data[v+1];

//                  indext++;
//             }
            
//             //now copy pos_x to raw_x with offset

//             for(int v = 0; v < mem_size; v++)
//             {
//                 raw_x[mem_size * i + v] = pos_x[v];
//                 raw_y[mem_size * i + v] = pos_y[v];
//             }

//             //this->unwrapData(&pos_x, n_polymers,l_polymer,t_adj);
//         }// end of collecting and unwrapping force data
//     }
// }


// vector<string> * GraphHandler::getSimulationDirectories(string path)
// {
//     int counter = 0;
//     vector<string> * dirs = new vector<string>();
//     for (auto&p: fs::directory_iterator(path))
//     {
//         if(p.is_directory() && contains(p.path(),".gsd"))
//             dirs->emplace_back(p.path().string());
            
    

//     }

//     return dirs;

    
// }

// bool GraphHandler::contains(fs::path path, string file_type)
// {
//     for (auto& p: fs::directory_iterator(path))
//     {
//         if(p.path().extension() == file_type)
//         {
//             return true;
//         }
//     }
//     return false;
// }


// int GraphHandler::createDensityPlots(Bounds bound, int xrez, int yrez)
// {

// }