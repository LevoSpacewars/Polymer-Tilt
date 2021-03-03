#include "Compiler.h"

Compiler::Compiler(string path, bool singular){

    if (singular)
    {
        string polyProfileName = path + "/_simulation_parameters.txt";
        string polyDataName = path + "/trajectory.gsd";

        this->current_path = path;
        this->definePolymerProfile(&polyProfileName, &profile);
        std::filesystem::remove("profileData.txt"); // this needs to be cleaned up
        this->compileData(&polyDataName, interval);
        this->writeResults(path);
    }
    else
    {
         vector<string> * paths = getSimulationDirectories(path);

        for (int i = 0; i < paths->size(); i++)
        {
            cout << i<< ". " + paths->at(i)<<endl;
        }
        cout <<paths->size()<< ". ALL"<<endl;

        cout << "choose which file to compile data for: ";
        int index = -1;
        cin >> index;
        if (index != paths->size())
        {
            
            string polyProfileName = paths->at(index) + "/_simulation_parameters.txt";
            string polyDataName = paths->at(index) + "/trajectory.gsd";

            this->current_path = paths->at(index);
            this->definePolymerProfile(&polyProfileName, &profile);
            std::filesystem::remove("profileData.txt");
            this->compileData(&polyDataName, interval);
            this->writeResults(paths->at(index));
            
        }

        else
        {
            for (int i = 0; i < paths->size(); i++)
            {
                cout << i<< ". " + paths->at(i)<<endl;
                string polyProfileName = paths->at(i) + "/_simulation_parameters.txt";
                string polyDataName = paths->at(i) + "/trajectory.gsd";


                this->definePolymerProfile(&polyProfileName, &profile);
                this->compileData(&polyDataName, interval);
                this->writeResults(paths->at(i));

            }
        }
    }
    
   
    
}

int Compiler::definePolymerProfile(string* parameter_file_location, PolymerProfile * p)
{
    ifstream inFile(*parameter_file_location);

    string x;
    if(!inFile){
        cout << "parameter file not found" <<endl;
        exit(1);
    }




    while ((getline(inFile,x))){
        if (x.find("sheerForceRange") != string::npos){
            float t [2] = {0,0};
            t[0] = atof( x.substr( x.find("=") + 1 , x.find(",") ).c_str() );
            t[1] = atof( x.substr( x.find(",") +1).c_str() );
            p->sheerForceRange[0] = t[0];
            p->sheerForceRange[1] = t[1];
        }
        else if (x.find("df")       !=string::npos){
            int v = 0;
            v = stoi( x.substr(x.find("=")+1).c_str());
            p->df = v;
        }
        else if (x.find("length")   !=string::npos){
            int v = 0;
            v = stoi( x.substr(x.find("=")+1).c_str());
            p->length = v;
        }
        else if (x.find("lines")    !=string::npos){
            int v = 0;
            v = stoi( x.substr(x.find("=")+1).c_str());
            p->lines = v;
        }
        else if (x.find("K")        !=string::npos){
            float v = 0;
            v = stof( x.substr(x.find("=")+1).c_str());
            p->k_amplitude = v;
        }
        else if (x.find("pull")     !=string::npos){
            float v = 0;
            v = stof( x.substr(x.find("=")+1).c_str());
            p->tension = v;
        }
        else if (x.find("amplitude")!=string::npos){
            float v = 0;
            v = stof( x.substr(x.find("=")+1).c_str());
            p->p_amplitude = v;
        }
        else if (x.find("gamma")    !=string::npos){
            float v = 0;
            v = stof( x.substr(x.find("=")+1).c_str());
            p->gamma = v;
        }
        else if (x.find("kbT")    !=string::npos){
            float v = 0;
            v = stof( x.substr(x.find("=")+1).c_str());
            p->kbT = v;
        }
        else if (x.find("dt")    !=string::npos){
            float v = 0;
            v = stof( x.substr(x.find("=")+1).c_str());
            p->dt =v ;

        }
        else if (x.find("probePeriod")    !=string::npos){
            int v = 0;
            v = stoi( x.substr(x.find("=")+1).c_str());
            p->sampleRate = v;
        }
        else if (x.find("runLength")    !=string::npos){
            int v = 0;
            v = stoi( x.substr(x.find("=")+1).c_str());
            p->runLength = v;
        }
        else if (x.find("BoxDimx")    !=string::npos){
            float v = 0;
            v = stof( x.substr(x.find("=")+1).c_str());
            p->boxdimx = v;

        }
        else if (x.find("BoxDimy")    !=string::npos){
            float v = 0;
            v = stof( x.substr(x.find("=")+1).c_str());
            p->boxdimy = v;
        }
        else if(x.find("Direction") !=string::npos)
        {
            p->read_direciton = x.substr(x.find("=")+1);
        }
    }

    inFile.close();



    return -1;
}
void Compiler::trackParticle(float** data, int index, int time_steps, int nParticles,float theta)
{
    float* x = *data;
    ofstream writeFile;
    string filename = this->current_path+"/" + "particle_"+ to_string(index) + ":"+to_string(theta) + ".txt";
    writeFile.open(filename,std::ios_base::trunc);
    for (int i = 0; i < time_steps -1; i++)
    {
        writeFile << x[i*nParticles + index]<< ",";
    }
    writeFile << x[(time_steps-1) * nParticles + index];
    writeFile.close();
    
    
}
int Compiler::compileData(string *filename, float interval)
{
    cout<<"compiling data"<<endl;
    runLength = (int) (this->profile.runLength / this->profile.sampleRate);
    truncateFiles();
    
    int l_polymer = this->profile.length;
    int n_polymers = this->profile.lines; 

    int df = (int)(this->profile.df);
    float conv  = (this->profile.sheerForceRange[1] - this->profile.sheerForceRange[0])/df;
    int t_step     = 0;
    int line     = 0;
    int n_runs   = (int)df;
    int c_run    = 0;
    int lstep    = 0;
    int t_adj   =0;
    float prev    = 0;
    float prevt   = 0;
    float L       = this->profile.boxdimx;
    int a       = ((int) (this->pi * l_polymer));
    int b       = ((int) (this->pf * l_polymer));
    float currentt = 0;
    int totalRunLength = runLength * df;

    int gsd_open_error = gsd_open(&this->handler,filename->c_str(), GSD_OPEN_READONLY);
    

    auto e = gsd_find_chunk(&this->handler,0,"particles/position");


    for (int i = 0; i< df;i++)
    {
        //Begin by allocating the arrays needed for position data;
        int current_run = i;
        float current_force = i * conv + this->profile.sheerForceRange[0];
        float theta = current_force/this->profile.tension;

        float * pos_x;
        float * pos_y;
        float * pos_xr;
        float * pos_yr;

        assert(l_polymer*n_polymers == e->N);


        int memblock = l_polymer*n_polymers * ((int)((1-interval)*runLength));
        pos_x = new float[memblock];
        pos_y = new float[memblock];
        pos_xr = new float[memblock];
        pos_yr = new float[memblock];


        cout<<i<<endl;
        //Next Sort through the data and "unwrap" the polymers from periodic boundary
        // Begin by 
        cout<<"creating raw data"<<endl;
        
        cout<<"done"<<endl;
        
        
        for(int j = (int)(runLength*interval); j < runLength; j++)
        {
            t_step = i*runLength + j;
            t_adj = j - (int)(runLength*interval);
            auto chunk_entry = gsd_find_chunk(&this->handler,t_step,"particles/position"); //retrives the chunk information from a time step from the gsd file
            float * raw_data = new float[e->N * e->M];
            int errorch = gsd_read_chunk(&this->handler,raw_data, chunk_entry); // retrives data from chunk
            
            if(errorch != 0){
                
            cout<<"read not valid"<<endl;
            cout<<"time-step:"<<t_step<<endl;
            cout << i*runLength<< "," << j<<endl;
            exit(1);

            } // if read is successfull

            int indext = 0;
            int base_offset = t_adj * (l_polymer*n_polymers);
            for (int v = 0; v < e->N*e->M; v=v+3) // 0,1,2 | ,3,4,5 |,6,7,8 ... // need to implement a proper starting index func 
            {
                 pos_x[indext + base_offset] = (float) raw_data[v];
                 pos_y[indext + base_offset] = (float) raw_data[v+1];
                 pos_xr[indext + base_offset] = (float) raw_data[v];
                 pos_yr[indext + base_offset] = (float) raw_data[v+1];

                 indext++;
            }
            
            // finally, correct for the boundary condition
            // kindof annoying to write out, so I just encaposlated into seperate function
            // boxdimy >> l_polymer, therefore no unwrap needed for y^hat

            this->unwrapData(&pos_x, n_polymers,l_polymer,t_adj);
            
            delete raw_data;
        }// end of collecting and unwrapping force data

        //Begin processing the data

        
        int adj_run = (int)((1-interval) * runLength);
        this->writeData("position_uw:" + to_string(theta),&pos_x,&pos_y,memblock,l_polymer*n_polymers);//debug only
        this->writeData("position_nw:" + to_string(theta),&pos_xr,&pos_yr,memblock,l_polymer*n_polymers);


        HeatMapParameters param;
        param.rezx = 100;
        param.rezy = 100;
        param.height = this->profile.length/3;
        param.width = 2*this->profile.lines;
        param.x = -this->profile.lines/2;
        param.y = 0;
        writeHeatMap(&pos_x,&pos_y, n_polymers*l_polymer,adj_run,i*conv,false,param,"sdf");
        

        //cout<<"tracking particles 0,100,200"<<endl;
        //trackParticle(&pos_x,0,adj_run,n_polymers*l_polymer,theta);
        //trackParticle(&pos_x,100,adj_run,n_polymers*l_polymer,theta);
        //trackParticle(&pos_x,200,adj_run,n_polymers*l_polymer,theta);

        cout<<"calculatig average position x"<<endl;
        float * avg_x = calcAveragePosition(&pos_x, n_polymers, l_polymer, adj_run);

        cout<<"calculating avg pos y"<<endl;
        float * avg_y = calcAveragePosition(&pos_y, n_polymers, l_polymer, adj_run);

        
        //float * avg_dx2 = calcAverageDxsqr(&pos_x, n_polymers, l_polymer, adj_run);
        cout<<"writePorfileoutput"<<endl;
        writePolymerSystem(&avg_x, &avg_y, n_polymers, l_polymer,current_path);

        cout<<"exportDensityFunction_avg"<<endl;
        exportDensityFunction_raw(&pos_x, &pos_y, n_polymers, l_polymer, adj_run,theta,current_path);
        exportDensityFunction_avg(&avg_x,&avg_y, n_polymers, l_polymer, theta, current_path);
        cout<<"calcAverageDx"<<endl;
        float * dx = calcAverageDx(&avg_x,n_polymers,l_polymer);

        cout<<"calcAverageLength"<<endl;
        float * length = calcAverageLength(&avg_x,&avg_y,n_polymers,l_polymer);

        cout<<"calcSystemoutput"<<endl; 
        float * output = calcSystemOutput(&dx,&length,0);



        this->output.emplace_back(output[0]);
        this->uoutput.emplace_back(output[1]);
        
        this->dx.emplace_back(dx[0]);
        this->udx.emplace_back(dx[1]);

        this->length.emplace_back(length[0]);
        this->ulength.emplace_back(length[1]);

        

        delete pos_x;
        delete pos_y;
        delete avg_x;
        delete avg_y;
        delete dx;
        delete length;
        delete output;
        delete pos_xr;
        delete pos_yr;







    }// end of force Iteration
    gsd_close(&this->handler);
    return -1;
}

void Compiler::writeData(string filename, float** xb, float**yb,int size,int blocksize)
{
    float * x = *xb;
    float * y = *yb;

    ofstream writeFile;
    filename = this->current_path+"/" + filename + ".txt";
    writeFile.open(filename,std::ios_base::trunc);
   
    for(int i =0; i < size;i++)
    {
        if(i%blocksize==0 & i !=0)
        {
            writeFile<<"\n";
        }
        writeFile << x[i]<<","<<y[i]<<",";
    } 

    writeFile.close();

    
    
    
}

bool Compiler::truncateFiles()
{
    ofstream writeFile;

    writeFile.open(this->current_path + "/DensityData_avg.txt",std::ios_base::trunc);
    writeFile.close();

    writeFile.open(this->current_path + "/data.txt",std::ios_base::trunc);
    writeFile.close();

    writeFile.open(this->current_path + "/heatmaps.txt",std::ios_base::trunc);
    writeFile.close();

    writeFile.open(this->current_path + "/profileData.txt",std::ios_base::trunc);
    writeFile.close();
    

}

bool Compiler::exportDensityFunction_avg(float** xa, float ** ya, int p_n, int p_length,float force_value, string path)
{

    float width = this->profile.boxdimx;
    float interval = width/p_n;
    
    


    float* x = *xa;
    float* y = *ya;
    float* sortOrder = new float[p_n];
    float* original= new float[p_n];
    int* offsetOrder = new int[p_n];
    
    for (int i = 0; i < p_n;i++)
    {
        sortOrder[i] = x[i*p_length+1];
        original[i] = sortOrder[i];
    }
    bool notdone = true;
    float temp = sortOrder[0];
    while(notdone == true){
        bool reordered = false;
        for (int i = 0; i < p_n - 1;i++)
        {
            if (sortOrder[i] > sortOrder[i+1]){
                temp = sortOrder[i];
                sortOrder[i] = sortOrder[i+1];
                sortOrder[i+1] = temp;
                reordered = true;
            }
        }
        if (reordered == false)
        {
            notdone = false;  
        }
    }
    
    for (int i = 0; i < p_n;i++)
    {
        cout<<sortOrder[i]<<",";
    }
    cout<<endl;
    for (int i = 0; i < p_n;i++)
    {
        cout<<original[i]<<",";
    }
    cout<<endl; 
    
    for (int i = 0; i < p_n;i++)
    {
        for(int j = 0; j < p_n;j++)
        {
            if(sortOrder[i] == original[j])
            {
                offsetOrder[i] = j;
            }
        }
    }
    for (int i = 0; i < p_n;i++)
    {
        cout<<offsetOrder[i]<<",";
    }
    cout<<endl;


    
    int unc_offset = p_n*p_length;
    ofstream writeFile;
    float conv = this->profile.boxdimx/p_n;
    writeFile.open(this->current_path + "/DensityData_avg:" + to_string(force_value) + ".txt",std::ios_base::trunc);

    writeFile << "parameters (p_n,p_l,force,Theta):" + to_string(p_n) + "," + to_string(p_length) + "," +to_string(force_value) + "," + to_string(0) <<endl;
    writeFile <<"x,ux,y"<<endl;
    int iter = 0;
    for (int i = 0; i < p_n; i++)
    {

        for(int j = offsetOrder[i]*p_length; j < (offsetOrder[i]+1)*p_length; j++)
        {
            writeFile<< x[j] - (i)*conv<< "," << x[j + unc_offset]<<","<< y[j] <<endl;
            iter++;
        }
        
        

    }
    cout<<iter<<endl;
    
    //cout<<endl;

  
    


    writeFile.close();
    
    return true;
}
bool Compiler::exportDensityFunction_raw(float** xa, float ** ya, int p_n, int p_length,int time_length,float force_value, string path)//M SAFE
{
    
 
    float width = this->profile.boxdimx;
    float interval = width/p_n;
    
    


    float* x = *xa;
    float* y = *ya;
    float* sortOrder = new float[p_n];
    float* original= new float[p_n];
    int* offsetOrder = new int[p_n];
    
    for (int i = 0; i < p_n;i++)
    {
        sortOrder[i] = x[i*p_length+1];
        original[i] = sortOrder[i];
    }
    bool notdone = true;
    float temp = sortOrder[0];
    while(notdone == true){
        bool reordered = false;
        for (int i = 0; i < p_n - 1;i++)
        {
            if (sortOrder[i] > sortOrder[i+1]){
                temp = sortOrder[i];
                sortOrder[i] = sortOrder[i+1];
                sortOrder[i+1] = temp;
                reordered = true;
            }
        }
        if (reordered == false)
        {
            notdone = false;  
        }
    }
    
    for (int i = 0; i < p_n;i++)
    {
        cout<<sortOrder[i]<<",";
    }
    cout<<endl;
    for (int i = 0; i < p_n;i++)
    {
        cout<<original[i]<<",";
    }
    cout<<endl; 
    
    for (int i = 0; i < p_n;i++)
    {
        for(int j = 0; j < p_n;j++)
        {
            if(sortOrder[i] == original[j])
            {
                offsetOrder[i] = j;
            }
        }
    }
    for (int i = 0; i < p_n;i++)
    {
        cout<<offsetOrder[i]<<",";
    }
    cout<<endl;
 
 
 
    int unc_offset = p_n*p_length*time_length/2;
    
   
    cout<<"writing raw"<<endl;
    ofstream writeFile;
    ofstream comFile;
    ofstream denFile;
    float conv = this->profile.boxdimx/p_n;
    writeFile.open(this->current_path + "/DensityData_raw:" + to_string(force_value) + ".txt",std::ios_base::trunc);
    comFile.open(this->current_path + "/COM:" + to_string(force_value) + ".txt",std::ios_base::trunc);
    denFile.open(this->current_path + "/DensityMap_raw:" + to_string(force_value) + ".txt",std::ios_base::trunc);
    writeFile << "parameters (p_n,p_l,force,Theta):" + to_string(p_n) + "," + to_string(p_length) + "," +to_string(force_value) + "," + to_string(0) <<endl;
    writeFile <<"x,y"<<endl;
    ofstream debugFile;
    debugFile.open(this->current_path + "/DebugDensityData_raw:" + to_string(force_value) + ".txt",std::ios_base::trunc);
        

    for (int k = 0; k<time_length;k++)
    {
        float com = 0;
        float avg = 0;
        float avg2= 0;
        int da = 0;
        int offset = k * p_length*p_n;
        for (int j = offsetOrder[0]*p_length; j < (offsetOrder[0]+1)*p_length; j++)
        {
            avg += x[j+offset];
        }
        avg= avg/(p_length);

        for (int i = 0; i < p_n; i++)
        {
            avg2=0;
            for (int j = offsetOrder[i]*p_length; j < (offsetOrder[i]+1)*p_length; j++)
            {
                avg2 += x[j+offset];
            }
            avg2 = avg2/(p_length);
            com += avg2;
            da = avg - avg2;
            for(int j = offsetOrder[i]*p_length; j < (offsetOrder[i]+1)*p_length; j++)
            {
                

                float calc = ((x[j+offset]) - avg) - i*conv;
                if (calc > 4 | calc < -5)
                {
                    debugFile << "T:" + to_string(k) + " P_n:" + to_string(i) + " " << x[j+offset] << " " + to_string(calc)<<endl;
                }
                denFile << ((x[j+offset]) - avg) << "," << y[j+offset]<<endl;
                writeFile<< calc<< "," << y[j+offset] <<endl;
            }
            

        }
        com = com/p_n;
        comFile << com <<","; 
    }
    


    comFile.close();
    debugFile.close();
    writeFile.close();

    
    return true;
}

bool Compiler::writeProfileOutput(float** xa, float ** ya, int p_n, int p_length,float force_value, string path) //memory safe
{
    
    float* x = *xa;
    float* y = *ya;
    int unc_offset = p_n*p_length;
    ofstream writeFile;
    writeFile.open(this->current_path + "/profileData.txt",std::ios_base::app);
    for (int i = 0; i < p_n; i++)
    {
        writeFile<< "Polymer," << i<<endl;
        writeFile<< "ForceValue," << force_value<<endl;

        for(int j = i*p_length; j < (i+1)*p_length; j++)
        {
            writeFile<< x[j] << "," << x[j+unc_offset] << "," << y[j] << "," << y[j + unc_offset]<<endl;
        }

    }

    writeFile.close();

    
    return true;

}

bool Compiler::writePolymerSystem(float** xa, float ** ya, int p_n, int p_length,string path)
{

    float* x = *xa;
    float* y = *ya;
    int unc_offset = p_n*p_length;
    ofstream writeFile;
    writeFile.open(this->current_path + "/polymersystem.txt",std::ios_base::trunc);
    for (int i = 0; i < p_n; i++)
    {
        for(int j = i*p_length; j < (i+1)*p_length; j++)
        {
            writeFile<< x[j] << "," << y[j] <<endl;
        }

    }

    writeFile.close();

    
    return true;


    return true;
}

bool Compiler::writeHeatMap(float** xs, float** ys, int time_steps, int nParticles, float force_value, bool normalized, HeatMapParameters param, string path) //memory safe
{
    float *x = *xs;
    float *y = *ys;

    
    //initialize spacing arrays
    float * tablex = new float[param.rezx];
    float * tabley = new float[param.rezy];

    //initalized output
    float * heatmap = new float[param.rezx * param.rezy];

    bool done = false;

    std::cout<<"writing HeatMap"<<endl;
    
    

    int offset = 0;

    float xconv = param.width/param.rezx;
    float yconv = param.height/param.rezy;
    for (int i = 0; i < param.rezx*param.rezy; i++) //zero heatmap
    {
        heatmap[i]=0;
    }

    //assign values for spacing on spacing arrays
    for (int i = 0; i < param.rezx; i++)
    {
        tablex[i] = xconv * (i+1) + param.x;
        
    }
    cout<<endl;

    for (int i = 0; i < param.rezy; i++)
    {
        tabley[i] = yconv * (i+1) + param.y;
        cout<<tabley[i]<<",";
    }
    cout<<endl;



    //begin sort
    for (int i = 0; i < time_steps;i++)
    {
        offset = nParticles * i;
        for (int j = 0; j < nParticles; j++)
        {
            int a = 0;
            int b = 0;

           
            bool donea = false;
            bool doneb = false;
            for (int k = 0; k < param.rezx;k++)
            {
                if (x[offset + j] <= tablex[k])
                {
                    a = k;
                    donea = true;
                    break;
                }
            }
            for (int k = 0; k < param.rezy;k++)
            {
                if (y[offset + j] <= tabley[k])
                {
                    b = k;
                    doneb=true;
                    break;
                }
            }

           
            
                

                
            
            if (donea == false | doneb == false)
            {
                cout<<"sort didn't work"<<endl;
                cout<<x[offset + j]<<","<<y[offset+j]<<endl;
                cout<<donea<<","<<doneb<<endl;
                exit(1);
            }
            int h_index = param.rezx* (b) + a;
            heatmap[h_index] +=1;
        }
        
    }

    //normalize?

    if(normalized)
    {
        int N = 0;
        float largest = 0;
        for (int i =1; i < param.rezy*param.rezx;i++)
        {
            
            largest = (heatmap[i] >= largest) * heatmap[i] + (heatmap[i] < largest) * largest;
        }

        float weight = 255.0/largest;

        for (int i =1; i < param.rezy*param.rezx;i++)
        {
            heatmap[i]= heatmap[i]*weight;
        }
    }

    ofstream myfile;
    myfile.open (this->current_path + "/heatmap"+ to_string(force_value/this->profile.tension) + ".txt",ios::trunc);
    myfile << "parameters (x,y,w,h,rx,ry):" + to_string(param.x) + "," + to_string(param.y) + "," +to_string(param.width) + "," + to_string(param.height) + "," + to_string(param.rezx) + "," + to_string(param.rezy) <<endl<<"{";

    for (int i = 0; i < param.rezy;i++)
    {
        for(int j = 0; j < param.rezx -1;j++)
        {
            myfile<< heatmap[i* param.rezx + j] << ",";
        }
        myfile<< heatmap[i* param.rezx + param.rezx] <<endl;
    }
    myfile.close();

    delete heatmap;
    delete tabley;
    delete tablex;



}

void Compiler::unwrapData(float ** data, int n_polymers, int l_polymer,int step){//M SAFE


    int base_offset = step*n_polymers*l_polymer;
    int base_offset_prev = (step-1)*n_polymers*l_polymer * int(bool(base_offset));
    int offset = 0;
    float* px = *data;

    float dl = 0;
    float diff = 0;

    float L = this->profile.boxdimx;

    for (int k = 0; k < n_polymers; k++)
    {
        int poffset = k* l_polymer;
        float dlt = px[base_offset + poffset] - px[base_offset_prev + poffset];
        if(abs(dlt/L) >= 0.5)
        {
            px[base_offset + poffset] -= L*(int(dlt/L));
            diff = dlt/L - int(dlt/L);
            if (abs(diff) > 0.7)
            {
                px[base_offset + poffset] -= copysign(L,dlt);
            }
        }
    
    }


    for (int i = 1 + base_offset; i < n_polymers*l_polymer + base_offset;i++)
    {
        dl = px[i] - px[i-1];
        if(abs(dl/L) >= 0.5)
        {
            px[i] = px[i] - L* (int(dl/L));
            diff = dl/L - int(dl/L);
            if (abs(diff) > 0.7)
            {
                px[i] -=copysign(L,dl);
            }
        }    
    }
    for (int k = 1; k < n_polymers; k++) // iterates over the number of polymers
    {
        offset = base_offset + l_polymer * k;
        
        for(int i = offset; i < offset + l_polymer; i++)
        {
            dl = px[i] - px[i-1];
            if(abs(dl/L) >= 0.5)
            {
                px[i] -= L* (int(dl/L));
                diff = dl/L - int(dl/L);
                if (abs(diff) > 0.7)
                {
                    px[i] -=copysign(L,dl);
                }
            }
        }

    }


}



float* Compiler::calcAveragePosition(float ** data, int n_polymers, int l_polymer,int sampleLength)//M SAFE
{
    float* r_data = *data;
    float *avg_unc = new float[(2 * n_polymers*l_polymer)]; // array where [0,N/2) is average positional data, and [N/2,N) the uncertainty on that
    // this shouldn't be needed anymore since I am now working with small N values. However, I am not sure how objects and memory work, so I am going to keep this clean
    
    int offset = 0;
    int index = 0;
    int unc_index = n_polymers*l_polymer;
    
    //quick test

    for (int i = 0; i < 2*n_polymers*l_polymer;i++)
    {
        avg_unc[i] = 0;
    }

    

    for (int polymer = 0; polymer < n_polymers; polymer++)
    {
        for(int particle  = 0; particle < l_polymer; particle++)
        {
            float sum = 0;
            float square = 0;

            for (int t_step = 0; t_step < sampleLength; t_step++)
            {
                offset = polymer * l_polymer + particle;
                index = t_step * l_polymer * n_polymers + offset;
                sum += r_data[index];
                square += r_data[index]*r_data[index];
            }

            float avg = sum / sampleLength;

            float unc = square - sum * sum / sampleLength;
            unc = unc / sampleLength;
            unc = pow(unc, 0.5);

            avg_unc[offset] = avg;
            avg_unc[unc_index + offset] = unc;
            

            
        }
    }


    

    return avg_unc;

}

float* Compiler::calcAverageDx(float ** avg_unc_x, int n_polymer, int l_polymer)
{
    float * system_dx = new float[2];
    float * data = *avg_unc_x;
    float tilt_temp = 0;
    int a = (int)(pi*l_polymer);
    int b = (int)(pf*l_polymer);
    
    
    float sum = 0;
    float square = 0;
    
    for ( int polymer = 0; polymer < n_polymer; polymer++)
    {
       

        sum += data[b] - data[a];
        // cout<< data[b]<<endl;
        square += pow(data[b] - data[a],2);
        a += l_polymer;
        b += l_polymer;
        assert(a != b);
        
    }
    
    
    system_dx[0] = sum/n_polymer;
    // cout<<sum<<endl;
    float unc = square - sum * sum / n_polymer;
    unc = unc / n_polymer;
    unc = pow(unc, 0.5);
    system_dx[1] = unc;


    return system_dx;
}

float* Compiler::calcAverageLength(float ** avg_unc_x, float ** avg_unc_y, int n_polymer, int l_polymer)
{
    float *system_length = new float[2]; 

    float *x = *avg_unc_x;
    float *y = *avg_unc_y;

    int a = (int)(l_polymer*pi);
    int b = (int)(l_polymer*pf);


    float sum = 0;
    float square = 0;
    for (int polymer = 0; polymer < n_polymer; polymer++)
    {
        

        float dx = x[b] - x[a];
        float dy = y[b] - y[a];
        cout << a << " " << b << " " << n_polymer * l_polymer << endl;
        float l  = pow(dx*dx+ dy*dy, 0.5);

        sum += l;
        square += l*l;
        a += l_polymer;
        b += l_polymer;
    }

    system_length[0] = sum/n_polymer;
    float unc = square - sum * sum / n_polymer;
    unc = unc / n_polymer;
    unc = pow(unc, 0.5);
    system_length[1] = unc;

    return system_length;
}

float* Compiler::calcAverageDxsqr(float ** data, int n_polymers, int l_polymer,int sampleLength)
{
    float* r_data = *data;
    float *avg_unc = new float[(n_polymers*l_polymer)]; // array where [0,N/2) is average positional data, and [N/2,N) the uncertainty on that
    // this shouldn't be needed anymore since I am now working with small N values. However, I am not sure how objects and memory work, so I am going to keep this clean
    
    int offset = 0;
    int index = 0;
    int pindex = 0;
    int unc_index = n_polymers*l_polymer;
    
    //quick test

    for (int i = 0; i < 2*n_polymers*l_polymer;i++)
    {
        avg_unc[i] = 0;
    }

    

    for (int polymer = 0; polymer < n_polymers; polymer++)
    {
        for(int particle  = 0; particle < l_polymer; particle++)
        {
            float sum = 0;
            float square = 0;

            for (int t_step = 1; t_step < sampleLength; t_step++)
            {
                offset = polymer * l_polymer + particle;
                index = t_step * l_polymer * n_polymers + offset;
                pindex = (t_step-1) * l_polymer * n_polymers + offset;// build off of this here
                sum += pow(r_data[index] - r_data[pindex],2);
            }

            float avg = sum / sampleLength;
            avg_unc[offset] = avg;
        }
    }


    

    return avg_unc;
}

float* Compiler::calcSystemOutput(float** sysdx, float ** syslength, float sheerTension)
{
    float * dx = *sysdx;
    float * length = *syslength;

    float * output = new float[2];

    output[0] = dx[0]/length[0];


    cout<<length[1]/length[0]<< ","<<dx[1]/dx[0]<<endl;
    output[1] = output[0] * pow( pow(0, 2) + pow(dx[1] / dx[0], 2), 0.5);

    return output;
}

float * Compiler::calcSheerTension()
{
    float * st = new float[this->profile.df];

    for ( int i = 0; i < this->profile.df; i++)
    {
        st[i] = this->profile.sheerForceRange[0] + i*(this->profile.sheerForceRange[1] -this->profile.sheerForceRange[0] ) / this->profile.df;
        st[i] = st[i] / this->profile.tension;
    }


    return st;

}

int Compiler::writeResults(string path)
{
    if (this->profile.read_direciton == "forward")
    {
        ofstream writeFile;
        writeFile.open(this->current_path + "/data.txt", ios::trunc);
        writeFile << "dx,";
        for (int i = 0; i < this->dx.size()-1;i++){
            writeFile << this->dx.at(i) << ",";
        }
        writeFile << this->dx.at(this->dx.size()-1)<<endl;

        writeFile << "udx,";
        for (int i = 0; i < this->dx.size()-1;i++){
            writeFile << this->udx.at(i) << ",";
        }
        writeFile << this->udx.at(this->dx.size()-1)<<endl;

        writeFile << "length,";
        for (int i = 0; i < this->dx.size()-1;i++){
            writeFile << this->length.at(i) << ",";
        }
        writeFile << this->length.at(this->dx.size()-1)<<endl;

        writeFile << "ulength,";
        for (int i = 0; i < this->dx.size()-1;i++){
            writeFile << this->ulength.at(i) << ",";
        }
        writeFile << this->ulength.at(this->dx.size()-1)<<endl;

        writeFile << "output,";
        for (int i = 0; i < this->dx.size()-1;i++){
            writeFile << this->output.at(i) << ",";
        }
        writeFile << this->output.at(this->dx.size()-1)<<endl;

        writeFile << "uoutput,";
        for (int i = 0; i < this->dx.size()-1;i++){
            writeFile << this->uoutput.at(i) << ",";
        }
        writeFile << this->uoutput.at(this->dx.size()-1)<<endl;

       




        writeFile.close();
    }

    else
    {
        ofstream writeFile;
        writeFile.open(this->current_path + "/data.txt", ios::trunc);
        writeFile << "dx,";
        for (int i = this->dx.size() -1; i > 0 ;i--){
            writeFile << this->dx.at(i) << ",";
        }
        writeFile << this->dx.at(0)<<endl;

        writeFile << "udx,";
        for (int i = this->udx.size() -1; i > 0 ;i--){
            writeFile << this->udx.at(i) << ",";
        }
        writeFile << this->udx.at(0)<<endl;


        writeFile << "length,";
        for (int i = this->length.size() -1; i > 0 ;i--){
            writeFile << this->length.at(i) << ",";
        }
        writeFile << this->length.at(0)<<endl;

        writeFile << "ulength,";
        for (int i = this->ulength.size() -1; i > 0 ;i--){
            writeFile << this->ulength.at(i) << ",";
        }
        writeFile << this->ulength.at(0)<<endl;

        writeFile << "output,";
        for (int i = this->output.size() -1; i > 0 ;i--){
            writeFile << this->output.at(i) << ",";
        }
        writeFile << this->output.at(0)<<endl;

        writeFile << "uoutput,";
        for (int i = this->uoutput.size() -1; i > 0 ;i--){
            writeFile << this->uoutput.at(i) << ",";
        }
        writeFile << this->uoutput.at(0)<<endl;

    }
    


    return 0;

}
vector<string> * Compiler::getSimulationDirectories(string path)
{
    int counter = 0;
    vector<string> * dirs = new vector<string>();
    for (auto&p: fs::directory_iterator(path))
    {
        if(p.is_directory() && contains(p.path(),".gsd"))
            dirs->emplace_back(p.path().string());
            
    

    }

    return dirs;

    
}

bool Compiler::contains(fs::path path, string file_type)
{
    for (auto& p: fs::directory_iterator(path))
    {
        if(p.path().extension() == file_type)
        {
            return true;
        }
    }
    return false;
}
