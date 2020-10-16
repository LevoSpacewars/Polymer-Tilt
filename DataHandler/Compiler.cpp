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

int Compiler::compileData(string *filename, float interval)
{
    cout<<"compiling data"<<endl;
    runLength = (int) (this->profile.runLength / this->profile.sampleRate);
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

    int gsd_open_error = gsd_open(&this->handler,filename->c_str(), GSD_OPEN_READONLY);

    auto e = gsd_find_chunk(&this->handler,0,"particles/position");


    for (int i = 0; i< df;i++)
    {
        //Begin by allocating the arrays needed for position data;
        int current_run = i;
        float current_force = i * conv + this->profile.sheerForceRange[0];

        float * pos_x;
        float * pos_y;
        float * pos_xr;
        float * pos_yr;

        assert(l_polymer*n_polymers == e->N);

        // pos_x = (float*)malloc(l_polymer*n_polymers * ((int)(1-interval))*runLength * sizeof(float)); //DEPRICATED
        // pos_y = (float*)malloc(l_polymer*n_polymers * ((int)(1-interval))*runLength * sizeof(float));
        int memblock = l_polymer*n_polymers * ((int)((1-interval)*runLength));
        pos_x = new float[memblock];
        pos_y = new float[memblock];
        pos_xr = new float[memblock];
        pos_yr = new float[memblock];

        
        cout<<i<<endl;
        //Next Sort through the data and "unwrap" the polymers from periodic boundary
        // Begin by 
        float * raw_data = new float[e->N * 3];
        
        for(int j = (int)(runLength*interval); j < runLength-1; j++)
        {
            t_step = i*runLength + j;
            cout<<t_step<<endl;
            cout<< runLength<<endl;
            t_adj = j - (int)(runLength*interval);
            
            auto chunk_entry = gsd_find_chunk(&this->handler,t_step,"particles/position"); //retrives the chunk information from a time step from the gsd file
            //cout<<"assigning raw data"<<endl;
            //cout<< e->N * 3 * sizeof(float) <<endl;
             // (number of particles) by (dimensions)
            // cout<<"failed?"<<endl;
            int errorch = gsd_read_chunk(&this->handler,raw_data, chunk_entry); // retrives data from chunk
            assert(errorch == 0); // if read is successfull

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
        }// end of collecting and unwrapping force data

        //Begin processing the data
        
        int adj_run = (int)((1-interval) * runLength);
        HeatMapParameters param;
        param.rezx = 100;
        param.rezy = 100;
        param.height = this->profile.length /10;
        param.width = this->profile.lines;
        param.x = -this->profile.lines/2;
        param.y = 0;
        writeHeatMap(&pos_xr,&pos_yr, n_polymers*l_polymer,adj_run,i*conv,true,param,"sdf");

        float * avg_x = calcAveragePosition(&pos_x, n_polymers, l_polymer, adj_run);
        float * avg_y = calcAveragePosition(&pos_y, n_polymers, l_polymer, adj_run);
        writeProfileOutput(&avg_x, &avg_y, n_polymers, l_polymer, current_force ,current_path);

        exportDensityFunction_avg(&avg_x, &avg_y, n_polymers, l_polymer, current_force ,current_path);
        

        float * dx = calcAverageDx(&avg_x,n_polymers,l_polymer);
        cout<<"here"<<endl;
        float * length = calcAverageLength(&avg_x,&avg_y,n_polymers,l_polymer);
        cout<<"here"<<endl;
        float * output = calcSystemOutput(&dx,&length,0);
        cout<<"here"<<endl; 
        this->output.emplace_back(output[0]);
        this->uoutput.emplace_back(output[1]);
        
        this->dx.emplace_back(dx[0]);
        this->udx.emplace_back(dx[1]);

        this->length.emplace_back(length[0]);
        this->ulength.emplace_back(length[1]);

        

        delete(pos_x);
        delete(pos_y);
        delete(avg_x);
        delete(avg_y);
        delete(dx);
        delete(length);
        delete(output);







    }// end of force Iteration

    try {
        std::filesystem::remove(current_path + "/profileData.txt");
        std::filesystem::copy("profileData.txt", current_path+"/profileData.txt");
    } catch (std::filesystem::filesystem_error& e) {
        std::cout << e.what() << '\n';
    }

    return -1;
}

bool Compiler::exportDensityFunction_avg(float** xa, float ** ya, int p_n, int p_length,float force_value, string path)
{
    float* x = *xa;
    float* y = *ya;
    int unc_offset = p_n*p_length;
    ofstream writeFile;
    float conv = this->profile.boxdimx/p_n;
    writeFile.open(this->current_path + "/DensityData_avg.txt",std::ios_base::app);
    for (int i = 0; i < p_n; i++)
    {
        writeFile<< "Polymer," << i<<endl;
        writeFile<< "ForceValue," << force_value<<endl;

        for(int j = i*p_length; j < (i+1)*p_length; j++)
        {
            writeFile<< x[j] - i*conv << "," << x[j+unc_offset] << "," << y[j] << "," << y[j + unc_offset]<<endl;
        }

    }

    writeFile.close();

    
    return true;
}
bool Compiler::exportDensityFunction_raw(float** xa, float ** ya, int p_n, int p_length,float force_value, string path)
{
    float* x = *xa;
    float* y = *ya;
    int unc_offset = p_n*p_length;
    ofstream writeFile;
    float conv = this->profile.boxdimx/p_n;
    writeFile.open("DensityData_raw.txt",std::ios_base::app);
    for (int i = 0; i < p_n; i++)
    {
        writeFile<< "Polymer," << i<<endl;
        writeFile<< "ForceValue," << force_value<<endl;

        for(int j = i*p_length; j < (i+1)*p_length; j++)
        {
            writeFile<< x[j] - i* conv<< "," << x[j+unc_offset] << "," << y[j] << "," << y[j + unc_offset]<<endl;
        }

    }

    writeFile.close();

    
    return true;
}

bool Compiler::writeProfileOutput(float** xa, float ** ya, int p_n, int p_length,float force_value, string path)
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



bool Compiler::writeHeatMap(float** xs, float** ys, int time_steps, int nParticles, float force_value, bool normalized, HeatMapParameters param, string path)
{
    float *x = *xs;
    float *y = *ys;

    

    float * tablex = new float[param.rezx];
    float * tabley = new float[param.rezx];
    float * heatmap = new float[param.rezx * param.rezy ];

    bool done = false;

    std::cout<<"writing HeatMap"<<endl;
    
    

    int offset = 0;

    float xconv = param.width/param.rezx;
    float yconv = param.height/param.rezy;
    for (int i = 0; i < param.rezx*param.rezy; i++)
    {
        heatmap[i]=0;
    }
    for (int i = 0; i < param.rezx; i++)
    {
        
        tablex[i] = xconv * i + param.x;
        cout <<tablex[i]<< " ";
    }
    cout<<endl;
    for (int i = 0; i < param.rezy; i++)
    {
        tabley[i] = yconv * i + param.y;
        
    }

    std::cout<<"done init"<<endl;
    for (int i = 0; i < time_steps;i++)
    {
        offset = nParticles * i;
        for (int j = 0; j< nParticles; j++)
        {
            int a = 0;
            int b = 0;

           

            int left = 0;
            int right = param.rezx;
            int m = 0;
            int iter = 0;
            done = false;
            while (left <= right && done == false)
            {
                
                
                m = int((left + right)/2);
                if (tablex[m] <= x[offset + j])
                {
                    
                    
                    if (tablex[m+1] >= x[offset + j])
                    {
                        
                        
                        a = m;
                        done = true;
                    }
                    else
                    {
                        left = m+1;
                    }
                }
                    
                else if(tablex[m] >= x[offset + j]) 
                {
                    
                    if (tablex[m-1] <= x[offset + j])
                    {
                        
                        a = m-1;
                        done = true;
                    }
                    else
                    {
                        right = m-1;
                    }
                    
                }
            }
            



            

            left = 0;
            done = false;
            right = param.rezy;

             while (left <= right && done == false)
            {
                
                
                m = int((left + right)/2);
                if (tabley[m] <= y[offset + j])
                {
                    
                    
                    if (tabley[m+1] >= y[offset + j])
                    {
                        
                        
                        b = m;
                        done = true;
                    }
                    else
                    {
                        left = m+1;
                    }
                }
                    
                else if(tabley[m] >= y[offset + j]) 
                {
                    
                    if (tabley[m-1] <= y[offset + j])
                    {
                        
                        b = m-1;
                        done = true;
                    }
                    else
                    {
                        right = m-1;
                    }
                    
                }
            }


           
            
                

                
            

            int h_index = param.rezx* (b) + a;
            heatmap[h_index] +=1;
        }
        cout<< float(i)/float(time_steps)<<endl;
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

    //print the values in a readable way

    for (int i = 1; i < param.rezy;i++)
    {
        for(int j = 0; j < param.rezx-1;j++)
        {
            cout<< heatmap[i* param.rezx + j] << ",";
        }
        cout<<endl;
    }

    ofstream myfile;
    myfile.open (this->current_path + "/heatmaps.txt",ios::app);
    myfile << "parameters (x,y,w,h,rx,ry):" + to_string(param.x) + "," + to_string(param.y) + "," +to_string(param.width) + "," + to_string(param.height) + "," + to_string(param.rezx) + "," + to_string(param.rezy) <<endl<<"{";

    for (int i = 1; i < param.rezy;i++)
    {
        for(int j = 0; j < param.rezx-1;j++)
        {
            myfile<< heatmap[i* param.rezx + j] << ",";
        }
        cout<<endl;
    }
    myfile.close();

    delete heatmap;
    delete tabley;
    delete tablex;



}

void Compiler::unwrapData(float ** data, int n_polymers, int l_polymer,int step){


    int base_offset = step*n_polymers*l_polymer;
    float* px = *data;
    float currentt = 0;
    float prevt   = 0;
    float dl = 0;
    float L = this->profile.boxdimx;
    for (int k = 0; k < n_polymers; k++) // iterates over the number of polymers
            {
               int line = k;
               int offset = base_offset + k * l_polymer;
                for (int particle = 1; particle < l_polymer; particle++) //iterates over the number of particles within a polymer
                {

                    // here I am compensating for the periodic box, like pac-man: if a particle leaves the right side it will come out the left side, +- L to fix this
                    currentt = px[offset + particle];
                    prevt = px[base_offset + (particle-1)];
                    dl  = currentt - prevt;
                    

                    if (  dl < -L/2 ){
                        px[offset + particle] = currentt + L;
                    }
                    else if ( dl > L/2 ){
                        px[offset + particle] = currentt - L;
                    }

                }

            }


}

float* Compiler::calcAveragePosition(float ** data, int n_polymers, int l_polymer,int sampleLength)
{
    float* r_data = *data;
    float *avg_unc = new float[(2 * n_polymers*l_polymer)]; // array where [0,N/2) is average positional data, and [N/2,N) the uncertainty on that
    // this shouldn't be needed anymore since I am now working with small N values. However, I am not sure how objects and memory work, so I am going to keep this clean
    
    int offset = 0;
    int index = 0;
    int unc_index = n_polymers*l_polymer;
    


    

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
        cout<< data[b]<<endl;
        square += pow(data[b] - data[a],2);
        a += l_polymer;
        b += l_polymer;
        assert(a != b);
        
    }
    
    
    system_dx[0] = sum/n_polymer;
    cout<<sum<<endl;
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


float* Compiler::calcSystemOutput(float** sysdx, float ** syslength, float sheerTension)
{
    float * dx = *sysdx;
    float * length = *syslength;

    float * output = new float[2];

    output[0] = dx[0]/length[0];


    output[1] = output[0] * pow( pow(length[1] / length[0], 2) + pow(dx[1] / dx[0], 2), 0.5);

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
        writeFile << this->output.at(this->dx.size()-1)<<endl;

       




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
    try {
        std::filesystem::remove(path + "/data.txt");
        std::filesystem::copy("data.txt", path+"/data.txt");
        //std::filesystem::remove("data.txt");
    } catch (std::filesystem::filesystem_error& e) {
        std::cout << e.what() << '\n';
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
