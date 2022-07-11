#include <iostream>
#include <fstream>
//#include "threadpool.hpp"

using namespace std;

void apply_mask(string file, int start_line, int n_lines, string mask, string& output){
    
    //validate
    if(start_line%2 != 0 || n_lines%2!=0){
        throw runtime_error("Thread cannot read odd number of lines ("
        +to_string(n_lines)+")or start at odd position ("+to_string(start_line)+")");
    }

    //reach start line
     ifstream input_stream(file);
    std::string line;
    for(int i=0; i<start_line;i++){
        getline(input_stream, line);
    }

    unsigned int abundance = 0;
    unsigned int kmer_id = start_line/2 + 1;
    int raw_mer_size = mask.size();

    for(int i=0; i<n_lines;i++){
        getline(input_stream, line);
        if(line[0] == '>')
        {
            //Read abundance
            abundance = std::stoi(line.substr(1));
            //std::cout << ">Abundance:" << abundance << std::endl;
        }
        else
        {
            if(line.size() != raw_mer_size){
                throw std::runtime_error("Malformed read, length: " + std::to_string(line.size())
                 + ". Expected length: " + std::to_string(raw_mer_size));
            }
            //Apply mask to read
            string masked_line("");
            for(unsigned int i = 0; i<raw_mer_size;i++){
                if(mask[i]=='1')
                {
                    masked_line.push_back(line[i]);
                }
            }
            //std::cout << masked_line << std::endl;

            //Write masked read to new file *abundance times
            for(unsigned int i = 1; i<= abundance; i++)
            {
                output.append(">K-mer "+ to_string(kmer_id) + " Abundance:" + 
                to_string(abundance)+ " Copy:" + to_string(i) + "\n");
                output.append(masked_line +"\n");
            }
            kmer_id++;
        }
    }
    input_stream.close();
}

int main(int argc, char* argv[])
{
    std::string k_mer_file = "mer_counts.txt";
    std::string mask_file = "mask.txt";
    std::string output_file = "masked_mer_counts.txt";

    //get path to mask file
    if(argc==2){
        mask_file = argv[1];
    } else if(argc >2){
        throw std::runtime_error("To many arguments, only path to mask file accepted.");
    }

    //read mask
    std::string mask = "";
    std::ifstream maskstream(mask_file);
    if(!maskstream.good()){
        throw std::runtime_error("mask file could not be opened!");
    }
    std::getline(maskstream,mask);
    maskstream.close();
    
    //check mask
    if(mask == ""){
        throw std::runtime_error("mask empty!");
    }
    int raw_mer_size = 0;
    int masked_mer_size = 0;
    for(char c : mask){
        if(c=='0'){
            raw_mer_size++;
        } else if(c=='1'){
            raw_mer_size++;
            masked_mer_size++;
        } else{
            throw std::runtime_error("invalid characters in mask. Only 0 and 1 allowed.");
        }
    }
    //Statistics
    std::cout << "Used Mask: " << mask << std::endl;
    std::cout << "Transform " << raw_mer_size << "-mers into " << masked_mer_size << "-mers" << std::endl;


    //Read k-mers file
    std::string line;

    std::ifstream infile(k_mer_file);
    if(!infile.good()){
        throw std::runtime_error("infile could not be opened!");
    }
    std::ofstream outfile(output_file);
    if(!outfile.good()){
        throw std::runtime_error("outfile could not be opened!");
    }
    
    //get number of lines in file
    int nr_lines = 0;
    while(std::getline(infile,line)){
        nr_lines++;
    }
    std::cout << "Number of lines: " << nr_lines << std::endl;

    //reset stream
    infile.close();

    

    //Initialize threads
    size_t nr_cores = 4;
    int step_size = ((nr_lines / 2)/nr_cores)*2;
    //ThreadPool threadpool(nr_cores);

    //execute threads
    string output("");
    string partial_output[nr_cores];
    for(size_t i = 0; i< nr_cores; i++){
        int start_line = step_size*i;
        
        if(i<nr_cores-1){
            apply_mask(k_mer_file,start_line,step_size,mask, partial_output[i]);
        } else{
            apply_mask(k_mer_file,start_line,nr_lines-(step_size*i),mask, partial_output[i]);
        }

        output.append(partial_output[i]);
    }    
    
    //write output
    outfile << output << endl;

    outfile.close();
    return 0;
}