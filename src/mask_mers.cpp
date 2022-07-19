    #include <iostream>
    #include <fstream>
    #include <thread>
    #include <vector>
    #include "command_line_parser.hpp"

    using namespace std;

    void apply_mask(string file, int start_line, int n_lines, string mask, string* output){

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

        //variables
        unsigned int abundance = 0;
        unsigned int kmer_id = start_line/2 + 1;
        int raw_mer_size = mask.size();

        for(int i=0; i<n_lines;i++){
            getline(input_stream, line);

            if(line[0] == '>') {
                 //Read abundance
                abundance = std::stoi(line.substr(1));
            }
            else{
                if(line.size() != raw_mer_size){
                    throw std::runtime_error("Malformed read, length: " + std::to_string(line.size())
                    + ". Expected length: " + std::to_string(raw_mer_size));
                }

                //Apply mask to read
                string masked_line("");
                for(unsigned int i = 0; i<raw_mer_size;i++){
                    if(mask[i]=='1'){
                        masked_line.push_back(line[i]);
                    }
                }

                //Write masked read to new file *abundance times
                for(unsigned int i = 1; i<= abundance; i++){
                    (*output).append(">K-mer "+ to_string(kmer_id) + " Abundance:" + 
                    to_string(abundance)+ " Copy:" + to_string(i) + "\n");
                    (*output).append(masked_line +"\n");
                }
                kmer_id++;
            }
        }
        input_stream.close();
    }

    void checkMask(string mask){
        
        if(mask == ""){
            throw std::runtime_error("mask empty!");
        }
        size_t raw_mer_size = 0;
        size_t masked_mer_size = 0;
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

        std::cerr << "Used Mask: " << mask << std::endl;
        std::cerr << "Transform " << raw_mer_size << "-mers into " << masked_mer_size << "-mers" << std::endl;
    }

    enum Mode {w,r};

    void checkFile(string path,Mode mode){

        if(mode==r){
            std::ifstream stream(path);
            if(!stream.good()){
                throw std::runtime_error(path + " could not be opened!");
            }
            stream.close(); 
        } else{
            std::ofstream stream(path);
            if(!stream.good()){
                throw std::runtime_error(path + " could not be opened!");
            }
            stream.close();
        }       
    }

    int main(int argc, char* argv[])
    {
        //Initialise command line parser
        CommandLineParser parser;
        parser.addAuthorInformation("MaskJelly,  author: Hartmut Haentze \nMaskJelly applies a mask of 'care' and 'dont care' positions to a set of k-mers read by Jellyfish");
        
        parser.addRequiredArgument('m',"Mask consisting of care (1) and don't care (0) positions. Needs to be saved as file.");
        parser.addRequiredArgument('i', "Input file");
        parser.addOptionalArgument('o',"","Output file");
        parser.addOptionalArgument('t',"1","number of threads to be used");

        //Parse command line
        parser.parse(argc,argv);
        string mask_file = parser.getArgument('m');
        string k_mer_file = parser.getArgument('i');
        string output_file = parser.getArgument('o');
        int nr_cores = stoi(parser.getArgument('t'));

        //test files
        checkFile(k_mer_file,r);
        checkFile(mask_file,r);
        if(output_file!="") checkFile(output_file,w);

        //read mask
        std::string mask = "";
        std::ifstream maskstream(mask_file);
        std::getline(maskstream,mask);
        maskstream.close();
        checkMask(mask);
       
        //Read k-mers file
        std::string line;
        std::ifstream infile(k_mer_file);

        //get number of lines in file
        int nr_lines = 0;
        while(std::getline(infile,line)){
            nr_lines++;
        }
        std::cerr << "Processing "<< nr_lines/2 <<" different k-mers" << std::endl;
        infile.close();

        //Initialize threads
        std::vector<std::thread> thread_vector;
        string partial_output[nr_cores];
        int step_size = ((nr_lines / 2)/nr_cores)*2;
        
        //execute threads
        for(size_t i = 0; i< nr_cores; i++){
            int f_nr_lines = i==nr_cores-1 ? nr_lines-(step_size*i) : step_size;
            thread_vector.push_back(
                thread(apply_mask,k_mer_file,step_size*i,f_nr_lines,mask,&(partial_output[i]))
              );
        }

        //wait for threads to finish
        for(auto &t : thread_vector){
            t.join();
        }
        
        //write output
        string output("");
        for(size_t i = 0; i< nr_cores; i++){
            output.append(partial_output[i]);
        }
        
        if(output_file!=""){
            std::ofstream outfile(output_file);
            outfile << output << endl;
            outfile.close();
        }
        else{
            cout << output <<endl;
        }

        return 0;
    }