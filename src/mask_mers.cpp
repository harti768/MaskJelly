    #include <iostream>
    #include <fstream>
    #include <thread>
    #include <vector>
    #include <mutex>
    #include "command_line_parser.hpp"

    using namespace std;
    
    struct deleteHelper{
        unsigned int kmers = 0;
        unsigned int abundance = 0;
        
    };

    //global variables
    bool print_names = false;
    bool write_cout = true;
    mutex output_mutex;
    ofstream outfile;

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

    void writeFile(string output){
        {
            lock_guard<mutex> lock_output(output_mutex);
            if(write_cout)  cout << output;
            else    outfile << output;   
        }
    }

    void apply_mask(string file, unsigned long long int start_line, unsigned long long int n_lines, unsigned int limit, string mask, deleteHelper* dHelper){

        //validate
        if(start_line%2 != 0 || n_lines%2!=0){
            throw runtime_error("Thread cannot read odd number of lines ("
            +to_string(n_lines)+")or start at odd position ("+to_string(start_line)+")");
        }

        //reach start line
        ifstream input_stream(file);
        std::string line;
        for(unsigned long long int i=0; i<start_line;i++){
            getline(input_stream, line);
        }

        //variables
        unsigned int abundance = 0;
        unsigned long long int kmer_id = start_line/2 + 1;
        size_t raw_mer_size = mask.size();
        string output("");

        for(unsigned long long int i=0; i<n_lines;i++){
            getline(input_stream, line);

            if(line[0] == '>') {
                 //Read abundance
                abundance = std::stoi(line.substr(1));

                //Ignore k-mer if abundance higher than limit
                if(abundance>limit){
                    dHelper->kmers++;
                    dHelper->abundance += abundance;

                    getline(input_stream, line);
                    i++;
                    continue;
                }
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
                    if(print_names){
                        output.append(">K:"+ to_string(kmer_id) + " A:" + 
                        to_string(i)+ "/" + to_string(abundance) + "\n");
                    } else{
                        output.append(">\n");
                    }
                    output.append(masked_line +"\n");
                }
                kmer_id++;
            }

            //Find stepsize with best trade-of for multi threading
            if(i%1000000==0){
                writeFile(output);
                output.clear();
            }
        }
        writeFile(output);
        input_stream.close();
    }

    int main(int argc, char* argv[])
    {
        //Initialise command line parser
        CommandLineParser parser;
        parser.addAuthorInformation("MaskJelly,  author: Hartmut Haentze \nMaskJelly applies a mask of 'care' and 'dont care' positions to a set of k-mers read by Jellyfish");
        
        parser.addRequiredArgument('m',"Mask consisting of care (1) and don't care (0) positions. Needs to be saved as file");
        parser.addRequiredArgument('i', "Input file");
        parser.addOptionalArgument('o',"","Output file");
        parser.addOptionalArgument('t',"1","Number of threads to be used");
        parser.addOptionalArgument('l',"1000000","Limit highest abundance of k-mers. Higher values will be neglected");
        parser.addFlagArgument('n',"Print names/ids of k-mers in output file. Increases file size");

        //Parse command line
        parser.parse(argc,argv);
        string mask_file = parser.getArgument('m');
        string k_mer_file = parser.getArgument('i');
        string output_file = parser.getArgument('o');
        int nr_cores = stoi(parser.getArgument('t'));
        int limit = stoi(parser.getArgument('l'));
        print_names = parser.getFlag('n');

        //test files
        checkFile(k_mer_file,r);
        checkFile(mask_file,r);
        if(output_file!=""){
            checkFile(output_file,w);
            outfile = ofstream(output_file);
            write_cout = false;
        }

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
        unsigned long long int nr_lines = 0;
        while(std::getline(infile,line)){
            nr_lines++;
        }
        std::cerr << "Processing "<< nr_lines/2 <<" different k-mers" << std::endl;
        infile.close();

        //Initialize threads
        std::vector<std::thread> thread_vector;
        deleteHelper dHelpers[nr_cores];
        unsigned long long int step_size = ((nr_lines / 2)/nr_cores)*2;
        
        //execute threads
        for(int i = 0; i< nr_cores; i++){
            unsigned long long int f_nr_lines = i==nr_cores-1 ? nr_lines-(step_size*i) : step_size;
            deleteHelper dHelper;
            thread_vector.push_back(
                thread(apply_mask,k_mer_file,step_size*i,f_nr_lines,limit,mask,&(dHelpers[i]))
              );
        }

        //wait for threads to finish
        for(auto &t : thread_vector){
            t.join();
        }
        
        //check how many k-mers were deleted
        unsigned int d_kmers = 0;
        unsigned long long int d_abundance = 0;
        for(auto& dHelper : dHelpers){
            d_kmers += dHelper.kmers;
            d_abundance += dHelper.abundance;
        }

        cerr << d_kmers << " different k-mers were deleted due to having an abundance higher than " << limit
        << ". Total abundance of deleted k-mers: " << d_abundance << endl;    

        return 0;
    }