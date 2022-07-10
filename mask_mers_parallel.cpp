#include <iostream>
#include <fstream>


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
    std::string masked_line = "";

    std::ifstream infile(k_mer_file);
    if(!infile.good()){
        throw std::runtime_error("infile could not be opened!");
    }
    std::ofstream outfile(output_file);
    if(!outfile.good()){
        throw std::runtime_error("outfile could not be opened!");
    }

    unsigned int abundance = 0;
    unsigned int kmer_id = 1;

    while(std::getline(infile, line)){
      if(line.size()==0) continue;
      //std::cout << "Size: " << line.size() << std::endl;
        //header
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
            masked_line="";
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
                outfile <<">K-mer "<< kmer_id << " Abundance:" << abundance << " Copy:" << i << std::endl;
                outfile << masked_line << std::endl;
            }
            kmer_id++;
        }
        //std::cout << "0:" << line[0] << " 1:" << line[1] << " 2:" << line[2] << std::endl;
        //std::cout << "1:" << line[1] << " 2:" << line[2] << " 3:" << line[3] << std::endl;
    }
    infile.close();
    outfile.close();
    return 0;
}