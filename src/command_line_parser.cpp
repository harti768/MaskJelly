    #include "command_line_parser.hpp"
    #include <iostream>

    using namespace std;

    CommandLineParser::CommandLineParser(){};

    void CommandLineParser::addAuthorInformation(string info){
        this->author_info = info;
    }

    void CommandLineParser::addRequiredArgument(char name, std::string description){

        checkUnique(name);
        Argument argument;
        argument.description = description;
        argument.required = true;
        this->arguments[name]=argument;
    }

    void CommandLineParser::addOptionalArgument(char name, std::string default_value, std::string description){

        checkUnique(name);
        Argument argument;
        argument.value = default_value;
        argument.description = description;
        argument.required = false;
        this->arguments[name]=argument;
    }

    void CommandLineParser::addFlagArgument(char name, std::string description){
        
        checkUnique(name);
        Flag flag;
        flag.value = false;
        flag.description = description;
        this->flags[name]= flag;
    }

    void CommandLineParser::checkUnique(char name){

        if(this->arguments.find(name) != this->arguments.end()
        || this->flags.find(name) != this->flags.end()){

            string s(1,name);
            throw runtime_error("Argument/Flag -"+s+" already used!");
        }
    }

    void CommandLineParser::parse(int argc, char* argv[]){
        
        //parse command line
        for(int i=1; i<argc; i++){
            char * arg = argv[i];
            if(arg[0]=='-'){
                char name = arg[1];

                //Case flag
                if(this->flags.find(name) != this->flags.end()){
                    this->flags[name].value = true;
                }

                //Case argument
                else if(this->arguments.find(name) != this->arguments.end()){
                    if(argc > i+1 && argv[i+1][0] != '-'){
                        this->arguments[name].value = argv[i+1];
                        i++;
                    } else{
                        string s(1,name);
                        throw runtime_error("No value for argument -"+s+" found");
                    }
                } 

                //default
                else {
                    info();
                    string s(1,name);
                    throw runtime_error ("Argument -"+s+" unkown");
                }
            }
        }

        //check whether all required arguments were used
        for(auto arg : this->arguments){
            if(arg.second.required && arg.second.value == ""){
                string s(1,arg.first);
                throw runtime_error("Mandatory argument -"+s+" not specified"); 
            }
        }
    }

    string CommandLineParser::getArgument(char name){
        
        if(this->arguments.find(name)==this->arguments.end()){
            string s(1,name);
            throw runtime_error("Argument -"+s+" not known");
        }
        return this->arguments.at(name).value;
    }

    bool CommandLineParser::getFlag(char name){

        if(this->flags.find(name)==this->flags.end()){
            string s(1,name);
            throw runtime_error("Argument -"+s+" not known");
        }
        return this->flags.at(name).value;
    }

    void CommandLineParser::info(){

        cerr << "\n" << this->author_info << "\n" << endl;
        cerr << "Required arguments:" << endl;
        for(auto arg : this->arguments){
            if(arg.second.required){
                cerr << "-" << arg.first << ":\t" << arg.second.description << endl;
            }
        }

        cerr << "Optional arguments:" << endl;
        for(auto arg : this->arguments){
            if(!arg.second.required){
                cerr << "-" << arg.first << ":\t" << arg.second.description <<  endl;
            }
        }

        if(!this->flags.empty()){
            cerr << "Flags:" << endl;
            for(auto flag : this->flags){
                cerr << "-" << flag.first << ":\t" << flag.second.description << endl;
            }
        }
    }