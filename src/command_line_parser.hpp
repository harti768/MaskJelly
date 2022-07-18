#ifndef COMMAND_LINE_PARSER_HPP
#define COMMAND_LINE_PARSER_HPP

#include <string>
#include <map>

class CommandLineParser{

struct Argument {
    std::string description;
    std::string value;
    bool required;
};

struct Flag {
    std::string description;
    bool value;
};

public:
    CommandLineParser();
    void addRequiredArgument(char name, std::string description);
    void addOptionalArgument(char name, std::string default_value, std::string description);
    void addFlagArgument(char name, std::string description);
    void parse(int argc, char* argv[]);
    void info();
    std::string getArgument(char c);
    bool getFlag(char c);
    void addAuthorInformation(std::string info);

private:
    std::map<char,Argument> arguments;
    std::map<char,Flag> flags;
    std::string author_info;

    void checkUnique(char name);    
};


#endif //COMMAND_LINE_PARSER_HPP