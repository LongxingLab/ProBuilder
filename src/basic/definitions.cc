#include "basic/global.hh"

namespace basic {
    
std::map<std::string, std::string> BASIC_AA3to1 = { 
        {"ALA","A"}, {"CYS","C"}, {"ASP","D"}, {"GLU","E"}, {"PHE","F"}, {"GLY","G"},
        {"HIS","H"}, {"ILE","I"}, {"LYS","K"}, {"LEU","L"}, {"MET","M"}, {"ASN","N"}, {"PRO","P"},
        {"GLN","Q"}, {"ARG","R"}, {"SER","S"}, {"THR","T"}, {"VAL","V"}, {"TRP","W"}, {"TYR","Y"},
        {"AFT", "X"},{"BLA","X"}
    };

std::map<char, std::string> BASIC_AA1to3 = { 
        { 'A',"ALA" }, {'C',"CYS"}, {'D',"ASP"}, {'E',"GLU"}, {'F',"PHE"}, {'G',"GLY"},
        {'H',"HIS"}, {'I',"ILE"}, {'K',"LYS"}, {'L',"LEU"}, {'M',"MET"}, {'N',"ASN"}, {'P',"PRO"},
        {'Q',"GLN"}, {'R',"ARG"}, {'S',"SER"}, {'T',"THR"}, {'V',"VAL"}, {'W',"TRP"}, {'Y',"TYR"},
        {'X', "BLA"},{'X',"AFT"}
    };
}
