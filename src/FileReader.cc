#include "FileReader.hh"
#include <iostream>
//#include <fstream>
//#include <map>
#include <sstream>
//#include <string>
#include <algorithm>

//register a new parameter with name key and initial int value
void FileReader::registerIntParameter(const std::string &key, int init)
{
   integers_[key]=init;
}

//register a new parameter with name key and initial double value
void FileReader::registerRealParameter(const std::string &key, real init)
{
   reals_[key]=init;
}

//register a new parameter with name key and initial string value
void FileReader::registerStringParameter(const std::string &key, const std::string &init)
{
    strings_[key]=init;
}

//set a value for the key string with value in
void FileReader::setParameter(const std::string &key, const std::string &in)
{
   //std::map <std::string, std::string>::iterator it = strings_.find(key);
   //integers_[key]=in;
   // (*it).second = in;
   strings_[key]=in;
}

//set a value for the key string with value in
void FileReader::setParameter(const std::string &key, real in)
{
   reals_[key]=in;
}

//set a value for the key string with value in
void FileReader::setParameter(const std::string &key, int in)
{
   integers_[key]=in;
   //strings_[key]=in;
}

//try to read all registered parameters from file name
bool FileReader::readFile(const std::string &name)
{
    std::ifstream reader;
    std::string temp;
    std::string pure_text;
    std::string key;
    bool found;
    reader.open(name);

    while(std::getline(reader, temp)){
    found=false;

    temp.erase(std::find( temp.begin(), temp.end(),'#'), std::find( temp.begin(), temp.end(),'\n'));
    std::stringstream ss(temp);
    ss >> key;

    pure_text = temp;
    pure_text.erase(remove_if(pure_text.begin(), pure_text.end(), isspace), pure_text.end());

    if(pure_text.empty()) found=true;

    if (found==false)
        for(std::map<std::string, int>::iterator a = integers_.begin(); a!=integers_.end(); a++)
            if((*a).first == key && found==false){
                ss >> (*a).second;
                found=true;
                 std::cout <<(*a).first<<"\t =  "<< (*a).second<<"\n";
            }

    if (found==false)
        for(std::map<std::string, real>::iterator a = reals_.begin(); a!=reals_.end(); a++)
            if((*a).first == key && found==false){
                ss >> (*a).second;
                found=true;
                std::cout <<(*a).first<<"\t =  "<< (*a).second<<"\n";
            }

    if (found==false)
        for(std::map<std::string, std::string>::iterator a = strings_.begin(); a!=strings_.end(); a++)
            if((*a).first == key && found==false){
                ss >> (*a).second;
                found=true;
                std::cout <<(*a).first<<"\t =  "<< (*a).second<<"\n";
            }

    if (found==false)   std::cout<<"parameter "<< key <<" has not registered!\n";
    }

   return true;
}

//print out all parameters to std:out
void FileReader::printParameters() const
{
   for(std::map<std::string, int>::const_iterator it=integers_.begin(); it!=integers_.end(); it++ )
        std::cout << (*it).first << '\t' << (*it).second << '\n';

   for(std::map <std::string, real>::const_iterator it=reals_.begin(); it!=reals_.end(); it++ )
        std::cout << (*it).first << '\t' << (*it).second << '\n';

   for(std::map < std::string, std::string>::const_iterator it=strings_.begin(); it!=strings_.end(); it++ )
        std::cout << (*it).first << '\t' << (*it).second << '\n';
}


