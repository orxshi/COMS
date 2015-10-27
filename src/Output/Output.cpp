#include "Output.h"

string createOutputDir()
{
    string folder = "../out/";
    
    struct stat sb;
    if ( !(stat(folder.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)) )
    {
        mkdir (folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }
    
    string outputDir;
    time_t t = time(0);
    struct tm* now = localtime(&t);
    
    string month = std::to_string(now->tm_mon + 1);
    string day = std::to_string(now->tm_mday);
    string year = std::to_string(now->tm_year + 1900);
    string hour = std::to_string(now->tm_hour);
    string min = std::to_string(now->tm_min);
    
    outputDir = folder;
    string us = "_";    
    outputDir.append (month);
    outputDir.append (us);
    outputDir.append (day);
    outputDir.append (us);
    outputDir.append (year);
    outputDir.append (us);
    outputDir.append (hour);
    outputDir.append (us);
    outputDir.append (min);
    
    mkdir (outputDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    
    return outputDir;
}