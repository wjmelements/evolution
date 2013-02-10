#include "include/Environment.h"
#include "string.h"
#include <iostream>
#include <fstream>
int main()
{
    std::ofstream out("hn10b.txt");
    out<<"Trial\tSave\tTotal\tw/Chlorophyll A\tw/Chlorophyll B\tw/ATP Synthase\tw/Topoisomerase\tw/Primase\tw/Ligase\tw/Helicase\tw/Polymerase\tw/CEP164"<<std::endl;
    for(unsigned int trial = 1;trial <= 250;trial++)
        for(unsigned int save = 1;save <= 10;save++)
        {
            std::string filename("data/hn10b/cyano");///
            char* trialStr = new char[11];
            sprintf(trialStr,"%u",trial);
            filename+=trialStr;
            filename+='-';
            char* saveStr = new char[11];
            sprintf(saveStr,"%u",save);
            filename+=saveStr;
            filename+=".dat";
            delete[]trialStr;
            delete[]saveStr;
            Environment instance(filename);
            out<<trial<<"\t"<<save;
            instance.printMicroscopeLine(out);
        }
    out.close();
    return 0;
}
