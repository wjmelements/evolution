#include <iostream>
#include <time.h>
#include <fstream>
#include "stdlib.h"
#include "include/Environment.h"
#include "string.h"

Data randomData(bool endOfChromosome=true)
{
    Data temp;
    temp.codons=rand();
    temp.codons-=temp.codons&15;//removes last 4 bits
    temp.codons<<=28;
    temp.codons|=rand();
    if(endOfChromosome)
    {
        do
            temp.numValidCodons = rand()&15;
            while(temp.numValidCodons > 10 || temp.numValidCodons == 0);
        unsigned long long bitMask=0;
        for(int i=0;i<temp.numValidCodons;i++)
        {
            bitMask<<=6;
            bitMask^=63;
        }
        temp.codons&=bitMask;
    }
    else
        temp.numValidCodons = 11;
    return temp;
}

int main()
{
    //init
    srand(time(NULL));
    for(unsigned int j = 1;j <= 500;j++)//Trials
    {
        /*
        Environment one(9,9,0);
        one.setLighting(4,4,4,.75);
        DNA* original = new DNA(randomData(false));
        for(unsigned int i = 0;i<15000;i++)
            *original<<randomData(false);
        *original<<randomData(true);
        */

        Data protochlorophyllideReductase;
        protochlorophyllideReductase.codons = 0x53FDDAF31ULL;  //Met Trp Pro Glu Gly Stop \Bla
        protochlorophyllideReductase.numValidCodons = 11;
        Data atpSynthase;
        atpSynthase.codons = 0x1C5C31C6571ULL; //Met Lys His Asp Val Glu Stop \Bla
        atpSynthase.numValidCodons = 11;
        DNA* original = new DNA(randomData(false));
        DNA* last = original;
        for(unsigned int i = 1;i<9;i++)
        {
            last = last->add(randomData(false));
        }
        last = last->add(protochlorophyllideReductase);
        for(unsigned int i = 10;i<1499;i++)
        {
            last = last->add(randomData(false));
        }
        last = last->add(randomData(true));
        for(unsigned int i = 0;i<9;i++)
        {
            last = last->add(randomData(false));
        }
        last = last->add(atpSynthase);
        for(unsigned int i = 10;i<1499;i++)
        {
            last = last->add(randomData(false));
        }
        last = last->add(randomData(true));
        for(unsigned int i = 2;i<10;i++)
        {
            for(unsigned int j = 0;j<1499;j++)
            {
                last=last->add(randomData(false));
            }
            last = last->add(randomData(true));
        }
        Organism* first = new Organism(0,original);
        Environment one(9,9,0);
        one.setLighting(4,4,11,.7);
        if(!one.insert(first,4,4))
        {
            delete first;
            std::cout<<"Initiation error"<<std::endl;
        }
        if(original->getHaploidNumber() != 10) //If haploid number isn't what is intended then it is an error
        {
            first->printDNA();
            exit(original->getHaploidNumber());
        }
        for(unsigned int i = 1; i <= 20; i++) //saves per trial
        {
            for(unsigned int generation = 0; generation<500; generation++) //generations per save
            {
                one.simulateGeneration();
                std::cout<<generation+500*i<<std::endl;
                one.printAnalysis();
                one.stir();///consider reimplementing similar to deterioration
                ///consider adding log
            }
            //output
            std::string filename("data/acc10a/cyano");
            char* jchar= new char[11];
            sprintf(jchar,"%u",j);
            filename+=jchar;
            filename+='-';
            char* ichar = new char[11];
            sprintf(ichar,"%u",i);
            filename+=ichar;
            filename+=".dat";
            one.saveToFile(filename);
            delete[] jchar;
            delete[] ichar;
        }

    }
    return 0;
}
