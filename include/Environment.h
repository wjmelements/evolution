#ifndef ENVIRONMENT
#define ENVIRONMENT
#include "Organism.h"
#include <fstream>
#include "string.h"

struct Location
{
    Organism* organism;

    unsigned int sunlight;
    unsigned long long protonGradient;
    unsigned long long glucose;
    unsigned long long glucose6Phosphate;
    unsigned long long fructose16Bisphosphate;
    unsigned long long G3P;
    unsigned long long DHAP;
    unsigned long long one3Bisphosphoglycerate;
    unsigned long long NADH;
    unsigned long long NADPH;
    unsigned long long threePhosphoglycerate;
    unsigned long long pyruvate;
    unsigned long long alphaKetoglutarate;
    unsigned long long succinylCoA;
    unsigned long long succinate;
    unsigned long long fumerate;
    unsigned long long FADH2;
};//Local

struct Atmosphere
{
    unsigned long long oxygen;
};//Ambient

class Environment
{
    public:
        Environment(std::string& filename);
        Environment(const unsigned int rows,const unsigned int columns,unsigned long long generationNumber);
        ~Environment();
        unsigned int getRows();
        unsigned int getCols();
        Location getLocation(unsigned int row,unsigned int col);
        Organism* getOrganism(unsigned int row,unsigned int col);
        bool insert(Organism* organism,unsigned int row,unsigned int col);//returns false if fail
        bool place(Organism* organism,unsigned int row,unsigned int col);//places *organism in Location adjacent to row,col; returns false if fail
        Organism* removeOrganism(unsigned int row,unsigned int col);//must call delete if you don't use the organism
        void setLighting(unsigned int row,unsigned int col,unsigned int intensity,long double decay);///consider reimplementing as circle
        unsigned long long countOrganisms();
        Organism* reproduce(unsigned int row,unsigned int col,long double turn);
        void stir();//moves all organisms randomly

        void simulateGeneration();
        void simulateUntilGeneration(unsigned long long number);
        void simulateGenerations(unsigned long long number);
        //IO
        void print();
        void saveToFile(std::string& filename);
        void printAnalysis();
        void printMicroscopeLine(std::ofstream& out);
    protected:
    private:
        unsigned int rows;
        unsigned int cols;//columns
        unsigned long long numLocations;//calculated
        Location** grid;
        Atmosphere atmosphere;
        unsigned long long generationNumber;
        void naturalSelection();
};

#endif // ENVIRONMENT
