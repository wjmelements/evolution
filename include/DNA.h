#ifndef DNA_H
#define DNA_H
struct Data
{
    unsigned long long codons : 60;
    unsigned long long numValidCodons : 4;   //0 to 15; if < 11, then end of chromosome; else, not end of chromosome
};

class DNA
{
    public:
        DNA();
        DNA(Data data,DNA* prev=0,DNA* next=0);
        virtual ~DNA();
        //List
        DNA* add(Data data);//returns last strand
        DNA* add(DNA* newDNA);//can be single unit or beginning of new strand; returns last strand
        bool hasNext();
        unsigned long long getCodonCount();
        unsigned long long getCodonCountForChromosome();
        unsigned int getHaploidNumber(bool startAtBeginning=true);
        DNA* last();
        DNA* lastInChromosome();
        DNA* nextChromosome();
        DNA* getNext();
        DNA* replicate(double accuracy = 1,DNA* andAddTo = 0);//must be deleted if not used or will cause memory leak
        //Individual
        bool isEndOfChromosome();
        Data getData(int index=0);
        friend DNA& operator<<(DNA& dna, Data data);//adds to end of referenced DNA
        friend DNA& operator>>(Data data,DNA& dna);
        //IO
        void print(bool isRNA=false);
        void printAll(bool isRNA=false,bool startFromBeginning=true);

        friend class Organism;
    protected:
    private:
        DNA*next;
        DNA*prev;
        Data data;
};

#endif // DNA_H
