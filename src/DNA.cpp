#include "../include/DNA.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>

DNA::DNA(Data myData, DNA* previous, DNA* following)
{
    data=myData;
    prev=previous;
    next=following;
}
DNA::DNA()
{
    data.codons=0;
    data.numValidCodons=0;
    prev=0;
    next=0;
}

DNA::~DNA()
{
    if(next)
        delete next;
}

DNA* DNA::add(Data newData)
{
    if(next)
        return next->add(newData);
    else
    {
        next = new DNA(newData,this);
        return next;
    }
}

DNA* DNA::add(DNA* newDNA)
{
    if(next)
    {
        return next->add(newDNA);
    }
    else
    {
        next=newDNA;
        newDNA->prev=this;
        return newDNA;
    }
}

bool DNA::hasNext()
{
    return next;
}

bool DNA::isEndOfChromosome()
{
    return !(next && data.numValidCodons>10);
}

DNA* DNA::getNext()
{
    return next;
}

DNA* DNA::last()
{
    if(next)
        return next->last();
    return prev->next;
}

DNA* DNA::lastInChromosome()
{
    if(isEndOfChromosome())
        return this;
    return next->lastInChromosome();
}

DNA* DNA::nextChromosome()
{
    if(isEndOfChromosome())
        return next;
    return next->nextChromosome();
}

unsigned int DNA::getHaploidNumber(bool startAtBeginning)
{
    if(startAtBeginning&&prev)
    {
        return prev->getHaploidNumber(true);
    }
    if(next)
        if(data.numValidCodons < 11)
            return 1+(next->getHaploidNumber(false));
        else
            return next->getHaploidNumber(false);
    else
        return 1;
}

unsigned long long DNA::getCodonCount()
{
    return(data.numValidCodons>9?10:data.numValidCodons)+(next?next->getCodonCount():0);
}

unsigned long long DNA::getCodonCountForChromosome()
{
    return(data.numValidCodons>9?10:data.numValidCodons)+(isEndOfChromosome()?0:next->getCodonCountForChromosome());
}

Data DNA::getData(int index)
{
    if(index>0)
        return next->getData(index-1);
    if(index<0)
        return prev->getData(index+1);
    return data;
}

void DNA::print(bool isRNA)
{
    using std::cout;
    using std::endl;
    Data temp = data;
    if(temp.numValidCodons>10)
        temp.numValidCodons=10;
    while(temp.numValidCodons)
    {
        for(int j=0;j<3;j++)
        {
            switch(temp.codons&3)
            {
                case 0:   cout<<'A';
                    break;
                case 1:   if(isRNA)
                                    cout<<'U';
                                else
                                    cout<<'T';
                    break;
                case 2:   cout<<'G';
                    break;
                case 3:   cout<<'C';
                    break;
                default:  cout<<(temp.codons&3);//debug
            }//print codon
            temp.codons >>= 2;

        }
        temp.numValidCodons--;
        cout<<' ';
    }
    if(data.numValidCodons<11||!next)
        cout<<'&';
    cout<<endl;
    if(temp.numValidCodons)
        cout<<"Error: Not all bits read"<<endl;//debug
}

void DNA::printAll(bool isRNA,bool startFromBeginning) {
    if(startFromBeginning&&prev)
        prev->printAll(isRNA,true);
    else
    {
        print(isRNA);
        if(next)
            next->printAll(isRNA,false);
    }
}

DNA& operator<< (DNA& dna, Data data)
{
    bool notEndOfChromosome=dna.data.numValidCodons>10;
    if(data.numValidCodons>10)
        data.numValidCodons=10;
    if(notEndOfChromosome)
        dna.data.numValidCodons=10;
    Data overflow;
    int numValidCodons=data.numValidCodons+dna.data.numValidCodons-10;
    if(numValidCodons<0)
    {
        overflow.numValidCodons=0;
        overflow.codons=0;
    }
    else
    {
        overflow.numValidCodons=numValidCodons;
        overflow.codons=data.codons>>6*(data.numValidCodons-overflow.numValidCodons);
    }
    dna.data.codons|=(data.codons<<(dna.data.numValidCodons*6));
    numValidCodons+=10;
    if(numValidCodons>10)
        dna.data.numValidCodons=11;
    else
        dna.data.numValidCodons=numValidCodons;
    if(overflow.numValidCodons)
    {
        if(dna.next&&notEndOfChromosome)
            overflow>>*dna.next;
        else
        {
            DNA*remainder=new DNA(overflow,&dna,dna.next);
            if(dna.next)
                dna.next->prev=remainder;
            dna.next=remainder;
        }
    }
    return dna;
}//DNA<<Data

DNA& operator>> (Data data,DNA& dna)
{
    bool notEndOfChromosome = dna.data.numValidCodons>10;
    if(data.numValidCodons>10)
        data.numValidCodons=10;
    if(notEndOfChromosome)
        dna.data.numValidCodons=10;
    Data overflow;
    int numValidCodons=data.numValidCodons+dna.data.numValidCodons-10;
    if(numValidCodons<0)
    {
        overflow.codons=0;
        overflow.numValidCodons=0;
    }
    else
    {
        overflow.numValidCodons=numValidCodons;
        overflow.codons=dna.data.codons>>(6*(dna.data.numValidCodons-overflow.numValidCodons));
    }
    dna.data.codons<<=(data.numValidCodons*6);
    dna.data.codons|=data.codons;
    numValidCodons+=10;
    if(numValidCodons>10)
        numValidCodons=11;
    dna.data.numValidCodons=numValidCodons;
    if(overflow.numValidCodons)
    {
        if(dna.next&&notEndOfChromosome)
            overflow>>*dna.next;
        else
        {

            DNA*remainder=new DNA(overflow,&dna,dna.next);
            if(dna.next)
                dna.next->prev=remainder;
            dna.next=remainder;
        }
    }
    return dna;
}//Data>>DNA

DNA* DNA::replicate(double accuracy,DNA* andAddTo)
{
    unsigned int dataCount = 0;//in nucleotides, resets every 30 to zero
    DNA* terminal = lastInChromosome();
    DNA* nextChrom = terminal->next;
    unsigned int terminalCount = terminal->getData().numValidCodons*3;
    DNA* currentParent = this;
    Data bogus;
    bogus.codons = 0;
    bogus.numValidCodons = 0;
    DNA* replicated = new DNA(bogus,andAddTo);
    DNA* currentChild = replicated;
    unsigned int replicatedCount = 0;//in nucleotides, resets every 30 to zero as flag for next in linked list
    Data replicatedData=bogus;
    while(currentParent != nextChrom && (currentParent!=terminal || dataCount<terminalCount))
    {
        if((double)rand()/RAND_MAX <= accuracy)
        {
            replicatedData.codons |= (((currentParent->getData().codons>>(2*dataCount))&3)<<(2*replicatedCount));
            replicatedCount++;
            if(replicatedCount==30)
            {
                replicatedData.numValidCodons = 11;
                currentChild->data = replicatedData;
                replicatedData = bogus;
                currentChild->next = new DNA(replicatedData,currentChild);
                currentChild = currentChild->next;
                replicatedCount = 0;
            }
            dataCount++;
            if(dataCount==30)
            {
                currentParent = currentParent->next;
                dataCount=0;
            }
        }//perfect nucleotide replication
        else
            switch(rand()&15)
            {
                case 0: //Deletion
                        dataCount++;
                        if(dataCount==30)
                        {
                            currentParent = currentParent->next;
                            dataCount=0;
                        }
                        break;
                case 1: //Substitution
                case 2: //Substitution
                case 3: //Substitution
                case 4: //Substitution
                case 5: //Substitution
                case 6: //Substitution
                case 7: //Substitution
                case 8: //Substitution
                case 9: //Substitution
                case 10://Substitution
                case 11://Substitution
                case 12://Substitution
                case 13://Substitution
                case 14://Substitution
                        dataCount++;
                        if(dataCount==30)
                        {
                            currentParent = currentParent->next;
                            dataCount=0;
                        }
                case 15://Insertion
                        replicatedData.codons |= (((unsigned long long)(rand() & 3)) << (replicatedCount*2));
                        replicatedCount++;
                        if(replicatedCount==30)
                        {
                            replicatedData.numValidCodons = 11;
                            currentChild->data = replicatedData;
                            replicatedData = bogus;
                            currentChild->next = new DNA(replicatedData,currentChild);
                            currentChild = currentChild->next;
                            replicatedCount = 0;
                        }
                        break;
            }//replication error
    }//replicate chromosome
    replicatedData.numValidCodons = replicatedCount/3;
    if(replicatedData.numValidCodons == 0 && andAddTo && currentChild->prev == andAddTo)
        return nextChrom?nextChrom->replicate(accuracy,currentChild):0;//no more chromosome :-(
    if(replicatedData.numValidCodons == 0 && currentChild->prev != andAddTo)
    {
        currentChild = currentChild->prev;
        if(currentChild->data.numValidCodons==11)
            currentChild->data.numValidCodons=10;
    }
    else
        currentChild->data = replicatedData;
    currentChild->next = nextChrom?nextChrom->replicate(accuracy,currentChild):0;
    return replicated;
}
