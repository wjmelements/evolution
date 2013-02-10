#include "../include/Organism.h"
#include <iostream>
#include <cstdlib>
#include <fstream>
/*
Organism::Organism()
{
    prev=0;
    next=0;
    dna=0;
    rna=0;
    traits.cost=0;
    traits.isNotSequenced=1;
    traits.hasHexokinase=0;
    traits.hasPhosphofructokinase=0;
    traits.hasAldolase=0;
    traits.hasTriosephosphateIsomerase=0;
    traits.hasTriosephosphateDehydrogenase=0;
    traits.hasPhosphoglycerateKinase=0;
    traits.hasPyruvateKinase=0;
    traits.hasIsocitrateDehydrogenase=0;
    traits.hasAlphaKetoglutarateDehydrogenase=0;
    traits.hasSuccinylCoASynthetase=0;
    traits.hasSdhA=0;
    traits.hasMalateDehydrogenase=0;
    traits.hasATPSynthase=0;
    traits.hasComplexI=0;
    traits.hasChain1=0;
    traits.hasChain2=0;
    traits.hasChain3=0;
    traits.hasChain4=0;
    traits.hasChain4L=0;
    traits.hasChain5=0;
    traits.hasChain6=0;
    traits.hasSdhB=0;
    traits.hasSdhC=0;
    traits.hasSdhD=0;
    traits.hasProtochlorophyllideReductase=0;
    traits.hasChlorophyllBSynthase=0;
    traits.hasRubisco=0;
    traits.hasRubiscoSmall=0;
    traits.hasRubiscoLarge=0;
    traits.hasPEPCarboxylase=0;
    traits.hasCEP164=0;
    traits.hasDNALigase=0;
    traits.hasDNAPolymerase=0;
    traits.hasDNAPrimase=0;
    traits.hasHelicase=0;
    traits.hasTopoisomerase=0;
    traits.accuracy=0.5;
}

Organism::Organism(DNA* data,Organism* previous,Organism* following)
{
    prev=previous;
    next=following;
    dna=data;
    rna=0;
    traits.cost=0;
    traits.isNotSequenced=1;
    traits.hasHexokinase=0;
    traits.hasPhosphofructokinase=0;
    traits.hasAldolase=0;
    traits.hasTriosephosphateIsomerase=0;
    traits.hasTriosephosphateDehydrogenase=0;
    traits.hasPhosphoglycerateKinase=0;
    traits.hasPyruvateKinase=0;
    traits.hasIsocitrateDehydrogenase=0;
    traits.hasAlphaKetoglutarateDehydrogenase=0;
    traits.hasSuccinylCoASynthetase=0;
    traits.hasSdhA=0;
    traits.hasMalateDehydrogenase=0;
    traits.hasATPSynthase=0;
    traits.hasComplexI=0;
    traits.hasChain1=0;
    traits.hasChain2=0;
    traits.hasChain3=0;
    traits.hasChain4=0;
    traits.hasChain4L=0;
    traits.hasChain5=0;
    traits.hasChain6=0;
    traits.hasSdhB=0;
    traits.hasSdhC=0;
    traits.hasSdhD=0;
    traits.hasProtochlorophyllideReductase=0;
    traits.hasChlorophyllBSynthase=0;
    traits.hasRubisco=0;
    traits.hasRubiscoSmall=0;
    traits.hasRubiscoLarge=0;
    traits.hasPEPCarboxylase=0;
    traits.hasCEP164=0;
    traits.hasDNALigase=0;
    traits.hasDNAPolymerase=0;
    traits.hasDNAPrimase=0;
    traits.hasHelicase=0;
    traits.hasTopoisomerase=0;
    traits.accuracy=0.5;
}
*/

Organism::Organism(long double turn, DNA* data)
{
    dna=data;
    rna=0;
    traits.cost=0;
    traits.isNotSequenced=1;
    traits.hasHexokinase=0;
    traits.hasPhosphofructokinase=0;
    traits.hasAldolase=0;
    traits.hasTriosephosphateIsomerase=0;
    traits.hasTriosephosphateDehydrogenase=0;
    traits.hasPhosphoglycerateKinase=0;
    traits.hasPyruvateKinase=0;
    traits.hasIsocitrateDehydrogenase=0;
    traits.hasAlphaKetoglutarateDehydrogenase=0;
    traits.hasSuccinylCoASynthetase=0;
    traits.hasSdhA=0;
    traits.hasMalateDehydrogenase=0;
    traits.hasATPSynthase=0;
    traits.hasComplexI=0;
    traits.hasChain1=0;
    traits.hasChain2=0;
    traits.hasChain3=0;
    traits.hasChain4=0;
    traits.hasChain4L=0;
    traits.hasChain5=0;
    traits.hasChain6=0;
    traits.hasSdhB=0;
    traits.hasSdhC=0;
    traits.hasSdhD=0;
    traits.hasProtochlorophyllideReductase=0;
    traits.hasChlorophyllBSynthase=0;
    traits.hasOEC=0;
    traits.hasOECEP1=0;
    traits.hasOECEP2=0;
    traits.hasOECEP3=0;
    traits.hasRubisco=0;
    traits.hasRubiscoSmall=0;
    traits.hasRubiscoLarge=0;
    traits.hasPEPCarboxylase=0;
    traits.hasCEP164=0;
    traits.hasSensoryRhodopsinII=0;
    traits.hasDNALigase=0;
    traits.hasDNAPolymerase=0;
    traits.hasDNAPrimase=0;
    traits.hasHelicase=0;
    traits.hasTopoisomerase=0;
    traits.accuracy=0.99;
    traits.ATP=0;
    sequence();
    traits.turn = turn + 1e-6L * traits.cost;
}

Organism::Organism(std::ifstream& in)
{
    char input;
    in >> traits.turn;
    in >> traits.ATP;
    in >> input;
    dna = 0;
    DNA* last = 0;//for efficiency
    while(input == '^')
    {
        unsigned long long compressed;
        in >> compressed;
        Data data;
        data.codons = compressed >> 4;
        data.numValidCodons = compressed&15;
        if(dna)
        {
            if(last)
            {
                last = last->add(new DNA(data));
            }
            else
                std::cout<<"ERROR: last is "<<last<<" while dna is "<<dna<<std::endl;
        }
        else
        {
            dna = new DNA(data);
            last = dna;
        }
        in >> input;
    }
    rna=0;
    traits.cost=0;
    traits.isNotSequenced=1;
    traits.hasHexokinase=0;
    traits.hasPhosphofructokinase=0;
    traits.hasAldolase=0;
    traits.hasTriosephosphateIsomerase=0;
    traits.hasTriosephosphateDehydrogenase=0;
    traits.hasPhosphoglycerateKinase=0;
    traits.hasPyruvateKinase=0;
    traits.hasIsocitrateDehydrogenase=0;
    traits.hasAlphaKetoglutarateDehydrogenase=0;
    traits.hasSuccinylCoASynthetase=0;
    traits.hasSdhA=0;
    traits.hasMalateDehydrogenase=0;
    traits.hasATPSynthase=0;
    traits.hasComplexI=0;
    traits.hasChain1=0;
    traits.hasChain2=0;
    traits.hasChain3=0;
    traits.hasChain4=0;
    traits.hasChain4L=0;
    traits.hasChain5=0;
    traits.hasChain6=0;
    traits.hasSdhB=0;
    traits.hasSdhC=0;
    traits.hasSdhD=0;
    traits.hasProtochlorophyllideReductase=0;
    traits.hasChlorophyllBSynthase=0;
    traits.hasOEC=0;
    traits.hasOECEP1=0;
    traits.hasOECEP2=0;
    traits.hasOECEP3=0;
    traits.hasRubisco=0;
    traits.hasRubiscoSmall=0;
    traits.hasRubiscoLarge=0;
    traits.hasPEPCarboxylase=0;
    traits.hasCEP164=0;
    traits.hasSensoryRhodopsinII=0;
    traits.hasDNALigase=0;
    traits.hasDNAPolymerase=0;
    traits.hasDNAPrimase=0;
    traits.hasHelicase=0;
    traits.hasTopoisomerase=0;
    traits.accuracy=0.99;
    sequence();
}

void Organism::saveToFile(std::ofstream& out)
{
    using std::endl;
    DNA* current = dna;
    out << traits.turn << endl;
    out << traits.ATP <<endl;
    while(current)
    {
        out << '^'<<endl;
        unsigned long long data = current->data.numValidCodons | current->data.codons << 4;
        out << data << endl;
        current = current->next;
    }
    out << ')'<<endl;
}

Organism::~Organism()
{
    delete dna;
    delete rna;
}

/*
Organism* Organism::get(int index)
{
    if(index>0)
        return next->get(index-1);
    if(index<0)
        return prev->get(index+1);
    return this;
}
*/
/*
Organism Organism::remove(int index)
{
    if(index>0)
        return next->remove(index-1);
    if(index<0)
        return prev->remove(index+1);
    if(prev)
    {
        prev->next=next;
        if(next)
            next->prev=prev;
        next=0;
        prev=0;
        return *this;
    }//if not first organism
    else
    {
        if(next)
        {   //copies first organism, replaces first organism with copy of second, removes second, returns copy of first
            DNA* dnaToReturn = dna;
            if(traits.isNotSequenced)
                sequence();
            Traits traitsToReturn = traits;
            dna=next->dna;
            if(next->traits.isNotSequenced)
                next->sequence();
            traits=next->traits;
            next->remove();
            Organism toReturn(dnaToReturn);
            toReturn.traits = traitsToReturn;
            return toReturn;
        }
        else
        {
            std::cerr<<"Fatal error: cannot remove the only organism"<<std::endl;
            exit(-1);
        }//Where do we go now?
    }
}

unsigned long long Organism::count(bool countFromBeginning)
{
    if(countFromBeginning)
        if(prev)
            return prev->count(true);
        else
            return 1+next->count(false);
    else
        if(next)
            return 1+next->count(false);
        else
            return 1;
}
*/
void Organism::printDNA()
{
    if(dna)
        dna->printAll(false);
}

void Organism::printRNA()
{
    if(rna)
        rna->printAll(true);
}

void Organism::printTraits()
{
    using namespace std;

    cout<<"Haploid Number: "<<traits.haploidNumber<<endl;
    cout<<"Cost: "<<traits.cost<<endl;//in ATP
    //Glycolysis
    if(traits.hasHexokinase)
        cout<<"Hexokinase"<<endl;
    if(traits.hasPhosphofructokinase)
        cout<<"Phosphofructokinase"<<endl;
    if(traits.hasAldolase)
        cout<<"Aldolase"<<endl;
    if(traits.hasTriosephosphateIsomerase)
        cout<<"Triosephosphate Isomerase"<<endl;
    if(traits.hasTriosephosphateDehydrogenase)
        cout<<"Triosephosphate Dehydrogenase"<<endl;
    if(traits.hasPhosphoglycerateKinase)
        cout<<"Phosphoglycerate Kinase"<<endl;
    if(traits.hasPyruvateKinase)
        cout<<"Pyruvate Kinase"<<endl;
    //Citric Acid Cycle
    if(traits.hasIsocitrateDehydrogenase)
        cout<<"Isocitrate Dehydrogenase"<<endl;
    if(traits.hasAlphaKetoglutarateDehydrogenase)
        cout<<"Alpha Ketoglutarate Dehydrogenase"<<endl;
    if(traits.hasSuccinylCoASynthetase)
        cout<<"Succinyl CoA Synthetase"<<endl;
    if(traits.hasSdhA)
        cout<<"Succinate Dehydrogenase subunit A"<<endl;
    if(traits.hasMalateDehydrogenase)
        cout<<"Malate Dehydrogenase"<<endl;
    //Chemiosmosis
    if(traits.hasATPSynthase)
        cout<<"ATP Synthase"<<endl;
    //Electron Transport Chain
    if(traits.hasComplexI)
        cout<<"Complex I"<<endl;
    if(traits.hasChain1)
        cout<<"Complex I chain 1"<<endl;
    if(traits.hasChain2)
        cout<<"Complex I chain 2"<<endl;
    if(traits.hasChain3)
        cout<<"Complex I chain 3"<<endl;
    if(traits.hasChain4)
        cout<<"Complex I chain 4"<<endl;
    if(traits.hasChain4L)
        cout<<"Complex I chain 4L"<<endl;
    if(traits.hasChain5)
        cout<<"Complex I chain 5"<<endl;
    if(traits.hasChain6)
        cout<<"Complex I chain 6"<<endl;
    if(traits.hasSdhB)
        cout<<"Succinate Dehydrogenase subunit B"<<endl;
    if(traits.hasSdhC)
        cout<<"Succinate Dehydrogenase subunit C"<<endl;
    if(traits.hasSdhD)
        cout<<"Succinate Dehydrogenase subunit D"<<endl;
    //Photosynthesis
    if(traits.hasProtochlorophyllideReductase)
        cout<<"Protochlorophyllide Reductase"<<endl;
    if(traits.hasChlorophyllBSynthase)
        cout<<"Chlorophyll B Synthase"<<endl;
    if(traits.hasOEC)
        cout<<"OEC"<<endl;
    if(traits.hasOECEP1)
        cout<<"OECEP1"<<endl;
    if(traits.hasOECEP2)
        cout<<"OECEP2"<<endl;
    if(traits.hasOECEP3)
        cout<<"OECEP3"<<endl;
    if(traits.hasRubisco)
        cout<<"RuBisCO"<<endl;
    if(traits.hasRubiscoSmall)
        cout<<"RuBisCO Small"<<endl;
    if(traits.hasRubiscoLarge)
        cout<<"RuBisCO Large"<<endl;
    if(traits.hasPEPCarboxylase)
        cout<<"PEP Carboxylase"<<endl;
    //Motion
    if(traits.hasCEP164)
        cout<<"CEP164"<<endl;
    if(traits.hasSensoryRhodopsinII)
        cout<<"SensoryRhodopsin II"<<endl;
}

void Organism::print()
{
    using std::cout;
    using std::endl;
    if(traits.isNotSequenced)
        sequence();
    cout<<"DNA:"<<endl;
    printDNA();
    cout<<"RNA"<<endl;
    printRNA();
    printTraits();
}

DNA* Organism::replicateDNA(double accuracy)
{
    return dna->replicate(accuracy);
}

void Organism::sequence()
{
    DNA*temp = dna;
    traits.haploidNumber=0;
    DNA* mRNA=0;
    DNA* lastInRNA=0;
    while(temp)
    {
        Data data = temp->getData();
        //transcription
        data.codons^= 0x555555555555555ULL;
        if(data.numValidCodons>10)
            data.numValidCodons=10;
        while(data.numValidCodons)
        {
            if(mRNA)
                switch(data.codons&63)
                {
                    case 1:    //UAA
                    case 9:    //UGA
                    case 33:    //UGA
                        translate(mRNA);//translation
                        if(rna)
                        {
                            if(lastInRNA)
                                lastInRNA = lastInRNA->add(mRNA);
                            else
                                std::cout<<"ERROR: lastInRNA not initialized"<<std::endl;
                        }
                        else
                        {
                            rna = mRNA;
                            lastInRNA = rna;
                        }
                        mRNA = 0;
                        break;
                    default:
                        Data codon;
                        codon.numValidCodons = 1;
                        codon.codons=data.codons&63;
                        *mRNA<<codon;
                }//apply codon
            else
                if((data.codons&63) == 36)//AUG
                {

                    Data aug;
                    aug.codons = 0;
                    aug.numValidCodons = 0;
                    mRNA = new DNA(aug);
                }
            data.codons>>=6;
            data.numValidCodons--;
        }
        if(temp->isEndOfChromosome())
        {
            traits.haploidNumber++;
            if(mRNA)
            {
                translate(mRNA);//translation
                if(rna)
                {
                    if(lastInRNA)
                        lastInRNA = lastInRNA->add(mRNA);
                    else
                        std::cout<<"ERROR: lastInRNA not initialized"<<std::endl;
                }
                else
                {
                    rna = mRNA;
                    lastInRNA = rna;
                }
                mRNA=0;
            }
        }
        temp=temp->next;
    }
    if(traits.hasDNALigase)
    {
        traits.accuracy += 0.0000550;
    }
    if(traits.hasDNAPolymerase)
    {
        traits.accuracy += 0.0099000;
    }
    if(traits.hasDNAPrimase)
    {
        traits.accuracy += 0.0000009;
    }
    if(traits.hasHelicase)
        traits.accuracy += 0.0000400;
    if(traits.hasTopoisomerase)
        traits.accuracy += 0.0000040;

    unsigned int position;
    if(traits.cost)
    {
        position=8;
        while(!(traits.cost>>position))
            position--;
    }
    else
        position=0;
    traits.step = 1.00L/(9-position);
    traits.isNotSequenced=false;
    if(traits.hasRubiscoLarge&&traits.hasRubiscoSmall)
        traits.hasRubisco=1;
    if(traits.hasOECEP1&&traits.hasOECEP2&&traits.hasOECEP3)
        traits.hasOEC=1;
}

Traits* Organism::getTraits()
{
    return &traits;
}

/*
    translate
        applies result of mRNA to Organism
    Start and stop codons should not be included in mRNA
    In this simplified model, 128 amino acids are reduced to 1. The number of amino acids are rounded up from the result of division.
    The necessary amino acids were determined randomly using a random number generator from 0 to 63, skipping the stops
*/
void Organism::translate(DNA* mRNA)
{
    traits.cost+=(mRNA->getCodonCount());
//    mRNA->printAll(true);
    switch(mRNA->getData().codons&63)
    {
        case 0:
        case 32:
                    switch(mRNA->getData().codons>>6&63)
                    {
                        case 19:
                        case 51:
                                    switch(mRNA->getData().codons>>12&63)
                                    {
                                        case 18:
                                        case 50:
                                                    switch(mRNA->getData().codons>>18&63)
                                                    {
                                                        case 6:
                                                        case 22:
                                                        case 38:
                                                        case 54:
                                                                    switch(mRNA->getData().codons>>24&63)
                                                                    {
                                                                        case 2:
                                                                        case 34:
                                                                                    if(mRNA->getCodonCount()==5)
                                                                                        traits.hasATPSynthase=1;
                                                                                break;//Glu
                                                                    }
                                                                break;//Val
                                                    }
                                                break;//Asp
                                    }
                                break;//His
                    }
                break;//Lys
        case 3:
        case 35:
                    switch(mRNA->getData().codons>>6&63)
                    {
                        case 19:
                        case 51:
                                    switch(mRNA->getData().codons>>12&63)
                                    {
                                        case 14:
                                        case 30:
                                        case 46:
                                        case 62:
                                                    switch(mRNA->getData().codons>>18&63)
                                                    {
                                                        case 41:
                                                                    if(mRNA->getCodonCount()==4)
                                                                        traits.hasRubiscoLarge=1;
                                                                break;//Trp
                                                    }
                                                break;//Ala
                                    }
                                break;//His
                    }
                break;//Gln
        case 4:
        case 20:
        case 52:
                    switch(mRNA->getData().codons>>6&63)
                    {
                        case 4:
                        case 20:
                        case 52:
                                    switch(mRNA->getData().codons>>12&63)
                                    {
                                        case 3:
                                        case 35:
                                                    switch(mRNA->getData().codons>>18&63)
                                                    {
                                                        case 5:
                                                        case 7:
                                                        case 23:
                                                        case 37:
                                                        case 39:
                                                        case 55:
                                                                    switch(mRNA->getData().codons>>24&63)
                                                                    {
                                                                        case 14:
                                                                        case 30:
                                                                        case 46:
                                                                        case 62:
                                                                                    if(mRNA->getCodonCount()==5)
                                                                                        traits.hasChain5=1;
                                                                                break;//Ala
                                                                    }
                                                                break;//Leu
                                                    }
                                                break;//Gln
                                    }
                                break;//Ile
                        case 8:
                        case 11:
                        case 27:
                        case 40:
                        case 43:
                        case 49:
                                    switch(mRNA->getData().codons>>12&63)
                                    {
                                        case 13:
                                        case 24:
                                        case 29:
                                        case 45:
                                        case 56:
                                        case 61:
                                                    if(mRNA->getCodonCount()==3)
                                                        traits.hasTriosephosphateDehydrogenase=1;
                                                break;//Ser
                                    }
                                break;//Arg
                        case 13:
                        case 24:
                        case 29:
                        case 45:
                        case 56:
                        case 61:
                                    switch(mRNA->getData().codons>>12&63)
                                    {
                                        case 19:
                                        case 51:
                                                    switch(mRNA->getData().codons>>18&63)
                                                    {
                                                        case 8:
                                                        case 11:
                                                        case 27:
                                                        case 40:
                                                        case 43:
                                                        case 49:
                                                                    if(mRNA->getCodonCount()==4)
                                                                        traits.hasIsocitrateDehydrogenase=1;
                                                                break;//Arg
                                                    }
                                                break;//His
                                    }
                                break;//Ser
                    }
                break;//Ile
        case 5:
        case 7:
        case 23:
        case 37:
        case 39:
        case 55:
                    switch(mRNA->getData().codons>>6&63)
                    {
                        case 2:
                        case 34:
                                    switch(mRNA->getData().codons>>12&63)
                                    {
                                        case 25:
                                        case 57:
                                                    switch(mRNA->getData().codons>>18&63)
                                                    {
                                                        case 13:
                                                        case 24:
                                                        case 29:
                                                        case 45:
                                                        case 56:
                                                        case 61:
                                                                    switch(mRNA->getData().codons>>24&63)
                                                                    {
                                                                        case 8:
                                                                        case 11:
                                                                        case 27:
                                                                        case 40:
                                                                        case 43:
                                                                        case 59:
                                                                                    switch(mRNA->getData().codons>>30&63)
                                                                                    {
                                                                                        case 12:
                                                                                        case 28:
                                                                                        case 44:
                                                                                        case 60:
                                                                                                    switch(mRNA->getData().codons>>36&63)
                                                                                                    {
                                                                                                        case 19:
                                                                                                        case 51:
                                                                                                                    switch(mRNA->getData().codons>>42&63)
                                                                                                                    {
                                                                                                                        case 8:
                                                                                                                        case 11:
                                                                                                                        case 27:
                                                                                                                        case 40:
                                                                                                                        case 43:
                                                                                                                        case 59:
                                                                                                                                    if(mRNA->getCodonCount()==8)
                                                                                                                                        traits.hasAlphaKetoglutarateDehydrogenase=1;
                                                                                                                                break;//Arg
                                                                                                                    }
                                                                                                                break;//His
                                                                                                    }
                                                                                                break;//Thr
                                                                                    }
                                                                                break;//Arg
                                                                    }
                                                                break;//Ser
                                                    }
                                                break;//Lys
                                    }
                                break;//Glu
                        case 12:
                        case 28:
                        case 44:
                        case 60:
                                    switch(mRNA->getData().codons>>12&63)
                                    {
                                        case 6:
                                        case 22:
                                        case 38:
                                        case 54:
                                                    switch(mRNA->getData().codons>>18&63)
                                                    {
                                                        case 10:
                                                        case 26:
                                                        case 42:
                                                        case 58:
                                                                    switch(mRNA->getData().codons>>24&63)
                                                                    {
                                                                        case 15:
                                                                        case 31:
                                                                        case 47:
                                                                        case 63:
                                                                                    switch(mRNA->getData().codons>>30&63)
                                                                                    {
                                                                                        case 36:
                                                                                                    switch(mRNA->getData().codons>>36&63)
                                                                                                    {
                                                                                                        case 18:
                                                                                                        case 50:
                                                                                                                    switch(mRNA->getData().codons)
                                                                                                                    {
                                                                                                                        case 13:
                                                                                                                        case 24:
                                                                                                                                    if(mRNA->getCodonCount()==7)
                                                                                                                                        traits.hasHexokinase=1;
                                                                                                                                    break;//Ser
                                                                                                                    }
                                                                                                                break;//Asp
                                                                                                    }
                                                                                                break;//Met
                                                                                    }
                                                                                break;//Pro
                                                                    }
                                                                break;//Gly
                                                    }
                                                break;//Val
                                    }
                                break;//Thr
                        case 13:
                        case 24:
                        case 29:
                        case 45:
                        case 56:
                        case 61:
                                    if(mRNA->getCodonCount()==2)
                                        traits.hasPyruvateKinase=1;
                                break;//Ser
                        case 15:
                        case 31:
                        case 47:
                        case 63:
                                    if(mRNA->getCodonCount()==2)
                                        traits.hasOECEP3=1;
                                    else switch(mRNA->getData().codons>>12&63)
                                    {
                                        case 14:
                                        case 30:
                                        case 46:
                                        case 62:
                                                    if(mRNA->getCodonCount()==3)
                                                        traits.hasMalateDehydrogenase=1;
                                                break;//Ala
                                    }
                                break;//Pro
                    }
                break;//Leu
        case 6:
        case 22:
        case 38:
        case 54:
                if(mRNA->getCodonCount()==1)
                    traits.hasChain4L=1;
                else
                    switch(mRNA->getData().codons>>6&63)
                    {
                        case 8:
                        case 11:
                        case 27:
                        case 40:
                        case 43:
                        case 59:
                                    switch(mRNA->getData().codons>>12&63)
                                    {
                                        case 25:
                                        case 57:
                                                    switch(mRNA->getData().codons>>18&63)
                                                    {
                                                        case 15:
                                                        case 31:
                                                        case 47:
                                                        case 63:
                                                                    switch(mRNA->getData().codons>>24&63)
                                                                    {
                                                                        case 10:
                                                                        case 26:
                                                                        case 42:
                                                                        case 58:
                                                                                    if(mRNA->getCodonCount()==5)
                                                                                        traits.hasChlorophyllBSynthase=1;
                                                                                break;//Gly
                                                                    }
                                                                break;//Pro
                                                    }
                                                break;//Cys
                                    }
                                break;//Arg
                        case 15:
                        case 31:
                        case 47:
                        case 63:
                                    switch(mRNA->getData().codons>>12&63)
                                    {
                                        case 15:
                                        case 31:
                                        case 47:
                                        case 63:
                                                    switch(mRNA->getData().codons>>18&63)
                                                    {
                                                        case 5:
                                                        case 7:
                                                        case 23:
                                                        case 37:
                                                        case 39:
                                                        case 55:
                                                                    switch(mRNA->getData().codons>>24&63)
                                                                    {
                                                                        case 19:
                                                                        case 51:
                                                                                    switch(mRNA->getData().codons>>30&63)
                                                                                    {
                                                                                        case 13:
                                                                                        case 24:
                                                                                        case 29:
                                                                                        case 45:
                                                                                        case 56:
                                                                                        case 61:
                                                                                                    switch(mRNA->getData().codons>>36&63)
                                                                                                    {
                                                                                                        case 12:
                                                                                                        case 28:
                                                                                                        case 44:
                                                                                                        case 60:
                                                                                                                    switch(mRNA->getData().codons>>42&63)
                                                                                                                    {
                                                                                                                        case 16:
                                                                                                                        case 48:
                                                                                                                                    if(mRNA->getCodonCount()==8)
                                                                                                                                        traits.hasDNAPolymerase=1;
                                                                                                                                break;//Asn
                                                                                                                    }
                                                                                                                break;//Thr
                                                                                                    }
                                                                                                break;//Ser
                                                                                    }
                                                                                break;//His
                                                                    }
                                                                break;//Leu
                                                    }
                                                break;//Pro
                                    }
                                break;//Pro
                        case 18:
                        case 50:
                                    switch(mRNA->getData().codons>>12&63)
                                    {
                                        case 10:
                                        case 26:
                                        case 42:
                                        case 58:
                                                    switch(mRNA->getData().codons>>18&63)
                                                    {
                                                        case 2:
                                                        case 34:
                                                                    switch(mRNA->getData().codons>>24&63)
                                                                    {
                                                                        case 5:
                                                                        case 7:
                                                                        case 23:
                                                                        case 37:
                                                                        case 39:
                                                                        case 55:
                                                                                    switch(mRNA->getData().codons>>30&63)
                                                                                    {
                                                                                        case 8:
                                                                                        case 11:
                                                                                        case 27:
                                                                                        case 40:
                                                                                        case 43:
                                                                                        case 59:
                                                                                                    if(mRNA->getCodonCount()==6)
                                                                                                        traits.hasTopoisomerase=1;
                                                                                                break;//Arg
                                                                                    }
                                                                                break;//Leu
                                                                    }
                                                                break;//Glu
                                                    }
                                                break;//Gly
                                    }
                                break;//Asp
                    }
                break;//Val
        case 8:
        case 11:
        case 27:
        case 40:
        case 43:
        case 59:
                if(mRNA->getCodonCount()==1)
                    traits.hasSdhD=1;
                else
                    switch(mRNA->getData().codons>>6&63)
                    {
                        case 4:
                        case 20:
                        case 52:
                                if(mRNA->getCodonCount()==2)
                                    traits.hasSdhC=1;
                                else
                                    switch(mRNA->getData().codons>>12&63)
                                    {
                                        case 10:
                                        case 26:
                                        case 42:
                                        case 58:
                                                    if(mRNA->getCodonCount()==3)
                                                        traits.hasSuccinylCoASynthetase=1;
                                                break;//Gly
                                    }
                                break;//Ile
                        case 8:
                        case 11:
                        case 27:
                        case 40:
                        case 43:
                        case 49:
                                    switch(mRNA->getData().codons>>12&63)
                                    {
                                        case 12:
                                        case 28:
                                        case 44:
                                        case 60:
                                                    switch(mRNA->getData().codons>>18&63)
                                                    {
                                                        case 25:
                                                        case 57:
                                                                    if(mRNA->getCodonCount()==4)
                                                                        traits.hasPhosphoglycerateKinase=1;
                                                                break;//Lys
                                                    }
                                                break;//Thr
                                        case 18:
                                        case 50:
                                                    switch(mRNA->getData().codons>>18&63)
                                                    {
                                                        case 13:
                                                        case 24:
                                                        case 29:
                                                        case 45:
                                                        case 56:
                                                        case 61:
                                                                    switch(mRNA->getData().codons>>24&63)
                                                                    {
                                                                        case 5:
                                                                        case 7:
                                                                        case 23:
                                                                        case 37:
                                                                        case 39:
                                                                        case 55:
                                                                                    if(mRNA->getCodonCount()==5)
                                                                                        traits.hasDNAPrimase=1;
                                                                                break;//Leu
                                                                    }
                                                                break;//Ser
                                                    }
                                                break;//Asp
                                    }
                                break;//Arg
                        case 12:
                        case 28:
                        case 44:
                        case 60:
                                    switch(mRNA->getData().codons>>12&63)
                                    {
                                        case 8:
                                        case 11:
                                        case 27:
                                        case 40:
                                        case 43:
                                        case 59:
                                                    switch(mRNA->getData().codons>>18&63)
                                                    {
                                                        case 13:
                                                        case 24:
                                                        case 29:
                                                        case 45:
                                                        case 56:
                                                        case 61:
                                                                    switch(mRNA->getData().codons>>24&63)
                                                                    {
                                                                        case 16:
                                                                        case 48:
                                                                                    switch(mRNA->getData().codons>>30&63)
                                                                                    {
                                                                                        case 15:
                                                                                        case 31:
                                                                                        case 47:
                                                                                        case 63:
                                                                                                    if(mRNA->getCodonCount()==6)
                                                                                                        traits.hasHelicase=1;
                                                                                                break;//Pro
                                                                                    }
                                                                                break;//Asn
                                                                    }
                                                                break;//Ser
                                                    }
                                                break;//Arg
                                    }
                                break;//Thr
                        case 25:
                        case 57:
                                    switch(mRNA->getData().codons>>12&63)
                                    {
                                        case 13:
                                        case 24:
                                        case 29:
                                        case 45:
                                        case 56:
                                        case 61:
                                                    if(mRNA->getCodonCount()==3)
                                                        traits.hasPhosphofructokinase=1;
                                                break;//Ser
                                    }
                                break;//Cys
                    }
                break;//Arg
        case 10:
        case 26:
        case 42:
        case 58:
                    switch(mRNA->getData().codons>>6&63)
                    {
                        case 14:
                        case 30:
                        case 46:
                        case 62:
                                    if(mRNA->getCodonCount()==2)
                                        traits.hasChain6=1;
                    }
                break;//Gly
        case 12:
        case 28:
        case 44:
        case 60:
                    switch(mRNA->getData().codons>>6&63)
                    {
                        case 0:
                        case 32:
                                    if(mRNA->getCodonCount()==2)
                                        traits.hasSensoryRhodopsinII=1;
                                break;//Lys
                        case 6:
                        case 22:
                        case 38:
                        case 54:
                                    switch(mRNA->getData().codons>>12&63)
                                    {
                                        case 5:
                                        case 7:
                                        case 23:
                                        case 37:
                                        case 39:
                                        case 55:
                                                    if(mRNA->getCodonCount()==2)
                                                        traits.hasAldolase=1;
                                            break;//Leu
                                    }
                                break;//Val
                        case 25:
                        case 57:
                                    switch(mRNA->getData().codons>>12&63)
                                    {
                                        case 36:
                                                    switch(mRNA->getData().codons>>18&63)
                                                    {
                                                        case 8:
                                                        case 11:
                                                        case 27:
                                                        case 40:
                                                        case 43:
                                                        case 59:
                                                                    switch(mRNA->getData().codons>>24&63)
                                                                    {
                                                                        case 3:
                                                                        case 35:
                                                                                    switch(mRNA->getData().codons>>30&63)
                                                                                    {
                                                                                        case 19:
                                                                                        case 51:
                                                                                                    switch(mRNA->getData().codons>>36&63)
                                                                                                    {
                                                                                                        case 3:
                                                                                                        case 35:
                                                                                                                    switch(mRNA->getData().codons>>42&63)
                                                                                                                    {
                                                                                                                        case 25:
                                                                                                                        case 57:
                                                                                                                                    switch(mRNA->getData().codons>>48&63)
                                                                                                                                    {
                                                                                                                                        case 6:
                                                                                                                                        case 22:
                                                                                                                                        case 38:
                                                                                                                                        case 58:
                                                                                                                                                    switch(mRNA->getData().codons>>54&63)
                                                                                                                                                    {
                                                                                                                                                        case 19:
                                                                                                                                                        case 51:
                                                                                                                                                                    switch(mRNA->getNext()->getData().codons&63)
                                                                                                                                                                    {
                                                                                                                                                                        case 4:
                                                                                                                                                                        case 20:
                                                                                                                                                                        case 52:
                                                                                                                                                                                    switch(mRNA->getNext()->getData().codons>>6&63)
                                                                                                                                                                                    {
                                                                                                                                                                                        case 14:
                                                                                                                                                                                        case 30:
                                                                                                                                                                                        case 46:
                                                                                                                                                                                        case 62:
                                                                                                                                                                                                    if(mRNA->getCodonCount())
                                                                                                                                                                                                        traits.hasCEP164=1;
                                                                                                                                                                                                break;//Ala
                                                                                                                                                                                    }
                                                                                                                                                                                break;//Ile
                                                                                                                                                                    }
                                                                                                                                                                break;//His
                                                                                                                                                    }
                                                                                                                                                break;//Val
                                                                                                                                    }
                                                                                                                                break;//Cys
                                                                                                                    }
                                                                                                                break;//Gln
                                                                                                    }
                                                                                                break;//His
                                                                                    }
                                                                                break;//Gln
                                                                    }
                                                                break;//Arg
                                                    }
                                                break;//Met
                                    }
                                break;//Cys
                    }
                break;//Thr
        case 13:
        case 24:
        case 29:
        case 45:
        case 56:
        case 61:
                    switch(mRNA->getData().codons>>6&63)
                    {
                        case 6:
                        case 22:
                        case 38:
                        case 54:
                                    if(mRNA->getCodonCount()==2)
                                        traits.hasRubiscoSmall=1;
                                break;//Val
                        case 8:
                        case 11:
                        case 27:
                        case 40:
                        case 43:
                        case 59:
                                    switch(mRNA->getData().codons>>12&63)
                                    {
                                        case 15:
                                        case 31:
                                        case 47:
                                        case 63:
                                                    switch(mRNA->getData().codons>>18&63)
                                                    {
                                                        case 8:
                                                        case 11:
                                                        case 27:
                                                        case 40:
                                                        case 43:
                                                        case 59:
                                                                    switch(mRNA->getData().codons>>24&63)
                                                                    {
                                                                        case 41:
                                                                                    if(mRNA->getCodonCount()==5)
                                                                                        traits.hasSdhA=1;
                                                                                break;//Trp
                                                                    }
                                                                break;//Arg
                                                    }
                                                break;//Pro
                                        case 16:
                                        case 48:
                                                    if(mRNA->getCodonCount()==3)
                                                        traits.hasChain1=1;
                                                break;//Asn
                                    }
                                break;//Arg
                        case 25:
                        case 57:
                                    if(mRNA->getCodonCount()==2)
                                        traits.hasSdhB=1;
                                break;//Cys
                        case 41:
                                    if(mRNA->getCodonCount()==2)
                                        traits.hasOECEP2=1;
                                break;//Trp
                    }
                break;//Ser
        case 15:
        case 31:
        case 47:
        case 63:
                if(mRNA->getCodonCount()==1)
                    traits.hasChain3=1;
                else
                    switch(mRNA->getData().codons>>6&63)
                    {
                        case 15:
                        case 31:
                        case 47:
                        case 63:
                                    switch(mRNA->getData().codons>>12&63)
                                    {
                                        case 15:
                                        case 31:
                                        case 47:
                                        case 63:
                                                    switch(mRNA->getData().codons>>18&63)
                                                    {
                                                        case 18:
                                                        case 50:
                                                                    switch(mRNA->getData().codons>>24&63)
                                                                    {
                                                                        case 6:
                                                                        case 22:
                                                                        case 38:
                                                                        case 54:
                                                                                    switch(mRNA->getData().codons>>30&63)
                                                                                    {
                                                                                        case 4:
                                                                                        case 20:
                                                                                        case 52:
                                                                                                    switch(mRNA->getData().codons>>24&63)
                                                                                                    {
                                                                                                        case 6:
                                                                                                        case 22:
                                                                                                        case 38:
                                                                                                        case 54:
                                                                                                                    if(mRNA->getCodonCount()==7)
                                                                                                                        traits.hasPEPCarboxylase=1;
                                                                                                                break;//Val
                                                                                                    }
                                                                                                break;//Ile
                                                                                    }
                                                                                break;//Val
                                                                    }
                                                                break;//Asp
                                                    }
                                                break;//Pro
                                    }
                                break;//Pro
                    }
                break;//Pro
        case 16:
        case 48:
                    switch(mRNA->getData().codons>>6&63)
                    {
                        case 13:
                        case 24:
                        case 29:
                        case 45:
                        case 56:
                        case 61:
                                    switch(mRNA->getData().codons>>12&63)
                                    {
                                        case 8:
                                        case 11:
                                        case 27:
                                        case 40:
                                        case 43:
                                        case 59:
                                                    if(mRNA->getCodonCount()==3)
                                                        traits.hasChain2=1;
                                                break;//Arg
                                    }
                                break;//Ser
                    }
                break;//Asn
        case 17:
        case 49:
                    switch(mRNA->getData().codons>>6&63)
                    {
                        case 41:
                                    switch(mRNA->getData().codons>>12&63)
                                    {
                                        case 10:
                                        case 26:
                                        case 42:
                                        case 58:
                                                    switch(mRNA->getData().codons>>18&63)
                                                    {
                                                        case 4:
                                                        case 20:
                                                        case 52:
                                                                    switch(mRNA->getData().codons>>24&63)
                                                                    {
                                                                        case 8:
                                                                        case 11:
                                                                        case 27:
                                                                        case 40:
                                                                        case 43:
                                                                        case 59:
                                                                                    switch(mRNA->getData().codons>>30&63)
                                                                                    {
                                                                                        case 4:
                                                                                        case 20:
                                                                                        case 52:
                                                                                                    switch(mRNA->getData().codons>>36&63)
                                                                                                    {
                                                                                                        case 5:
                                                                                                        case 7:
                                                                                                        case 23:
                                                                                                        case 37:
                                                                                                        case 39:
                                                                                                        case 55:
                                                                                                                    switch(mRNA->getData().codons>>42&63)
                                                                                                                    {
                                                                                                                        case 17:
                                                                                                                        case 49:
                                                                                                                                    if(mRNA->getCodonCount()==8)
                                                                                                                                        traits.hasDNALigase=1;
                                                                                                                                break;//Tyr
                                                                                                                    }
                                                                                                                break;//Leu
                                                                                                    }
                                                                                                break;//Ile
                                                                                    }
                                                                                break;//Arg
                                                                    }
                                                                    break;//Ile
                                                    }
                                                break;//Gly
                                    }
                                break;//Trp
                    }
                break;//Tyr
        case 21:
        case 53:
                    switch(mRNA->getData().codons>>6&63)
                    {
                        case 13:
                        case 24:
                        case 29:
                        case 45:
                        case 56:
                        case 61:
                                    switch(mRNA->getData().codons>>12&63)
                                    {
                                        case 10:
                                        case 26:
                                        case 42:
                                        case 58:
                                                    if(mRNA->getCodonCount()==3)
                                                        traits.hasTriosephosphateIsomerase =1;
                                                break;//Gly
                                    }
                                break;//Ser
                    }
                break;//Phe
        case 25:
        case 57:
                    switch(mRNA->getData().codons>>6&63)
                    {
                        case 13:
                        case 24:
                        case 29:
                        case 45:
                        case 56:
                        case 61:
                                    switch(mRNA->getData().codons>>12&63)
                                    {
                                        case 21:
                                        case 53:
                                                    if(mRNA->getCodonCount() == 3)
                                                        traits.hasOECEP1 = 1;
                                                break;//Phe
                                    }
                                break;//Ser
                    }
                break;//Cys
        case 36:
                    switch(mRNA->getData().codons>>6&63)
                    {
                        case 5:
                        case 7:
                        case 23:
                        case 37:
                        case 39:
                        case 55:
                                    switch(mRNA->getData().codons>>12&63)
                                    {
                                        case 12:
                                        case 28:
                                        case 44:
                                        case 60:
                                                    switch(mRNA->getData().codons>>18&63)
                                                    {
                                                        case 10:
                                                        case 26:
                                                        case 42:
                                                        case 58:
                                                                    if(mRNA->getCodonCount()==4)
                                                                        traits.hasChain4=1;
                                                                break;//Gly
                                                    }
                                                break;//Thr
                                    }
                                break;//Leu
                    }
                break;//Met
        case 41:
                    switch(mRNA->getData().codons>>6&63)
                    {
                        case 15:
                        case 31:
                        case 47:
                        case 63:
                                    switch(mRNA->getData().codons>>12&63)
                                    {
                                        case 2:
                                        case 34:
                                                    switch(mRNA->getData().codons>>18&63)
                                                    {
                                                        case 10:
                                                        case 26:
                                                        case 42:
                                                        case 58:
                                                                    if(mRNA->getCodonCount()==4)
                                                                        traits.hasProtochlorophyllideReductase=1;
                                                                break;//Gly
                                                    }
                                                break;//Glu
                                    }
                                break;//Pro
                    }
                break;//Trp
    }//<Met>
}//translation
