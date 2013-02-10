#include "../include/Environment.h"
#include "stdio.h"
#include "stdlib.h"
#include <iostream>
#include <fstream>
#include "string.h"
Environment::Environment(const unsigned int rows,const unsigned int columns,unsigned long long generationNumber)
{
    this->rows = rows;
    this->cols = columns;
    this->generationNumber = generationNumber;
    grid = new Location*[rows];
    for(unsigned int row = 0;row<rows;row++)
        grid[row]=new Location[columns];
    numLocations = rows*cols;
    atmosphere.oxygen = numLocations;
    for(unsigned int row = 0;row<rows;row++)
        for(unsigned int col = 0;col<cols;col++)
        {
            grid[row][col].organism = 0;
            grid[row][col].sunlight = 0;
            grid[row][col].protonGradient = 0;
            grid[row][col].glucose = 0;
            grid[row][col].glucose6Phosphate = 0;
            grid[row][col].fructose16Bisphosphate = 0;
            grid[row][col].G3P = 0;
            grid[row][col].DHAP = 0;
            grid[row][col].one3Bisphosphoglycerate = 0;
            grid[row][col].NADH = 0;
            grid[row][col].NADPH = 0;
            grid[row][col].threePhosphoglycerate = 0;
            grid[row][col].pyruvate = 0;
            grid[row][col].alphaKetoglutarate = 0;
            grid[row][col].succinylCoA = 0;
            grid[row][col].succinate = 0;
            grid[row][col].fumerate = 0;
            grid[row][col].FADH2 = 0;
        }
    std::cout<<"Environment created"<<std::endl;
}

Environment::~Environment()
{
    for(unsigned int row = 0;row<rows;row++)
    {
        for(unsigned int col = 0;col<cols;col++)
        {
            delete grid[row][col].organism;
            grid[row][col].organism = 0;
        }
        delete[] grid[row];
        grid[row] = 0;
    }
    delete[] grid;
    grid = 0;
    std::cout<<"Environment deleted"<<std::endl;
}

unsigned int Environment::getRows()
{
    return rows;
}

unsigned int Environment::getCols()
{
    return cols;
}

Location Environment::getLocation(unsigned int row,unsigned int col)
{
    return grid[row][col];
}

Organism* Environment::getOrganism(unsigned int row,unsigned int col)
{
    return grid[row][col].organism;
}

bool Environment::insert(Organism*organism,unsigned int row,unsigned int col)
{
    if(row>=rows||col>=cols||grid[row][col].organism)
            return false;
    grid[row][col].organism = organism;
    return true;
}

bool Environment::place(Organism*organism,unsigned int centerRow,unsigned int centerCol)
{
    unsigned int attempt = rand()&7;
    switch(attempt)
    {
        case 0: if(insert(organism,centerRow-1,centerCol-1))
                    return true;
        case 1: if(insert(organism,centerRow-1,centerCol))
                    return true;
        case 2: if(insert(organism,centerRow-1,centerCol+1))
                    return true;
        case 3: if(insert(organism,centerRow,centerCol-1))
                    return true;
        case 4: if(insert(organism,centerRow,centerCol+1))
                    return true;
        case 5: if(insert(organism,centerRow+1,centerCol-1))
                    return true;
        case 6: if(insert(organism,centerRow+1,centerCol))
                    return true;
        case 7: if(insert(organism,centerRow-1,centerCol+1))
                    return true;
                if(attempt == 0)
                    return false;
                if(insert(organism,centerRow-1,centerCol-1))
                    return true;
                if(attempt == 1)
                    return false;
                if(insert(organism,centerRow-1,centerCol))
                    return true;
                if(attempt == 2)
                    return false;
                if(insert(organism,centerRow-1,centerCol+1))
                    return true;
                if(attempt == 3)
                    return false;
                if(insert(organism,centerRow,centerCol-1))
                    return true;
                if(attempt == 4)
                    return false;
                if(insert(organism,centerRow,centerCol+1))
                    return true;
                if(attempt == 5)
                    return false;
                if(insert(organism,centerRow+1,centerCol-1))
                    return true;
                if(attempt == 6)
                    return false;
                if(insert(organism,centerRow+1,centerCol))
                    return true;
                return false;
        default:break;//Impossible
    }
    std::cout<<"Unreachable"<<std::endl;
    return false;//Appeases compiler
}

Organism* Environment::removeOrganism(unsigned int row, unsigned int col)
{
    Organism* toReturn = grid[row][col].organism;
    grid[row][col].organism=0;
    return toReturn;
}

void Environment::setLighting(unsigned int focusRow,unsigned int focusCol,unsigned int intensity,long double decay)
{
    for(unsigned int row = 0;row < rows; row++)
        for(unsigned int col = 0;col < cols; col++)
        {
            grid[row][col].sunlight = 0;
        }
    unsigned int distance = 0;
    while(intensity)
    {
        for(unsigned int row = (focusRow - distance)>rows?0 : focusRow - distance;row < focusRow + distance + 1;row++)
        {
            if(focusCol - distance < cols && row < rows)
                grid[row][focusCol - distance].sunlight = intensity;
            if(focusCol + distance < cols && row < rows)
                grid[row][focusCol + distance].sunlight = intensity;
        }
        for(unsigned int col = (focusCol - distance + 1)>cols?0 : focusCol - distance + 1; col < focusCol + distance;col++)
        {
            if(focusRow - distance < rows && col < cols)
                grid[focusRow - distance][col].sunlight = intensity;
            if(focusRow + distance < rows && col < cols)
                grid[focusRow + distance][col].sunlight = intensity;
        }
        intensity *= decay; //reduces and rounds down
        distance++;
    }

}

unsigned long long Environment::countOrganisms()
{
    unsigned long long count = 0;
    for(unsigned int row = 0;row<rows;row++)
        for(unsigned int col = 0;col<cols;col++)
            if(grid[row][col].organism)
                count++;
    return count;
}

Organism* Environment::reproduce(unsigned int row,unsigned int col,long double turn)
{
    if(grid[row][col].organism->getTraits()->isNotSequenced)
    {
        grid[row][col].organism->sequence();
    }
    return new Organism(turn,grid[row][col].organism->replicateDNA(grid[row][col].organism->getTraits()->accuracy));
}

void Environment::stir()
{
    for(unsigned int row = 0; row < rows; row++)
        for(unsigned int col = 0; col < cols; col++)
            if(grid[row][col].organism)
            {
                if(place(grid[row][col].organism,row,col))
                    grid[row][col].organism = 0;
            }
}

void Environment::simulateGeneration()
{
    generationNumber++;
    naturalSelection();
}

void Environment::simulateGenerations(unsigned long long number)
{
    generationNumber+=number;
    naturalSelection();
}

void Environment::simulateUntilGeneration(unsigned long long number)
{
    if(generationNumber < number)
    {
        generationNumber = number;
        naturalSelection();
    }
}

void Environment::print()
{
    using std::cout;
    using std::endl;
    cout<<"Oxygen:\t"<<atmosphere.oxygen<<endl;
    for(unsigned int row=0;row<rows;row++)
        for(unsigned int col=0;col<cols;col++)
        {
            cout<<'('<<row<<','<<col<<')'<<endl;
            if(grid[row][col].organism)
                grid[row][col].organism->print();
            if(grid[row][col].sunlight)
                cout<<"Sunlight:\t"<<grid[row][col].sunlight<<endl;
            if(grid[row][col].protonGradient)
                cout<<"Proton gradient:\t"<<grid[row][col].protonGradient<<endl;
            if(grid[row][col].glucose)
                cout<<"Glucose:\t"<<grid[row][col].glucose<<endl;
            if(grid[row][col].glucose6Phosphate)
                cout<<"Glucose 6-phosphate:\t"<<grid[row][col].glucose6Phosphate<<endl;
            if(grid[row][col].fructose16Bisphosphate)
                cout<<"Fructose 1,6-Bisphosphate:\t"<<grid[row][col].fructose16Bisphosphate<<endl;
            if(grid[row][col].G3P)
                cout<<"Glyceraldehide 3-phosphate:\t"<<grid[row][col].G3P<<endl;
            if(grid[row][col].DHAP)
                cout<<"Dihydroxyacetone phosphate:\t"<<grid[row][col].DHAP<<endl;
            if(grid[row][col].one3Bisphosphoglycerate)
                cout<<"1,3-Bisphosphoglycerate:\t"<<grid[row][col].one3Bisphosphoglycerate<<endl;
            if(grid[row][col].NADH)
                cout<<"NADH:\t"<<grid[row][col].NADH<<endl;
            if(grid[row][col].NADPH)
                cout<<"NADPH:\t"<<grid[row][col].NADPH<<endl;
            if(grid[row][col].pyruvate)
                cout<<"Pyruvate:\t"<<grid[row][col].pyruvate<<endl;
            if(grid[row][col].alphaKetoglutarate)
                cout<<"α-Ketoglutarate:\t"<<grid[row][col].alphaKetoglutarate<<endl;
            if(grid[row][col].succinylCoA)
                cout<<"SuccinylCoA:\t"<<grid[row][col].succinylCoA<<endl;
            if(grid[row][col].succinate)
                cout<<"Succinate:\t"<<grid[row][col].succinate<<endl;
            if(grid[row][col].fumerate)
                cout<<"Fumerate:\t"<<grid[row][col].fumerate<<endl;
            if(grid[row][col].FADH2)
                cout<<"FADH2:\t"<<grid[row][col].FADH2<<endl;
            cout<<endl;
        }
}

Environment::Environment(std::string& filename)
{
    using std::ifstream;
    ifstream in(filename.c_str());
    in>>rows;
    in>>cols;
    grid = new Location*[rows];
    in>>atmosphere.oxygen;
    for(unsigned int row = 0;row<rows;row++)
        grid[row] = new Location[cols];
    in>>generationNumber;
    for(unsigned int row = 0;row<rows;row++)
        for(unsigned int col = 0;col<cols;col++)
        {
            char input;
            in >> input;
            if(input == '@')
            {
                grid[row][col].organism = new Organism(in);
                ///load organism
            }
            else
                grid[row][col].organism = 0;
            in>>grid[row][col].sunlight;
            in>>grid[row][col].protonGradient;
            in>>grid[row][col].glucose;
            in>>grid[row][col].glucose6Phosphate;
            in>>grid[row][col].fructose16Bisphosphate;
            in>>grid[row][col].G3P;
            in>>grid[row][col].DHAP;
            in>>grid[row][col].one3Bisphosphoglycerate;
            in>>grid[row][col].NADH;
            in>>grid[row][col].NADPH;
            in>>grid[row][col].threePhosphoglycerate;
            in>>grid[row][col].pyruvate;
            in>>grid[row][col].alphaKetoglutarate;
            in>>grid[row][col].succinylCoA;
            in>>grid[row][col].succinate;
            in>>grid[row][col].fumerate;
            in>>grid[row][col].FADH2;
        }
    in.close();
    std::cout<<"Environment loaded"<<std::endl;
}

void Environment::saveToFile(std::string& filename)
{
    using std::ofstream;
    using std::endl;
    ofstream out(filename.c_str());
    out<<rows<<endl;
    out<<cols<<endl;
    out<<atmosphere.oxygen<<endl;
    out<<generationNumber<<endl;
    for(unsigned int row = 0;row<rows;row++)
        for(unsigned int col = 0;col<cols;col++)
        {
            if(grid[row][col].organism)
            {
                out<<'@'<<endl;
                grid[row][col].organism->saveToFile(out);
            }
            else
            {
                out<<'/'<<endl;
            }
            out<<grid[row][col].sunlight<<endl;
            out<<grid[row][col].protonGradient<<endl;
            out<<grid[row][col].glucose<<endl;
            out<<grid[row][col].glucose6Phosphate<<endl;
            out<<grid[row][col].fructose16Bisphosphate<<endl;
            out<<grid[row][col].G3P<<endl;
            out<<grid[row][col].DHAP<<endl;
            out<<grid[row][col].one3Bisphosphoglycerate<<endl;
            out<<grid[row][col].NADH<<endl;
            out<<grid[row][col].NADPH<<endl;
            out<<grid[row][col].threePhosphoglycerate<<endl;
            out<<grid[row][col].pyruvate<<endl;
            out<<grid[row][col].alphaKetoglutarate<<endl;
            out<<grid[row][col].succinylCoA<<endl;
            out<<grid[row][col].succinate<<endl;
            out<<grid[row][col].fumerate<<endl;
            out<<grid[row][col].FADH2<<endl;
        }
    out.close();
    std::cout<<"Environment saved"<<std::endl;
}

void Environment::printAnalysis()
{
    unsigned int total = 0;
    unsigned int countProto = 0;
    for(unsigned int row = 0;row<rows;row++)
        for(unsigned int col = 0;col<cols;col++)
            if(grid[row][col].organism)
            {
                total++;
                if(grid[row][col].organism->getTraits()->hasProtochlorophyllideReductase)
                    countProto++;
            }
    std::cout<<atmosphere.oxygen<<' '<<total<<' '<<countProto<<std::endl;
}

void Environment::printMicroscopeLine(std::ofstream& out)
{
    unsigned long long totalCount = 0;
    unsigned long long chlorophyllACount = 0;
    unsigned long long chlorophyllBCount = 0;
    unsigned long long atpSynthaseCount = 0;
    unsigned long long topoIsomeraseCount = 0;
    unsigned long long primaseCount = 0;
    unsigned long long ligaseCount = 0;
    unsigned long long helicaseCount = 0;
    unsigned long long polymeraseCount = 0;
    unsigned long long cep164Count = 0;
    for(unsigned int row = 0;row<rows;row++)
        for(unsigned int col = 0;col<cols;col++)
        {
            if(grid[row][col].organism)
            {
                totalCount++;
                if(grid[row][col].organism->getTraits()->hasProtochlorophyllideReductase)
                    chlorophyllACount++;
                if(grid[row][col].organism->getTraits()->hasChlorophyllBSynthase)
                    chlorophyllBCount++;
                if(grid[row][col].organism->getTraits()->hasATPSynthase)
                    atpSynthaseCount++;
                if(grid[row][col].organism->getTraits()->hasTopoisomerase)
                    topoIsomeraseCount++;
                if(grid[row][col].organism->getTraits()->hasDNAPrimase)
                    primaseCount++;
                if(grid[row][col].organism->getTraits()->hasDNALigase)
                    ligaseCount++;
                if(grid[row][col].organism->getTraits()->hasHelicase)
                    helicaseCount++;
                if(grid[row][col].organism->getTraits()->hasDNAPolymerase)
                    polymeraseCount++;
                if(grid[row][col].organism->getTraits()->hasCEP164)
                    cep164Count++;
            }
        }
    out<<"\t"<<totalCount<<"\t"<<chlorophyllACount<<"\t"<<chlorophyllBCount<<"\t"<<atpSynthaseCount<<"\t"<<topoIsomeraseCount<<"\t"<<primaseCount<<"\t"<<
        ligaseCount<<"\t"<<helicaseCount<<"\t"<<polymeraseCount<<"\t"<<cep164Count<<std::endl;
}

/*
    Organisms survive if they do not die and if they gather energy somehow, either by consumption of other organisms or by gathering resources
*/
int const MAX_ATP = 40000; //prevents overflow into the negatives
void Environment::naturalSelection()
{
    long double lowestOrganism = generationNumber;
    Location* lowest = 0;//for direct access
    unsigned int lowestRow = 0;
    unsigned int lowestCol = 0;
    for(unsigned int row=0;row<rows;row++)
        for(unsigned int col=0;col<cols;col++)
            if(grid[row][col].organism && grid[row][col].organism->getTraits()->turn < lowestOrganism)
            {
                lowest = &grid[row][col];
                lowestOrganism = lowest->organism->getTraits()->turn;
                lowestRow = row;
                lowestCol = col;
            }
    while(lowestOrganism<generationNumber)
    {
        ///organism at location acts
        //Photosynthesis
        if(lowest->organism->getTraits()->hasProtochlorophyllideReductase)
        {
        //    std::cout<<"Photosynthesis"<<std::endl;
            if(lowest->organism->getTraits()->ATP < MAX_ATP)
            {
                if(lowest->organism->getTraits()->hasChlorophyllBSynthase)
                {
                    if(lowest->organism->getTraits()->hasOEC)
                    {
                        atmosphere.oxygen += (2 * lowest->sunlight);
                        lowest->protonGradient += lowest->sunlight;
                        lowest->NADPH += (lowest->sunlight * 2/3);
                    }
                    else
                    {
                        lowest->protonGradient += lowest->sunlight / 2;
                    }
                }
                else
                {
                    if(lowest->organism->getTraits()->hasOEC)
                    {
                        atmosphere.oxygen += (2 * lowest->sunlight);
                        lowest->protonGradient += lowest->sunlight * 2/3;
                        lowest->NADPH += (lowest->sunlight * 4/9);
                    }
                    else
                    {
                        lowest->protonGradient += lowest->sunlight / 3;
                    }
                }
            }
        }
        if(lowest->organism->getTraits()->hasRubisco)
        {
            std::cout<<"RuBisCO"<<std::endl;
            if(lowest->organism->getTraits()->hasPEPCarboxylase)
            {
                if(lowest->organism->getTraits()->ATP>17 && lowest->NADPH>11)
                {
                    lowest->glucose+=2;
                    lowest->organism->getTraits()->ATP-=18;
                    lowest->NADPH-=12;
                }
                else
                {
                    if(lowest->organism->getTraits()->ATP>8 && lowest->NADPH > 5)
                    {
                        lowest->glucose++;
                        lowest->organism->getTraits()->ATP-=9;
                        lowest->NADPH-=6;
                    }
                }
            }
            else
            {
                if(lowest->organism->getTraits()->ATP>8 && lowest->NADPH > 5)
                {
                    lowest->glucose++;
                    lowest->organism->getTraits()->ATP-=9;
                    lowest->NADPH-=6;
                }
            }
        }
        //Glycolysis
        if(lowest->organism->getTraits()->hasHexokinase)
        {
            lowest->glucose6Phosphate += lowest->glucose;
            lowest->organism->getTraits()->ATP -= lowest->glucose;
            lowest->glucose = 0;
        }
        if(lowest->organism->getTraits()->hasPhosphofructokinase)
        {
            lowest->fructose16Bisphosphate += lowest->glucose6Phosphate;
            lowest->organism->getTraits()->ATP-=lowest->glucose6Phosphate;
            lowest->glucose6Phosphate = 0;
        }
        if(lowest->organism->getTraits()->hasAldolase)
        {
            lowest->G3P += lowest->fructose16Bisphosphate;
            lowest->DHAP += lowest->fructose16Bisphosphate;
            lowest->fructose16Bisphosphate = 0;
        }
        if(lowest->organism->getTraits()->hasTriosephosphateIsomerase)
        {
            lowest->G3P += lowest->DHAP;
            lowest->DHAP = 0;
        }
        if(lowest->organism->getTraits()->hasTriosephosphateDehydrogenase)
        {
            lowest->one3Bisphosphoglycerate += lowest->G3P;
            lowest->NADH += 2 * lowest->G3P;
            lowest->G3P = 0;
        }
        if(lowest->organism->getTraits()->hasPhosphoglycerateKinase)
        {
            lowest->threePhosphoglycerate += lowest->one3Bisphosphoglycerate;
            lowest->organism->getTraits()->ATP += 2 * lowest->one3Bisphosphoglycerate;
            lowest->one3Bisphosphoglycerate = 0;
        }
        if(lowest->organism->getTraits()->hasPyruvateKinase)
        {
            lowest->pyruvate += lowest->threePhosphoglycerate;

            lowest->organism->getTraits()->ATP += 2 * lowest->threePhosphoglycerate;
            lowest->threePhosphoglycerate = 0;
        }
        //Citric Acid Cycle
        if(lowest->organism->getTraits()->hasIsocitrateDehydrogenase)
        {
            lowest->alphaKetoglutarate += lowest->pyruvate;
            lowest->NADH += lowest->pyruvate;
            lowest->pyruvate = 0;
        }
        if(lowest->organism->getTraits()->hasAlphaKetoglutarateDehydrogenase)
        {
            lowest->succinylCoA += lowest->alphaKetoglutarate;
            lowest->NADH += lowest->alphaKetoglutarate;
            lowest->alphaKetoglutarate = 0;
        }
        if(lowest->organism->getTraits()->hasSuccinylCoASynthetase)
        {
            lowest->succinate += lowest->succinylCoA;
            lowest->organism->getTraits()->ATP +=  lowest->succinylCoA;
            lowest->succinylCoA = 0;
        }
        if(lowest->organism->getTraits()->hasSdhA)
        {
            lowest->fumerate += lowest->succinate;
            lowest->FADH2 += lowest->succinate;
            lowest->succinate = 0;
        }
        if(lowest->organism->getTraits()->hasMalateDehydrogenase)
        {
            lowest->NADH += lowest->fumerate;
            lowest->fumerate = 0;
        }
        //Electron Transport Chain
        if(lowest->organism->getTraits()->hasComplexI)
        {
            lowest->protonGradient += 3 * lowest->NADH;
            lowest->NADH = 0;
        }
        if(lowest->organism->getTraits()->hasSdhB && lowest->organism->getTraits()->hasSdhC && lowest->organism->getTraits()->hasSdhD)
        {
            lowest->protonGradient += 2 * lowest->FADH2;
            lowest->FADH2 = 0;
        }
        //Chemiosmosis
        if(lowest->organism->getTraits()->hasATPSynthase)
        {
            lowest->organism->getTraits()->ATP += lowest->protonGradient;
            lowest->protonGradient = 0;
        }
        //Movement
        if(lowest->organism->getTraits()->hasCEP164)
        {
            std::cout<<"Movement"<<std::endl;
            if(lowest->organism->getTraits()->hasSensoryRhodopsinII)
            {
                unsigned int mostSunlightRow = lowestRow;
                unsigned int mostSunlightCol = lowestCol;
                for(unsigned int row = lowestRow-1;row<lowestRow+2;row++)
                    for(unsigned int col = lowestCol-1;col<lowestCol+2;col++)
                    {
                        if(grid[row][col].sunlight > grid[mostSunlightRow][mostSunlightCol].sunlight)
                        {
                            mostSunlightRow = row;
                            mostSunlightCol = col;
                        }
                    }
                if(mostSunlightRow != lowestRow && mostSunlightCol != lowestCol)
                {
                    if(insert(lowest->organism,mostSunlightRow,mostSunlightCol))
                    {
                        lowest->organism->getTraits()->ATP -= 2;
                        lowest->organism = 0;
                        lowest = &grid[mostSunlightRow][mostSunlightCol];
                        lowestRow = mostSunlightRow;
                        lowestCol = mostSunlightCol;
                    }
                }
            }
            else
            {
                switch(rand()&7)
                {
                    case 0: if(insert(lowest->organism,lowestRow-1,lowestCol-1))
                            {
                                lowest->organism->getTraits()->ATP -= 2;
                                lowest->organism = 0;
                                lowestRow--;
                                lowestCol--;
                            }
                            break;
                    case 1: if(insert(lowest->organism,lowestRow,lowestCol-1))
                            {
                                lowest->organism->getTraits()->ATP -= 2;
                                lowest->organism = 0;
                                lowestCol--;
                            }
                            break;
                    case 2: if(insert(lowest->organism,lowestRow+1,lowestCol-1))
                            {
                                lowest->organism->getTraits()->ATP -= 2;
                                lowest->organism = 0;
                                lowestRow++;
                                lowestCol--;
                            }
                            break;
                    case 3: if(insert(lowest->organism,lowestRow-1,lowestCol))
                            {
                                lowest->organism->getTraits()->ATP -= 2;
                                lowest->organism = 0;
                                lowestRow--;
                            }
                            break;
                    case 4: if(insert(lowest->organism,lowestRow+1,lowestCol))
                            {
                                lowest->organism->getTraits()->ATP -= 2;
                                lowest->organism = 0;
                                lowestRow++;
                            }
                            break;
                    case 5: if(insert(lowest->organism,lowestRow-1,lowestCol+1))
                            {
                                lowest->organism->getTraits()->ATP -= 2;
                                lowest->organism = 0;
                                lowestRow--;
                                lowestCol++;
                            }
                            break;
                    case 6: if(insert(lowest->organism,lowestRow,lowestCol+1))
                            {
                                lowest->organism->getTraits()->ATP -= 2;
                                lowest->organism = 0;
                                lowestCol++;
                            }
                            break;
                    case 7: if(insert(lowest->organism,lowestRow+1,lowestCol+1))
                            {
                                lowest->organism->getTraits()->ATP -= 2;
                                lowest->organism = 0;
                                lowestRow++;
                                lowestCol++;
                            }
                            break;
                }//random motion
                lowest = &grid[lowestRow][lowestCol];
            }//no phototaxis
        }//if hasCEP164

        //Deterioration
        if((unsigned long long)(lowest->organism->getTraits()->turn + lowest->organism->getTraits()->step) > (unsigned long long) lowest->organism->getTraits()->turn)
        {
            unsigned int localOxygen = atmosphere.oxygen/(numLocations);
            lowest->organism->getTraits()->ATP -= localOxygen;
            atmosphere.oxygen -= localOxygen;
        }//time-based as opposed to turn-based
        //Death
        if(lowest->organism->getTraits()->ATP <= 0)
        {
            std::cout<<"Death"<<std::endl;
            lowestOrganism += 500;
            delete lowest->organism;
            lowest->organism = 0;
        }
        else
        {
            //calculates turn
            lowest->organism->getTraits()->turn += lowest->organism->getTraits()->step;
            lowestOrganism += lowest->organism->getTraits()->step;
            //Reproduction
            if(lowest->organism->getTraits()->ATP >= 5)
            {
                Organism* child = reproduce(lowestRow,lowestCol,lowest->organism->getTraits()->turn);
                if(place(child,lowestRow,lowestCol))
                {
                    lowest->organism->getTraits()->ATP -= 5;
                }
                else
                {
                    delete child;
                    child = 0;
                }
            }
        }
        //Finds next organism
        bool extinction = true;
        for(unsigned int row=0;row<rows;row++)
            for(unsigned int col=0;col<cols;col++)
                if(grid[row][col].organism)
                {
                    extinction = false;
                    if(grid[row][col].organism->getTraits()->turn<lowestOrganism)
                    {
                        lowest = &grid[row][col];
                        lowestOrganism = lowest->organism->getTraits()->turn;
                        lowestRow = row;
                        lowestCol = col;
                    }
                }
        if(extinction)
        {
            std::cout<<"Extinction"<<std::endl;
            lowestOrganism = generationNumber;
        }
    }//natural selection time loop
}
