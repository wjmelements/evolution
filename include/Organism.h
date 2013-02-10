#ifndef ORGANISM
#define ORGANISM
#include "DNA.h"
#include <fstream>

typedef unsigned long long ull;

struct Traits
{
    ull isNotSequenced : 1;                        //true if the genes have been evaluated
    /*  Energy production notes:
                -the functions of several enzymbes have been combined for simplicity, including:
                    -Phosphoglucoisomerase and Phosphofrucokinase
                    -Enolase and Pyruvatekinase
                    -Citrate synthase, aconitase, and isocitrate dehydrogenase
                    -Fumerase and Malate dehydrogenase
                    -Complex I (1.6.5.3) was simplified to the seven core components of its NADH dehydrogenase subunit
    */

    //Glycolysis
    ull hasHexokinase : 1;                      //Glucose and ATP to Glucose 6-phosphate and ADP
    //917 http://www.uniprot.org/uniprot/P19367
    //8 Leu Thr Val Gly Pro Met Asp Ser
    ull hasPhosphofructokinase : 1;             //Glucose 6-phosphate and ATP to Fructose 1,6-bisphosphate and ADP
    //337 http://www.uniprot.org/uniprot/Q99ZD0
    //3 Arg Cys Ser
    ull hasAldolase : 1;                        //Fructose 1,6-bisphosphate to Glyceraldehyde 3-phosphate and Dihydroxyacetone phosphate
    //358 http://www.uniprot.org/uniprot/P0AB71
    //3 Thr Val Leu
    ull hasTriosephosphateIsomerase : 1;        //Dihydroxyacetone phosphate to Glyceraldehyde 3-phosphate
    //286 http://www.uniprot.org/uniprot/P60174
    //3 Phe Ser Gly
    ull hasTriosephosphateDehydrogenase : 1;    //Glyceraldehyde 3-phosphate and 2NAD+ to 1,3-Bisphosphoglycerate, 2 NADH
    //334 http://www.uniprot.org/uniprot/P04406
    //3 Ile Arg Ser
    ull hasPhosphoglycerateKinase : 1;            //1,3-Bisphosphoglycerate and 2 ADP to 3-Phosphoglycerate and 2 ATP
    //416 http://www.uniprot.org/uniprot/P07205
    //4 Arg Arg Thr Lys
    ull hasPyruvateKinase : 1;                  //Phosphoenolpyruvate and 2 ADP to Pyruvate and 2 ATP
    //216 http://www.uniprot.org/uniprot/Q9PU23
    //2 Leu Ser

    //Citric acid cycle
    ull hasIsocitrateDehydrogenase : 1;         //Oxaloacetate, Acetyl COA and NAD+ to α-Ketoglutarate and NADH
    //427 http://www.uniprot.org/uniprot/Q9ZH99
    //4 Ile Ser His Arg
    ull hasAlphaKetoglutarateDehydrogenase: 1;  //α-Ketoglutarate and NAD+ to NADH and Succinyl-CoA
    //1023 http://www.uniprot.org/uniprot/Q02218
    //8 Leu Glu Lys Ser Arg Thr His Arg
    ull hasSuccinylCoASynthetase : 1;             //Succinyl-CoA and ADP to Succinate and ATP
    //346 http://www.uniprot.org/uniprot/P53597
    //3 Arg Ile Gly
    ull hasSdhA : 1;                        //Succinate and FAD to Fumerate and FADH2
    //588 http://www.uniprot.org/uniprot/P0AC41
    //5 Ser Arg Pro Arg Trp
    ull hasMalateDehydrogenase : 1;             //Fumerate and NAD+ to NADH and Oxaloacetate
    //338 http://www.uniprot.org/uniprot/P40926
    //3 Leu Pro Ala

    //Chemiosmosis
    ull hasATPSynthase : 1;                     //proton to ATP
    //529 http://www.uniprot.org/uniprot/P06576
    //5 Lys His Asp Val Glu

    //Electron Transport Chain
    ull hasComplexI : 1;                        //NADH to 3 protons; requires all seven chains
        ull hasChain1 : 1;
        //318 uniprot
        //3 Ser Arg Asn
        ull hasChain2 : 1;
        //347 uniprot
        //3 Asn Ser Arg
        ull hasChain3 : 1;
        //115 uniprot
        //1 Pro
        ull hasChain4 : 1;
        //459 uniprot
        //4 Met Leu Thr Gly
        ull hasChain4L : 1;
        //63 uniprot
        //1 Val
        ull hasChain5 : 1;
        //603 uniprot
        //5 Ile Ile Gln Leu Ala
        ull hasChain6 : 1;
        //174 uniprot
        //2 Gly Ala
    ull hasSdhB : 1;                        //FADH2 to 2 protons, given ShdC and ShdD
    //238 http://www.uniprot.org/uniprot/P07014
    //2 Ser Cys
    ull hasSdhC : 1;                        //FADH2 to 2 protons, given ShdB and ShdD
    //129 http://www.uniprot.org/uniprot/P69054
    //2 Arg Ile
    ull hasSdhD : 1;                        //FADH2 to 2 protons, given ShdB and ShdC
    //115 http://www.uniprot.org/uniprot/P0AC45
    //1 Arg

    //Photosynthesis
    ull hasProtochlorophyllideReductase : 1;    //makes chlorophyll a
    //398 http://www.uniprot.org/uniprot/Q9SDT1
    //4 Trp Pro Glu Gly
    ull hasChlorophyllBSynthase : 1;            //chlorophyll b given chlorophyll a
    //541 http://www.uniprot.org/uniprot/Q8S7E1
    //5 Val Arg Cys Pro Gly
    ull hasOEC : 1;                             //Photosystem II, allows extra electron for electron transfer to NADPH
        ull hasOECEP1 : 1;
        //291 P12853
        //3 Cys Ser Phe
        ull hasOECEP2 : 1;
        //245 P11471
        //2 Ser Trp
        ull hasOECEP3 : 1;
        //199 P12852
        //2 Leu Pro
    ull hasRubisco : 1;                         //ATP and 6CO2 to Glucose, inefficient; requires RubiscoSmall and RubiscoLarge
        ull hasRubiscoSmall : 1;                    //Component of RuBisCO
        //180 http://www.uniprot.org/uniprot/P69249
        //2 Ser Val
        ull hasRubiscoLarge : 1;                    //Component of RuBisCO
        //477 http://www.uniprot.org/uniprot/P00876
        //4 Gln His Ala Trp
    ull hasPEPCarboxylase : 1;                  //makes Rubisco efficient
    //883 http://www.uniprot.org/uniprot/Q8Z307
    //7 Pro Pro Pro Asp Val Ile Val

    //Motion
    ull hasCEP164 : 1;                      //enables motion
    //1460 http://www.uniprot.org/uniprot/Q9UPV0
    //12 Thr Cys Met Arg Gln His Gln Cys Val His Ile Ala

    //Taxis and Kinesis
    ull hasSensoryRhodopsinII : 1;  //Positive phototaxis
    //239 http://www.uniprot.org/uniprot/P42196
    //2 Thr Lys

    //Replication
    ull hasTopoisomerase : 1;
    //765 uniprot
    //6 Val Asp Gly Glu Leu Arg
    ull hasDNAPrimase : 1;
    //581 uniprot
    //5 Arg Arg Asp Ser Leu
    ull hasDNALigase : 1;
    //919 uniprot
    //8 Tyr Trp Gly Ile Arg Ile Leu Tyr
    ull hasHelicase : 1;
    //720 http://www.uniprot.org/uniprot/P03018
    //6 Arg Thr Arg Ser Asn Pro
    ull hasDNAPolymerase : 1;
    //928 http://www.uniprot.org/uniprot/P00582
    //8 Val Pro Pro Leu His Ser Thr Asn

    //Calculated
    long double accuracy;
    unsigned int cost;                                   //in ATP
    unsigned int haploidNumber;
    //Natural Selection
    long double turn;
    long double step;
    int ATP;//must be signed
};

//http://www.ncbi.nlm.nih.gov/pmc/articles/PMC248711/pdf/jbacter00372-0173.pdf
//http://www.ncbi.nlm.nih.gov/pmc/articles/PMC202969/

class Organism
{
    public:
        Organism(long double turn, DNA* dna=0);
        Organism(std::ifstream& in);
        virtual ~Organism();

        //IO
        void printDNA();
        void printRNA();
        void printTraits();
        void print();
        void saveToFile(std::ofstream& file);
        DNA* replicateDNA(double accuracy);
        void sequence();
        Traits* getTraits();
    protected:
    private:
        //hereditary
        DNA* dna;//all DNA
        DNA* rna;//all mRNA
        Traits traits;

        void translate(DNA* mRNA);//translation

/*            //the minimum temperature in which the organism can function
            uint minTemp;
            //the maximum temperature in which the organism can function
            uint maxTemp;
            //determines whether or not the organism can reproduce
            bool canReproduce;
            //determines capability of survival
            uint size;
            //determines energy cost of survival
            uint resourcesConsumed;
            //can gain energy by consuming light
            bool canPhotosynthesize;
            //can move with direction towards light
            bool hasPhototaxis;
            //moves randomly to optomize light
            bool hasPhotokinesis;
*/
};

#endif // ORGANISM
