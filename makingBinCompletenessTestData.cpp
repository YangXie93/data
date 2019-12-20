#include <iostream>
#include <vector>
#include <list>
#include <iterator>
#include <algorithm>
#include <Rcpp.h>
#include <time.h>
#include <random>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <numeric>
#include <string>
#include <cmath>

//[[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

//###################### bestimmung der Anzahl und laenge der Contigs ####################################################

std::vector<int> randomContigs(int minContigLength,int meanContigLength,int covering,std::string distr = "normal",int seed = 0){
    if(covering <= 0){
        return std::vector<int> {};
    }
    if(covering <= minContigLength){
        return std::vector<int> (1,covering);
    }
    std::default_random_engine generator;
    generator.seed(seed);
    std::normal_distribution<double> normDist(meanContigLength,sqrt(meanContigLength));
    std::poisson_distribution<int> poisDist(meanContigLength);
    std::geometric_distribution<int> expDist;
    std::uniform_int_distribution<int> uniformDist(minContigLength,meanContigLength*2);
    std::vector<int> res;
    int tmp = 0;
    while(covering > 0){
        if(distr == "poisson" || distr == "normal" || distr == "exponential" || distr == "uniform"){
            if(distr == "poisson"){
                tmp = poisDist(generator);
                poisDist.reset();
            }
            if(distr == "normal"){
                tmp = (int) normDist(generator);
                normDist.reset();
            }
            if(distr == "exponential"){
                
                tmp = minContigLength + (expDist(generator));
                expDist.reset();
            }
            if(distr == "uniform"){
                tmp = uniformDist(generator);
                uniformDist.reset();
            }
        }
        else{
            Rcerr << "die gewählte verteilung steht nicht zur wahl." << std::endl << "verfügbar sind: \"normal\", \"poisson\",\"exponential\" und \"uniform\" " << std::endl;
            Rcerr << "Das Programm wird jetzt mit einer Normal verteilung fortgesetzt." << std::endl;
            tmp = (int) normDist(generator);
            normDist.reset();
        }
        if(tmp >= covering){
            if(covering >= minContigLength){
                res.push_back(covering);
            }
            covering = 0;
        }
        else{
            if(tmp > 0){
                covering -= tmp;
                res.push_back(tmp);
            }
        }
    }
    return res;
}

//############ Setzung der Abstaende zwischen den Contigs ##########

//[[Rcpp::export]]
std::vector<int> randomSpaces(int numSp,int free,int seed){
    if(free <= 0){
        return std::vector<int> (numSp,0);
    }
    std::vector<int > spaces;
    spaces.reserve(numSp);

    std::default_random_engine generator;
    generator.seed(seed);

    std::vector<int> t;
    t.reserve(numSp);
    double sum = 0.0;
    int flex;
    int i;
    int sum2 = 0;
    int* elem = t.data();

    for(i = 0;i < numSp;i++){
        flex = generator() % 10000;
        t.push_back(flex);
        sum += flex;
    }

    for(i = 0;i < numSp;i++){
        *(elem +i) =  (int)(free*(*(elem +i)/sum));
        sum2 += *(elem +i);
    }
    if(sum2 > free){
        while(sum2 > free){
            i = generator() % numSp;
            if(*(elem +i) > 0){
                (*(elem +i))--;
                sum2--;
            }
            i++;
        }
    }
    if(sum2 < free){
        while(sum2 < free){
            i = generator() % numSp;
            (*(elem +i))++;
            sum2++;
        }
    }
    std::vector<int>::iterator it;
    while((int) t.size() > 0){
         it =  t.begin() +(generator() % ((int) t.size()));
         spaces.push_back( *(it));
         t.erase(it);
    }
    return spaces;
}

int whichToSmall;

//################# Auswahl von welchem Chromosom eines Genoms wie viele Basen fuer die Contigerstellung genommen werden ####################################

//[[Rcpp::export]]
std::vector<int> fromWhichHowMany(int minContigLength, int meanContigLength,int totalLength,std::vector<int> lengths,int needed,int seed,bool debugInfo = false){
    
    std::default_random_engine generator;
    generator.seed(seed);
    
    std::normal_distribution<double> distr (meanContigLength,sqrt(meanContigLength));
    
    std::vector<int> res (lengths.size(),0);
    int share;
    int at;
    int tmp = needed;
    int drawnNr;
    
    std::vector<int> chromsGreaterMin;
    std::vector<int> accession;
    std::vector<int> relFreq;
    relFreq.reserve(totalLength*2/minContigLength);
    std::vector<int>::iterator rF; 

    for(int i = 0;i < (int) lengths.size();i++){
        if(lengths[i] > minContigLength){
            chromsGreaterMin.push_back(lengths[i]);
            accession.push_back(i);
        }
        else{
            totalLength -= lengths[i];
        }
    }
    
    if(totalLength < needed){
        needed = totalLength;
        if(debugInfo)
            Rcout << whichToSmall << ": die angefragte Menge ist nicht erreichbar\n";
    }
    
    
    
    for(int i = 0;i < chromsGreaterMin.size();i++){
        share = round(chromsGreaterMin[i]/minContigLength);
        for(int j = 0;j < share;j++ ){
            relFreq.push_back(i);
        }
        
    }
    
    bool noMore = false;
    int count;
    while(needed >= 0 && !noMore){
        at = generator() %relFreq.size();
        rF = next(relFreq.begin(),at);
        
        count = 0;
        while(chromsGreaterMin[(*rF)] < minContigLength && count < relFreq.size()){
            rF++;
            count++;
            if(rF == relFreq.end()){
                rF = relFreq.begin();
            }
            if(count == relFreq.size()){
                noMore = true;
            }
        }
        if(!noMore){
            drawnNr = round(distr(generator));
            
            if(drawnNr > chromsGreaterMin[(*rF)]){
                drawnNr = chromsGreaterMin[(*rF)];
            }
            chromsGreaterMin[(*rF)] -= drawnNr;
            
            res[accession[(*rF)]] += drawnNr;

            relFreq.erase(rF);
            needed -= drawnNr;
        }
    }
    if(debugInfo){
        Rcout << "gewollte Anzahl: " << tmp << " gezogene Anzahl: " << tmp - needed << std::endl << std::endl;
    }
    
    return res;
}

int minL;

bool isBiggerMinL(int i){
    return(i > minL);
}


//############################ Erstellen von Contigs zufaelliger Laenge von einem zufaelligem Genom, kontaminiert mit Contigs zufaelliger Laenge eines anderen zufaelligen Genoms #############

//[[Rcpp::export]]

std::list<std::list<std::list<std::vector<int> > > > mkContigs(std::list<std::vector<int> >& lengths,std::list<std::vector<int> > &IDs,std::vector<int>& lengthSums,int minContigLength,int meanContigLength,int number,std::vector<double> comp,std::vector<double> cont,int seed = 0,std::string distr = "normal",bool debugInfo = false){
    whichToSmall = 0;
    minL = minContigLength;
    
    std::default_random_engine generator;
    generator.seed(seed);
    
    std::list<std::list<std::list<std::vector<int> > > > res;
    std::vector<int>::iterator totLen;
    std::vector<int> baseNrs;
    std::vector<int> indicies;
    std::vector<int> contigs;
    std::vector<int> spaces;
    std::vector<int>::iterator co;
    std::vector<int>::iterator sp;
    std::vector<int> starts;
    std::vector<int> ends;
    std::vector<int> chromBaseNrs;
    std::vector<int>::iterator which;
    std::vector<int>::iterator max;
    std::list<std::vector<int> >::iterator tester;
    
    bool swtch;
    int index;
    int l;
    int n;
    int j;
    int count;
    int accuContigs;
    int at = 0;
    int contigSum;
    double partCovered;
    bool contIsNull;
    bool justZero = false;
    int testTmp;
    
    //--------- Prüfen von comp und cont ---------------
    
    if(cont[1] > comp[0]){
        cont[1] = comp[0];
    }
    if(cont[0] == cont[1]){
        cont[0] = cont[1] -0.01;
    }
    if(comp[0] == comp[1]){
        comp[0] = comp[1] -0.01;
    }
    
    //--------------------------------------------------
    
    for(int i = 0; i < number; i++){
        
        whichToSmall++;
        count = 0;
        contIsNull = false;
        
        //------------------------ Aussuchen des bin angebenden Genoms und der completness ----------------------
        
        
        partCovered = ((comp[0] *100) + (generator() % (int)(round((comp[1]-comp[0]) *100))))/100.0;
        
        testTmp = (generator() % lengthSums.size());
        which = next(lengthSums.begin(),testTmp);
        tester = next(lengths.begin(),testTmp);
        

        while(((*which)* partCovered) < minContigLength && count < (int)lengthSums.size()){
            which++;
            count++;
            tester++;

            if(which != lengthSums.end() && (find_if((*tester).begin(),(*tester).end(),isBiggerMinL)) == (*tester).end()){
                which++;
                tester++;
                count++;
            }
            
            if(which == lengthSums.end()){
                which = lengthSums.begin();
                tester = lengths.begin();
            }
        }
        
        if(count == (int)lengthSums.size() && (*which) < minContigLength){
            Rcerr << "Die mindest Länge ist zu groß für den Datensatz" << std::endl;
            return res;
        }
        totLen = which;
        baseNrs.push_back((int) ((*totLen) *partCovered));
        
        //-------------------------------------------------------------------------------------------------------
        
        //------------------------- Aussuchen des genoms zur contamination --------------------------------------
        
        double contPart = ((cont[0] *100) + (generator() % ((int)(cont[1]*100) - (int)(cont[0] *100))))/100.0;
        if((*totLen) *contPart <= minContigLength){
            if(contPart != 0){
                baseNrs.push_back(minContigLength);
            }
            else{
                contIsNull = true;
            }
        }
        else{
            baseNrs.push_back((*totLen) *contPart);
        }
        if(!contIsNull){
            count = 0;
            testTmp = (generator() % lengthSums.size());
            which = next(lengthSums.begin(),testTmp);
            tester = next(lengths.begin(),testTmp);
            max = lengthSums.begin();
            while(((*which) < *prev(baseNrs.end()) && count < (int)lengthSums.size()) || which == totLen){
                if(which == lengthSums.end()){
                    which = lengthSums.begin();
                    tester = lengths.begin();
                }
                else{
                    
                    if(which > max && which != totLen){
                        max = which;
                    }
                    which++;
                    count++;
                    tester++;

                    if(which != lengthSums.end() && (find_if((*tester).begin(),(*tester).end(),isBiggerMinL)) == (*tester).end()){
                        which++;
                        tester++;
                        count++;
                    }
                }
            }
            if(((*which) < *prev(baseNrs.end()) && count == (int)lengthSums.size()) || which == lengthSums.end() || which == totLen){
                baseNrs.pop_back();
            }
        }
        //--------------------------------------------------------------------------------------------------------
        
        index = distance(lengthSums.begin(),which);
        indicies.push_back(distance(lengthSums.begin(),totLen));
        indicies.push_back(index);
        
        std::list<std::list<std::vector<int> > > r;
        
        for(n = 0; n < (int) baseNrs.size();n++){       // für completeness und contamination 
            
            std::list<std::vector<int> > re;
            accuContigs = 0;
            if(debugInfo){
                int accuChromBaseNrs = accumulate(chromBaseNrs.begin(),chromBaseNrs.end(),0);
                if(n == 0){
                    Rcout << i+1 << std::endl << std::endl;
                    Rcout << "completeness: " << partCovered  << std::endl << std::endl;
                    Rcout << "Gesamtlaenge: " << (*totLen) << std::endl;
                }
                else{
                    Rcout << "contamination: " << contPart << std::endl << std::endl; 
                }
            }
            
            chromBaseNrs = fromWhichHowMany(minContigLength,meanContigLength,(*next(lengthSums.begin(),indicies[n])),(*next(lengths.begin(),indicies[n])),baseNrs[n],seed+i,debugInfo);
            justZero = (accumulate(chromBaseNrs.begin(),chromBaseNrs.end(),0) == 0);
            
            if(justZero && debugInfo){
                Rcout << *((*next(IDs.begin(),indicies[n])).begin()) << " " << i << " " << indicies[n] << " " << contPart << std::endl;
            }
            if(!justZero){
                
                for(j = 0; j < (int) chromBaseNrs.size();j++){      // für alle chromosomen
                    contigs = randomContigs(minContigLength,meanContigLength,chromBaseNrs[j],distr,seed+j+n+1);
                    if(debugInfo){
                        Rcout << " von: " << (*next(IDs.begin(),indicies[n]))[j] << " mit Laenge: " << (*next(lengths.begin(),indicies[n]))[j] << " nimm: " << chromBaseNrs[j] << std::endl;
                        Rcout << "Anzahl der Contigs " << contigs.size() << std::endl;
                    }
                    if( contigs.size() > 0){
                        
                        contigSum = accumulate(contigs.begin(),contigs.end(),0);
                        
                        spaces = randomSpaces((int)contigs.size() +1,(*next(lengths.begin(),indicies[n]))[j] -contigSum,seed);
                        co = contigs.begin();
                        sp = spaces.begin();
                        starts.reserve(contigs.size());
                        ends.reserve(contigs.size());
                        
                        swtch = true;
                        at = 0;
                        for(int m = 0;m < (int)contigs.size();m++){     // für alle Contigs
                            at += (*sp);
                            starts.push_back(at+1);
                            at += (*co);
                            sp++;
                            co++;
                            ends.push_back(at);
                        }           // ende alle Contigs
                        
                        at += (*sp);
                        
                        re.push_back(std::vector<int> {(*next(IDs.begin(),indicies[n]))[j]});
                        re.push_back(starts);
                        re.push_back(ends);
    
                        accuContigs += contigSum;
                        
                        starts.clear();
                        ends.clear();
                    }
                }       // ende alle Chromosomen
                re.push_back(std::vector<int> {accuContigs});
                re.push_back(std::vector<int> {(*totLen)});
                r.push_back(re);
            }
            if(debugInfo)
                Rcout << std::endl << std::endl;
        }       // ende completness und contamination
        if(debugInfo)
            Rcout << "-------------" << std::endl << std::endl;
        if(!justZero){    
            res.push_back(r);
        }
        baseNrs.clear();
        indicies.clear();
    }
    return res;
}

//##################### wende die mkContigs function auf ein einzelnes Genom an ########################################

//[[Rcpp::export]]

std::list<std::list<std::vector<int> > > singleGenomeMkContigs(std::vector<int> lengths,std::vector<int> IDs,double comp,int minContigLength,int meanContigLength,int seed = 0,std::string distr = "normal",bool debugInfo = false){
    
    std::list<std::list<std::vector<int> > > res;
    
    int lengthSum = accumulate(lengths.begin(),lengths.end(),0);
    
    std::list<std::vector<int> > wrapLengths;
    wrapLengths.push_back(lengths);
    
    std::vector<int> wrapLengthSum;
    wrapLengthSum.push_back(lengthSum);
    
    std::list<std::vector<int> > wrapIDs;
    wrapIDs.push_back(IDs);
    
    std::list<std::list<std::list<std::vector<int> > > > tmp = mkContigs(wrapLengths,wrapIDs,wrapLengthSum,minContigLength,meanContigLength,1,{comp,comp+0.01},{1,0},seed,distr,debugInfo);
    
    if(tmp.size() > 0){
        res = (*(tmp.begin()));
    }
    
    return res;
}

//############ vergleiche die Vollen PFAM counts aller vorhandenen Genome mit dem eines zu untersuchendem Bins ##########

//[[Rcpp::export]]

List comparePfamCounts(List refList, std::vector<int> orgVec, std::vector<int> query){
    
    std::vector<int> missing (orgVec.size(),0);
    std::vector<int> exceeding (orgVec.size(),0);
    std::vector<int> total (orgVec.size(),0);
    std::vector<bool> isIn (orgVec.size(),true);
    
    std::vector<int>::iterator missingIt = missing.begin();
    std::vector<int>::iterator exceedingIt = exceeding.begin();
    std::vector<int>::iterator totalIt = total.begin();
    std::vector<bool>::iterator isInIt = isIn.begin();
    
    std::vector<int>::iterator queryIt = query.begin();
    
    List::iterator orgListIt;
    std::vector<int>::iterator pfamCountVecIt;
    std::vector<int>::iterator orgListVecIt;
    
    std::vector<int>::iterator orgIt;
    std::vector<int> countVec;
    std::vector<int> orgV;
    
    
    int dictSize = *std::max_element(orgVec.begin(),orgVec.end());
    std::vector<int> dict (dictSize,-1);
    
    for(int t = 0;t < orgVec.size();t++){
        dict[orgVec[t]] = t;
    }
    
    int i = 0;
    
    
    for(List::iterator refListIt = refList.begin(); refListIt != refList.end(); refListIt++){
        
        List refListList = as<List>((*refListIt));
        List::iterator orgL_u_countV = refListList.begin();
        
        if(refListList.size() == 0){
            if(*next(queryIt,i) != 0){
                for(int j =0 ;j < exceeding.size();j++){
                    *next(exceedingIt,j) += *next(queryIt,i);
                }
            }
        }
        else{
            countVec = as<std::vector<int> >((*orgL_u_countV));
            orgL_u_countV++;
            List orgL = as<List>(*orgL_u_countV);
            
            pfamCountVecIt = countVec.begin();
            orgListIt = orgL.begin();
            int tmp;
            
            for(int j = 0; j < countVec.size();j++){
                
                orgV = as<std::vector<int> >((*orgListIt));
                
                for(int n = 0; n < orgV.size();n++){
                    
                    tmp = (*pfamCountVecIt) - *next(queryIt,i);
                    
                    int dist = dict[orgV[n]];
                    
                    
                    if(tmp < 0){
                        *next(exceedingIt,dist) += -1* tmp;
                    }
                    else{
                        *next(missingIt,dist) += tmp;
                    }
                    *next(totalIt,dist) += (*pfamCountVecIt);
                }
                
                
                pfamCountVecIt++;
                orgListIt++;
            }
            
        }
        
        
        i++; 
    }
    
    List res = List::create(Named("total") = total,Named("missing") = missing,Named("exceeding") = exceeding);
    
    return res;
}
