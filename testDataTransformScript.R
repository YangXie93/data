
library(Biostrings)

Rcpp::sourceCpp('~/CompletnessTestData/makingBinCompletenessTestData.cpp') # hier den Pfad zu der Datei makingBinCompletenessTestData.cpp angeben
Rcpp::sourceCpp('~/sequenceToFasta.cpp') # hier den Pfad zu der Datei sequenceToFasta.cpp angeben

################### in diesem skript wird jede Datei aus dem Ordner einmal mit einer anderen zufaellig contaminiert

path = "~/testDaten/" # hier der Pfad zu deinen heruntergeladenen Genomen
run = 10 ############################## run ist zur unterscheidung falls zufaelliger weise ein paar zweimal genommen werden sollte (also bei jedem durchlauf aendern)

gens = dir(path)
datNms = gsub(".fna|.fasta","",gens)
gens = gens[grep(".*\\.fna|.*\\.fasta",gens)]

minCon = 50000
meanCon = 100000

gens = paste0(path,gens)

for(i in 1:length(gens)){
    
    rest = gens[-i]
    restNms = datNms[-i]
    
    comp = readDNAStringSet(gens[i])
    compLens = width(comp)
    compLen = sum(compLens)
    compNms = c(1:length(comp))
    
    index = sample(c(1:length(rest)),1)
    
    cont = rest[index]
    cont = readDNAStringSet(cont)
    contLens = width(cont)
    contLen = sum(contLens)
    contNms = c((length(comp) +1):(length(comp)+ length(cont)))
     
    lens = list(compLens,contLens)
    len = c(compLen,contLen)
    nms = list(compNms,contNms)
    
    contigs = mkContigs(lens,nms,len,minCon,meanCon,1,comp = c(0.7,0.8),cont = c(0.01,0.5),seed = run + i)
    
    compl = contigs[[1]][[1]][[length(contigs[[1]][[1]])-1]]/contigs[[1]][[1]][[length(contigs[[1]][[1]])]]
    if(length(contigs[[1]]) > 1){
        conta = contigs[[1]][[2]][[length(contigs[[1]][[2]])-1]]/contigs[[1]][[2]][[length(contigs[[1]][[2]])]]
        fastaName = paste0("testBins/",datNms[i],"_X_",restNms[index],"comp:",compl,"_cont:",conta,".",run,".fasta")
    }
    else{
        fastaName = paste0("testBins/",datNms[i],"_X_",restNms[index],"comp:",compl,".",run,".fasta")
    }
    
    compCont = list(comp,cont)
    
    x = 0
    for(j in 1:length(contigs[[1]])){
        if(contigs[[1]][[j]][[1]] > length(comp)){
            y = 2
            x = length(comp)
        }
        else{
            y = 1
            x = 0
        }
        for(n in seq(1,length(contigs[[1]][[j]])-2,3)){
            starts = contigs[[1]][[j]][[n+1]]
            widths = contigs[[1]][[j]][[n+2]] - starts
            nametag = names(compCont[[y]])[contigs[[1]][[j]][[n]] -x]
            sequenceToFastaConts(starts,widths,toString(compCont[[y]][contigs[[1]][[j]][[n]] -x]),fastaName,nametag)
        }
    }
}





