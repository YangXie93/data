
# create fasta file with random sequences using different alphabets

alphabet = c("A","T","G","C")

alphabet = c("A","T","G")

alphabet = c("A","T","G","C","Y")


file = "~/Z.fasta"
n = FALSE

for(i in 1:sample(100:150,1)){
    
    nr = sample(10000:20000,1)
    if(i > 1){
        n = TRUE
    }
    write(paste0(">",i,"\n",paste(sample(alphabet,nr,replace = TRUE),collapse = ""),"\n"),file,append = n)
}





