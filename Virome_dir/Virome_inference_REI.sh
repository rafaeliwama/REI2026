## requires: GetVirome.sh, Accession_list.txt, bbmap


## Settinng up pipeline requirements

mkdir virome_databases
wget -P virome_databases https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec
wget -P virome_databases https://github.com/biocore/sortmerna/releases/download/v4.3.4/database.tar.gz
tar -xvzf virome_databases/database.tar.gz
wget -P virome_databases https://kaiju-idx.s3.eu-central-1.amazonaws.com/2024/kaiju_db_viruses_2024-08-15.tgz
tar -xvzf virome_databases/kaiju_db_viruses_2024-08-15.tgz



## Virome inferece for paired-ended transcriptomes

./GetVirome.sh Accession_list_paired.txt

## generate kaiju tables
for file in *.kaiju; do kaiju-addTaxonNames -r superkingdom,phylum,class,order,family,genus,species -t virome_databases/nodes.dmp -n virome_databases/names.dmp -i $file -o $file.names; done

# SraRunInfo.csv is a csv containing information for all SRAs used in the study 
python3 KaijuToDF.py virome_databases/nodes.dmp virome_databases/names.dmp REI_SraRunInfo.csv genus


# generate counts of viral reads per sample
for f in *.kaiju.names; do
    base="${f%.kaiju.names}"
    count=$(grep -c '^C' "$f")
    echo "$count,$base" >> Kaiju_total_reads.txt
done
