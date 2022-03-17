# Build database
../Trinotate_softwares/Trinotate-Trinotate-v3.2.2/admin/Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate_Greg

## BlastDB
../Trinotate_softwares/ncbi-blast-2.7.1+/bin/makeblastdb -in uniprot_sprot.pep -dbtype prot

## Pfam
gunzip Pfam-A.hmm.gz
../Trinotate_softwares/hmmer-3.2/src/hmmpress Pfam-A.hmm
