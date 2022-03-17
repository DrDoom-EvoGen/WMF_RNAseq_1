TRINOTATE=./Trinotate_softwares/Trinotate-Trinotate-v3.2.2


${TRINOTATE}/auto/autoTrinotate.pl --Trinotate_sqlite ./Trinotate_Greg.sqlite \
    --transcripts ./Trinity.fasta --gene_to_trans_map ./TrinotateNamesFinale.txt  \
    --conf ./conf.txt --CPU 10
