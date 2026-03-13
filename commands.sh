conda install pandas matplotlib mygene -c bioconda -y

python perform_analysis.py --features-file-k562 cv_k562__extracted_features_for_GO_enrichment --features-file-hepg2 cv_hepg2__extracted_features_for_GO_enrichment --features-file-sknsh cv_sknsh__extracted_features_for_GO_enrichment --rna-seq-file-k562 K562/RNA-seq/ENCFF928NYA.tsv --rna-seq-file-hepg2 HepG2/RNA-seq/ENCFF863QWG.tsv 