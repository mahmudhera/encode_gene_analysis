# encode_gene_analysis
```bash
conda install pandas matplotlib mygene -c bioconda -y

python perform_analysis.py --features-file cv_k562__extracted_features_for_GO_enrichment --rna-seq-file-k562 K562/RNA-seq/ENCFF928NYA.tsv --rna-seq-file-hepg2 HepG2/RNA-seq/ENCFF863QWG.tsv --output-file k562_features_to_TPMs_in_k562_and_hepg2


```