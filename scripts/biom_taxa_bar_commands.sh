#!/bin/bash

# TO CHANGE TAXONOMIC LEVEL GO INTO SCRIPT AND CHANGE 

./biom_taxa_bar.py -i ../data/drosophila/5_5_filter_samples_from_otu_table.biom -m ../data/drosophila/drosophila_map.txt -o ../data/drosophila/taxa_bar_chart.pdf
./biom_taxa_bar.py -i ../data/hominid/5_5_filter_samples_from_otu_table.biom -m ../data/hominid/map.txt -o ../data/hominid/taxa_bar_chart.pdf
./biom_taxa_bar.py -i ../data/mosquito/5_5_filter_samples_from_otu_table.biom -m ../data/mosquito/mosquito_map.txt -o ../data/mosquito/taxa_bar_chart.pdf
./biom_taxa_bar.py -i ../data/nasonia/5_5_filter_samples_from_otu_table.biom -m ../data/nasonia/nasonia_map.txt -o ../data/nasonia/taxa_bar_chart.pdf
./biom_taxa_bar.py -i ../data/peromyscus/5_5_filter_samples_from_otu_table.biom -m ../data/peromyscus/map.txt -o ../data/peromyscus/taxa_bar_chart.pdf
