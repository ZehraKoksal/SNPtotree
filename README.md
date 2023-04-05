# SNPtotree - A software to resolve the phylogeny of haploid SNPs

Zehra Köksal¹, Leonor Gusmão², Claus Børsting¹, Vania Pereira¹

¹Section of Forensic Genetics, Department of Forensic Medicine, Faculty of Health and Medical Sciences,
University of Copenhagen, Denmark

²DNA Diagnostic Laboratory (LDD), State University of Rio de Janeiro (UERJ), Brazil

### 1) About SNPtotree
SNPtotree determines the hierarchical order of biallelic haploid variants - variants on haploid markers, which do not undergo recombination. Even for sequencing data with high percentage of missing information, SNPtotree reliably generates a phylogenetic tree without error-prone manual sorting. This is unique to SNPtotree, when compared with alternative methods, like maximum likelihood (ML) based approaches.

SNPtotree allows the creation of phylogenetic trees of variants based on genetically similar and dissimilar sequences with low or high amount of missing data. The algorithm conducts comparisons between the allelic states (ancestral "A" or derived "D") of all variant pairs to infer their relationships and generate a phylogenetic tree. The tree is more accurate and complete, the more accurate and complete the input sequencing data. Variants that predict contradictory pairwise relationships or ambiguous positions in the tree are ignored

### 2) Installation
Operating system: linux

Type in the shell (linux):
```
git clone https://github.com/ZehraKoksal/SNPtotree.git
cd SNPtotree/
SNPtotree.py -h
```

### 3) Running SNPtotree
The software is simple and requires only one input file and creates one main output file, the phylogenetic tree. Three additional output files are optional and give background information on your genetic variants.

#### a) Command
```
python SNPtotree.py /path_to_input_file.csv/ /path_to_output_folder/ 
```
optional:
```
-contradictory_variants /path_to_input_data_folder/ 
-ambiguous_variants /path_to_input_data_folder/
-metadata_individuals /path_to_input_data_folder/  
```
#### b) Input file
The user is required to provide the path to the input file in the _.csv_ format. The input file contains the ancestral **A** or derived **D** allelic state or missing information **X** for each **_polymorphic variant_** in a tab-separated format. The rows present variants, the columns the individuals.
The header row should present the individuals' labels and the first (index) column the variant names.

<img src="/Images/inputfile_snptotree.png" alt="Input file style" width="580"/>

The allelic states ancestral and derived of the most used model organisms are reported in public repositories.

Human Y-chromosome: [Y-DNA Haplogroup Tree 2019-2020](https://isogg.org/tree/)
Human mt DNA:

#### c) Output file(s)



### 4) Additional information and Contact
More information on the software are available in our publication: XX (link)

For reporting bugs, comments or questions, you are welcome to contact zehra.koksal@sund.ku.dk.

### 5) Referencing

Please cite:



