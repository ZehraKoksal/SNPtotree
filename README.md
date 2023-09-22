# SNPtotree - A software to resolve the phylogeny of SNPs on non-recombining DNA

Zehra Köksal¹, Leonor Gusmão², Claus Børsting¹, Vania Pereira¹

¹Section of Forensic Genetics, Department of Forensic Medicine, Faculty of Health and Medical Sciences,
University of Copenhagen, Denmark

²DNA Diagnostic Laboratory (LDD), State University of Rio de Janeiro (UERJ), Brazil

### 1) About SNPtotree
SNPtotree determines the hierarchical order of biallelic variants on non-recombining DNA. Even for sequencing data with high percentage of missing information, SNPtotree is able to generate a phylogenetic tree without error-prone manual sorting. This is unique to SNPtotree, when compared with alternative methods, like maximum likelihood (ML) based approaches.

SNPtotree allows the creation of phylogenetic trees of variants based on sequences from very high to very low similarity and amounts of missing data. The algorithm conducts comparisons between the allelic states (ancestral "A" or derived "D") of all variant pairs to infer their relationships and generate a phylogenetic tree. The tree is more accurate and complete, the more accurate and complete the input sequencing data. Variants that predict contradictory pairwise relationships or ambiguous positions in the tree are ignored.

For more in-depth explanations on the algorithms used and information on the tool, please see: [[Click here!](https://www.mdpi.com/2073-4425/14/10/1837)

### 2) Installation
Operating system: Linux

Type in the shell:
```
git clone https://github.com/ZehraKoksal/SNPtotree.git
cd SNPtotree/
python SNPtotree.py -h
```


### 3) Algorithm and Commands
SNPtotree uses a simple algorithm, which conducts pairwise comparisons of the allelic states of all variants. These are then combined to the final hierarchical tree order.

The required information in the shell are the paths to the _input_ and _main output_ files, which are explained below:

```
python SNPtotree.py path_to_input_file.csv/ path_to_output_folder/ 
```
Additionally, there are three optional output files that are explained below:
```
python SNPtotree.py path_to_input_file.csv/ path_to_output_folder/ -contradictory_variants -ambiguous_variants -metadata_individuals
```


#### a) Input file
The user is required to provide the path to the input file in the _.csv_ format. The input file contains the ancestral **A** or derived **D** allelic state or missing information **X** for each **_polymorphic variant_** in a tab-separated format. The rows present variants, the columns individuals.
The header row should present the individuals' labels and the first (index) column the variant names.

<img src="/Images/inputfile_snptotree.png" alt="Input file style" width="700"/>

The allelic states "ancestral" and "derived" of the most used model organisms are reported in public repositories. For novel SNPs or for not well investigated organisms without already reported relevant SNP information, the ancestral and derived allelic states have to be identified by the user. The ancestral allele is found in a common ancestor of the group of analysed individuals. Thus, it is helpful to conduct sequence alignments to a common ancestor rather than an arbitrarily selected reference genome, e.g. GRCh38 for humans.


#### b) Main Output Files: Phylogenetic Tree

SNPtotree generates the phylogenetic tree in two file formats: in a tab-separated csv file and a traditional phyloxml file.

The csv output file is to be read from left to right. Downstream variants are located in the cells below and to the right. In this example variants b, d, e and f are downstream of variant a, and variant c is downstream of variant b.
Variants in parallel branches (sister clades) within a clade are presented in a column: Variants g and h are parallel to each other. Not separable variants based on the available data are presented in one cell divided by a comma, like variants i and j.

<img src="/Images/output_phyltree.png" alt="Input file style" width="700"/>

The phyloxml output file contains annotated branches and nodes and can be viewed in phylogenetic tree visualization tools that support phyloxml format, e.g. the interactive Tree Of Life (iTOL). This tree contains the same information as the csv output file.

<img src="/Images/output_phyloxml.png" alt="Input file style" width="900"/>
<img src="/Images/output2_phyloxml.png" alt="Input file style" width="900"/>

#### c) Additional Output File: Certainty Value File

In a separate csv file, statistical support values for each variant present in the phylogenetic tree are given. The support or certainty values are the fraction of variants in the tree that support that variant's position based on their informative (=upstream/basal/rootward, downstream/terminal, parallel) pairwise relationships out of all remaining variants in the tree.

<img src="/Images/certainty_values_example.png" alt="Input file style" width="400"/>

In this example, the tree location of variant i is supported by informative pairwise relationships of 8 of the remaining 9 variants (8/9 = 0.88888). The ninth variant has no informative pairwise relationship to variant i. Variants with contradictory relationships are not present in the final tree, since they have been filtered out during the tree generation process.

#### d) Optional Output Files:

**metadata_individuals**


In the metadata output file, the individuals presented in each row correspond to the respective row of the phylogenetic tree (tree layer). The variants in each tree layer were observed in the sequences of the respective metadata output row. In this example, variant a was found in all individuals 1 to 10, and variant b was only found in individuals 2, 3 and 4. Variants that could not be separated into different branches were represented in one tree layer (like variants i and j). In this case, the sequences corresponding to this tree layer (individuals 11 and 12) were each found in at least one of the variants (i and j).

<img src="/Images/output_phyltree_metadata.png" alt="Input file style" width="700"/>


**contradictory_variants**

For certain variants - including those resulting from sequencing errors, recurrent mutations or backmutations - variants with contradictory pairwise hierarchical relationships are ignored for the tree generation, but saved as "contradictory variants". 

**ambiguous_variants**

Based on the pairwise relationships, the final hierarchical order of the variants is established.
For some variants, an explicit position in the tree could not be determined. These variants have ambiguous positions in the tree.


### 4) Additional Information and Contact
More information on the software are available in our publication: [Paper submitted]

For reporting bugs, comments or questions, you are welcome to contact zehra.koksal@sund.ku.dk.

### 5) Referencing

Please cite: [Paper submitted]



