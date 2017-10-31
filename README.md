# iSeeCE

iSeeCE is a bioinformatic pipeline that allows identifying recurrent gene conversion across lineages.

# Flowchart of iSeeCE
<div style="display:table-cell; vertical-align:middle; text-align:center">
<img src=images/Flowchart.png width=800 height=500></img>
</div>

# Installation
iSeeCE is run in Linux environment.

Many scripts included in this product require packages from
* [BioRuby](http://bioruby.org)
* [BioPerl](http://bioperl.org)

This product also employs several other computational tools. Please ensure that you have them installed.
* [OrthoMCL](http://orthomcl.org/orthomcl/) 
* [MAFFT](http://mafft.cbrc.jp/alignment/software/)
* [FastTree](http://darlinglab.org/blog/2015/03/23/not-so-fast-fasttree.html)
* [RAxML](https://sco.h-its.org/exelixis/software.html) (optional)
* [Mauve](http://darlinglab.org/mauve/mauve.html) (optional)

To see whether required software and libraries have been installed, please type 
```bash
bash check_requirements.sh
```

# Usage
iSeeCE will identify recurrent concerted evolution (CE) or gene conversion (GC) across species based on gene phylogeny and synteny. In brief, paralogs from the same species that form monophyly and have syntenic (positional) orthologs across species will be identified as converted genes.

The input files should be genomes in the format of [Genbank](https://www.ncbi.nlm.nih.gov/genbank/). Please make sure that the nucleotide sequences of the genome are included in each Genbank-formatted input file.

For the full usage, please type 
```bash
bash iSeeCE.sh -h
```

As follows are some examples:
```bash
#
bash iseeCE.sh --indir genebank_dir --outdir new_dir
# No. of converted genes >= 8, each conversion event supported by bootstrap >= 0.9
bash iseeCE.sh --indir genebank_dir --outdir new_dir --gc_count_min 8 --bootstrap 0.9
```

# Notes
By default, syntenic orthologs will be identified by using only best reciprocal BLAST hits. The results of converted genes will be output to the file **FastTree/GC.result**.

When "--is_mauve" is specified, syntenic orthologs will additionally be identified by Mauve whose results will be output to the file **FastTree/GC.mauve.result**. Also, it is necessary to specify the path to Mauve.jar and progressiveMauve using --mauve_jar and --mauve, respectively. So please make sure that Mauve.jar has been installed if you want to perform the analysis with Mauve. Note that '--is_mauve' is optional. When '--is_mauve' is **NOT** specified, the analysis will be performed based on the best reciprocal BLAST hits.

# Copyright and Licence:
Please see [CC0](https://creativecommons.org/share-your-work/public-domain/cc0/) for licensing information.

You are welcome to cite our paper coming soon for iSeeCE if you feel it is useful.

You may need to cite corresponding papers for the use of OrthoMCL, MAFFT, FastTree, RAxML, Mauve, BioRuby and BioPerl.

# Acknowledgements
Thanks for Michael Imelfort and Kranti Konganti for sharing the free script gff2fasta.pl for sequence processing.


