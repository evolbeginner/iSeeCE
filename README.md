# iSeeGC

iSeeGC is a bioinformatic pipeline that allows identifying recurrent gene conversion across lineages. It is written by [Sishuo Wang](http://www.researcherid.com/rid/F-8081-2015) (sishuowang{at}hotmail.com) from University of British Columbia.

# Installation
iSeeGC is run in Linux environment.

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
To the full usage, please type 
```bash
bash iSeeGC.sh -h
```

iSeeGC will identify recurrent gene conversion (GC) across species based on gene phylogeny and synteny. In brief, paralogs from the same species that form monophyly and have syntenic (positional) orthologs across species will be identified as converted genes.

The results of converted genes will be output to the file **FastTree/GC.result**.

By default, syntenic orthologs will be identified by using only best reciprocal BLAST hits. When "--mauve" is specified, syntenic orthologs will additionally be identified by Mauve.

# Caveats
![alt text](https://upload.wikimedia.org/wikipedia/commons/thumb/2/2f/Google_2015_logo.svg/220px-Google_2015_logo.svg.png)

# Copyright and Licence:
Please see [CC0](https://creativecommons.org/share-your-work/public-domain/cc0/) for licensing information.

You are welcome to cite our paper coming soon for iSeeGC if you feel it is useful.

You may need to cite corresponding papers for the use of OrthoMCL, MAFFT, FastTree, RAxML, Mauve, BioRuby and BioPerl.

# Acknowledgements
Thanks for Michael Imelfort and Kranti Konganti for sharing the free script gff2fasta.pl for sequence processing.


