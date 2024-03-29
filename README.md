# bpptools
A collection of R functions assisting phylogenetic analyses in BPP software.

Among others, it covers:
1. parsing of input files (sequence, classification and control files) for BPP (Yang & Rannala 2010, Flouri et al. 2018)
2. importing outputs (posterior sample)
3. converting outputs into file formats used, e.g., by FigTree (https://github.com/rambaut/figtree/releases) or Tracer (Rambaut et al. 2018)
4. calculation of genealogical divergence index (Jackson et al. 2017, Leaché et al. 2019)
5. plotting of trees annotated either with posterior probabilities or parameters of multispecies coalescent model (Rannala & Yang 2003)
6. plotting of prior and posterior densities

Two example scripts are provided. The script 'bpptools_examples.R' shows a couple of typical operations that can be performed by bpptools. The script 'bpptools_pipeline.R' provides an overview of bpptools functionality as it starts with preparation of input files for A11 analysis and proceeds by processing of A11 outputs, which are also used to prepare inputs for A00 analysis and then it concludes with processing of A00 outputs. The scripts partly re-create results already presented in standard output of BPP, but they offer much more flexibility and publication quality visualisation.

The main toolkit is in the file func/bpptools.R. Apart from these functions, the scripts also employ functions using posterior sample of species trees to find the maximum clade credibility tree ('find_mcct', also at https://github.com/onmikula/find_mcct) and maximum credibility species delimitation ('find_spdelim'). Included is also a function 'read_loci' for import of multi-locus sequence data. The functions in bpptools.R depend on R packages ape (Paradis & Schliep 2019) and invgamma (Kahle & Stamey 2017). The other functions depend also on phangorn (Schliep 2011)

For more details about the functions and their arguments, see function headers in the file bpptools.R and comments in the scripts.
The example data are taken for the description of Harenna mouse (*Mus harennensis*) from southern Ethiopia (Krásová et al. 2022).

**References:**
- Flouri T, Jiao X, Rannala B, Yang Z (2018) Species tree inference with BPP using genomic sequences and the multispecies coalescent. Mol Biol Evol 35: 2585-2593.
- Jackson N, Carstens B, Morales A, O’Meara BC (2017) Species delimitation with gene flow. Syst Biol 66: 799-812.
- Kahle D, Stamey J (2017) invgamma: The Inverse Gamma Distribution. R package version 1.1. https://CRAN.R-project.org/package=invgamma
- Krásová J, Mikula O, Lavrenchenko LA, Šumbera R, Meheretu Y, Bryja J (2002) A new rodent species of the genus Mus (Rodentia: Muridae) confirms the biogeographical uniqueness of the isolated forests of southern Ethiopia. Org Divers Evol https://doi.org/10.1007/s13127-022-00539-x
- Leaché AD, Zhu T, Rannala B, Yang Z (2019) The spectre of too many species. Syst Biol 68: 168-181.
- Paradis E, Schliep K (2019) ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. Bioinformatics 35, 526–528. 
- R Core Team (2022) R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org
- Rambaut A, Drummond AJ, Xie D, Baele G, Suchard MA (2018) Posterior summarization in Bayesian phylogenetics using Tracer 1.7. Syst Biol 67: 901–904.
- Rannala B, Yang Z (2003) Bayes estimation of species divergence times and ancestral population sizes using DNA sequences from multiple loci. Genetics 164: 1645–1656.
- Schliep KP (2011) phangorn: phylogenetic analysis in R. Bioinformatics, 27: 592-593.
- Yang Z, Rannala B (2010) Bayesian species delimitation using multilocus sequence data. Proc Natl Acad Sci U S A 107: 9264–9269.
