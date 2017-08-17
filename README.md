
# transmod - Transcriptogram and Modularity

R package that implements a method to specify list of most relevant genes among differential expression profiles based on gene network knowledge. This is done through [transcriptogram](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-15-1181) of gene expression profiles and analysis of differentially expressed modules.

Supporting Institutions:

* [Federal Institute of Pará (IFPA)](http://www.ifpa.edu.br/) - Funding
* [Evandro Chagas Institute (MS/SVS/IEC)](http://www.iec.gov.br/)
* [Vale Institute of Technology (ITV)](http://www.vale.com/itv/)
* [French National Centre for Scientific Research (CNRS)](http://www.cnrs.fr/)

**Main reference:** DIAS JÚNIOR, J.F.S.; ALVES, R.; COMMES, T. [A module-based approach for evaluating differential genome-wide expression profiles](https://doi.org/10.1109/BRACIS.2016.069). In the 5th Brazilian Conference on Intelligent System ([BRACIS](http://cin.ufpe.br/%7Ebracis2016/)). Recife, PE, Brazil, October 9-12, 2016.

## Package installation

This package was developed and tested only on version 3.3.1 of R platform.

Since *transmod* is still under development, it is not yet available on [CRAN repository](https://cran.r-project.org/).

It is necessary to install other packages for the full operation of the *transmod*:

``` R
# devtools
install.packages("devtools")

# GeneSelector
source("https://bioconductor.org/biocLite.R")
biocLite("GeneSelector")
```

Installing the *transmod* package:


``` R
library(devtools)
install_github("joseflaviojr/transmod")
```

## Example usage

Example of selection and analysis of differentially expressed genes.

The *transmod* package contains sample data, which will be used here:

``` R
library(transmod)

seriation <- seriation_dias2016
expression <- expression_gse48173
```

The `seriation_dias2016` is a list of human genes ordered according to a network. Reference: [https://doi.org/10.1109/BRACIS.2016.069](https://doi.org/10.1109/BRACIS.2016.069). For more details, run `?seriation_dias2016`.

The `expression_gse48173` is a RNA-Seq data corresponding to [experimental study about leukemia](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48173). The table contains 72 human samples: 43 Acute Myeloid Leukemia (AML), 12 Acute Lymphoblastic Leukemia (ALL) and 17 Healthy (HEA). The samples cover 21865 genes:

``` no-highlight
             GSM1185603_HEA GSM1185604_HEA ...
WASH7P              1.28108        0.76803
OR4F5               0.00000        0.00000
LOC100133331        0.61782        0.58245
LOC100288069        3.95892        3.08883
NCRNA00115          1.50934        1.08628
LOC643837           2.03448        1.45142
...
```

Arranging the expression table according to the seriation:

``` R
seriation <- intersect(seriation,rownames(expression))
expression <- expression[seriation,]
```

Calculating the [transcriptogram](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-15-1181) of each sample (expression profile); because this method is pure R, its execution can be time consuming:

``` R
tgram <- transcriptogram(expression)
```

Calculating the differentiation level of each gene between AML and healthy samples:

``` R
diffg <- differentiate(tgram, 1:17, 30:72)
```

The `head(diffg)` returns `1.928160 2.060090 2.038235 2.027537 2.030800 2.028469`.

The differentiation series is inspected in order to detect modules:

``` R
modules <- modularize(diffg)
modules_summary <- summarize_modules(diffg, modules)
```

The variable `modules_summary` contains the summary of each module detected:

``` no-highlight
  module begin end size     mean      max      min
1      1     1  91   91 2.163497 3.028608 1.089711
2      2    92  92    1 1.090482 1.090482 1.090482
3      3    93 270  178 3.339707 5.585274 1.088448
4      4   271 276    6 2.377144 2.428238 2.337269
5      5   277 277    1 2.343614 2.343614 2.343614
6      6   278 392  115 5.962280 7.465978 2.181652
...
```

To view the levels of differentiation for each gene and the detected modules, run:

``` R
palette(c("black","gray"))
plot(diffg, col=modules, type="h", main="Differentially Expressed Modules", xlab="Genes", ylab="Differentiation Level")
```

Result:

![Image - Differentially Expressed Modules](https://github.com/joseflaviojr/transmod/raw/master/proj/example_modules.png)

Selecting the 100 most relevant genes among the modules:

``` R
selection_index <- select_from_modules(diffg, modules, select=100)
selection <- seriation[selection_index]
cat(selection, sep="\n")
```

Result:

``` no-highlight
COPE
LYRM2
HDHD2
BCL7B
CSNK1E
REL
POP5
EFHC1
KIAA0408
AURKAIP1
...
```

The selected genes can be submitted to enrichment tools to find scientific data related. Using the list above, for example, in the online tool [Enrichr](http://amp.pharm.mssm.edu/Enrichr/), it is obtained more significantly from databases/ontologies:

* OMIM Disease = Leukemia, Cataract
* Jensen DISEASES = Acute Promyelocytic Leukemia (AML subtype), Corneal disease
* MSigDB Computational = [MODULE_13](http://robotics.stanford.edu/~erans/cancer/modules/module_13) (related to leukemia and B lymphoma)
* LINCS L1000 Chem Pert up = CPC011 HT29 6H-idarubicin hcl-10.0 (the *idarubicin* is related to the [leukemia treatment](http://emedicine.medscape.com/article/2004793-overview))
* dbGaP = Keratoconus (There are studies that relate leukemia and cataract: [1](https://www.ncbi.nlm.nih.gov/pubmed/24116693), [2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3861856/), [3](http://www.nature.com/eye/journal/v18/n7/full/6701308a.html?foxtrotcallback=true), [4](https://www.ncbi.nlm.nih.gov/pubmed/12579166))
* GO Biological Process = positive regulation of transcription from RNA polymerase II promoter (this class contains 10 term members directly related to leukemia: B-cell lymphoma/leukemia 11A (BCL11A), B-cell lymphoma/leukemia 11B (BCL11B), Hepatic leukemia factor (HLF), Friend leukemia integration 1 transcription factor (FLI1), Pre-B-cell leukemia transcription factor 2 (PBX2), Pre-B-cell leukemia transcription factor 3 (PBX3), T-cell leukemia homeobox protein 1 (TLX1), T-cell acute lymphocytic leukemia protein 1 (TAL1), Leukemia inhibitory factor (LIF) and T-cell leukemia homeobox protein 2 (TLX2))
* GO Cellular Component = spindle pole centrosome (There are studies that relate leukemia and centrosome: [1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5411748/), [2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4113111/), [3](http://onlinelibrary.wiley.com/store/10.1111/j.1365-2141.2009.07772.x/asset/j.1365-2141.2009.07772.x.pdf;jsessionid=76406C06B17CD3A74442763FB5709C6E.f02t01?v=1&t=j6fi9k9v&s=72c1748f26b43786413fdc073bb635f809116634), [4](https://www.ncbi.nlm.nih.gov/pubmed/11309836))
* GO Molecular Function = RNA binding ([1](http://cancerdiscovery.aacrjournals.org/content/7/6/OF16), [2](http://www.lrjournal.com/article/S0145-2126(17)30015-2/pdf), [3](http://www.bloodjournal.org/content/128/22/739?sso-checked=true), [4](https://www.sciencedaily.com/releases/2016/03/160314211412.htm), [5](http://www.nature.com/ni/journal/v11/n8/full/ni.1901.html))
* Human Phenotype Ontology = Pectus excavatum ([1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4267483/), [2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3401115/)), Zonular cataract, Congenital cataract
* Jensen TISSUES = Corpus callosum (hemorrhagic complications are common in patients with leukemia: [1](http://www.japi.org/october_2015/13_cr_acute_myeloid_leukemia.pdf), [2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4895777/), [3](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3099428/)), Bone marrow

## Contributors

* José Flávio de Souza Dias Júnior (Researcher/Coordinator) - joseflaviojr@gmail.com
* Ronnie Alves (Researcher) - alvesrco@gmail.com
* Thérèse Commes (Researcher) - therese.commes@gmail.com
* Andréa do Socorro Bolhosa Sarmento (Scientific Initiation/Scholarship Student) - andreassarmento@yahoo.com.br
* Guilherme Lopes Sousa (Scientific Initiation/Volunteer Student) - guilherme.lops@hotmail.com
