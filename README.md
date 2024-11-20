# R code and data to the article: Class-focused variable importance via multi forests

Authors: Roman Hornung<sup>1,2,*</sup> and Alexander Hapfelmeier<sup>3</sup>

1. Institute for Medical Information Processing, Biometry and Epidemiology, LMU Munich, Marchioninistr. 15, Munich, 81377, Germany, ORCID: 0000-0002-6036-1495.
2. Munich Center for Machine Learning (MCML), Munich, Germany.
3. Institute of AI and Informatics in Medicine, TUM School of Medicine and Health, Technical University of Munich, Ismaninger Str. 22, Munich, 81675, Germany, ORCID: 0000-0001-6765-6352.

\* For questions, please contact: hornung@ibe.med.uni-muenchen.de

---

## Program and Platform

- **Program**: R, versions 4.1.2 and 4.3.2.
- The raw results of the simulation and the benchmark study were obtained on a Linux cluster, and the evaluation of the raw results to produce the final results (i.e., the figures and tables) was performed on Windows 11.
- Below is the output of the R command `sessionInfo()` on the Linux machine and on the Windows machine. 
  The output specifies which R packages and versions of those packages were used to generate the raw results 
  and to evaluate them.

### sessionInfo() on the Linux cluster

```R
> sessionInfo()
R version 4.1.2 (2021-11-01)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: SUSE Linux Enterprise Server 15 SP1

Matrix products: default
BLAS/LAPACK: /dss/dsshome1/lrz/sys/spack/release/22.2.1/opt/haswell/intel-mkl/
2020.4.304-gcc-3e7v2iy/compilers_and_libraries_2020.4.304/linux/mkl/lib/
intel64_lin/libmkl_gf_lp64.so

Random number generation:
 RNG:     Mersenne-Twister
 Normal:  Inversion
 Sample:  Rounding

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
[1] nnet_7.3-17            diversityForestD_0.5.0 diversityForestC_0.5.0
[4] diversityForestB_0.5.0 diversityForestA_0.5.0 ranger_0.14.1

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.8.3     magrittr_2.0.3   ggpubr_0.6.0     tidyselect_1.2.0
 [5] munsell_0.5.0    colorspace_2.0-3 lattice_0.20-45  R6_2.5.1
 [9] rlang_1.1.2      rstatix_0.7.2    carData_3.0-5    fansi_1.0.3
[13] car_3.1-2        dplyr_1.1.4      grid_4.1.2       broom_1.0.5
[17] gtable_0.3.0     utf8_1.2.2       cli_3.6.2        abind_1.4-5
[21] tibble_3.2.1     lifecycle_1.0.4  ggsignif_0.6.4   Matrix_1.6-4
[25] tidyr_1.3.0      purrr_1.0.2      ggplot2_3.4.4    vctrs_0.6.5
[29] glue_1.6.2       compiler_4.1.2   pillar_1.9.0     backports_1.4.1
[33] generics_0.1.2   scales_1.2.0     pkgconfig_2.0.3
```

NOTE: The packages `diversityForestA`, `diversityForestB`, `diversityForestC`,
and `diversityForestD` are not available on the CRAN repository. Instead, 
these packages are included in this electronic appendix (see section "General 
Information and Contents of this Electronic Appendix" below for details) for 
the sake of reproducibility. These packages are preliminary versions of the 
`diversityForest` R package used in the simulation and the benchmark study.

### sessionInfo() on the Windows machine

```R
> sessionInfo()
R version 4.3.2 (2023-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default


locale:
[1] LC_COLLATE=German_Germany.utf8  LC_CTYPE=German_Germany.utf8    LC_MONETARY=German_Germany.utf8
[4] LC_NUMERIC=C                    LC_TIME=German_Germany.utf8    

time zone: Europe/Berlin
tzcode source: internal

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ranger_0.16.0         diversityForest_0.5.0 stringdist_0.9.12     pmlbr_0.2.1          
 [5] OpenML_1.12           patchwork_1.2.0       scales_1.3.0          rstatix_0.7.2        
 [9] cowplot_1.1.2         RColorBrewer_1.1-3    gridExtra_2.3         ggplot2_3.4.4        
[13] tidyr_1.3.0           forcats_1.0.0         stringr_1.5.1         xtable_1.8-4         
[17] dplyr_1.1.4          

loaded via a namespace (and not attached):
 [1] utf8_1.2.4        generics_0.1.3    lattice_0.21-9    stringi_1.8.3     digest_0.6.34    
 [6] magrittr_2.0.3    BBmisc_1.13       fastmap_1.1.1     Matrix_1.6-1.1    jsonlite_1.8.8   
[11] backports_1.4.1   httr_1.4.7        purrr_1.0.2       fansi_1.0.6       XML_3.99-0.16.1  
[16] abind_1.4-5       cli_3.6.2         rlang_1.1.3       munsell_0.5.0     cachem_1.0.8     
[21] withr_3.0.0       parallel_4.3.2    tools_4.3.2       ggsignif_0.6.4    memoise_2.0.1    
[26] checkmate_2.3.1   colorspace_2.1-0  ggpubr_0.6.0      broom_1.0.5       curl_5.2.0       
[31] vctrs_0.6.5       R6_2.5.1          lifecycle_1.0.4   car_3.1-2         pkgconfig_2.0.3  
[36] pillar_1.9.0      gtable_0.3.4      Rcpp_1.0.12       glue_1.7.0        data.table_1.15.4
[41] tibble_3.2.1      tidyselect_1.2.1  rstudioapi_0.15.0 carData_3.0-5     compiler_4.3.2
```

---

## General Information and Contents of this Electronic Appendix

### Preliminary Remark
Readers who are not interested in the detailed 
  contents of this electronic appendix, but only in the evaluation of 
  the results or the full reproduction of the results, may skip to 
  the sections "Evaluation of the Results" or "Full Reproduction 
  of the Results", respectively.

### Contents
- **simulation**: This subfolder contains the R scripts `simulation.R`,
    `simulation_functions.R`, and `simulation_evaluation.R` as well
    as the subfolder `intermediate_results`.
    
    The R script `simulation.R` can be used to run the simulation study 
    using parallel computing on the Linux cluster, producing the 
    raw results (see the section "Full Reproduction of the Results" 
    below for details).

    The R script `simulation_functions.R` is sourced by `simulation.R` 
    and contains all of the functions required for the simulation study.

    The R script `simulation_evaluation.R` evaluates the raw results of 
    the simulation study to produce all figures and tables associated with 
    the simulation study.

    The subfolder `intermediate_results` contains the intermediate results 
    of the simulation study.
- **benchmark_study**: This subfolder contains the R scripts `benchmark_study.R`, 
    `benchmark_study_functions.R` and `benchmark_study_evaluation.R` as 
     well as the subfolders `data` and `intermediate_results`.

    The R script `benchmark_study.R` can be used to run the benchmark study 
    in parallel on the Linux cluster, producing the raw results. 

    The R script `benchmark_study_functions.R` is sourced by `benchmark_study.R` 
    and contains all the functions needed for the benchmark study.

    The R script `benchmark_study_evaluation.R` evaluates the raw results 
    of the benchmark study to produce all figures and tables associated with 
    the benchmark study.
    
    The subfolder `intermediate_results` contains the intermediate results 
    of the benchmark study.

    The subfolder `data` contains the following files: 
    1. The R script `download_data.R` downloads, selects, and preprocesses 
         the datasets used in the benchmark study.
    2. The Rda file `datainfo.Rda` was generated by `download_data.R` and 
         contains a data.frame with meta information about the datasets, such 
         as sample sizes and numbers of covariates.
    3. The RData file `classifTasks.infos.RData` was also generated by 
         `download_data.R` and contains the list of all OpenML tasks that 
         met our filtering criteria.
    4. The text file `OpenML-CC18_table.txt` contains meta information about 
         all 72 datasets included in the OpenML-CC18 machine learning benchmark 
         suite.
    5. The subfolder `datasets` contains the pre-processed versions of all 121 datasets 
         included in the benchmark study in the form of Rda files.
- **diversityForest_versions**: This subfolder contains four subfolders 
    `diversityForestA`, `diversityForestB`, `diversityForestC`, and `diversityForestD`. 

    These subfolders contain different versions of the R package `diversityForest` 
    that were used in the simulation study and the benchmark study:
    1. `diversityForestA` implements the version of multi forests with variants 
         "Gini" and "Squared" ("wsquared_wgini"),
    2. `diversityForestB` implements the version of multi forests with variants 
         "Assign Classes" and "Squared" ("wsquared_wogini"),
    3. `diversityForestC` implements the version of multi forests with variants 
         "Gini" and "Non-Squared" ("wosquared_wgini"), and
    4. `diversityForestD` implements the version of multi forests with variants 
         "Assign Classes" and "Non-Squared" ("wosquared_wogini").
- **figures**: This subfolder contains all figures shown in the main paper  
    and in the supplementary material as eps and pdf files, respectively.
- **tables**: This subfolder contains the tables shown in the main paper and 
    in the supplementary material as tex files. Note that when the tables were included 
    in the main paper and in the supplementary material, the tex code of these files was 
    slightly modified for visual reasons (without changing the values in the tables).
    Note further that Table 1 is not included here because this table is related
    to the design of the simulation study and thus does not show empirical or simulation
    results.

---

## Evaluation of the Results

For the evaluation of the results it is not necessary to re-perform the analyses:

  The R scripts `simulation_evaluation.R` and `benchmark_study_evaluation.R` contained 
  in the subfolders `simulation` and `benchmark_study`, respectively, produce all 
  results shown in the main paper and in the supplementary material without the need
  of re-performing the analyses. These R scripts read in Rda files (stored in the 
  subfolders `simulation/intermediate_results` and 
  `benchmark_study/intermediate_results`) that contain the raw results.

---

## Full Reproduction of the Results
- As a first step, the folder `MultiForests_code_and_data` this README is contained
  in has to be placed in the home directory (`~/`) of a Linux machine.
  Alternatively, the folder can be placed in a directory of your choice, with the paths 
  in line 3 of the R scripts `simulation/simulation.R` and `benchmark_study/benchmark_study.R` 
  being changed accordingly.
- An MPI environment is required.
- The R scripts `simulation/simulation.R` and `benchmark_study/benchmark_study.R` 
  perform the simulation study and the benchmark study, respectively.
  These R scripts require the RMPISNOW shell script from the R package `snow`.
  Therefore, before executing these scripts you need to install the RMPISNOW shell script 
  from the installed `snow` R package or `inst` directory of the package sources
  of the `snow` R package in an appropriate location, preferably
  on your path. 
  See http://homepage.divms.uiowa.edu/~luke/R/cluster/cluster.html (last accessed: 
  20th November 2024) for more details.
  Subsequently, you need to create two sh files, each for a different of the
  above R scripts. 

  The following is the content of an example sh file `simulation.sh`:

  ```bash
    #!/bin/bash
    #SBATCH -o /myoutfiledirectory/myjob.%j.%N.out
    #SBATCH -D /myhomedirectory
    #SBATCH -J simulation
    #SBATCH --get-user-env 
    #SBATCH --clusters=myclustername
    #SBATCH --partition=mypartitionname
    #SBATCH --qos=mypartitionname
    #SBATCH --nodes=??
    #SBATCH --ntasks-per-node=??
    #SBATCH --mail-type=end
    #SBATCH --mail-user=my@mail.de
    #SBATCH --export=NONE
    #SBATCH --time=??:??:??  

    module load slurm_setup
    module load r/4.1.2-gcc11-mkl
    module load openmpi

    mpiexec -n $SLURM_NTASKS /pathtoRMPISNOW/RMPISNOW \
      < ./MultiForests_code_and_data/simulation/simulation.R
  ```

  The above sh file of course has to be adjusted to be useable (e.g., the "?"s have
  to replaced by actual numbers, the directories have to be adjusted and
  you need to specify your e-mail address; an e-mail will be sent to this address
  once the job is finished).

  Note that it is possible to use other parallelization techniques (e.g., the parallel R 
  package) than RMPISNOW to reproduce the results. This is because we use a specific seed 
  for each line in the scenariogrid data frames created by the `simulation.R` and 
  `benchmark_study.R` scripts. Each line in these data frames correspond to one iteration 
  in the simulation study and the benchmark study, respectively (see the corresponding 
  files for details). This makes the reproducibility independent of the specific type 
  of parallelization. However, to use a different type of parallelization than RMPINOW, 
  it is necessary to modify the  `simulation.R` and `benchmark_study.R` scripts accordingly.
