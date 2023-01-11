# Divergent cranial suture fusion at the marsupial-placental boundary

__Authors:__
[Heather E White](mailto:heather.white.17@ucl.ac.uk), 
Abigail S Tucker,
Anjali Goswami

__To cite the paper__: 

>White HW, Tucker AS, Goswami A. Divergent cranial suture fusion at the marsupial-placental boundary. Journal details needed. 2023.

Available at: https://github.com/HeatherEWhite/mammal_cranial_development

If using any of this code or data please cite the paper above and this repo

__To cite this repo__: 

> White HW, Tucker AS, Goswami A. Divergent cranial suture fusion at the marsupial-placental boundary. Journal details needed. 2023. Github: https://github.com/HeatherEWhite/mammal_suture_fusion plus the Zenodo DOI: https://doi.org/10.5281/zenodo.4037220

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4037220.svg)](https://doi.org/10.5281/zenodo.4037220)

UPDATE ZENODO DOI'S ABOVE


## Data

In this folder you will find data stored in the `Data` folder used for running the scripts stored in the `Code_for_analyses` folder.

All .csv documents have single tabs 

The data used here are provided in the `Data` folder
1. `Specimen_info` - full dataset specimen details (n=165), including: continuous and discrete age, species, centroid size, raw closure scores (.csv)
2. `Specimen_info_adults` - details for the adult only specimens (n=22) as above (.csv)
3. `TSCS_adults_CS_pgls` - adult total suture closure scores and centroid size for PGLS of closure and skull size correlation (.csv)
4. `average_adult_score_Krogman_rank` - average closure score for each Krogman region combining all adult specimens and rank order (.csv)
5. `average_total_suture_closure_adults_marsupials_vs_placentals_Krogman` - average closure by Krogman region for marsupials and placentals (.csv)
6. `my_mammal_tree` - phylogeny containing all 22 species within the dataset, trimmed from Upham et al. (2019) (.nexus)
7. `rank_order_adults_Krogman_region` - order of fusion for each species using Krogman regions (.csv)
8. `rank_order_sutures` - order of suture fusion for each species (.csv)
9. `rank_order_sutures_boxplot` - order of fusion for each species formatted to produce boxplot (.csv)
10. `raw_closure_status` - raw closure of all specimens (n=165) as scored from 3D surface meshes (.csv)
11. `suture_closure_average_by_age` - total suture closure score for each species at each discrete age category (.csv)
12. `suture_closure_score_per_specimen_per_Krogman_region` - suture closure score for each specimen for each Krogman region (.csv)
13. `suture_closure_score_species_average_per_suture` - suture closure score for each suture combining all specimens for one species (.csv)
14. `suture_closure_scores_adults_dev_origin` - suture closure scores for each suture (adult specimens) and associated developmental origins (.csv)
15. `suture_closure_scores_all_specimens_dev_origin` - suture closure scores for each (all specimens) and associated developmental origins (.csv)
16. `suture_closure_status_adults` - raw closure for adult specimens (n=22) as scored from 3D surface meshes (.csv)
17. `total_suture_closure_scores_adults` - total suture closure scores for adults only (combining all sutures for each specimen) (.csv)
18. `total_suture_closure_scores_all` - total suture closure scores for all specimens (combining all sutures for each specimen) (.csv)
19. `tree_taxa_names` - taxa names matching those within the phylogeny for the full dataset (.csv)

## Analysis
In this repository you will find raw data (.csv, .Rdata, .nexus) and scripts for analyses (scripts supplied as .R files)

 :file_folder:
* **Data**

All data is outlined above in the 'Data' section. The `Data` folder includes .csv files and .nexus files for analysis

 :file_folder:
* **Code_for_analyses**

`ancestral_estimation.R`

`fusion_vs_skull_size.R`

`rank_order_analysis.R`

`suture_closure_scores.R`



## License
This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/HeatherEWhite/suture_complexity_metrics/blob/master/LICENSE) file for details

## Session Info
For reproducibility purposes, here is the output of `devtools::session_info()` used to perform the analyses in the publication. 

```{r}
─ Session info ──────────────────────────────────────────────────────────────────────────────────
 setting  value                       
 version  R version 3.6.1 (2019-07-05)
 os       macOS Catalina 10.15.7      
 system   x86_64, darwin15.6.0        
 ui       RStudio                     
 language (EN)                        
 collate  en_GB.UTF-8                 
 ctype    en_GB.UTF-8                 
 tz       Europe/London               
 date     2023-01-11                  

─ Packages ──────────────────────────────────────────────────────────────────────────────────────
 package     * version date       lib source        
 assertthat    0.2.1   2019-03-21 [1] CRAN (R 3.6.0)
 cachem        1.0.5   2021-05-15 [1] CRAN (R 3.6.2)
 callr         3.7.0   2021-04-20 [1] CRAN (R 3.6.2)
 cli           3.0.1   2021-07-17 [1] CRAN (R 3.6.2)
 crayon        1.4.1   2021-02-08 [1] CRAN (R 3.6.1)
 desc          1.2.0   2018-05-01 [1] CRAN (R 3.6.0)
 devtools      2.3.2   2020-09-18 [1] CRAN (R 3.6.2)
 ellipsis      0.3.2   2021-04-29 [1] CRAN (R 3.6.2)
 fastmap       1.1.0   2021-01-25 [1] CRAN (R 3.6.2)
 fs            1.5.0   2020-07-31 [1] CRAN (R 3.6.2)
 glue          1.4.2   2020-08-27 [1] CRAN (R 3.6.2)
 lifecycle     1.0.1   2021-09-24 [1] CRAN (R 3.6.2)
 magrittr      2.0.3   2022-03-30 [1] CRAN (R 3.6.1)
 memoise       2.0.0   2021-01-26 [1] CRAN (R 3.6.2)
 pkgbuild      1.2.0   2020-12-15 [1] CRAN (R 3.6.2)
 pkgload       1.1.0   2020-05-29 [1] CRAN (R 3.6.2)
 prettyunits   1.1.1   2020-01-24 [1] CRAN (R 3.6.0)
 processx      3.5.2   2021-04-30 [1] CRAN (R 3.6.2)
 ps            1.6.0   2021-02-28 [1] CRAN (R 3.6.2)
 purrr         0.3.4   2020-04-17 [1] CRAN (R 3.6.2)
 R6            2.5.1   2021-08-19 [1] CRAN (R 3.6.2)
 remotes       2.4.0   2021-06-02 [1] CRAN (R 3.6.2)
 rlang         1.0.6   2022-09-24 [1] CRAN (R 3.6.1)
 rprojroot     2.0.2   2020-11-15 [1] CRAN (R 3.6.2)
 rstudioapi    0.13    2020-11-12 [1] CRAN (R 3.6.2)
 sessioninfo   1.1.1   2018-11-05 [1] CRAN (R 3.6.0)
 testthat      3.0.1   2020-12-17 [1] CRAN (R 3.6.2)
 usethis       2.0.0   2020-12-10 [1] CRAN (R 3.6.2)
 withr         2.4.2   2021-04-18 [1] CRAN (R 3.6.2)
```
