# sand_dollar_response
## Intro
This repository contains the data and results supporting the manuscript Mongiardino Koch (2023) _Sci. Rep._

The R script _sand_dollar_phylogeny.R_ can reproduce all steps of the analysis. It relies on R functions from packages _ape_ (Paradis & Schliep 2019), _phangorn_ (Schliep 2011), _tidyverse_ (Wickham et al. 2019), _combinat_ (Chasalow 2012), and _MetBrewer_ (Mills 2012). Analyses also relied on MAFFT v7.505 (Katoh & Standley 2013) for multiple sequence alignment, IQ-TREE v1.6.12 (Nguyen et al. 2015) for model selection and maximum likelihood inference, and TNT v1.5 (Goloboff & Catalano 2016) for morphological inference under maximum parsimony.

# References
Chasalow S. 2012. combinat: combinatorics utilities. R package version 0.0-8. <https://CRAN.R-project.org/package=combinat>.
Goloboff P.A. & Catalano S.A. 2016. TNT version 1.5, including a full implementation of phylogenetic morphometrics. Cladistics 32: 221-238.
Katoh K. & Standley D.M. 2013. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Molecular Biology and Evolution, 30: 772-780.
Mills B.R. 2022. MetBrewer: Color Palettes Inspired by Works at the Metropolitan Museum of Art. R package version 0.2.0. <https://CRAN.R-project.org/package=MetBrewer>.
Nguyen L.-T., Schmidt H.A., von Haeseler A., Minh B.Q. 2015. IQ-TREE: A fast and effective stochastic algorithm for estimating maximum likelihood phylogenies. Molecular Biology and Evolution, 32: 268-274.
Paradis E. & Schliep K. 2019. ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. Bioinformatics, 35: 526-528.
Schliep K. 2011. phangorn: phylogenetic analysis in R. Bioinformatics, 27: 592-593.
Wickham H., Averick M., Bryan J., Chang W., McGowan L.D., et al. (2019). Welcome to the tidyverse. Journal of Open Source Software, 4: 1686.
