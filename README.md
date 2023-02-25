# kCRFC

This is the R code to implement the simulations and real data analysis of the following paper:
> Hyunsung Kim and Yaeji Lim (2022+). Functional clustering on a sphere via Riemannian functional principal components, *under review.*

## Description

- **functions.R** : R functions to implement the proposed **kCRFC** and other functions are included.
- **simulation_S2.R** : Simulation code for the curves on a 2-dimensional sphere
- **simulation_S3.R** : Simulation code for the curves on a 3-dimensional sphere
- **bird_bigration_clust.R** : Real data clustering of Egyptian vultures in the Middle East and East Africa (2-dimensional sphere)
- **medfly_clust.R** : Real data clustering of compositional medfly data (3-dimensional sphere)
- **clust_validation.R** : Clustering validation from real data clustering of compositional medfly data (3-dimensional sphere)
- **data/** : Preprocessed bird migration data is contained.
- **figure/** : Figures used on the paper are contained.



## Data source

- **Bird migration data** : You can download at the [**movebank**](https://www.movebank.org/cms/webapp?gwt_fragment=page=studies,path=study9651291) repository.

- **Fruit fly behaviors data** : Available on the supplement material of 

  > Chiou, J. M., & Müller, H. G. (2014). Linear manifold modelling of multivariate functional data. *Journal of the Royal Statistical Society: Series B: Statistical Methodology*, 605-626.
