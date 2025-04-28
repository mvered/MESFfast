# MESFfast: Faster Moran Eigenvector Spatial Filtering (MESF)

Quickly performs Moran eigenvector spatial filtering (MESF). MESF involves
removing/reducing spatial autocorrelation present in the residuals of
generalized linear models (GLMs). MESF searches eigenvectors of a transformed
spatial weights matrix (a doubly centered matrix). One-by-one, eigenvectors
are added to the right-hand side of the GLM and Moran's I statistic (a measure
of spatial autocorrelation) is calculated for the residuals of the newly fitted
model. The eigenvector which produces the lowest Moran's I is then permanently
added to the model and the search process repeated to choose additional
eigenvectors to add to the model. This process continues until the p-value
for Moran's I exceeds the specified stopping threshold alpha.

This implementation optimizes MESF for bigger data applications in several ways:
1. Reduces the set of eigenvectors over which to conduct a brute force search
to only the n_eigs largest eigenvectors of the transformed spatial weights
matrix. This speeds up the process by only computing the smaller number of
eigenvectors we are going to consider (using the RSpectra package) and by
requiring us to test fewer eigenvectors in the model re-fitting stage.
2. Runs the search for eigenvectors to add to the model in parallel instead
of sequentially.
3. Calculates Moran's I in a faster and more memory efficient way using the
moranfast() function implemented in this package. This function calculates
Moran's I on the fly rather than by calculating and storing in memory a new
distance matrix for every new set of residuals produced, and uses Rcpp/C++
to speed up calculation. (Credit to Matt Cooper, moranfast package for this).
4. Creates the transformed spatial weights matrix (the eigenvectors
of which will be considered for inclusion in the model) using faster matrix
computations from Rfast package which runs in parallel in C++.

# Installation & Use


```
# installation
library(devtools)
install_github('mvered/MESFfast')


# example use
# for more detailed explanation, see vignette
library(MESFfast)
data(columbus_sp)
neighbors <- dnearneigh(columbus_sp, d1 = 0.0001, d2 = 75) 
idw_weights <- nb2listwdist(neighbors, columbus_sp, type = "idw", zero.policy=TRUE, style="W")
find_mesf <- mesf(crime ~ home_value + income,
                   data = columbus_sp,
                   family = gaussian,
                   listw = idw_weights, verbose = T, alpha = .05,
                   n_eigs = 10, lat = columbus_sp$lat, lon = columbus_sp$lon)
final_model <- glm(crime ~ home_value + income 
                   + find_mesf$vectors[,1] + find_mesf$vectors[,2] + find_mesf$vectors[,3],
                   data = columbus_sp)
moranfast(final_model$residuals, lat=columbus_sp$lat, lon=columbus_sp$lon)
```

