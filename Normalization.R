
#Description and effect of normalization techniques
#Min–Max normalization, Z-score normalization, log2 normalization, upper quartile and
#whitening methods were included in this step
#musa article



#   In MuSA tool, radiomic and genomic features can be normalized accordingly to five most common nor-
#   malization methods. Min–Max normalization, Z-score normalization, log2 normalization, upper quartile and
# whitening methods were included in this step. Min–max normalization is one of the most common method
# to normalize data. For every feature, the minimum value of a feature gets transformed into a 0, the maximum
# value gets transformed into a 1, and every other value is rescaled to lie within a range of 0 to 1^36
# .Regarding stand-
#   ardized z-score normalization, each feature was normalized as z = (x-−
#                                                                      x)/s, where x, −
# x and s are the feature, the
# meanand the standard deviation respectively ^37 . Features can also be standardized as Log-transformation (base
#                                                                                                           2), a constant value a = b – min(x) where b is 1 and x is the feature, was added to the data for handling negative
# values ^38 . The upper quartile normalization divides each read count by the 75th percentile of the read counts in
# its sample ^39 ; lastly, whitening normalization technique from the principle component analysis (PCA), is based
# on a linear transformation that converts a vector of random variables with a known covariance matrix into a
# set of new variables whose covariance is the identity matrix, meaning that they are uncorrelated and each have
# variance equal to one ^40 .*/
  
#quantile-quantile
#upper quantile
#z-score
#min-max
#log2
#trimmed