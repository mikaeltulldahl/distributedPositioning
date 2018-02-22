# DistributedPositioning

## TODO improve difficulty:
* generalize to 3d
* increase time between each communication
* non-zero-mean noise

## TODO improve accuracy estimate:
* use number of knownNodes and their accuracy to estimate own accuracy
* if knownNodes <= 2 then cost should be 0 and can't be used for estimating accuracy
* stddev(posEst) could be a good estimate of accuracy if it is close to correct position
* switch from scalar 'accuracy' to covariance matrix for each direction

## TODO misc:
* find why system does not converge in beginning
* improve search algortihm in findNode

