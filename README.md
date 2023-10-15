# Efficient and Scalable Clustering Algorithms for Optimizing the Average Silhouette Width

### Motivation
<p align="justify"> 
The Average Silhouette Width (ASW) is a popular internal clustering validation index used to measure clustering quality and determine an optimal number of clusters. There have been only a few proposals in the literature, such as the Optimum Silhouette (OSil) algoritm and the PAMSil algorithm, using the ASW as an objective function in cost-based clustering. OSil and PAMSil are computationally expensive with $\mathcal{O}(qkN^3)$ time complexity, where $q$ is the number of iterations needed for convergence, $k$ is the number of clusters, and $N$ is the dataset size. Moreover, for large datasets, $q$ can be really large. Consequently, both PAMSil and OSil are only suitable for clustering small datasets. </p>
<p align="justify"> 
We have proposed a new algorithm called Efficient Optimum Silhouette (effOSil), which performs the exact OSil algorithm but $\mathcal{O}(N)$ times faster than the original OSil algorithm at the cost of storing $\mathcal{O}(N)$ additional values. </p>

</p>

### Description

<p align="justify"> 
This R-package (EfficientOASW v.0.0.0.9000) implements the computationally expensive OSil algorithm and provides the $\mathcal{O}(N)$ times faster implementation of the exact OSil algorithm (effOSil). 
</p>

### Future releases
A new scalable approximation algorithm of OSil (scalOSil) for clustering and PAMSil will be implemented in the next releases.

### References

[Batool, F. and Hennig, C., 2021. Clustering with the average silhouette width. Computational Statistics & Data Analysis, 158, p.107190.](https://www.sciencedirect.com/science/article/abs/pii/S0167947321000244)

[Van der Laan, M., Pollard, K. and Bryan, J., 2003. A new partitioning around medoids algorithm. Journal of Statistical Computation and Simulation, 73(8), pp.575-584.](https://www.tandfonline.com/doi/abs/10.1080/0094965031000136012)

[Rousseeuw, P.J., 1987. Silhouettes: a graphical aid to the interpretation and validation of cluster analysis. Journal of computational and applied mathematics, 20, pp.53-65.](https://www.sciencedirect.com/science/article/pii/0377042787901257)



