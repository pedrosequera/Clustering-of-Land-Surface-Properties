# Clustering of land surface properties
 Creation of new Land-Cover-Land-Use Urban land classes from satellite sensor data using Cluster Large Applications (CLARA) in R 

Creation of new urban Land Cover Land Use (LCLU) subclasses result by use of the Clustering for Large Applications (CLARA) algorithm on the broadband albedos in R. The CLARA algorithm extends the k-medoids approach for a large number of objects by clustering a sample from the dataset and then assigning all objects in the dataset to these
clusters. 
The k-medoids approach is a partitioning algorithm that clusters a dataset of n objects into k clusters known a priori by choosing a data point (or medoid) for
each cluster at each iteration. Medoids for each cluster are calculated by finding a
point within each cluster that minimizes the distance between itself and the other
points of that cluster (Kaufman and Rousseeuw 1990).
The optimum number of subclasses was selected by using the average silhouette
criterion (Rousseeuw 1987).
The silhouette width is a measure of dissimilarity
between one point and all other points of the same cluster. Observations with a
large silhouette (almost 1) are very well clustered, whereas a small silhouette
(around 0) means that the observation lies between two clusters, and observations
with a negative silhouette are probably in the wrong cluster. Hence, the optimum
number of new urban classes is the one with the highest average silhouette.
