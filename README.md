# wine_report_data_analysis
Code and report for the exam of Data Analysis. 1st year of MS in Data Science for Management.
To see the pictures run the code

 Casella Bruno 1000014143
UNIVERSITÀ DEGLI STUDI DI CATANIA
DIPARTIMENTO DI ECONOMIA E IMPRESA
CORSO DI LAUREA IN DATA SCIENCE FOR MANAGEMENT
  Report on “wine” dataset Library(ContaminatedMixt)
Data Analysis ACADEMIC YEAR 2019 – 2020
  
“WINE” DATASET
This report is about the “wine” dataset, contained in the R library “ContaminatedMixt”. This dataset is from the UCI machine learning repository and it is available at http://archive.ics.uci.edu/ml/datasets/Wine These data are the result of a chemical analysis of wines grown in the same region in Italy but derived from three different cultivars.
The analysis determined the quantities of 13 constituents found in each of the three types of wine: Barolo, Grignolino, Barbera.
This dataset contains 178 rows, each corresponding to a different cultivar of wine produced in Piedmont (Italy) and 14 columns.
The first column is the type of wine (Type), a factor variable with the following levels: Barolo, Grignolino, Barbera. The variables measured on the three types of wines are the following:
• Alcohol
• Malic
• Ash
• Alcalinity
• Magnesium • Phenols
• Flavanoids
• Nonflavanoids
• Proanthocyanins • Color
• Hue
• Dilution
• Proline
All variables but the label class are continuous.
The original dataset comprises 27 variables. Here a subset of 14 variables only has been included.
I attached the dataset to the R search path, and I displayed its structure.
 
 UNIVARIATE ANALYSIS
1. TYPE
“Type” is a nominal variable. It assumes different values for identifying 3 different types of wine: Barolo, Grignolino, Barbera.
As we can see, this dataset is not perfectly balanced with respect to these three categories, but this is not a problem. Indeed it should not cause any significant performance degradation, because the imbalance is very low. Almost 40% of values is in Grignolino, 33% in Barolo and 27% in Barbera, so the extent is similar.
2. ALCOHOL
 
“Alcohol” is a numeric and continuous variable.
 This variable assumes 126 values, from 11.03 to 14.83.
The mode assumes 2 different values (bimodal distribution): 12.37 and 13.05. The modal values are found using the names() function. I passed as argument the frequencyAlcohol; I used the square brackets for saying that when the argument of these brackets it’s the maximum, returns this value.
The coefficient variation is about 0.06%.
The distribution is platykurtic, so it has shorter tails than a normal distribution. The function used is kurtosis() from the EnvStats library.
 
By default the parameter excess is equal to TRUE, so in the first case the result is less than 0. Indeed, for a normal distribution the coefficient of kurtosis is equal to 3, and the coefficient of excess kurtosis is 0. Distributions with kurtosis less than 3 (excess kurtosis less than 0) are called platykurtic and they have shorter tails than a normal distribution. Distributions with kurtosis greater than 3 (excess kurtosis greater than 0) are called leptokurtic, and they have heavier tails than a normal distribution.
The distribution is left skewed. The function skew() is from the labstatR library. It assumes positive values for right skewness and negative values for left skewness.
    Above I present the boxplot of Alcohol. There isn’t any outlier (a data point that differs significantly from other observations) in the distribution.

We can find the best model to fit this distribution. Note that in this case of a single number giving the number of cells for the histogram in the parameter breaks, is only a suggestion to R.
The best way to model this kind of distribution is using mixture models.
The first model is a gamma mixture with number of distribution K equal to 2.
  
  
 Now I’m going to try to fit the distribution with K = 3.
 
  
 Comparing the Akaike information criterion (AIC) and the Schwarz information criterion (SBC), said also Bayesian information criterion (BIC), of the two models, the best model to fit the Alcohol distribution is the gamma mixture model with K = 2. Indeed, the models with the lowest AIC and SBC are preferred.
In the picture below there is a comparison between the 2 indexes.
 
3. MALIC
“Malic” is a numeric and continuous variable.
 This variable assumes 133 values, from 0.74 to 5.8. The mode is equal to 1.73.
The coefficient of variation is about 0.48%.
The distribution is leptokurtic because it has excess kurtosis greater than 0 and kurtosis greater than 3. It has heavier tails than a normal distribution.
 
 The distribution is right skewed because the function skew() returns positive values for right skewed distribution and negative values for left skewed distribution.
  The boxplot representation of Malic says that there are some outliers in the distribution. Their values are 5.80, 5.51 and 5.65.
Now we can find the best model to fit this distribution.
 
 Now I show different models to fit the Malic distribution
   
       Comparing AIC and SBC of the models, we can see which of them fits better.

We can see that the Inverse Gaussian distribution fits better the distribution. However we can try to fit this distribution using the gamma mixture model with K = 2.
  
 As we can see, both AIC and SBC have an improvement in this case, with respect to the Inverse Gaussian distribution.
Moreover, now I try with K = 3.
 
  
 Comparison between the models with K = 2 and K = 3
 
We can see that with K = 3 there is another improvement, so it is the best model to fit the Malic distribution.
4. ASH
“Ash” is a numeric and continuous variable.
 This variable assumes 178 values, from 1.36 to 3.23.
The mode assumes 2 different values (bimodal distribution): 2.28 and 2.3.
The coefficient of variation is about 0.12%.
The distribution is leptokurtic because it has excess kurtosis greater than 0 and kurtosis greater than 3. It has heavier tails than a normal distribution.
 
 The distribution is left skewed because the function skew() returns positive values for right skewed distribution and negative values for left skewed distribution.
   The boxplot representation of Ash says that there are some outliers in the distribution. Their values are 3.22, 1.36 and 3.23.
We can find the best model to fit this distribution.
 
 Different models to fit the Ash distribution:
 
        Comparing AIC and SBC of the models, we can see that the Logistic distribution fits better the distribution.
Now I try to fit the distribution using the gamma mixture model with K=2.

  
 AndnowItrywithK=3.
 
  
  In the last picture we can see that the model with K = 2 has a better BIC, but an AIC worse than K = 3. However, neither of them has better values with respect to the values of the Logistic distribution, although a better

graphical representation with K = 3. So the Logistic distribution is the best model to fit the Ash distribution.
5. ALCALINITY
“Alcalinity” is a numeric and continuous variable.
 This variable assumes 178 values, from 10.6 to 30. The modal value is 20.
The coefficient of variation is about 0.17%.
The distribution is leptokurtic because it has excess kurtosis greater than 0 and kurtosis greater than 3. It has heavier tails than a normal distribution.
The distribution is left skewed because the function skew() returns positive values for right skewed distribution and negative values for left skewed distribution.
Boxplot of Alcalinity. We can see that there are some outliers, of values 10.6, 30.0, 28.5 (x2).
  
  We can find the best model to fit this distribution.
 
 Different models to fit the Alcalinity distribution:
   
       
Comparing AIC and SBC of the models, we can see that the Logistic distribution fits better the distribution.
Now I try to fit the distribution using the gamma mixture model with K=2.
  
 I try with K = 3.
 
  
  As we can see, gamma model with K = 2 has better indexes with respect to K = 3. However, Logistic distribution has better values, so it is the best model to fit Alcalinity distribution.

6. MAGNESIUM
“Magnesium” is numeric and discrete variable.
 This variable assumes 178 values, from 70 to 162. The mode assumes the value 88.
The coefficient of variation is about 0.14%
The distribution is leptokurtic because it has excess kurtosis greater than 0 and kurtosis greater than 3. It has heavier tails than a normal distribution.
The distribution is positive skewed because the function skew() returns positive values for right skewed distribution and negative values for left skewed distribution.
  
  The boxplot representation of Magnesium says that there are some outliers, and their values are 151, 139, 136 and 162.
We can find the best model to fit this distributioin.
  
Different models to fit the Alcalinity distribution:
       
   Comparing AIC and SBC of the models, we can see that the Logistic distribution fits better the distribution.
Now I try to fit the distribution using the gamma mixture model with K=2.

  
 Logistic distribution has better values of AIC and SBC, so it is the best model to fit Magnesium distribution.
7. PHENOLS
“Malic” is a numeric and continuous variable.

 This variable assumes 178 values, from 0.98 to 3.88. The mode is equal to 2.2.
The coefficient of variation is about 0.27%.
The distribution is platykurtic because it has excess kurtosis lower than 0 and kurtosis lower than 3. It has shorter tails than a normal distribution.
The distribution is right skewed because the function skew() returns positive values for right skewed distribution and negative values for left skewed distribution.
  
The boxplot representation of Phenols says that there aren’t outliers in the distribution.
Now we can find the best model to fit this distribution.
    
Now I show different models to fit the Phenols distribution.
       
  Comparing AIC and SBC of the models, we can see which of them fits better.
 We can see that the Skew Power Exp. Type 2 distribution fits better the distribution. However we can try to fit this distribution using the gamma mixture model with K = 2.
 
  
As we can see, both AIC and SBC have an improvement in this case, with respect to the Skew Power Exp. Type 2 distribution. Moreover, now I try with K = 3.
  
  
Comparison between the models with K = 2 and K = 3.
 WecanseethatwithK=2isbetterthanthemodelwithK=3,soit is the best model to fit the Phenols distribution.
8. FLAVANOIDS
“Flavanoids” is a numeric and continuous variable.
 This variable assumes 178 values, from 0.34 to 5.08. The mode is equal to 2.65.
The coefficient of variation is about 0.49%.

 The distribution is platykurtic because it has excess kurtosis lower than 0 and kurtosis lower than 3. It has shorter tails than a normal distribution.
The distribution is right skewed because the function skew() returns positive values for right skewed distribution and negative values for left skewed distribution.
   The boxplot representation of Phenols says that there aren’t outliers in the distribution.

Now we can find the best model to fit this distribution.
  Now I show different models to fit the Phenols distribution.
   
      Comparing AIC and SBC of the models, we can see which of them fits better.
 
We can see that the Skew Power Exp. Type 2 distribution fits better the distribution, because has the lowest BIC value. The Weibull instead, has the lowest SBC. However we can try to fit this distribution using the gamma mixture model with K = 2.
  
 As we can see, both AIC and SBC have an improvement in this case, with respect to the Skew Power Exp. Type 2 distribution. Moreover, now I try with K = 3.

   
 Comparison between the models with K = 2 and K = 3.
 
WecanseethatwithK=3isbetterthanthemodelwithK=2,soit is the best model to fit the Flavanoids distribution.
9. NONFLAVANOID
“Nonflavanoid” is a numeric and continuous variable.
 This variable assumes 178 values, from 0.13 to 0.66.
The mode assumes 2 different values (bimodal distribution): 0.26 and 0.43.
The coefficient of variation is about 0.34%.
The distribution is platykurtic because it has excess kurtosis lower than 0 and kurtosis lower than 3. It has shorter tails than a normal distribution.
The distribution is right skewed because the function skew() returns positive values for right skewed distribution and negative values for left skewed distribution.
  
The boxplot representation of Phenols says that there aren’t outliers in the distribution.
  Now we can find the best model to fit this distribution.
  
Now I show different models to fit the Nonflavanoid distribution.
       
  Comparing AIC and SBC of the models, we can see which of them fits better.
 We can see that the Skew Power Exp. Type 2 distribution fits better the distribution. However we can try to fit this distribution using the gamma mixture model with K = 2.
 
  
As we can see, there isn’t any improvement in terms of AIC and SBC in this case, with respect to the Skew Power Exp. Type 2 distribution.
So, we can say that Skew Power Exp. Type 2 is the best model to fit the Nonflavanoid distribution.
10.PROANTHOCYANINS
“Proanthocyanins” is a numeric and continuous variable.
 This variable assumes 178 values, from 0.41 to 3.58. The mode is equal to 1.35.
The coefficient of variation is about 0.36%.
The distribution is leptokurtic because it has excess kurtosis greater than 0 and kurtosis greater than 3. It has heavier tails than a normal distribution.
 
 The distribution is right skewed because the function skew() returns positive values for right skewed distribution and negative values for left skewed distribution.
The boxplot representation of Phenols says that there are some outliers in the distribution and their value is 3.28 and 3.58.
  Now we can find the best model to fit this distribution.
 
 Now I show different models to fit the Phenols distribution.
    
    Comparing AIC and SBC of the models, we can see which of them fits better.
 We can see that the Gamma distribution fits better the distribution.

11.COLOR
Histogram of the Color distribution
 12.HUE
Histogram of the Hue distribution.
 13.DILUTION

Histogram of the Dilution distribution.
 14.PROLINE
Histogram of the Proline distribution.
 
UNIVARIATE ANALYSIS FOR EACH TYPE
The dataset contains 3 different types of wines grown in the same region in Italy but derived from three different cultivars. So, we can create three new datasets for each Type and then observe histograms of variable’s distribution.
   
  
   
  
  
   
  
  The 39.9% of the data are stored in Grignolino, 33,15% in Barolo and 26.95% in Barbera. We know that using the gamma mixture models using K = 2 or K = 3 can be useful to fit the distributions better. However, this histograms representation seem good, although it is better to use gamma mixture for minimizing AIC and SBC.

MULTIVARIATE ANALYSIS: PRINCIPAL COMPONENT ANALYSIS
Principal Component Analysis (PCA) involves the process by which principal components are computed, and their role in understanding the data. PCA is an unsupervised approach, which means that it is performed on a set of variables X1, X2, ..., Xn with no associated response Y. PCA reduces the dimensionality of the data set, allowing most of the variability to be explained using fewer variables. PCA is commonly used as one step in a series of analyses. We can use PCA to reduce the number of variables and avoid multicollinearity, or when we have too many predictors relative to the number of observations.
  
 I excluded from the correlation matrix the Type variable because it is categorical
We can see that there is a high positive correlation between Phenols and Flavanoids (0.8645635).
Now I use the function prcomp(), that performs PCA on the given data matrix. By default, prcomp() centers the variables to have mean zero. By using the option scale = TRUE, we scale the variables to have standard deviation one. In this way the variables will be comparable.
 
    
    
The “sdev” component contains the standard deviations of the principal components (i.e., the square roots of the eigenvalues of the covariance/correlation matrix).
The center and scale components correspond to the means and standard deviations of the original variables.
The “x” matrix provides the principal component scores.
The “rotation” matrix provides the principal component loading. Each column of pr.out$rotation contains the corresponding PC loading vector.
By default, eigenvectors in R point in the negative direction. We can adjust this with a simple change, a sign flip.
  
 The first principal component (PC1) has positive values for all but Malic, Ash, Alcalinity, Nonflavanoid and Color.
The second principal component (PC2) has negative values for all but Alcalinity, Flavanoids, Hue and Dilution.
The third principal component (PC2) has high negative contributes for Ash and Alcalinity.
Using the function fviz_pca_ind (from the factoextra package) we can obtain a graph for individual.
  
Now we can plot the first two principal components using biplot, a type of exploratory graph allowing information on both the n sample units and the d variables of a data matrix X to be displayed on a bidimensional Cartesian plane spanned by (usually the fist) two PCs. Principal components are displayed as axes. Sample units are displayed as points (PC scores). Original variables are displayed as arrows (PC loadings). The angle between variables reflects the correlation matrix: angle close to 0 mean high positive correlation, angle close to 90° mean absence of correlation, and angle of 180° mean perfect negative correlation.
  For example we can see that there is a high positive correlation between Nonflavanoid and Alcalinity, Ash and Color, Magnesium and Alcohol, Phenols and Proanthocyanins. Instead there is a negative correlation between Alcalinity and Nonflavanoid with respect to Phenols, Proanthocyanins, Flavanoids and Dilution, and between Malic and Hue and Dilution.
The parameter cex indicates the dimensions of the numbers and characters.

 The Kaiser’s rule suggests to retain as many principal components the ones with variance greater than 1 (for standardized data). So the first three principal components can represent the phenomena.
   
   
MULTIVARIATE ANALYSIS: CLUSTER ANALYSIS
Cluster Analysis (CA), simply said clustering, is one of the most important statistical methods for discovering knowledge in multidimensionality data. The goal of CA is to identify patterns (or groups cluster) of similar units with a data set X.
Before applying any clustering method on our data, it’s important to evaluate whether the data set contains meaningful clusters or not.
  We can see, from the graphical pairs representation that the data presents a grouping propensity, conversely to a random dataset generated from the wine dataset, presented below:
  
   
 It can be seen that the standardized wine data contain 3 clusters in the PC space. Note that I have used the 3 levels of the available nominal variable Type to color the points. Instead, it can be seen, in the PC space too, that the standardized randomly generated uniform data do not contain meaningful clusters.
For evaluating the clustering tendency, we can use a statistical method called Hopkins statistic, that takes values in [0,1].
The value H of this method close to 0 indicates clustered data; H close to 0.5 indicates uniformly distributed data (no meaningful clusters).
 
It can be seen that the wine data set (df) is clusterable because his H value (0.2830363) is close enough to 0. However, the random_df dataset is not clusterable (H = 0.5062973).
There is also a visual method for confirming the cluster tendency of a dataset, and it is the visual assessment of cluster tendency (VAT). It computes the dissimilarity matrix (DM) between the units in the dataset using the Euclidean distance, and it reorders it (the DM) so that similar units are close to one another. This process creates an ordered dissimilarity matrix (ODM) that is displayed as an ordered dissimilarity image (ODI), which is the visual output of the VAT algorithm. For the visual assessment of clustering tendency, we start by computing the dissimilarity matrix between observations using the function dist() and then using the function fviz_dist() to display the dissimilarity matrix.
  
 The color level is proportional to the value of the dissimilarity between observations: red denotes high similarity (i.e. low dissimilarity); blue denotes low similarity (i.e. high dissimilarity). The dissimilarity matrix image confirms that there is a cluster structure in the standardized wine dataset but not in the random one. The VAT algorithm detects the clustering tendency in a visual form by counting the number of a square shaped red blocks along the diagonal in a VAT image.
Now, we can proceed with the cluster analysis.
The first clustering algorithm that I applied is the agglomerative hierarchical clustering, that is a “bottom-up” approach: each observation starts in its own cluster (leaf), and pairs of cluster are merged as one moves up the hierarchy. This process goes on until there is just one single big cluster (root). In this case the leaves are 178, and for every step, the algorithm merges pairs of cluster with the lowest dissimilarity according to a different given linkage method, and then it moves up the hierarchy, until there is just one single big cluster, with all the observations within, by building a sort of tree diagram called dendrogram.

I will calculate the dissimilarity matrix using two different distances: the Euclidean, that is the most used, and the Manhattan, that reduces the effects of outliers.
I start using the Euclidean distance, and respectively single, average and complete linkage method.
   
    
    
  The single linkage method has the so called “chain effect”. Because of the chain effect, the method can bring together in a single cluster even very distant units when there is a succession of intermediate points between them. Because of the chain effect, the single linkage method is sensitive to outliers.
Thanks to the chain effect, the single linkage method has the advantage to handle groups with shapes which are different from the hyper-spherical ones, which are instead guaranteed by the complete linkage method.
The complete linkage method is less susceptible to noise and outliers, but it can break large cluster and can favors hyper- spherical shapes.
The average linkage method is an intermediate approach between the single and complete link approaches.
There are (agglomerative) hierarchical methods that also use the data matrix, and not just the distance/dissimilarity matrix. We will consider the centroid linkage method and the Ward’s (minimum deviance) method.
Single, average and complete linkage method use the distance (minimum, average and maximum) between the units of clusters.

Instead, centroid and Ward’s linkage methods assume that a cluster is represented by its centroid (cluster center). So they don’t use a unit, but a point, and for this reason they are feasible for numeric dataset only; Ward’s method minimizes iteratively the sum of the squared Euclidean distances of units within a cluster for every cluster (WSS – Within deviance), and maximizes the sum of the squared average distances of units between centroids (BSS – Between deviance). So, the Ward’s method accepts as input only a Euclidean distance matrix (the centroid method is often used with Euclidean distances too).
  
    
The difference between the methods Ward.D and Ward.D2 is the distance matrix to be given as an input to hclust(). In Ward.D it is squared, in Ward.D2 it is not squared.
We can see in the picture representing the dendrogram of the centroid linkage method that it is not monotonic. So called inversions occurred. This happens when similarity increase during clustering. For example when there are 3 points, d1, d2 and d3: the distance between d1 and d2 is the lowest, but then the distance of the centroid of d1 and d2 with d3 is lower than the previous distance. In brief, it happens when the union of two cluster 1 and 2 to a third cluster is less than the distance between 1 and 2. This case is not easy to understand.
The cophenetic dissimilarity/distance of two units is a measure of how similar those two units have to be in order to be grouped into the same cluster. From a practical point of view, the cophenetic distance between two units is the height of the dendrogram where the two branches that include the two units merge into a single branch (height of the fusion).
The function daisy() (in the cluster package) provides a solution (Gower’s distance) for computing the distance matrix, in the situation where the data contain non numeric columns, as in our case of the factor variable Type. The Gower’s distance is good to account categorical variables, such Type.
  
   After linking the units in a dataset into a hierarchical cluster tree, we might want to assess that the distances (i.e. heights) in the tree reflect the original distances accurately.

One way to measure how well the cluster tree generated by the hclust() function reflects our data is to compute the correlation between the cophenetic distances and the original distances generated by the dist() function. If the clustering is valid, the linking of units in the cluster tree should have a strong correlation with the distances between units in the original distance matrix.
The closer the value of the correlation coefficient is to 1, the more accurately the clustering solution reflects our data. Values above 0.75 are felt to be good.
In this case the average linkage method and the Gower’s distance produces the highest value of this statistic: 0.8742858.
Determining the optimal number of clusters is a fundamental issue. Unfortunately, there is no definitive solution to this issue. The optimal number of cluster is somehow subjective and depends on the clustering method used.
There are two methods to determine the optimal number of clusters: direct methods like Elbow and Silhouette, and statistical methods like Gap Statistic.
Elbow method roughly measures the quality of a clustering by determining how compact clusters are in terms of within-cluster sum of squares (WSS), a classical measure of compactness or cohesion.
The average silhouette method roughly measures the quality of a clustering by determining how well each unit lies within its cluster. It assumes values from -1 (units badly matched, perfect assignment for the neighboring cluster) to 1 (perfect assignment). We can compute the mean of the silhouette values in the dataset for each value of K, that is the number of clusters, e.g. between 1 and 10, and then we select that K which maximizes the silhouette value, a high average silhouette width indicates good clustering.
The Gap statistic provides a statistical procedure to formalize the heuristic elbow method.
 
       
For most combinations, according to the silhouette method, the best value of K is 2.
For each K, the gap statistic compares the WSS value observed to its expectation under an appropriate null reference distribution (i.e. a distribution with no obvious clustering, typically a uniform distribution). The estimate of the optimal number of cluster will be the value that maximizes the gap statistic. This means that the clustering structure is far away from the random uniform distribution of points.
     
   These charts indicate that the optimal number of clusters is 1.
Now I evaluate the results of the bar charts for the preferred value of K for all indices and I can see that the optimal number of cluster is 2.
 
       Now I will apply partitioning clustering methods as the k-means and the k-medoid.
K-means clustering is very simple and fast algorithm. K-means can efficiently deal with very large data sets. However K-means requires the analyst to choose the appropriate number of clusters

K in advance. Moreover, the final results obtained is sensitive to the initial random selection of cluster centers. This is a problem because for every different run of the algorithm on the same dataset, you may choose a different set of initial centers. This may lead to different clustering results on different runs of the algorithm. Moreover, K-means is sensitive to outliers.
K-means selects K points as centroids from the dataset and assigns every unit in the dataset to the closest one, according to Euclidean distance between the observation and the cluster mean, to form clusters, then it computes the mean for each cluster, using it as centroid, and iterates. At every iteration, it minimizes the increasing of the WSS value. It cannot be applied to datasets with categorical variables. We can see the results for the silhouette method and the gap statistics of the K-means algorithm.
     
The K-medoids algorithm, also known as partitioning around medoids (PAM), where each cluster is represented by one of the points in the cluster. PAM is less sensitive to outliers compared to K-means.
The most common K-medoids clustering method is the PAM algorithm, that is based on the search for K representative medoids among the observations of the data set. After finding a set of K medoids, clusters are constructed by assigning each observation to the nearest medoid. Next, the objective function is computed
The goal is to find K representative units which minimize the sum of the dissimilarities of the observations to their closest representative unit (medoid).
It can be applied according to the Euclidean distance and the Manhattan distance to have more robust results for data containing outliers, and it is applicable also for categorical datasets.
     According to these results, the optimal number of clusters is 3.

CLUSTER VALIDATION
Now, the most frequent values of the optimal number of clusters K are 1, 2 and 3.
I can evaluate the goodness of these results using cluster validation, that can be categorized in internal and external. The internal can be also used to find the number of clusters and the appropriate clustering algorithm without any external data. The external consists in comparing the results of a cluster analysis to an externally known result, such as externally provided class labels; since we know the “true” cluster number in advance, this approach is mainly used for selecting the right clustering algorithm for a specific data set.
Internal validation measures reflect often the compactness, connectedness and separation of the cluster partitions. Compactness or cluster cohesion measures how close are the units within the same cluster. A lower within-cluster variation is an indicator of a good compactness (i.e. a good clustering). Separation measures how we—separated a cluster is from other clusters. We will use two indexes which include both these measures: the silhouette width to identify the optimal K and the Dunn index that is the ratio between the minimal inter-cluster separation and the maximal intra-cluster distance, and it should be maximized.
The stability measures are a special version of internal measures, which evaluate the stability of a clustering result by comparing
it with the clusters obtained after each column (variable) is removed, one at a time.
Cluster stability measures include the average proportion of non- overlap (APN), the average distance (AD), the average distance between means (ADM) and the figure of merit (FOM).
APN measures the average proportion of observations not placed in the same cluster by clustering based on the full data and clustering based on the data with a single column removed. The AD measures the average distance between observations placed in the same cluster under both cases. The ADM measures the average distance between cluster centers for observations placed in the same cluster under both cases. The FOM measures the average intra-cluster variance of the deleted column, where the clustering is based on the remaining (undeleted) columns. In all cases the average is taken over all the deleted columns, and all measures should be minimized.

We will use the R package clValid, which contains the clValid() function for computing the above mentioned validation statistics for hierarchical clustering, k-means and k-medoids with K = 2, 3 for every combination of linkage method and distance metric.
  
   
   
   
  
  
    
    
 We can note that the value of K = 3 in general, is chosen more times than K = 2.
The hierarchical clustering with 2 clusters retains the best results for cohesion/separation, in particular Dunn index and connectivity for almost all combinations, and in terms of stability, in particular APN and ADM it is the best in both single linkage method and Euclidean and Manhattan distance. The k- means algorithm with 3 clusters best results in almost all combination in terms of stability APN, AD and ADM and in all combinations in terms of cohesion/separation, in particular Silhouette index.
In the external validation the aim is to compare the identified clusters (by k-means, PAM or hierarchical clustering) to an external reference. It’s possible to quantify the agreement between partitioning clusters and external reference using either the corrected version of Rand index or Meila’s variation index (VI), which are both implemented in the R function cluster.stats() (included in the fpc package).
The corrected Rand index varies from -1 (no agreement) to 1 (perfect agreement).
We know that the wine dataset contains exactly 3 types of wine. The question is: the k-means clustering matches with the true structure of the data?
Let start by computing a cross-tabulation (confusion matrix) between K-means clusters and the reference Type labels:

 It can be seen that: all Barbera Type (n=48) has been classified in cluster 1. All Barolo Type (n=59) has been classified in cluster 2. A large number of Grignolino Type (n=65) has been classified in cluster 3 and some of them have been classified in cluster 1 and 2 (n=3 respectively).
It’s possible to quantify the agreement between Type and K- means clusters using either the corrected Rand index and Meila’s VI provided as follows:
 The same analysis can be computed for both PAM and hierarchical clustering.
 
 In hierarchical clustering, the choice of a determined number of clusters corresponds to a certain height in the dendrogram, that will be cut according to it. I will evaluate the dendrograms’ cut for values of K=2 and K=3.
  
      Now I will cut the dendrograms with K = 3.

     
    
MODEL BASED CLUSTERING
One disadvantage of traditional clustering methods, such as hierarchical and partitioning clustering algorithms, is that they are largely heuristic and not based on formal models. This means that formal inference is not possible.
An alternative is model-based clustering, which considers the sample observations x1, ..., xn as coming from a distribution that is mixture of two or more, say K, distribution, commonly one for each cluster.
Each component, k=1,...,K, is described by a probability density or mass function and has an associated probability density or mass function and has an associated probability or “weight” in the mixture.
Unlike hierarchical and partitioning clustering methods, the model-based clustering uses a soft assignment, where each data point has a probability of belonging to each cluster.
  
 We can take a look at the information contained in the model and we can see that the BIC has a log-likelihood of -2292.525, that it uses 158 parameters, that the BIC value is -5403.772 and that the ICL (Integrated Completed Likelihood) is -5404.735, and that there are 3 clusters.
The three best model configurations are
We can take a look at the results of BIC also from a graphical point of view. These are the BIC curves for the models in our family. We have a different curve for each number of clusters and for each parsimonious configuration we have a symbol. We can see that the maximum is the one with 3 clusters and the symbol of the model is VVE. So in this way we have also a graphical counterpart of the results that we obtained above.
 
  
   
I plot the data using the colors of the classification vector:
  I compare the external information of the Type variable with the result of the model based clustering with the external cluster validation, the confusion matrix, and then with another external cluster validation, the adjusted Rand index:
 
Now I plot the objects showing the clustering:
   
  
