<style type="text/css">
h1.title {
  text-align: center;
}
h3.subtitle {
  text-align: center;
}
h4.author {
  text-align: center;
}
</style>
## Theory

### Data setup

Let *X* ∈ ℝ<sup>*n* × *p*</sup> be a column-centered data matrix with
*n* samples in the rows and *p* features in the columns. The covariance
matrix of the features is given by
$$C=\\frac{X^\\top X}{n-1} \\in\\mathbb{R}^{p\\times p},$$
and satisfies
*C*<sub>*i**j*</sub> = Cov(*x*<sub>*i*</sub>,*x*<sub>*j*</sub>),
where *x*<sub>*i*</sub> is the *i*-th column of X. It is important that
*X* is column-centered. Otherwise, *C* is *not* the covariance matrix of
the features.

### Principal axes and components

Since *C* is a covariance matrix, it is positive semi-definite and can
thus be diagonalized as
*C* = *V**L**V*<sup>⊤</sup>
where *L* is a diagonal matrix of non-negative eigenvalues
*λ*<sub>1</sub>, …, *λ*<sub>*n*</sub> and *V* is an orthonormal matrix,
i.e. *V*<sup>⊤</sup>*V* = 𝟙<sub>*p*</sub>, of eigenvectors, i.e. the
*i*-th column *v*<sub>*i*</sub> of *V* satisfies
*C**v*<sub>*i*</sub> = *λ*<sub>*i*</sub>*v*<sub>*i*</sub>. The vectors
*v*<sub>1</sub>, …, *v*<sub>*n*</sub> are called the *principal axes* or
*principal directions* of the data. Projections of the data onto the
principal axes are called *principal components*. Since the columns of
*V* have unit norm, the entry (*X**V*)<sub>*i**j*</sub> is the
projection of the *i*-th sample (row) of *X* onto the *j*-th principal
axis (column) of *V*. Thus, the *i*-th row of *X**V* is the projection
of the *i*-th sample on the principal axes and hence the coordinates of
this sample in the new coordinate system given by the principal axes.
The *j*-th column of *X**V* is the *j*-th principal component.

### Singular value decomposition

As for any matrix, there exits a *singular value decomposition* (SVD)
*X* = *U**D**V*<sup>⊤</sup>
where *U* is a *n* × *n* orthonormal matrix, *D* a *n* × *p* matrix with
only non-zero values on the main diagonal, and *V* is a *p* × *p*
orthonormal matrix. The entries on the main diagonal of *D* are called
singular values and are the square-roots of the eigenvalues of
*X*<sup>⊤</sup>*X* and *X**X*<sup>⊤</sup>. The columns of *U* (of *V*)
are the left-singular (right-singular) vectors of *X*, i.e. the
eigenvectors of *X**X*<sup>⊤</sup> (of *X*<sup>⊤</sup>*X*).

It follows from the definition of the SVD that
$$C=V\\frac{D^2}{n-1}V^\\top,$$
which implies that the right-singular vectors *V* are the principal
axes, the eigenvalues *λ*<sub>*i*</sub> and singular values
*d*<sub>*i*</sub> are related via
*λ*<sub>*i*</sub> = *d*<sub>*i*</sub><sup>2</sup>/(*n*−1) and the
principal components are given by *X**V* = *U**D*. One can achieve a
lower rank (rank *k*) approximation of *X* by multiplying the first *k*
principal components (*U**D*)<sub>*k*</sub> with the first *k* principal
axes *V*<sub>*k*</sub><sup>⊤</sup>.

If *n* &gt; *p*, the last *n* − *p* rows of *D* are zero and thus the
last *n* − *p* columns of *U* can be arbitrary and do not contribute to
the factorization. In that case, a reduced SVD can be achived where *D*
is a *p* × *p* and *U* a *n* × *p* matrix. Likewise, if *n* &lt; *p*, a
reduced SVD has *D* as *n* × *n* matrix and *V* as *p* × *n* matrix.

### Variance decomposition

Since *X* is column-centered, the sum of variances of all features is
given by
$$
\\mathbb{V} =\\frac{1}{n-1}\\sum\_{j=1}^p\\sum\_{i=1}^n(X\_{ij})^2
=\\mathop{\\mathrm{tr}}(C)
=\\sum\_{j=1}^p \\frac{d\_j^2}{n-1}
=\\sum\_{j=1}^p \\lambda\_j.
$$

The second equality holds because the trace of a square matrix is the
sum of its diagonal entries. The third equality results from the trace
being the sum of the eigenvalues and
{*d*<sub>1</sub>, …, *d*<sub>*p*</sub>} being the eigenvalues of
*X*<sup>⊤</sup>*X*. This shows that the variance of the data *X* is
decomposed by its singular values contained in *D*, and by the
eigenvalues of *X*<sup>⊤</sup>*X* contained in *L*.

The variance attributed to each principal component can be further
decomposed onto the features: A rank 1 approximation of *X* from the
first principal axis and component is given by the matrix *X̃* with
*X̃*<sub>*i**j*</sub> = (*U**D*)<sub>*i*1</sub>*V*<sub>*j**i*</sub> = *d*<sub>1</sub>*U*<sub>*i*1</sub>*V*<sub>*j*1</sub>.
The variance of feature *j* satisfies
$$\\mathrm{Var}\\left(\\tilde X\_{.j}\\right) = \\frac{d\_1^2}{n-1}(V\_{j1})^2.$$
Together with the fact that $\\sum\_{j=1}^p (V\_{j1})^2$, it implies
that the squared entries of the principal axis vector contain the
fraction of variance that is explained by the respective feature.

## Practical example

### The data

We analyze the `iris` dataset with principal component analysis (PCA)
and SVD using *R*’s `prcomp` and `svd` functions.

    # load data
    data <- as.matrix(iris[,1:4])

    # setup variables as in theory section
    X <- data
    n <- nrow(data)
    rownames(X) <- paste0("sample", 1:n)

    # center the features
    X <- scale(X, center = T, scale = F)

    # do PCA and SVD
    pca <- prcomp(X)
    sv <- svd(X)

    head(X)

    ##         Sepal.Length Sepal.Width Petal.Length Petal.Width
    ## sample1   -0.7433333  0.44266667       -2.358  -0.9993333
    ## sample2   -0.9433333 -0.05733333       -2.358  -0.9993333
    ## sample3   -1.1433333  0.14266667       -2.458  -0.9993333
    ## sample4   -1.2433333  0.04266667       -2.258  -0.9993333
    ## sample5   -0.8433333  0.54266667       -2.358  -0.9993333
    ## sample6   -0.4433333  0.84266667       -2.058  -0.7993333

### Results of `prcomp`

    # principal axes
    head(pca$rotation)

    ##                      PC1         PC2         PC3        PC4
    ## Sepal.Length  0.36138659 -0.65658877  0.58202985  0.3154872
    ## Sepal.Width  -0.08452251 -0.73016143 -0.59791083 -0.3197231
    ## Petal.Length  0.85667061  0.17337266 -0.07623608 -0.4798390
    ## Petal.Width   0.35828920  0.07548102 -0.54583143  0.7536574

    # projection of X onto principal axes
    head(pca$x)

    ##               PC1        PC2         PC3          PC4
    ## sample1 -2.684126 -0.3193972  0.02791483  0.002262437
    ## sample2 -2.714142  0.1770012  0.21046427  0.099026550
    ## sample3 -2.888991  0.1449494 -0.01790026  0.019968390
    ## sample4 -2.745343  0.3182990 -0.03155937 -0.075575817
    ## sample5 -2.728717 -0.3267545 -0.09007924 -0.061258593
    ## sample6 -2.280860 -0.7413304 -0.16867766 -0.024200858

### Compare `prcomp` and `svd`

    # don't compare object attributes
    all_equal <- function(...) all.equal(..., check.attributes = FALSE)

    # principal axes (principal directions) 
    all_equal(pca$rotation, sv$v)

    # projections of X onto principal axes
    all_equal(pca$x, X %*% pca$rotation)
    all_equal(pca$x , sv$u %*% diag(sv$d))

    # full data reconstruction
    all_equal(X, pca$x %*% t(pca$rotation))

    # variance decomposition
    all_equal(pca$sdev^2, ((sv$d)^2 /(n-1)))
    all_equal(sum(pca$sdev^2), sum(X^2)/(n-1))

    # eigenvalues
    lambda <- eigen((t(X) %*% X)/(n-1))$values
    all_equal(pca$sdev^2, lambda[lambda>0])

    ## [1] TRUE
    ## [1] TRUE
    ## [1] TRUE
    ## [1] TRUE
    ## [1] TRUE
    ## [1] TRUE
    ## [1] TRUE

### Low-rank approximation from principal components

    # number of principal components to use
    n_comp = 3
    # reconstruction of column-centered X
    X_approx <-  pca$x[,1:n_comp] %*% t(pca$rotation[,1:n_comp])

    # add original column means to approximate original data
    data_approx <-  scale(X_approx, center = -colMeans(data), scale = FALSE)

    # check approximation error
    all_equal(data_approx, data)

    ## [1] "Mean relative difference: 0.01562306"

### Loadings decompose variance by feature

    # 1D approximation of the data
    projections <- pca$x
    loadings <- pca$rotation
    X1 <- projections[,1] %*% t(loadings[,1])

    # feature-wise variances
    vars <- matrixStats::colVars(X1, useNames = T)

    # relative proportions of feature-wise variances 
    # are the squares of the loadings
    all_equal(sum(vars), (pca$sdev^2)[1])
    all_equal(vars / sum(vars), loadings[,1]^2)

    ## [1] TRUE
    ## [1] TRUE

## References

\[1\]
<https://stats.stackexchange.com/questions/134282/relationship-between-svd-and-pca-how-to-use-svd-to-perform-pca>
