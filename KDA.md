KDA in Ovals dataset
================
Roberto Rey
13/3/2020

# KDA

We load the data

``` r
# Load data
load("ovals.RData")
head(ovals)
```

    ##             x.1         x.2 labels
    ## 2963  1.7422628 -0.05399626      3
    ## 1320 -1.5463537  2.71793820      2
    ## 2155 -1.3282129 -2.53476352      3
    ## 1869 -2.0215933  0.88696110      2
    ## 1219 -1.0848312  2.75145649      2
    ## 736  -0.2747645 -2.19227121      1

Now we have to plot the data and split between a training sample and a
test sample.

``` r
plot(ovals$x.1,ovals$x.2,col = ovals$labels)
```

![](KDA_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

It is going to be hard because they do not look like separate groups
where the intersection between ovals occur, those points are probably
going to be assigned to wrong groups.

``` r
indx = 2000
train = ovals[1:indx,]
test = ovals[-c(1:indx),]
```

We separate our train dataset and we calculate our bandwidth matrix.

``` r
x = train[,1:2]
groups = train$labels

# Manual specification of bandwidths via ks::Hkda
Hs <- ks::Hkda(x = x, x.group = groups, bw = "plugin")
```

Now we perform kernel discriminant analyisis with the bandwidth obtained
by plugin.

``` r
kda <- ks::kda(x = x, x.group = groups, Hs = Hs)
```

We plot kernel discriminant analyisis

``` r
# Plot of classification regions
plot(kda, col = rainbow(3), lwd = 2, col.pt = 1, cont = seq(5, 85, by = 20),
     col.part = rainbow(3, alpha = 0.25), drawpoints = TRUE)
```

![](KDA_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
# Classification error
error = ks::compare(x.group = kda$x.group, est.group = kda$x.group.estimate)
success = 1-error$error
success
```

    ## [1] 0.731

``` r
lda <- MASS::lda(labels ~ x.1 + x.2, data = train)
pred_lda <- predict(lda, newdata = test)
mean(pred_lda$class != test$labels)
```

    ## [1] 0.664

So as we can see the precision is 0.663 in the LDA, it is smaller than
the KDA which is 0.731.

``` r
qda <- MASS::qda(labels ~ x.1 + x.2, data = train)
pred_qda <- predict(qda, newdata = test)
mean(pred_qda$class != test$labels)
```

    ## [1] 0.566

So as we can see the precision is 0.566 in the QDA, it is smaller than
the KDA which is 0.731 and it is also smaller than the LDA one, so QDA
seems to perform the worst in this dataset.
