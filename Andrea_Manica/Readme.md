If you want to run the practical on your computer, you will have to copy over
the data (with sftp or scp), and install two `R` packages (they are not on cran, but you can get
precompiled versions from my personal r universe repository). In RStudio, simply
type:

```r
install.packages("tidypopgen", 
                 repos = c("https://evolecolgroup.r-universe.dev",
                           "https://cloud.r-project.org"))
install.packages("admixtools", 
                 repos = c("https://evolecolgroup.r-universe.dev",
                           "https://cloud.r-project.org"))
```

If you want to explore more what `tidypopgen` can do, look at its [website](https://evolecolgroup.github.io/tidypopgen/). For `admixtools`, there
is also a dedicated [website](https://evolecolgroup.github.io/tidypopgen/).
 
