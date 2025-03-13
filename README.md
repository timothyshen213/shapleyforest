# `shapleyforest`: SHAPley Forest

UNDER DEVELOPMENT R Package

An interpretable Random-Forest based feature selection method through SHAPley values.

```
devtools::install_github("timothyshen213/shapleyforest")
```

### A note on downloading `WGCNA`:

First, download `BiocManager`
```
if (!require("BiocManager", quietly = TRUE))
       install.packages("BiocManager")
BiocManager::install(version = "3.20")
```

Before downloading `WGCNA`, install its dependencies.

````
BiocManager::install(c("GO.db", "preprocessCore", "impute") )
```

Finally, download `WGCNA`

```
install.packages("WGCNA")
```
