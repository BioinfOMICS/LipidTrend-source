## ----style, echo=FALSE, results='asis'----------------------------------------
BiocStyle::markdown(css.files = c('custom.css'))
knitr::opts_chunk$set(collapse=TRUE, comment="#>")

## ----install_Bioconductor, eval=FALSE-----------------------------------------
# if (!requireNamespace("BiocManager", quietly=TRUE))
#     install.packages("BiocManager")
# BiocManager::install()

## ----install_package, eval=FALSE----------------------------------------------
# if (!requireNamespace("BiocManager", quietly=TRUE))
#     install.packages("BiocManager")
# BiocManager::install("LipidTrend")

## ----load, message=FALSE------------------------------------------------------
library(LipidTrend)

## ----abundance----------------------------------------------------------------
# load example abundance data
data("abundance_2D")

# view abundance data
head(abundance_2D, 5)

## ----char_table---------------------------------------------------------------
# load example lipid characteristic table (2D)
data("char_table_2D")

# view lipid characteristic table
head(char_table_2D, 5)

## ----groupInfo----------------------------------------------------------------
# load example group information table
data("group_info")

# view group information table
group_info

## ----build_se-----------------------------------------------------------------
se_2D <- SummarizedExperiment::SummarizedExperiment(
    assays=list(abundance=abundance_2D),
    rowData=S4Vectors::DataFrame(char_table_2D),
    colData=S4Vectors::DataFrame(group_info))

## ----set_seed-----------------------------------------------------------------
set.seed(1234)

## ----LipidTrend_1D------------------------------------------------------------
# load example data
data("lipid_se_CL")

# quick look of SE structure
show(lipid_se_CL)

## ----countRegion1D------------------------------------------------------------
# run analyzeLipidRegion
res1D <- analyzeLipidRegion(
    lipid_se_CL, ref_group="sgCtrl", split_chain=TRUE, 
    chain_col="chain", radius=3, own_contri=0.5, test="t.test", 
    abund_weight=TRUE, permute_time=1000)

# view result summary 
show(res1D)

## ----Region1D_result----------------------------------------------------------
# view even chain result (first 5 lines)
head(even_chain_result(res1D), 5)

# view odd chain result (first 5 lines)
head(odd_chain_result(res1D), 5)

## ----plotRegion1D-------------------------------------------------------------
# plot result
plots <- plotRegion1D(res1D, p_cutoff=0.05, y_scale='identity')

# even chain result
plots$even_result
# odd chain result
plots$odd_result

## ----LipidTrend_2D------------------------------------------------------------
# load example data
data("lipid_se_2D")

# quick look of SE structure
show(lipid_se_2D)

## ----countRegion2D------------------------------------------------------------
# run analyzeLipidRegion
res2D <- analyzeLipidRegion(
    lipid_se_2D, ref_group="sgCtrl", split_chain=TRUE, 
    chain_col="Total.C", radius=3, own_contri=0.5, test="t.test", 
    abund_weight=TRUE, permute_time=1000)

# view result summary 
show(res2D)

## ----Region2D_result----------------------------------------------------------
# view even chain result (first 5 lines)
head(even_chain_result(res2D), 5)
# view odd chain result (first 5 lines)
head(odd_chain_result(res2D), 5)

## ----LipidTrend_2D_plot-------------------------------------------------------
# plot result
plot2D <- plotRegion2D(res2D, p_cutoff=0.05)

# even chain result
plot2D$even_result
# odd chain result
plot2D$odd_result

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

