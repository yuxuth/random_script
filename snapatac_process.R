library(SnapATAC)
library(tidyverse)
library(GenomicRanges)

## load the snap file directory

pdf_dir = './pdf'
if(! dir.exists(pdf_dir)) dir.create(pdf_dir)

file.list = 'Mouse-cortex_2nd.snap'
sample.list = c("Mouse-cortex")
x.sp = createSnap(file=file.list, sample=sample.list, do.par =T, num.cores = 24 )

# barcode selection qc 
pdf(file.path( pdf_dir, 'quliaty_score.pdf'), width =  8, height =  6)
plotBarcode(x.sp, col="blue", border="grey")
dev.off()

summarySnap(x.sp)


## filter cell by tge UMI and fragment 
x.sp <- filterCells(
    obj=x.sp,
    subset.names=c( "UMI"),
    low.thresholds=c(200), ## set 200 to test first 
    high.thresholds=c( Inf)
)
x.sp


x.sp = addBmatToSnap(
    obj=x.sp, 
    bin.size=5000,  do.par = T, 
    num.cores=24
)


## examine the promoter ration
library(Rsubread)
mm10_promoter = promoterRegions(annotation="mm10", upstream=1000L, downstream=1000L)

promoter.gr = GRanges(mm10_promoter[,2], IRanges(mm10_promoter[,3], mm10_promoter[,4]))
ov = findOverlaps(x.sp@feature, promoter.gr)
idy = queryHits(ov)


promoter_ratio = SnapATAC::rowSums(x.sp[,idy, mat="bmat"], mat="bmat") / SnapATAC::rowSums(x.sp, mat="bmat")

promoter = SnapATAC::rowSums(x.sp[,idy, mat="bmat"], mat="bmat")
row_sum = SnapATAC::rowSums(x.sp, mat="bmat")

pdf(file.path( pdf_dir, 'reads_in_promoter_ratio.pdf'), width =  8, height =  6)
plot(
    x=log(row_sum + 1,10), 
    y=promoter_ratio, 
    cex=0.5, 
    col="grey", 
    xlab="log(count)", 
    ylab="FIP Ratio",
    ylim=c(0,1 )
)
title( main = 'Reads in promoter ratio')

idx = which(promoter_ratio > 0.2 )

x.sp = x.sp[idx,]

x.sp

# system("wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz")

black_list = read.table("mm10.blacklist.bed.gz");
black_list.gr = GRanges(
    black_list[,1], 
    IRanges(black_list[,2], black_list[,3])
);
idy1 = queryHits(findOverlaps(x.sp@feature, black_list.gr));  ## feature is the 5k bin 
idy2 = grep("chrM|random", x.sp@feature); ## filter out some unwanted chr
idy = unique(c(idy1, idy2));
x.sp = x.sp[,-idy, mat="bmat"]
idy3 = grep("chrUn", x.sp@feature); ## filter out some unwanted chr
x.sp = x.sp[,-idy3, mat="bmat"]

x.sp

x.sp = makeBinary(x.sp, mat="bmat")

## bin coverage
plotBinCoverage(
    obj=x.sp,
    col="grey",
    border="grey",
    breaks=20,
)

bin.cov = log10(Matrix::colSums(x.sp@bmat)+1)
hist(
    bin.cov[bin.cov > 0], 
    xlab="log10(bin cov)", 
    main="log10(Bin Cov)", 
    col="lightblue", 
    xlim=c(0, 5)
)

bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95)  ## remove the top 5% feature as the 

idy = which(bin.cov <= bin.cutoff & bin.cov > 0)
x.sp = x.sp[, idy, mat="bmat"]
x.sp

idx = which(Matrix::rowSums(x.sp@bmat) > 100)  ## filter by bin coverage
x.sp = x.sp[idx,]
x.sp


## diffusionmaps

x.sp = runDiffusionMaps(
    obj=x.sp,
    input.mat="bmat", 
    num.eigs=50
)  ## include the diffusion map 

plotDimReductPW(
    obj=x.sp, 
    eigs.dims=1:50,
    point.size=0.3,
    point.color="grey",
    point.shape=19,
    point.alpha=0.6,
    down.sample=5000,
    pdf.file.name=NULL, 
    pdf.height=7, 
    pdf.width=7
)

x.sp = runKNN(
    obj=x.sp,
    eigs.dims=1:20,
    k=15
)

# x.sp=runCluster(
#     obj=x.sp,
#     tmp.folder=tempdir(),
#     louvain.lib="R-igraph",
#     seed.use=10
# )

# library("reticulate")  ## install the required python lib
# py_config()
# py_install("python-igraph")
# py_install("leidenalg", forge = TRUE)

library(leiden)
x.sp=runCluster(
    obj=x.sp,
    tmp.folder=tempdir(),
    louvain.lib="leiden",
    seed.use=10,
    resolution = 0.5
)

x.sp@metaData$cluster = x.sp@cluster;


x.sp = runViz(
    obj=x.sp, 
    tmp.folder=tempdir(),
    dims=2,
    eigs.dims=1:20, 
    method="tsne",
    seed.use=10
)

eigs.dims=1:20
umap_dime = umap::umap(x.sp@smat@dmat[,eigs.dims], method = 'umap-learn')
x.sp@umap = umap_dime$layout


plotViz(
    obj=x.sp,
    method="tsne",
    point.size=0.5,
    point.shape=19,
    point.alpha=0.8,
    point.color=x.sp@cluster,
    text.add=TRUE,
    text.size=1.2,
    text.color="black",
    text.halo.add=TRUE,
    text.halo.color="white",
    text.halo.width=0.2,
    # down.sample=10000,
    pdf.file.name=NULL,
    pdf.width=7,
    pdf.height=7,
    legend.add=T
)


## generate track


fn_keep_minimal <- function(input_snap){
    new_snap = newSnap()
    new_snap@des = input_snap@des
    new_snap@metaData = input_snap@metaData
    new_snap@file = input_snap@file
    new_snap@sample = input_snap@sample
    new_snap@barcode = input_snap@barcode
    new_snap
}



info_test = fn_keep_minimal(x.sp)

# rm(x.sp)
gc()

## test the barcode extract information 

clusters.sel = names(table(info_test@metaData$cluster))[which(table(info_test@metaData$cluster) > 40)]

peaks.ls = mclapply(seq(clusters.sel), function(i){
    print(clusters.sel[[i]])
    runMACS(
        obj=info_test[which(info_test@metaData$cluster %in% clusters.sel[[i]]),], 
        output.prefix= paste0("cell_type_track/mouse_cortex_", gsub(" |\\/", "_", clusters.sel)[i]),
        # path.to.snaptools="/apps/miniconda3/envs/py2/bin/snaptools",
        path.to.snaptools="/apps/miniconda3/bin/snaptools",
        path.to.macs="/apps/miniconda3/envs/py2/bin/macs2",
        gsize="mm", # mm, hs, etc
        buffer.size= 10000, 
        num.cores=1,
        macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",
        tmp.folder= './tmp_cell_types'
    )
}, mc.cores=13)


## convert the bdg to bw, which is ready for igv 

