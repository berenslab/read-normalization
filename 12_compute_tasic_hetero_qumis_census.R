suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(EDASeq))
suppressPackageStartupMessages(library(Matrix))
library(scRNAseq)
library(quminorm)
library(monocle)
library(scran)

print('reading data')
mtx_tpm = readMM('data//tasic//tasic_tpm_transposed.mtx')

print('make sce object')
sce = SingleCellExperiment(list(tpm=mtx_tpm))
cellnames = read.csv('data//tasic//tasic_tpm_cells.csv')$cells
genenames = read.csv('data//tasic//tasic_tpm_genes.csv')$gene_name
rownames(sce) = genenames
colnames(sce) = cellnames

print('compute qUMIs')
system.time(sce<-quminorm(sce,assayName="tpm",shape=2,mc.cores=16))
qumi=assay(sce,"qumi_poilog_2")

print('write out qUMIs')
writeMM(qumi,'data//tasic//qumi_transposed.mtx')


#monocle/census in batches
print('preparing batches')
batchsize_max = 1000
batch_id <- ceiling(seq(nrow(colData(sce)))/batchsize_max)
batches <- lapply(unique(batch_id), function(x) sce[, batch_id == x])

print('compute census counts')
print('number of batches:')
print(length(unique(batch_id)))

system.time(
for (i in 1:length(batches)) {
    print('working on batch')
    print(i)
    sce_batch = batches[[i]]
    cds_batch = convertTo(sce_batch,"monocle",assay.type='tpm')
    census_batch = relative2abs(cds_batch,cores=16,verbose=TRUE,return_all=TRUE,method="num_genes")
    census_counts_batch = census_batch$norm_cds
    if (i==1){
        results = census_counts_batch
    }        
    else {
        results = cbind(results,census_counts_batch)
    }        
}
)

print('write out census counts')
writeMM(results,'data//tasic//census_transposed.mtx')

