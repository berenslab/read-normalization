suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(EDASeq))
suppressPackageStartupMessages(library(Matrix))
library(scRNAseq)
library(quminorm)
library(monocle)
library(scran)
library(glue)



path = glue('data//tasic//runtime_exp//tpm_transposed_100000.mtx')

print(glue('reading data from {path}'))

mtx_tpm = readMM(path)

print('make sce object')
sce = SingleCellExperiment(list(tpm=mtx_tpm))

print('compute qUMIs')
qumi_runtime = system.time(sce<-quminorm(sce,assayName="tpm",shape=2,mc.cores=40))
print('done. qUMI time:')
print(qumi_runtime)


#monocle/census in batches
print('preparing batches')
batchsize_max = 1000
batch_id <- ceiling(seq(nrow(colData(sce)))/batchsize_max)
batches <- lapply(unique(batch_id), function(x) sce[, batch_id == x])

print('compute census counts')
print('number of batches:')
print(length(unique(batch_id)))

census_runtime = system.time(
for (i in 1:length(batches)) {
    print('working on batch')
    print(i)
    sce_batch = batches[[i]]
    cds_batch = convertTo(sce_batch,"monocle",assay.type='tpm')
    census_batch = relative2abs(cds_batch,cores=40,verbose=TRUE,return_all=TRUE,method="num_genes")
    census_counts_batch = census_batch$norm_cds
    if (i==1){
        results = census_counts_batch
    }        
    else {
        results = cbind(results,census_counts_batch)
    }        
}
)
print('done. census time:')
print(census_runtime)

