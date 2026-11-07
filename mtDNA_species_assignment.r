setwd("C:/Users/topalw/Desktop/ants/Behavioural/genetic_analyses")

library(pegas)
library(ape)
cols5m <- c("#E46B6B", "#527FBD", "#38C4B9","#55A85B" , "#E9B82F", "#A95AAE")
s <- read.table('1.data/samples_species_bamorder.txt')
b <- read.table('1.data/full_observed_nests.txt',sep='\t',h=T)
# read file
seqs <- read.dna("1.data/for_mega_full_aligned.fas", format = "fasta")
seqmat <- as.character(as.matrix(seqs))  # character matrix, rows=individuals, cols=sites
n_ind <- nrow(seqmat)

# count Ns per site
callable <- apply(seqmat, 2, function(col) sum(col %in% c("a", "c", "g", "t")))

# filter on 100% callability -- sites genotyped in every single sample
full_sites <- which(callable == n_ind)
cat(length(full_sites), "out of", ncol(seqmat), "sites fully callable\n")

seqs_filtered <- seqs[, full_sites]
h <- haplotype(seqs_filtered)
net <- haploNet(h)
ids <- rownames(seqs_filtered)

species <- rep("unknown", length(ids))
idx <- match(ids, s$V1)
species[!is.na(idx)] <- ifelse(s$V2[idx[!is.na(idx)]] == "lugubris", "lugubris", "paralugubris")
species <- factor(species, levels = c("paralugubris", "lugubris", "unknown"))

ind.hap <- with(
  stack(setNames(attr(h, "index"), rownames(h))),
  table(hap = ind, pop = species[values])
)

pdf('mtDNA_assignment.pdf',6,6)
plot(net, size = attr(net, "freq"), pie = ind.hap,
     bg = c(cols5m[1], cols5m[2], "grey70"),   # matches factor levels: paralugubris, lugubris, unknown
     labels = FALSE, show.mutation = 1, threshold = 0)
legend("topright", legend = levels(species), fill = c(cols5m[1], cols5m[2], "grey70"))
dev.off()

hap_samples <- lapply(attr(h, "index"), function(idx) ids[idx])
names(hap_samples) <- rownames(h)
hap_samples

ids2 <- ids[ids %in% b$Nest_ID]
out <- data.frame('id' = ids2, 
                  'sp' = ifelse(ids2 %in% hap_samples$I, 'paralugubris','lugubris'))
out$id <- gsub('\\(2\\)','',out$id)
write.table(out,'1.data/species_assignment_mtDNA.txt',quote=F,row.names=F)
