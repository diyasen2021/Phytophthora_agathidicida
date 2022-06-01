# First install
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(version = '3.14', ask= FALSE, "~/R/x86_64-pc-linux-gnu-library/4.0")
#BiocManager::install("karyoploteR", force = TRUE)

# Update karyoploteR
#BiocManager::valid()$out_of_date
#BiocManager::install()

library(karyoploteR)

# This makes linear plots for multiple genomes for read density, gc, repeats and allele frequency

custom.genome <- toGRanges("mygenome.txt")
gc.cont <- toGRanges("gc10000.bed")
repeats <- toGRanges("repeat.txt")
AF <- toGRanges("Allele_frequency-1.txt")
rxlr <- toGRanges("Paga_3770v2.rxlrs_forcircos.bed")
head(rxlr)

#png("plot.png", width = 1500, height = 700)
# Ideogram
kp <- plotKaryotype(genome = custom.genome, plot.type = 4, labels.plotter = NULL)
kpAddChromosomeNames(kp, cex=0.5, srt=45, yoffset = -2)
kp <- kpAddBaseNumbers(kp, tick.dist = 2e6, cex=0.5, units="Mb", tick.len = 3)

# rxlr
at <- autotrack(current.track = 1, total.tracks = 19)
kpDataBackground(kp, r0=at$r0, r1=at$r1)
kpAddLabels(kp, labels = "rxlr", r0=at$r0, r1=at$r1, cex=0.5, side="right")
kp<-kpPlotDensity(kp, data=rxlr, window.size = 1e4, col="#AAAAAA", r0=at$r0, r1=at$r1)

# GC
at <- autotrack(current.track = 2, total.tracks = 19)
kpDataBackground(kp, r0=at$r0, r1=at$r1)
kpAddLabels(kp, labels = "GC%", r0=at$r0, r1=at$r1, cex=0.5, side="right")
kpLines(kp,data=gc.cont, y=as.numeric(gc.cont$name), r0=at$r0, r1=at$r1, col="darkgreen")
kpAxis(kp, r0=at$r0, r1=at$r1, cex=0.3)

# Repeats
at <- autotrack(current.track = 3, total.tracks = 19)
kpDataBackground(kp, r0=at$r0, r1=at$r1)
kpAddLabels(kp, labels = "Repeats", r0=at$r0, r1=at$r1, cex=0.5, side="right")
kp<-kpPlotDensity(kp, data=repeats, window.size = 1e4, col="#AAAAAA", r0=at$r0, r1=at$r1)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density,  r0=at$r0, r1=at$r1, cex=0.3)
# Allele Frequency
at <- autotrack(current.track = 4, total.tracks = 19)
kpDataBackground(kp, r0=at$r0, r1=at$r1)
kpAddLabels(kp, labels = "AF", r0=at$r0, r1=at$r1, cex=0.5, side="right")
kp<-kpPoints(kp,data=AF, y=as.numeric(AF$V4), r0=at$r0, r1=at$r1,  col="#FFAAAA44")
kpAxis(kp, ymin=0, ymax=1, r0=at$r0, r1=at$r1, cex=0.3)

# Bam read density
myFiles<-list.files(path="../gatk", pattern="*coordsort.bam$", full.names = TRUE)
myFiles
sampleid <- read.table("samplename.txt", header = TRUE)
id <-sampleid[,1]
j=0
for (f in 1:length(myFiles)){
  j=f+3
  at <- autotrack(current.track = j, total.tracks = length(myFiles))
  kpDataBackground(kp, r0=at$r0, r1=at$r1)
  kpAddLabels(kp, labels = id[f], r0=at$r0, r1=at$r1, cex=0.5, side="right")
  kp<-kpPlotBAMDensity(kp, data=myFiles[f], ymax=50000, window.size=1e4, r0=at$r0, r1=at$r1, col="darkorange")
  kpAxis(kp, ymax=50000, r0=at$r0, r1=at$r1, cex=0.3)
  j=j+1
}


