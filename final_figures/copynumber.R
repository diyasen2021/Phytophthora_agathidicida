library(karyoploteR)

# Read in the genome file
custom.genome <- toGRanges("mygenome_nomt.txt") # remove mt
chr.names<- as.character(c(1:10))
#Add space for the legend on the right margin
pp <- getDefaultPlotParams(4)
pp$rightmargin <- 0.09
pp$leftmargin <- 0.04
pp$topmargin <- 10
pp$bottommargin <- 18

# Create tiff file
tiff("copynumber.tiff", units="in", width=7, height=5, res=300)

# Create the plot

# Ideogram
kp <- plotKaryotype(genome = custom.genome,  plot.params = pp, plot.type = 4, labels.plotter = NULL)
kpAddChromosomeNames(kp, chr.names = chr.names , cex=0.8, yoffset = -2)
kp <- kpAddBaseNumbers(kp, tick.dist = 2e6, cex=0.6, units="Mb", tick.len = 2)

# Tracks
sample <- read.table("samplename_no3128.txt", header = TRUE) # remove 3128
id <-sample[,2]
file <- sample[,1]
i=0
for (f in file){
  i=i+1
  gfile<-paste0(f, "_gain.txt")
  lfile<-paste0(f,"_loss.txt")
  gdf <-read.table(gfile, header = TRUE)
  ldf<-read.table(lfile, header = TRUE)
  gain<-toGRanges(gdf[,-1])
  loss<-toGRanges(ldf[,-1])
  #print (id[i])
  at <- autotrack(current.track = i, total.tracks = length(file))
  kpAddLabels(kp, labels = id[i], r0=at$r0, r1=at$r1, cex=0.8, side="right")
  kpPlotRegions(kp, data=gain, r0=at$r0, r1=at$r1, col="darkblue", border = NA)
  kpPlotRegions(kp, data=loss, r0=at$r0, r1=at$r1, col="darkred", border = NA)
}

dev.off()
