library(karyoploteR)

# This script will display copy number calls from CNVpytor output as linear plots across chromosomes 

# Read in the genome file
custom.genome <- toGRanges("mygenome.txt")

# Add space for the legend on the right margin
pp <- getDefaultPlotParams(4)
pp$rightmargin <- 0.15

# Create the ideogram
kp <- plotKaryotype(genome = custom.genome,  plot.params = pp, plot.type = 4, labels.plotter = NULL)
kpAddChromosomeNames(kp, cex=0.5, srt=45, yoffset = -2)
kp <- kpAddBaseNumbers(kp, tick.dist = 2e6, cex=0.5, units="Mb", tick.len = 3)

# Read in *cn files from CNVpytor
myFiles <- list.files(pattern="*_cn.txt")
sampleid <- read.table("samplename.txt", header = TRUE)
id <-sampleid[,1]
for (f in 1:length(myFiles)){
  df <-read.table(myFiles[f], header = TRUE)
  CN<-toGRanges(df[,-1])
  print (myFiles[f])
  print (id[f])
  at <- autotrack(current.track = f, total.tracks = length(myFiles))
  kpDataBackground(kp, r0=at$r0, r1=at$r1)
  kpAddLabels(kp, labels = id[f], r0=at$r0, r1=at$r1, cex=0.5, side="right")
  kpHeatmap(kp, data=CN, y=CN$CNV_level, ymin=min(CN$CNV_level), ymax = max(CN$CNV_level),  colors = c("red", "white", "blue"), r0=at$r0, r1=at$r1, cex=0.3)
}


#Add the legend
#Simulate a continuous legend with a tightly packed discrete legend with fake data
#temporarily remove all margins, so we can control with precision where to plot the legend
pars.old <- graphics::par(no.readonly = TRUE) #Save the old arams
par(mar=c(5.1, 4.1, 4.1, 5.1), xpd=TRUE)

lgd_ = rep(NA, 101) #The legend will be created with 101 color boxes
lgd_[c(101, 50, 1)] = c(0, 2, 4) #Add the legend labels to the top, center and bottom elements
legend(x = "topright",  #legend positioning (stadard base R legend, accepts named positions or numeric positions)
       legend = lgd_,
       y.intersp = 0.07,
       fill = colorRampPalette(colors = c("red", "white", "blue"))(101),
       border = NA, box.lwd = NA, box.col = transparent("white", 1),
       cex=0.4) #set the length of the legend

graphics::par(pars.old) #reset the original pars
