library(dplyr)

# Read in repeats file with column for length of repeats in no of bases
repeats<-read.table("repeat_withlength_oct2022.txt", header=TRUE)
str(repeats)

# First create the bins column (tried adding bins to dplyr code but was giving wrong sum values)
group = cut(repeats$end, breaks = seq(1, max(repeats$end), 100000))
repeats <-cbind(repeats, group)

# Group by chromosome and bins and sum over lengths variable
newdf<-repeats %>%
  group_by(chromosome, group) %>%
  mutate(sumlength = sum(length)) %>%
  distinct(chromosome, group, sumlength)

# Save output
write.table(newdf, "repeats_withlength_forcircos_oct2022.txt", sep="\t")

