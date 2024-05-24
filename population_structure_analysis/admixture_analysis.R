
# Install ggplot2 
if(!require(ggplot2)) {
  install.packages("ggplot2")
  library(ggplot2)
}


### Cross-Validation Errors###

# Create df with cross validation error 
all_samples_QC_MAC.cv_error <- data.frame(
  K = c(2, 3, 4, 5, 6, 7),
  Error = c(0.20099, 0.19048, 0.20664, 0.20969, 0.22221, 0.23231, 
            0.19440, 0.18434, 0.20011, 0.20302, 0.21492, 0.22479,
            0.21920, 0.20643, 0.21552, 0.23199, 0.23410, 0.24731, 
            0.21885, 0.20565, 0.21454, 0.23082, 0.23279, 0.24620),
  Samples = rep(c("all_samples_QC_MAC", "all_samples_QC_noMAC", "all_samples_QC_MAC_LD", "all_samples_QC_noMAC_LD"), each=6)
)

# Create the plot 
plot <- ggplot(all_samples_QC_MAC.cv_error, aes(x = K, y = Error, color=Samples, group=Samples)) +
  geom_line() +
  geom_point() +
  labs(x = "K", y = "Cross-validation Error", color="Sample") +
  ggtitle("CV Error all_samples_QC_MAC")
  theme_minimal()

# Display the plot
print(plot)



### Q Estimates ###

# all_samples_QC_MAC #
all_samples_QC_MAC.3.Q=read.table("all_samples_QC_MAC.3.Q")
sample_origin = c(rep("DRC",25),rep("Gambia",25),rep("Kenya",25), rep("Tanzania",14))

# Create the barplot and store bar midpoints
all_samples_QC_MAC.3.Q_barplot <- barplot(t(as.matrix(all_samples_QC_MAC.3.Q)), 
                                          col=rainbow(3), 
                                          xlab="all_samples_QC_MAC", 
                                          ylab="Ancestry", 
                                          border=NA, 
                                          names.arg=sample_origin, 
                                          las=2, 
                                          cex.names=0.5)

# add black lines to sample origin labels
#axis(1, at=all_samples_QC_MAC.3.Q_barplot, labels=sample_origin, las=2, cex.axis=0.5)

# all_samples_QC_MAC_LD #
all_samples_QC_MAC_LD.3.Q=read.table("all_samples_QC_MAC_LD.3.Q")
sample_origin = c(rep("DRC",25),rep("Gambia",25),rep("Kenya",25), rep("Tanzania",14))

# Create the barplot and store bar midpoints
all_samples_QC_MAC_LD.3.Q_barplot <- barplot(t(as.matrix(all_samples_QC_MAC_LD.3.Q)), 
                           col=rainbow(3), 
                           xlab="all_samples_QC_MAC_LD", 
                           ylab="Ancestry", 
                           border=NA, 
                           names.arg=sample_origin, 
                           las=2, 
                           cex.names=0.5)

# add black lines to sample origin labels
#axis(1, at=all_samples_QC_MAC_LD.3.Q_barplot, labels=sample_origin, las=2, cex.axis=0.5)


# all_samples_QC_noMAC #
all_samples_QC_noMAC.3.Q=read.table("all_samples_QC_noMAC.3.Q")
sample_origin = c(rep("DRC",25),rep("Gambia",25),rep("Kenya",25), rep("Tanzania",14))

# Create the barplot and store bar midpoints
all_samples_QC_noMAC.3.Q_barplot <- barplot(t(as.matrix(all_samples_QC_noMAC.3.Q)), 
                                          col=rainbow(3), 
                                          xlab="all_samples_QC_noMAC", 
                                          ylab="Ancestry", 
                                          border=NA, 
                                          names.arg=sample_origin, 
                                          las=2, 
                                          cex.names=0.5)

# add black lines to sample origin labels
#axis(1, at=all_samples_QC_noMAC.3.Q_barplot, labels=sample_origin, las=2, cex.axis=0.5)

# all_samples_QC_noMAC_LD #
all_samples_QC_noMAC_LD.3.Q=read.table("all_samples_QC_noMAC_LD.3.Q")
sample_origin = c(rep("DRC",25),rep("Gambia",25),rep("Kenya",25), rep("Tanzania",14))

# Create the barplot and store bar midpoints
all_samples_QC_noMAC_LD.3.Q_barplot <- barplot(t(as.matrix(all_samples_QC_noMAC_LD.3.Q)), 
                                             col=rainbow(3), 
                                             xlab="all_samples_QC_noMAC_LD", 
                                             ylab="Ancestry", 
                                             border=NA, 
                                             names.arg=sample_origin, 
                                             las=2, 
                                             cex.names=0.5)

# add black lines to sample origin labels
#axis(1, at=all_samples_QC_noMAC_LD.3.Q_barplot, labels=sample_origin, las=2, cex.axis=0.5)


