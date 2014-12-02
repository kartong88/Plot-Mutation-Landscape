
# Get some of the required libraries
require(stringr)


generate_plot<-function(mutlist, pct_yaxs_max = 5.5, filename=0){
  ## Perform analysis for the trinucleotide variants
  
  ### Check which of the mutation variant is missing and add that to the named vector
  Letters<-c("A", "T", "C", "G")
  conversion<-c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  trinuc.combi<-c()
  
  
  # Get all combinations of the conversions
  conv.combi<-c()
  for(i in 1:length(Letters)){
    for(j in 1:length(conversion)){
      for(k in 1:length(Letters)){
        combination<-paste(Letters[i], "[" , conversion[j], "]" , Letters[k], sep="")
        conv.combi<-c(conv.combi, combination)
      }
    }
  }
  
  
  
  conv.data<-table(toupper(mutlist$Somatic_mutation_type))
  
  # Set missing mutations as value of zero
  conv.data[conv.combi[!(conv.combi %in% names(conv.data))]] = 0
  
  
  
  names(conv.data)
  data.conv<-names(conv.data)
  
  
  conv.label<-str_extract(data.conv, ".>.")
  firstchar<-str_extract(data.conv, "^.")
  lastchar<-str_extract(data.conv, ".$")
  conv.label.all<-paste(conv.label, "_", firstchar, lastchar, sep="")
  conv.label.all.sorted<-order(conv.label.all)
  
  
  # Generate the trinucleotide label
  midrefchar<-str_extract(str_extract(data.conv, "(.)>"), "^.")
  trinuc.lab<-paste(firstchar, midrefchar, lastchar, sep="")
  trinuc.lab.sorted<-trinuc.lab[conv.label.all.sorted]
  
  
  
  # Labels information for subsequent barplot
  colour_array<-c("red","blue","green","yellow", "purple", "black")
  text_array=c("C>A","C>G","C>T","T>A","T>C","T>G")
  
  
  # Convert the data to percentage
  conv.data.norm<-conv.data/sum(conv.data)*100
  

  
  ### Barplot for Trinucleotide Mutation rate (Counts)
  #barplot(conv.data[conv.label.all.sorted], col=rep(colour_array,each=16), cex.names=0.45, las=3, names.arg=trinuc.lab.sorted, ylim=c(0,100                                                                                                                                                                   ))
  
  for(i in 1:6){
    # Size of 1.2 per bar.
    # Total size of 19.2 for 16 barplots
    left<-(i -1)*19.2 + 0.2 # to create a bit of white space
    right<-i* 19.2 - 0.2
    label_mid<-19.2/2 +(i-1)*19.2
    
    rect(left,95,right,96, col=colour_array[i], border=NA)
    text(x=label_mid, y=98, labels=text_array[i])
  }
  
  
  
  
  ### Barplot for Trinucleotide Mutation rate (Percentage)
  # Make the PDF plot of the graph
  if(!filename==0){
    pdf(filename, family="Arial")
  }
  barplot(conv.data.norm[conv.label.all.sorted], col=rep(colour_array,each=16), cex.names=0.3, las=3, names.arg=trinuc.lab.sorted,ylim=c(0,pct_yaxs_max), ylab="Percentage of Somatic Mutations (%)", xlab="Trinucleotide")
  
  for(i in 1:6){
    # Size of 1.2 per bar.
    # Total size of 19.2 for 16 barplots
    left<-(i -1)*19.2 + 0.2 # to create a bit of white space
    right<-i* 19.2 - 0.2
    label_mid<-19.2/2 +(i-1)*19.2
    
    #rect(left,5.1,right,5.2, col=colour_array[i], border=NA)
    rect(left,0.95*pct_yaxs_max,right,0.96*pct_yaxs_max, col=colour_array[i], border=NA)
    text(x=label_mid, y=0.98*pct_yaxs_max, labels=text_array[i])
  }
  
  if(!filename==0){
    dev.off()
  }


}
