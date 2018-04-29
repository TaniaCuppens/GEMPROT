#!/bin/R
 #Tania Cuppens
 #INSERM U1078, Ph.D.
 #Protein Plotting Script
 #Programming Language: R
 #Updated 09/15/2017

 #NOTE: All files should be referring to the same isoform of the protein. This is imperative for drawing the plot correctly.


 #Usage:
 ##Basic
 ###Rscript plotProtein.R -m psen1_mutation_file.txt -a psen1_architecture_file.txt -l 463


 ##Advanced
 ###Rscript plotProtein.R -m psen1_mutation_file.txt -a psen1_architecture_file.txt -l 463 -n Disease -t 25 -s yes -z yes -b 50 -c 100


 #will have to run the install for optparse if never installed before
 #install.packages("optparse")


library("optparse")


 option_list <- list(
     make_option(c('-m', '--mutations'), action='store', type='character', default='mutationFile.txt', help='This is the mutation file. It should be a tab-delimited file containing 5 columns (ProteinId, GeneName, ProteinPositionOfMutation, ReferenceAminoAcid, AlternateAminoAcid) NO HEADER FOR NEEDED FOR THIS FILE. (REQUIRED)'),
     make_option(c('-d', '--domain'), action='store', type='character', default='domainFile.txt', help='This is the protein domain file. It should be a tab-delimited file containing 3 columns (domain_name, start_site, end_site). This file NEEDS the header and it is the same as what was previously written. (REQUIRED)'),
     make_option(c('-p', '--pdf'), action='store', type='character', default='gene', help='This is the pdf id. (REQUIRED)'),
     make_option(c('-l', '--length'), action='store', type='numeric', default=100, help='protein length (REQUIRED)'),
     make_option(c('-n', '--name'), action='store', type='character', default='Test', help='Name of your query. Default is Test'),
     make_option(c('-t', '--ticksize'), action='store', type='numeric', default=10, help='Size of ticks on x-axis. Default is 10'),
     make_option(c('-s', '--showlabels'), action='store', type='character', default='no', help='Option to show labels. Default is no'),
     make_option(c('-z', '--zoom'), action='store', type='character', default='no', help='Option to zoom in. Default is no'),
	 make_option(c('-a', '--mode'), action='store', type='character', default='indiv', help='Option to define mode. Default is indiv')
 )
 opt <- parse_args(OptionParser(option_list = option_list))


 set = opt$set #global quality threshold (can vary)
 mutationFile <- opt$mutations
 proteinDomainFile <- opt$domain
 proteinLength <- opt$length
 nameOfYourQuery <- opt$name
 tickSize <- opt$ticksize
 showLabels <- opt$showlabels
 pdf <- opt$pdf
 zoomIn <- opt$zoom
 mode <- opt$mode
 if(zoomIn == "yes"){
 	     zoomStart <- opt$zoomstart
 	     zoomEnd <- opt$zoomend
 }


 ####################ANALYSIS####################
 #Read in the files
 var <- read.table(mutationFile, sep="\t", colClasses = "character")
 pd <- read.table(proteinDomainFile, sep="\t", header=TRUE)
 
 ############PLOTTING#############
 #x is the input data, y is rpt, z is rpa from HPRD
pdf(paste(pdf, "_protein_plot.pdf", sep=""), height=7.5, width=12)

 par(mar=c(0, 0, 0, 0))
 nameOfYourQuery <- as.character(var[1,2])
 



 xlimRegion <- c(0,(as.numeric(proteinLength)))
 	      if(zoomIn == "yes") {
           xlimRegion <- c(as.numeric(zoomStart), as.numeric(zoomEnd))
 	     }

plot((1:as.numeric(proteinLength)), rep(2.5, as.numeric(proteinLength)), type="l", lwd=10, main=paste("Amino Acid Changes in", " ", as.character(var[1,2]), " ", "(", as.character(var[1,1]), ")", sep=""), xlab="Amino Acid Position", ylab="", ylim=c(0,5), cex.lab=0.9, cex.main=1, yaxt="n", xlim=xlimRegion, xaxt="n", col="black",ann=FALSE, bty="n")
#limite mutation size
limt <- ceiling((0.162*proteinLength)/12)
#limite domain size
limd <- ((0.101*proteinLength)/12)

if(showLabels == "yes"){
			th1 <- 0
			th2 <- 0
			tm1 <- 0
			tm2 <- 0
			db2 <-(2.8)
			pt2 <-(2.75)
			db3 <-(2.2)
			pt3 <-(2.25)
	if (nrow(var)==1) {
	
	}	
	else {
 		 #Label mutations
 		 for(i in 2:nrow(var)){
			cmut <- as.numeric(var[i,3])
			csample <- var[i,1]
			#Plot mutations
			if (as.character(var[i,8])=="red"){
				segments(as.numeric(var[i,3]), 4, as.numeric(var[i,3]), 1, col="wheat3", lwd=0.8, lty=2)
				#abline(v=as.numeric(var[i,3]), col="darkred", lwd=1, lty=2)
			}
			else if (grepl("hap1",csample,)){
				db <-(2.8)
				pt <-(2.75)
				if (th1!=0) {
					
					if ((as.numeric(var[i,3])-as.numeric(var[(i-1),3]))<limt && (as.numeric(var[i,3])-as.numeric(var[(i-tm1),3]))<limt && tm1 <= 4 ) {
						if (as.character(var[i,6])=="frameshift") {
							nmut <- paste(as.character(var[(i-1),4]),as.character(var[(i-1),3]),as.character(var[(i-1),5]), sep="")
							db2 <-(db2+(0.05*(nchar(nmut))))
							pt2 <- (pt2+(0.05*(nchar(nmut))))
							text(as.character(var[i,3]), rep(db2, length(var[i,3])),paste("FS:",as.character(var[i,3]),"-",as.character(var[i,5]), sep=""), font=2, family="mono",col=as.character(var[i,8]), cex=1, srt=90, adj = 0)
							for (f in as.numeric(var[i,3]):as.numeric(var[i,5])){
								points(f, rep(pt2, length(f)), pch=as.numeric(var[i,7]), col=as.character(var[i,8]), cex=0.8)
							}
						th1 <- (th1+1)
						tm1 <- (tm1+1)	
						}
						else {
							nmut <- paste(as.character(var[(i-1),4]),as.character(var[(i-1),3]),as.character(var[(i-1),5]), sep="")
							db2 <-(db2+(0.09*(nchar(nmut))))
							pt2 <- (pt2+(0.09*(nchar(nmut))))
							points(var[i,3], rep(pt2, length(var[i,3])), pch=as.numeric(var[i,7]), col=as.character(var[i,8]), cex=0.8)
							text(as.character(var[i,3]), rep(db2, length(var[i,3])),paste(as.character(var[i,4]),as.character(var[i,3]),as.character(var[i,5]), sep=""), font=2, family="mono",col=as.character(var[i,8]), cex=1, srt=90, adj = 0)
							th1 <- (th1+1)
							tm1 <- (tm1+1)
						}
					}
					else if (tm1>=4){
						points(var[i,3], rep(pt2, length(var[i,3])), pch=as.numeric(var[i,7]), col=as.character(var[i,8]), cex=0.8)
						text(as.character(var[i,3]), rep(db2, length(var[i,3])),paste(".",".",".", sep=""), font=2, family="mono",col="red", cex=1, srt=90, adj = 0)
						tm1 <- 1
						th1 <- (th1+1)
						db2 <-(2.8)
						pt2 <-(2.75)
					}
					
					else if (as.character(var[i,6])=="frameshift") {
						text(as.character(var[i,3]), rep(db, length(var[i,3])),paste("FS:",as.character(var[i,3]),"-",as.character(var[i,5]), sep=""), font=2,family="mono",col=as.character(var[i,8]), cex=1, srt=90, adj = 0)	
						for (f in as.numeric(var[i,3]):as.numeric(var[i,5])){
							points(f, rep(pt, length(f)), pch=as.numeric(var[i,7]), col=as.character(var[i,8]), cex=0.8)
						}
						th1 <- (th1+1)
						db2 <-(2.8)
						pt2 <-(2.75)
					}
					
					else {
						points(var[i,3], rep(pt, length(var[i,3])), pch=as.numeric(var[i,7]), col=as.character(var[i,8]), cex=0.8)
						text(as.character(var[i,3]), rep(db, length(var[i,3])),paste(as.character(var[i,4]),as.character(var[i,3]),as.character(var[i,5]), sep=""), font=2, family="mono",col=as.character(var[i,8]), cex=1, srt=90, adj = 0)
						th1 <- (th1+1)
						tm1 <- 1
						db2 <-(2.8)
						pt2 <-(2.75)
					}
				}
				else {
					if (as.character(var[i,6])=="frameshift") {
						text(as.character(var[i,3]), rep(db, length(var[i,3])),paste("FS:",as.character(var[i,3]),"-",as.character(var[i,5]), sep=""), font=2, family="mono",col=as.character(var[i,8]), cex=1, srt=90, adj = 0)
						for (f in as.numeric(var[i,3]):as.numeric(var[i,5])){
							points(f, rep(pt2, length(f)), pch=as.numeric(var[i,7]), col=as.character(var[i,8]), cex=0.8)
						}
						th1 <- (th1+1)
						tm1 <- (tm1+1)
						db2 <-(2.8)
						pt2 <-(2.75)
					}
					else {
						points(var[i,3], rep(pt2, length(var[i,3])), pch=as.numeric(var[i,7]), col=as.character(var[i,8]), cex=0.8)
						text(as.character(var[i,3]), rep(db2, length(var[i,3])),paste(as.character(var[i,4]),as.character(var[i,3]),as.character(var[i,5]), sep=""), font=2, family="mono",col=as.character(var[i,8]), cex=1, srt=90, adj = 0)
						th1 <- (th1+1)
						db2 <-(2.8)
						pt2 <-(2.75)
					}
				}
		
			 }
			else {
				pt1 <- (2.25)
				db1 <- (2.2)
				#nmut2 <- paste(as.character(var[i,4]),as.character(var[i,3]),as.character(var[i,5]), sep="")
			
				#ll <- (0.08*(nchar(nmut2)))
				#db1 <- (pt1-ll)
				
				if (th2!=0) {
					if ((as.numeric(var[i,3])-as.numeric(var[(i-1),3]))<limt && (as.numeric(var[i,3])-as.numeric(var[(i-tm2),3]))<limt && tm2<=4)  {
						if (as.character(var[i,6])=="frameshift") {
							nmut3 <- paste(as.character(var[(i-1),4]),as.character(var[(i-1),3]),as.character(var[(i-1),5]), sep="")
							db3 <-(db3-(0.09*(nchar(nmut3))))
							pt3<- (pt3-(0.09*(nchar(nmut3))))
							text(as.numeric(var[i,3]), rep(db3, length(var[i,3])),paste("FS:",as.character(var[i,3]),"-",as.character(var[i,5]), sep=""), font=2,family="mono",col=as.character(var[i,8]), pos=2, offset=0, cex=1, srt=90, adj = 0)	
							for (f in as.numeric(var[i,3]):as.numeric(var[i,5])){
								points(f, rep(pt3, length(f)), pch=as.numeric(var[i,7]), col=as.character(var[i,8]), cex=0.8)
							}
							th2 <- (th2+1)
							tm2 <- (tm2+1)
						}	
						else {
							nmut3 <- paste(as.character(var[(i-1),4]),as.character(var[(i-1),3]),as.character(var[(i-1),5]), sep="")
							db3 <-(db3-(0.09*(nchar(nmut3))))
							pt3<- (pt3-(0.09*(nchar(nmut3))))
							points(var[i,3], rep(pt3, length(var[i,3])), pch=as.numeric(var[i,7]), col=as.character(var[i,8]), cex=0.8)
							text(as.numeric(var[i,3]), rep(db3, length(var[i,3])),paste(as.character(var[i,4]),as.character(var[i,3]),as.character(var[i,5]), sep=""), font=2,family="mono",col=as.character(var[i,8]), cex=1, pos=2, offset=0, srt=90, adj = 0)
							th2 <- (th2+1)
							tm2 <- (tm2+1)
							}
					}
					else if (tm2>=4){
						 points(var[i,3], rep(pt3, length(var[i,3])), pch=as.numeric(var[i,7]), col=as.character(var[i,8]), cex=0.8)
						 text(as.numeric(var[i,3]), rep(db3, length(var[i,3])),paste(".",".",".", sep=""), font=2, family="mono",col="red", cex=1, pos=2, offset=0, srt=90, adj = 0)
						 tm2 <- 1
						 th2 <- (th2+1)
						 pt3 <- (2.25)
						 db3 <- (2.2)
					 }
					else if (as.character(var[i,6])=="frameshift") {
						text(as.numeric(var[i,3]), rep(db1, length(var[i,3])),paste("FS:",as.character(var[i,3]),"-",as.character(var[i,5]), sep=""), font=2,family="mono",col=as.character(var[i,8]), pos=2, offset=0, cex=1, srt=90, adj = 0)	
						for (f in as.numeric(var[i,3]):as.numeric(var[i,5])){
							points(f, rep(pt1, length(f)), pch=as.numeric(var[i,7]), col=as.character(var[i,8]), cex=0.8)
						}
						th2 <- (th2+1)
						pt3 <- (2.25)
						db3 <- (2.2)
							
					}		
					else {	
						points(var[i,3], rep(pt1, length(var[i,3])), pch=as.numeric(var[i,7]), col=as.character(var[i,8]), cex=0.8)
						text(as.numeric(var[i,3]), rep(db1, length(var[i,3])),paste(as.character(var[i,4]),as.character(var[i,3]),as.character(var[i,5]), sep=""), font=2, family="mono",col=as.character(var[i,8]), pos=2 , offset=0, cex=1, srt=90, adj = 0)
						th2 <- (th2+1)
						tm2 <- 1
						pt3 <- (2.25)
						db3 <- (2.2)
					}
				}
				else {	
					if (as.character(var[i,6])=="frameshift") {
						text(as.numeric(var[i,3]), rep(db1, length(var[i,3])),paste("FS:",as.character(var[i,3]),"-",as.character(var[i,5]), sep=""), font=2,family="mono",col=as.character(var[i,8]), pos=2, offset=0, cex=1, srt=90, adj = 0)	
						for (f in as.numeric(var[i,3]):as.numeric(var[i,5])){
							points(f, rep(pt1, length(f)), pch=as.numeric(var[i,7]), col=as.character(var[i,8]), cex=0.8)
						}
						th2 <- (th2+1)
						pt3 <- (2.25)
						db3 <- (2.2)
					}
					else {
						points(var[i,3], rep(pt1, length(var[i,3])), pch=as.numeric(var[i,7]), col=as.character(var[i,8]), cex=0.8)
						text(as.numeric(var[i,3]), rep(db1, length(var[i,3])),paste(as.character(var[i,4]),as.character(var[i,3]),as.character(var[i,5]), sep=""), font=2, family="mono",col=as.character(var[i,8]), pos=2, offset=0, cex=1, srt=90, adj = 0)
						th2 <- (th2+1)
						pt3 <- (2.25)
						db3 <- (2.2)
					}
				}
				
		       	
			}
			
		}
	}
}

if (th1==0) {
	text(proteinLength/2, 2.8, labels="no changes", font=2, family="mono",col="blue2", cex=0.9, srt=0, adj = 0)
}

if (th2==0 && mode=="indiv") {
	text(proteinLength/2, 2.2, labels="no changes", font=2, family="mono",col="orangered2", cex=0.9, srt=0, adj = 0)
}

#ticks=seq(0,as.numeric(proteinLength), by=tickSize)
#axis(side = 1, at = ticks, las=3)

#labels
if (mode=="indiv"){
	text(-(limt*2.5), 2.6, labels="hap1", font=1, family="mono",col="blue2", cex=0.9, srt=0, adj = 0)
	text(-(limt*2.5), 2.4, labels="hap2", font=1, family="mono",col="orangered2", cex=0.9, srt=0, adj = 0)
}
ndomain_out <- 0
m <- matrix(, nrow = 0, ncol = 2)
m = rbind(m, c("domain","name")) 
l <-list()
#print (l)
#print (as.character(pd[1,1]))

if ( !(as.character(pd[1,2])=="No_domains")){
	for(i in 1:length(pd$start_site)){
      rect(as.numeric(pd$start_site[i]), 2.6, as.numeric(pd$end_site[i]), 2.4, col="palegreen3")
	}
	for(i in 1:nrow(pd)){
		domain_size <- (as.numeric(pd[i,4])-as.numeric(pd[i,3]))	
		#print("----------------")
		#print (i)
		#print ("domain_size")
		#print (domain_size)
		#print ("domain_name")
		#print (pd[i,1])
		#print ("size")
		#print (nchar(as.character(pd[i,1])))
		#print (limd*nchar(as.character(pd[i,1])))
		#print (l)
		if (domain_size<=ceiling(limd*nchar(as.character(pd[i,2])))){
			if( !(as.character(pd[i,2]) %in% m[,1]) ) { 
				#print (m[,1])
				ndomain_out <- (ndomain_out+1)
				#print (as.character(pd[i,2]))
				m = rbind(m, c(as.character(pd[i,2]),LETTERS[ndomain_out]) )
				text(median(c(as.numeric(pd[i,3]), as.numeric(pd[i,4]))), 2.5, LETTERS[ndomain_out],font=2, family="mono", cex=1)
				l[ndomain_out] <- paste(LETTERS[ndomain_out]," = ",as.character(pd[i,2]),sep="")
				
			}
			else {
				text(median(c(as.numeric(pd[i,3]), as.numeric(pd[i,4]))), 2.5, m[m[,1] == as.character(pd[i,2]), 2],font=2, family="mono", cex=0.9)
			}
			#l[ndomain_out] <- paste(LETTERS[ndomain_out]," = ",as.character(pd[i,1]),sep="")
			#text(median(c(as.numeric(pd[i,2]), as.numeric(pd[i,3]))), 2.5, LETTERS[ndomain_out],font=2, family="mono", cex=1) 
			
		}
		else {
			text(median(c(as.numeric(pd[i,3]), as.numeric(pd[i,4]))), 2.5, pd[i,2],font=2, family="mono", cex=0.9)
		}
	}
}

#print (m)
#print (l)
if (ndomain_out!=0){
	legend ("topright", legend = l , box.col="white", bg="white", cex=0.8)
}

legend ("topleft", legend = (as.character(pd[1,1])) , box.col="white", bg="white", cex=0.8)
legend("bottomright", c("mismatch", "synonymous", "stop", "elongation", "deletion", "insertion"), 
      pch = c(19, 1, 8, 13, 17, 6), box.col="white", bg="white", cex=0.8)
#legend("bottomright", legend =c("Protein Domain", nameOfYourQuery,"Haplotype1","Haplotype2"), fill = c("palegreen3", "black","blue2","indianred3"), box.col="white", bg="white", cex=0.7)
#dev.off()

invisible(dev.off())


