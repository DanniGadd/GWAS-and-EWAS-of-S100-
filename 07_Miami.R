
## Miami Plot ##
# Load in EWAS and GWAS summstats and assign them to "ewas" and "gwas" respectively
# Make sure your EWAS and GWAS summstats have the following columns:
# CHR, BP, P, SNP, Data
# CHR = Chromosome
# BP = Genomic coordinates
# SNP = rsID or cgID depending on GWAS or EWAS
# Data = Data type. "GWAS" in your GWAS summstats and "EWAS" in your EWAS summstats



miami <- function(x, chr="CHR", bp="BP", p="P", snp="SNP", data="Data", data.upper=NA, data.lower=NA,
                      col=c("firebrick1","firebrick4"), pfilt=1, chrlabs=NULL,
                      lowersuggestiveline=-log10(1e-5), lowergenomewideline=-log10(3.6e-8), 
					  uppersuggestiveline=-log10(1e-5), uppergenomewideline=-log10(5e-8), 
                      highlight=NULL, logp=TRUE, annotatePval = NULL, annotateTop = TRUE, ...) {

    # Not sure why, but package check will warn without this.
    CHR=BP=P=index=NULL
    
    # Check for sensible dataset
    ## Make sure you have chr, bp and p columns.
    if (!(chr %in% names(x))) stop(paste("Column", chr, "not found!"))
    if (!(bp %in% names(x))) stop(paste("Column", bp, "not found!"))
    if (!(p %in% names(x))) stop(paste("Column", p, "not found!"))
    ## warn if you don't have a snp column
    if (!(snp %in% names(x))) warning(paste("No SNP column found. OK unless you're trying to highlight."))
    ## make sure chr, bp, and p columns are numeric.
    if (!is.numeric(x[[chr]])) stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
    if (!is.numeric(x[[bp]])) stop(paste(bp, "column should be numeric."))
    if (!is.numeric(x[[p]])) stop(paste(p, "column should be numeric."))
    
    # Create a new data.frame with columns called CHR, BP, and P.
    # d=data.frame(CHR=x[[chr]], BP=x[[bp]], P=x[[p]], pos = NA, index = NA) # with millions of SNPs, create dataframe at once 
	                                                         #  rather than dynamically allocated(see line 72-73, and remove line 87 and line 91 )
    
    # If the input data frame has a SNP column, add it to the new data frame you're creating.
    if (!is.null(x[[snp]])) d = data.frame(CHR=x[[chr]], BP=x[[bp]], P=x[[p]], Data=x[[data]], pos = NA, index = NA ,SNP=x[[snp]], stringsAsFactors = FALSE) else 
	    d = data.frame(CHR=x[[chr]], BP=x[[bp]], P=x[[p]], Data=x[[data]], pos = NA, index = NA)
	    
	    
	    
	    
    
    # Set positions, ticks, and labels for plotting
    ## Sort and keep only values where is numeric.
    #d <- subset(d[order(d$CHR, d$BP), ], (P>0 & P<=1 & is.numeric(P)))
    #  d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))       ## unused, all three variables are numeric, line:63-65 
    d <- d[order(d$CHR, d$BP), ]
    d = d[which(d$P < pfilt), ]
    #d$logp <- ifelse(logp, yes=-log10(d$P), no=d$P)
    if (logp) {
        d$logp <- -log10(d$P)
        d[which(d$Data==data.lower), "logp"] <- d[which(d$Data==data.lower), "logp"] * -1
    } else {
        d$logp <- d$P
    }
   # d$pos=NA
    
    
    # Fixes the bug where one chromosome is missing by adding a sequential index column.
   # d$index=NA
   # ind = 0
   # for (i in unique(d$CHR)){
   #     ind = ind + 1
   #     d[d$CHR==i,]$index = ind
   # }
   d$index = rep.int(seq_along(unique(d$CHR)), times = tapply(d$SNP,d$CHR,length))  # replcace the for loop of line 92-96 to improve efficiency
    
    # This section sets up positions and 	. Ticks should be placed in the
    # middle of a chromosome. The a new pos column is added that keeps a running
    # sum of the positions of each successive chromsome. For example:
    # chr bp pos
    # 1   1  1
    # 1   2  2
    # 2   1  3
    # 2   2  4
    # 3   1  5
    nchr = length(unique(d$CHR))
    if (nchr==1) { ## For a single chromosome
        ## Uncomment the next two linex to plot single chr results in Mb
        #options(scipen=999)
	    #d$pos=d$BP/1e6
        d$pos=d$BP
      #  ticks=floor(length(d$pos))/2+1          ## unused, from code line: 169
        xlabel = paste('Chromosome',unique(d$CHR),'position')
      #  labs = ticks          ## unused, from code line: 169
    } else { ## For multiple chromosomes
        lastbase=0
        ticks=NULL
        for (i in unique(d$index)) {
            if (i==1) {
                d[d$index==i, ]$pos=d[d$index==i, ]$BP
            } else {
		## chromosome position maybe not start at 1, eg. 9999. So gaps may be produced. 
		lastbase = lastbase +max(d[d$index==(i-1),"BP"])   # replace line 128
		d[d$index == i,"BP"] = d[d$index == i,"BP"]-min(d[d$index==i,"BP"]) +1
		d[d$index == i, "pos"] = d[d$index == i,"BP"] + lastbase    # replace line 129
                # lastbase=lastbase+tail(subset(d,index==i-1)$BP, 1)
                # d[d$index==i, ]$pos=d[d$index==i, ]$BP+lastbase
		   
            }
            # Old way: assumes SNPs evenly distributed
            # ticks=c(ticks, d[d$index==i, ]$pos[floor(length(d[d$index==i, ]$pos)/2)+1])
            # New way: doesn't make that assumption
           # ticks = c(ticks, (min(d[d$index == i,]$pos) + max(d[d$index == i,]$pos))/2 + 1)  # see line 136, to reduce the burden of for loop 
        }
	ticks <-tapply(d$pos,d$index,quantile,probs=0.5)   # replace line 135
        xlabel = 'Chromosome'
        #labs = append(unique(d$CHR),'') ## I forgot what this was here for... if seems to work, remove.
        labs <- unique(d$CHR)
    }
    
    # Initialize plot
    xmax = max(d$pos) + 100000
    xmin = floor(max(d$pos) * -0.03)
    
    # The old way to initialize the plot
    # plot(NULL, xaxt='n', bty='n', xaxs='i', yaxs='i', xlim=c(xmin,xmax), ylim=c(ymin,ymax),
    #      xlab=xlabel, ylab=expression(-log[10](italic(p))), las=1, pch=19, cex=0.5, ...)

    
    # The new way to initialize the plot.
    ## See http://stackoverflow.com/q/23922130/654296
    ## First, define your default arguments
    def_args <- list(xaxt='n', bty='n', xaxs='i', yaxs='i', las=1, pch=19, cex=0.5,
                     xlim=c(xmin,xmax), ylim=c(floor(min(d$logp)),ceiling(max(d$logp))),
                     xlab=xlabel, ylab=expression(-log[10](italic(p))))
    ## Next, get a list of ... arguments
    #dotargs <- as.list(match.call())[-1L]
    dotargs <- list(...)
    ## And call the plot function passing NA, your ... arguments, and the default
    ## arguments that were not defined in the ... arguments.
    do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
    
    # If manually specifying chromosome labels, ensure a character vector and number of labels matches number chrs.
    if (!is.null(chrlabs)) {
        if (is.character(chrlabs)) {
            if (length(chrlabs)==length(labs)) {
                labs <- chrlabs
            } else {
                warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
            }
        } else {
            warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
        }
    }
    
    # Add an axis. 
    if (nchr==1) { #If single chromosome, ticks and labels automatic.
        axis(1, ...)
    } else { # if multiple chrs, use the ticks and labels you created above.
        axis(1, at=ticks, labels=labs, ...)
    }
    
    # Create a vector of alternatiting colors
    #col=rep(col, max(d$CHR))  # replaced by line 187
    col = rep_len(col, max(d$index))  ## mean this one?  the results are same

    # Add points to the plot
    if (nchr==1) {
        with(d, points(pos, logp, pch=19, cex=0.5, col=col[1], ...))
    } else {
        # if multiple chromosomes, need to alternate colors and increase the color index (icol) each chr.
        icol=1
        for (i in unique(d$index)) {
            #with(d[d$index==unique(d$index)[i], ], points(pos, logp, col=col[icol], pch=19, cex=0.5, ...))
	    points(d[d$index==i,"pos"], d[d$index==i,"logp"], col=col[icol], pch=19, cex=0.5, ...)
            icol=icol+1
        }
    }
    
    # Add suggestive and genomewide lines
    if (uppersuggestiveline){abline(h=uppersuggestiveline, col="blue")}
    if (uppergenomewideline){abline(h=uppergenomewideline, col="red")}
    if (lowersuggestiveline){abline(h=lowersuggestiveline*-1, col="blue")}
    if (lowergenomewideline){abline(h=lowergenomewideline*-1, col="red")}

    # Highlight snps from a character vector
    if (!is.null(highlight)) {
        if (any(!(highlight %in% d$SNP))) warning("You're trying to highlight SNPs that don't exist in your results.")
        d.highlight=d[which(d$SNP %in% highlight), ]
        with(d.highlight, points(pos, logp, col="green3", pch=19, cex=0.5, ...)) 
    }
    
    # Highlight top SNPs
    if (!is.null(annotatePval)) {
        # extract top SNPs at given p-val
        topHits = subset(d, P <= annotatePval)
        par(xpd = TRUE)
        # annotate these SNPs
        if (annotateTop == FALSE) {
            with(subset(d, P <= annotatePval), 
                 textxy(pos, -log10(P), offset = 0.625, labs = topHits$SNP, cex = 0.45), ...)
        }
        else {
            # could try alternative, annotate top SNP of each sig chr
            topHits <- topHits[order(topHits$P),]
            topSNPs <- NULL
            
            for (i in unique(topHits$CHR)) {
                
                chrSNPs <- topHits[topHits$CHR == i,]
                topSNPs <- rbind(topSNPs, chrSNPs[1,])
                
            }
            textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.625, labs = topSNPs$SNP, cex = 0.5, ...)
        }
    }  
    par(xpd = FALSE)
}


### LOAD AND PLOT MIAMI DATA 

EWAS <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/s100b_project/OSCA_EWAS_060921/OSCA/result_s100b.csv")
GWAS <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/s100b_project/PLINK_GWAS_060921/s100b_GWAS_output.txt", sep = ",")

ewas <- EWAS[c(1,3,8,2)]
gwas <- GWAS[c(1,2,13,14)]

ewas$Data <- "EWAS"
gwas$Data <- "GWAS"

names(ewas) <- c("CHR", "BP", "P", "SNP", "Data")
names(gwas) <- c("CHR", "BP", "P", "SNP", "Data")

miami_dat = rbind(gwas[,c("CHR", "BP", "SNP", "Data", "P")], ewas[,c("CHR", "BP", "SNP", "Data", "P")])

# Miami plot
png("/Cluster_Filespace/Marioni_Group/Danni/s100b_project/PLOTS/miami_correct_211021.png", width = 1900, height = 1000)
par(mar=c(5,6,4,1)+.1)
miami(miami_dat, main="", data.upper="GWAS", data.lower="EWAS", col=c("darkslategray", "darkseagreen2"),
	cex.main = 2, cex.lab = 2, cex.axis = 1.5, ylim = c(-10,15), lowersuggestiveline=-log10(1e-5), lowergenomewideline=-log10(3.6e-8), 
                      uppersuggestiveline=-log10(1e-5), uppergenomewideline=-log10(5e-8))
abline(h=0, col="black", lwd = 3)
dev.off()



# png("/Cluster_Filespace/Marioni_Group/Danni/s100b_project/PLOTS/miamiplot_V5.png", width = 1200, height = 800)
# par(mfrow=c(2,1))
# par(mar=c(0,5,3,3))
# manhattan(gwas,ylim=c(0,10), col = c("darkslategray", "darkseagreen2"), cex=2.2,cex.lab=2.5,font.lab=2,font.axis=2,cex.axis=1.6,las=2,font=4,xlab="",xaxt="n",
#     suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08))
# par(mar=c(5,5,3,3))
# manhattan(ewas,ylim=c(10,0), col = c("darkslategray", "darkseagreen2"), cex=2.2,cex.lab=2.5,font.lab=2,font.axis=2,cex.axis=1.6,las=2,font=4,xlab="",xaxt="n",
#   suggestiveline = -log10(1e-05), genomewideline = -log10(3.6e-08))
# dev.off()

