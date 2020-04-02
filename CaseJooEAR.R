library(data.table)
library(dplyr)
library(stringi)
library(rio)
library(easyPubMed)
library(XML)
library(tidyr)
library(ggplot2)
library(xml2)
library(quanteda)
library(reshape2)


#to include arguments, use Rscript --verbose CaseJooEAR.R search-terms.txt tsv-files.txt >& CaseJooEAR.Rout

memory.limit(size = 56000)
args <- commandArgs(trailingOnly = TRUE)

query.file <- args[1]
tsv.files <- args[2] #allow multiple tsv files

#get the query terms
#query.dt <- as.data.frame(fread(query.file,header=FALSE,sep=","))
query_terms <- unique(scan(query.file, what=character(), sep='\n'))
query_terms

input_files <- unique(scan(tsv.files, what=character(), sep='\n',comment.char = "#"))
input_files

# will want to allow for multiple tsv files to upload
#df <- import("/net/jarvik/vol1/home/erosen/UDN/UDN436557/005/UDN436557.AR.stdGATKstdfxn.tsv",format = 'tsv')
uniqGenes <- c()
aliases <- c()
gene.alias.mat <- c("GENE","ALIAS")
for(i in 1:length(input_files)){
#i<-1
  df <- import(input_files[i],format = 'tsv')
  tm1 <- df
  keep.cols <- cbind(tm1$gene,tm1$aliases)
  gene.alias.mat <- rbind(gene.alias.mat,keep.cols) 
  gene.alias.mat <- gene.alias.mat[-1,]
  gene.alias.mat <- as.data.table(gene.alias.mat)
  setnames(gene.alias.mat,c("GENE","ALIAS"))
  setkey(gene.alias.mat,"GENE")
  gene.alias.mat <- unique(gene.alias.mat,by="GENE")

#can I create a lookup table?
  gene.alias.mat[,ALIAS:=ifelse(is.na(ALIAS),GENE,ALIAS)]

#start creating the lookup table so that each alias is associated with its current gene symbol
  gene.alias.mat$ALIAS
  my.list <- strsplit(gene.alias.mat$ALIAS,',')
  lookup.table <- (do.call(rbind.data.frame,my.list))
  lookup.table <- cbind(gene.alias.mat$GENE, lookup.table)
# convert the wide table to a long table
  colnames(lookup.table) <- c()
  lookup.table <- sapply(lookup.table,as.character)
  num_cols <- length(lookup.table[1,])
  gene.symbols <- c(rep((lookup.table[,1]),num_cols))
  gene.alias <- as.data.table(melt(as.character(lookup.table),id.vars=1,value.name="ALIAS"))
  gene.alias[,ALIAS:=as.character(ALIAS)] #necessary so that they are not though of as factors

  lookup.table <- as.data.table(cbind(gene.symbols,gene.alias))
  setkey(lookup.table,"ALIAS")
  lookup.table <- unique(lookup.table, by="ALIAS")
  setkey(lookup.table,"gene.symbols")

  uniqGenes <- unique(c(lookup.table$gene.symbols,lookup.table$ALIAS))
  gcnt <- 0
  gcnt2 <- 0
  num_terms <- length(query_terms)
  for(i in 1:length(query_terms)){
    query_terms[i] <- paste0(query_terms[i],'[Title/Abstract])')
  }
  base_query <- paste(query_terms,collapse=' OR ')
#add in the opening brackets
  base_query <- paste0(paste(rep('(',length(query_terms)),collapse=""),base_query, collpase="")

  pubmedDf <- data.frame(matrix(ncol = 5, nrow = 0))
  for (g in uniqGenes){
    print(gcnt2)
    gcnt2 <- gcnt2 + 1

    my_query <- paste(base_query,' AND ',g,'[Title/Abstract]',sep = '')
    out.A <- batch_pubmed_download(pubmed_query_string = my_query, 
                                   format = "xml", 
                                   batch_size = 5000,
                                   dest_file_prefix = "project4", encoding = "ASCII") ##EAR will want to make this something the user can change easily

    
    
    if (!is.null(out.A[1])){
        a <- read_xml(out.A)
        b <- xml_find_all(a, "//PubmedBookArticle")
        for (i in seq_along(b)){
            # extracting the abstract, pub_year, pmid, title from xml data
            abstract <- ""
            abstract <- xml_text(xml_find_all(read_xml(as.character(b[[i]])),"//Abstract"))
if(length(abstract)==0){
  abstract <- NA
}
            pub_year <- xml_text(xml_find_all(read_xml(as.character(b[[i]])),"//ContributionDate//Year"))
if(length(pub_year)==0){
  pub_year <- NA
}
            pmid <- xml_text(xml_find_all(read_xml(as.character(b[[i]])),"//BookDocument//PMID"))
if(length(pmid)==0){
  pmid <- NA
}
            title <- xml_text(xml_find_all(read_xml(as.character(b[[i]])),"//BookDocument//ArticleTitle"))
if(length(title)==0){
  title <- NA
}
            
            # the authorlist include the authors and editors groups, so I just extract the authors group
            authorlist <- xml_find_all(read_xml(as.character(b[[i]])),"//AuthorList")
            for(j in seq_along(authorlist)){
                if(xml_attr(authorlist[[j]],"Type")=='authors'){
                    authors <- xml_find_all(read_xml(as.character(authorlist[[j]])), "//Author")
                    author_list <- c()
                    for(k in seq_along(authors)){
                        author_f <- xml_text(xml_find_all(read_xml(as.character(authors[[k]])), "//ForeName"))
                        author_l <- xml_text(xml_find_all(read_xml(as.character(authors[[k]])), "//LastName"))
                        author_list <- c(author_list , paste(author_f, author_l, sep = " "))
                    }
                }
            }
            # binding the new row to the dataframe pmid, title, pub_year, abstract, "gene"
            # Note: in the line below, I just added the "gene" as a string, we need to put a variable here
            pubmedDf <- rbind(pubmedDf,as.data.frame(t(c(pmid, title, pub_year, abstract, g))))
        }
    }
    
}
x <- c("pmid","title","year","abstract","gene")
colnames(pubmedDf) <- x
    



names(pubmedDf)
pubmedDf
output <- as.data.table(pubmedDf)
output <- output[,list(gene,pmid,year,title)]
setnames(output, "gene","ALIAS")
setkey(output,"ALIAS")
setkey(lookup.table,"ALIAS")
output <- lookup.table[output]

pubmedDT <- copy(output) #to be used in the stats part
setkeyv(output,c("gene.symbols","pmid"))
output <- unique(output,by=c("gene.symbols","pmid"))
output[,LINK:=paste0("https://www.ncbi.nlm.nih.gov/pubmed/",pmid)]
output <- output[,list(gene.symbols,ALIAS,year,LINK,title)]
setnames(output,"gene.symbols","SYMBOL")
write.csv(output,"pubmed-search.csv",quote=FALSE,row.names=FALSE)

num_genes=length(unique(output$SYMBOL))

q()
setnames(pubmedDT,"gene.symbols","gene")
pubmedDT[,ALIAS:=NULL]
summGenes <- as.data.frame(pubmedDT) %>% group_by(gene) %>% summarise(cnt = n_distinct(title)) #why not distinct by pmid?
#summGenes <- pubmedDf %>% group_by(gene) %>% summarise(cnt = n_distinct(title))
c <- tm1[which(tm1$gene %in% summGenes$gene[which(summGenes$cnt>0)]),]

c
# need to have a way to check for an empty table. 

similarityCalc <- function(abst){
    #will need to also update this to use an input file
    train <- data.frame(id = c(1), text = c(paste0(query_terms,collapse=' '),abst), Text1 = c('a','b'))
    #train <- data.frame(id = c(1), text = c('Flexion contracture Umbilical hernia Camptodactyly Talipes equinovarus Focal seizures',abst), Text1 = c('a','b'))
    # Tokenize abstracts.
    train.tokens <- tokens(as.character(train$text), what = "word", 
                           remove_numbers = TRUE, remove_punct = TRUE,
                           remove_symbols = TRUE, remove_hyphens = TRUE, remove_separators = TRUE)
    #Added by EAR
    train.tokens <- tokens_remove(train.tokens, "\\p{Z}", valuetype = "regex")
    
    # Lower case the tokens.
    train.tokens <- tokens_tolower(train.tokens)
    
    # Use quanteda's built-in stopword list for English.
    # NOTE - You should always inspect stopword lists for applicability to
    #        your problem/domain.
    train.tokens <- tokens_select(train.tokens, stopwords(), 
                                  selection = "remove")
    
    # Perform stemming on the tokens.
    train.tokens <- tokens_wordstem(train.tokens, language = "english")
    
    
    
    # Create our first bag-of-words model.
    train.tokens.dfm <- dfm(train.tokens, tolower = FALSE)
    
    
    
    
    # Transform to a matrix and inspect.
    train.tokens.matrix <- as.matrix(train.tokens.dfm)
    
    # Our function for calculating relative term frequency (TF)
    term.frequency <- function(row) {
        row / sum(row)
    }
    
    # Our function for calculating inverse document frequency (IDF)
    inverse.doc.freq <- function(col) {
        corpus.size <- length(col)
        doc.count <- length(which(col > 0))
        if (doc.count > 0){
            return(1+log10(corpus.size / doc.count))
        }
        else{
            return(1)
        }
    }
    
    # Our function for calculating TF-IDF.
    tf.idf <- function(x, idf) {
        x * idf
    }
    
    
    # First step, normalize all documents via TF.
    train.tokens.df <- apply(train.tokens.matrix, 1, term.frequency)
    
    
    
    # Second step, calculate the IDF vector that we will use - both
    # for training data and for test data!
    train.tokens.idf <- apply(train.tokens.matrix, 2, inverse.doc.freq)
    
    
    
    # Lastly, calculate TF-IDF for our training corpus.
    train.tokens.tfidf <-  apply(train.tokens.df, 2, tf.idf, idf = train.tokens.idf)
    
    
    
    # Transpose the matrix
    train.tokens.tfidf <- t(train.tokens.tfidf)
    a <- train.tokens.tfidf[1,]
    b <- train.tokens.tfidf[2,]
    
    theta <- ( sum(a*b) / ( sqrt(sum(a * a)) * sqrt(sum(b * b)) ) )
    return(theta)
}

simDf <- data.frame()
for(j in summGenes$gene){
#j <- 1
    my_query <- paste(base_query,' AND ',i,'[Title/Abstract]',sep = '')
    #my_query <- paste('((((( Flexion contracture[Title/Abstract]) OR Umbilical hernia[Title/Abstract]) OR Camptodactyly[Title/Abstract]) OR Talipes equinovarus[Title/Abstract]) OR Focal seizures[Title/Abstract])   AND ',j,'[Title/Abstract]',sep = '')
    #my_query <- paste('(((((Ichthyosis[Title/Abstract]) OR Palmoplantar hyperkeratosis[Title/Abstract]) OR Anhidrosis[Title/Abstract]) OR Erythroderma[Title/Abstract]) OR Ectropion[Title/Abstract]) AND ',j,'[Title/Abstract]',sep = '')
    first_PM_records <- get_pubmed_ids(my_query)   # submit the query to PubMed
    totalXML <- fetch_pubmed_data(first_PM_records)      # retrieve the first output record
    
    rootnode <- xmlRoot(xmlTreeParse(totalXML, useInternalNodes = TRUE))
    
    xml_cap <- capture.output(rootnode)
    val_df <- as.data.frame(xml_cap[which(grepl('AbstractText|PMID',xml_cap))])
    names(val_df)[1] <- 'abst'
    val_df_sep <- val_df %>% separate(abst, into = paste0('AbstractText', 1:4), sep = '[><]')
    val_df_sep$rowtype <- NA
    val_df_sep$AbstractText1 <- seq.int(nrow(val_df_sep))
    val_df_sep$rowtype[which(grepl('PMID',val_df_sep$AbstractText4))] <- val_df_sep$AbstractText1[which(grepl('PMID',val_df_sep$AbstractText4))]
    
    val_df_sep <- val_df_sep %>% fill(rowtype) 
    
    val_df_sep_agg <- aggregate(val_df_sep$AbstractText3, list(val_df_sep$rowtype), paste, collapse=" ")
    val_df_sep_agg_filtered <- val_df_sep_agg[which(nchar(val_df_sep_agg$x)>10),'x']
    
    simlist <- c()
    val_df_sep_agg_filtered
# It looks like I have non-English words
# for now, see if it works with the first 10 (hopefully they are ok
# the first one has an unknown character. Maybe that is the problem
#    val_df_sep_agg_filtered <- val_df_sep_agg_filtered[2:10]
# I'm getting the same error. Hmmm.
#    val_df_sep_agg_filtered
    for (i in val_df_sep_agg_filtered){

        simDf <- rbind(simDf,data.frame(t(c('X1' = j, 'X2' = similarityCalc(i)))))
    }
    Sys.sleep(2)
}

simDf
simDf <- simDf[-4,]

simDf$X2 <- as.numeric(as.character(simDf$X2))

meanSim <- simDf %>%
    group_by(X1) %>%
    summarise(count=n(),m = mean(X2))
meanSim
output2 <- as.data.table(meanSim)
write.csv(output2,"meanSim.csv",quote=FALSE, row.names=FALSE)
print("simDf\n")
q()
ggplot(meanSim,aes(m,count,label=X1))+
    geom_point(color = "blue", size = 3)+
    geom_point(data=meanSim[meanSim$X1=='KRT10',],
               pch=21, fill=NA, size=37, colour="red", stroke=2) +
    geom_text(aes(label=X1),hjust=-0.1, vjust=-0.1)+
    xlab('Average of similarity measures')+
    ylab('No. of documents retrieved')+
    theme_bw(base_size = 16)+
    scale_x_continuous(limits = c(0, 0.1))+
    scale_y_continuous(limits = c(0, 40))



