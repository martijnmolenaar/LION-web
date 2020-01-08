convertLipidNames <- function(name){
  
  options(stringsAsFactors = FALSE)
  
  sapply(name, function(input) {
    processed_input <- input
    if (processed_input==""){processed_input <- " "}
    
    #aliases 
    {
      processed_input <- gsub("DAG|diacylglycerol","DG",processed_input)
      processed_input <- gsub("TAG|triacylglycerol","TG",processed_input)
      
      processed_input <- gsub("GPA","PA",processed_input)
      processed_input <- gsub("GPEtn|phosphatidylethanolamine|Phosphatidylethanolamine","PE",processed_input)
      processed_input <- gsub("GPCho|phosphatidylcholine|Phosphatidylcholine","PC",processed_input)
      processed_input <- gsub("GPSer|phosphatidylserine|Phosphatidylserine","PS",processed_input)
      processed_input <- gsub("GPIns|phosphatidylinositol|Phosphatidylinositol","PI",processed_input)
      processed_input <- gsub("GPGro|phosphatidylglycerol|Phosphatidylglycerol","PG",processed_input)
      
      
      processed_input <- gsub("LBPA","BMP",processed_input)
      processed_input <- gsub("lyso|Lyso","L",processed_input)
      
      processed_input <- gsub("^FA","FFA",processed_input)
      processed_input <- gsub("ChE|CholE","CE",processed_input)
      
      processed_input <- gsub("\\Q(3'-sulfo)Galβ-Cer\\E|\\Q(3'-sulfo)LacCer\\E|\\Q(3'-sulfo)GalCer\\E|Sulfogalactosyl ceramide[ ]*|[Ss]+(Gal|Hex)Cer|Sulfatide",
                              "SHexCer",processed_input)
      
      
      processed_input <- gsub("LacCer","Hex2Cer",processed_input)
      processed_input <- gsub("Gl[u]*c[β-]*Cer|Gal[β-]*Cer|GluCer","HexCer",processed_input)
      
      processed_input <- gsub("^aPE","NAPE",processed_input)
      
    }    ### end aliases
    
    ### special cases:  'DG(0:0/18:1/16:0)' >> DG(18:1/16:0)
    if(grepl('DG',processed_input) & grepl('\\D0:0',processed_input)){
      processed_input <- gsub('/0:0',"",processed_input)
      processed_input <- gsub('\\(0:0/',"(",processed_input)
      
    }
    
    ### end special cases
    if(grepl("\\d+:\\d+;\\d+",processed_input) & !grepl("\\d+:\\d+;\\d+:\\d+",processed_input)){    ### 32:1;1 format >> 32:1
      processed_input <- gsub(";\\d+","",processed_input)
    }
    if(grepl("\\d+:\\d+[cΔ]\\d+",processed_input) & !grepl("\\d+:\\d+[cΔ]\\d+:\\d+",processed_input)){    ### 32:1c1 format >> 32:1
      processed_input <- gsub("[cΔ]\\d+","",processed_input)
    }
    ###############    
    if(grepl("\\d+:\\d+:\\d+",processed_input)){    ### 32:1;1 format >> 32:1
      #processed_input <- gsub(";\\d+","",processed_input)
      parts <- unlist(strsplit(x = paste("_",processed_input,"_",sep=""), split = "\\d+:\\d+:\\d+"))
      specialFAs <- unlist(regmatches(processed_input, gregexpr("\\d+:\\d+:\\d+", processed_input)))
      specialFAs <- unlist(regmatches(processed_input, gregexpr("\\d+:\\d+", processed_input)))
      processed_input[seq(from = 1, to = length(parts)*2, by = 2)] <- parts
      processed_input[seq(from = 2, to = length(specialFAs)*2, by = 2)] <- specialFAs
      processed_input <- gsub("^_|_$","",paste(processed_input, collapse = ""))
      
    }
    
    ###  >> remove (FA 16:0), not intensively tested
    processed_input <- gsub("\\(FA( )*", "(",processed_input)   
    
    ### 'CL(32:0)(34:1)' or 'TAG 32:0(FA 16:0)' >> sum FAs between (), not intensively tested
    FAs <- unlist(regmatches(processed_input, gregexpr("\\d+:\\d+", processed_input)))
    if(length(FAs)>0){
      FAs <- FAs[sapply(regmatches(FAs, gregexpr("\\d+", FAs)), function(i){  as.numeric(i)[1] >= as.numeric(i)[2]  })]  ## C should be higher than DB  
    }
    if((grepl("CL", processed_input ) & (length(FAs) == 2 | length(FAs) == 3)) |
       (grepl("TAG|TG", processed_input ) & (length(FAs) == 2)))     {
      parts <- unlist(regmatches(processed_input, gregexpr("[\\(]*\\d+:\\d+[\\)]*", processed_input)) )
      prefix <- unlist(strsplit(processed_input, "[\\(]*\\d+:\\d+[\\)]*"))[1]
      suffix <- unlist(strsplit(processed_input, "[\\(]*\\d+:\\d+[\\)]*"))[2]
      processed_input <- paste(prefix, "(",
                               sum(sapply(regmatches(
                                 parts, gregexpr("\\d+", parts)
                               ), function(i) {
                                 as.numeric(i)[1]
                               })), ## number of Cs
                               ":",
                               sum(sapply(regmatches(
                                 parts, gregexpr("\\d+", parts)
                               ), function(i) {
                                 as.numeric(i)[2]
                               })), ## number of DBs
                               ")", suffix, sep = "")
    }
    ### end 'CL(32:0)(34:1)'
    
    
    if(length(FAs)==2 & grepl("^SM",    processed_input)){                ### voor 'SM 18/16:0)'
      if(!grepl("[dt]\\d+",processed_input)){                             ### double check: no 'd'18...?
        processed_input <- gsub("SM 18","SM d18",processed_input)         ### add 'd'
        processed_input <- gsub("SM\\(18","SM(d18",processed_input)
      }
      
    }
    
    processed_input <- gsub("_\\d*\\.*\\d*/\\d*\\.*\\d*", "",processed_input)   ## remove amu/RT info
    
    processed_input <- gsub(";\\d$","",processed_input)   ### in case of SHexCer 30:2;1
    
    if (grepl("DHHexosylceramide|Hexosylceramide|DHSphingomyelin|Sphingomyelin|DHCeramide|Ceramide", processed_input)){   ## for notation 'DHSphingomyelin 16'
      if (grepl("C\\d+",processed_input)){   ## replace 'C20:3 to 20:3
        C_containing_notation <- unlist(regmatches(processed_input,  gregexpr("C\\d+", processed_input)))
        processed_input <- sub(C_containing_notation,  gsub('C',"",C_containing_notation), processed_input)
      }
      
      if (!grepl("\\d+:\\d+",processed_input)){   ## replace 'C16 to 16:0
        C_number <- regmatches(processed_input, regexpr("\\d+", processed_input)) 
        processed_input <- sub(C_number, paste(C_number,":0",sep=""), processed_input)
      }
      
      if (grepl("^DH",processed_input)){   ## DH >> d18:0, otherwise >> d18:1
        processed_input <- gsub("DH","",processed_input)
        FAs <- regmatches(processed_input, regexpr("\\d+.*\\d*", processed_input)) 
        processed_input <-  sub(FAs, paste("d18:0/",FAs,sep=""), processed_input)
      } else {
        #processed_input <- gsub("DH","",processed_input)
        FAs <- regmatches(processed_input, regexpr("\\d+.*\\d*", processed_input)) 
        processed_input <-  sub(FAs, paste("d18:1/",FAs,sep=""), processed_input)
      }
      
      processed_input <- gsub("Hexosylceramide","HexCer",processed_input)
      processed_input <- gsub("Sphingomyelin","SM",processed_input)
      processed_input <- gsub("Ceramide","Cer",processed_input)
    }
    
    processed_input <-
      gsub("_|;", "/", processed_input)                                ## _ and ; into /
    
    for(substitute in unlist(regmatches(processed_input, gregexpr("\\d+-\\d+", processed_input)) )){
      processed_input <- gsub(substitute, gsub("-","/",substitute), processed_input) 
    }   ## - into / if it's between numbers
    
    
    processed_input <-
      gsub("\\((\\d+[EZ]{1},*)+\\)", "", processed_input)   ## removing (5Z,8Z,11Z,14Z)-like info
    
    if (grepl("\\d+:\\d+[pe]{1}", processed_input)){             ### ether lipids: in case of 18:1p in stead of P-18:1
      
      ether_FAs <- unlist(regmatches(processed_input, gregexpr("\\d+:\\d+[pe]{1}", processed_input)))  ## extract 18:1
      
      for(i in 1:length(ether_FAs)){
        processed_input <-
          gsub(ether_FAs[i],              ## replace 18:1p for P-18:1
               paste(
                 ifelse(grepl("p", ether_FAs[i]), "P-", "O-"),
                 regmatches(ether_FAs[i], regexpr("\\d+:\\d+", ether_FAs[i])),
                 sep = ""
               ),
               processed_input)
      }
    }
    
    if (!grepl("\\(", processed_input)) {
      ## if they don't contain /
      ## process PC 18:1/20:4 format to PC(..)
      
      middlePart <-         ### ..O-18:1/16:0
        regmatches(processed_input,
                   regexpr("[dtAOP-]*(\\d+:\\d+/*)+", processed_input))
      otherParts <-         ## SM
        unlist(strsplit(processed_input, "[dtAOP-]*(\\d+:\\d+/*)+"))
      otherParts <- gsub(" ", "", otherParts)
      text <- ""
      
      for (i in 1:max(c(length(otherParts),   length(middlePart)))) {
        ## reconstruct
        if (!is.na(otherParts[i])) {
          text <- paste(text, otherParts[i], sep = "")
        }
        
        if (!is.na(middlePart[i])) {
          text <- paste(text, "(", middlePart[i], ")", sep = "")
        }
      }
      processed_input <- text
    }
    
    
    if(grepl("-DiHex",processed_input)){          # in case of d18:0/12:0-DiHex
      processed_input <- gsub("-DiHex","",processed_input)
      processed_input <- paste("Hex2Cer",processed_input,sep="")
    }
    
    if(grepl("-MonoHex",processed_input)){          # in case of d18:2/26:0-MonoHex
      processed_input <- gsub("-MonoHex","",processed_input)
      processed_input <- paste("HexCer",processed_input,sep="")
    }
    
    if(grepl("\\D+ [\\(|\\d+]", processed_input)){        ## in case of 'CL (70:2)'
      processed_input <- gsub("\\D+ [\\(|\\d+]", 
                              gsub(" ","", regmatches(processed_input, regexpr("\\D+ [\\(|\\d+]", processed_input)) ),
                              processed_input)
    }
    
    if (grepl("^L[PC].*\\((\\d+:\\d+/*)+", processed_input)) {
      ## if LPC format >> format to PC(xx:0:0)
      
      
      middlePart <-         ### ..18:1 + /0:0
        paste(regmatches(
          processed_input,
          #regexpr("(\\d+:\\d+/*)+", processed_input)
          regexpr("([OP-]*\\d+:\\d+[ep/]*)+", processed_input)
        ), sep = "")
      
      if(identical(middlePart,character(0))){middlePart <- ""}
      ### if no 0:0 yet, add this
      if(!grepl("\\D0:0|^0:0", middlePart)){middlePart <- paste(middlePart,"/0:0",sep="")}
      
      
      otherParts <-         ## LPC
        unlist(strsplit(processed_input, "([OP-]*\\d+:\\d+[ep/]*)+"))
      otherParts <- gsub("^L", "", otherParts)
      if(identical(otherParts,character(0))){otherParts <- ""}
      
      text <- ""
      
      for (i in 1:max(c(length(otherParts),   length(middlePart)))) {
        ## reconstruct
        if (!is.na(otherParts[i])) {
          text <- paste(text, otherParts[i], sep = "")
        }
        
        if (!is.na(middlePart[i])) {
          text <- paste(text, middlePart[i], sep = "")
        }
      }
      processed_input <- text
      
      
    }
    
    
    if(grepl('^P\\D\\(A-*',processed_input) ){      ########## ether lipids 20180605
      processed_input <- gsub('\\(A-*',"(O-",processed_input)    ##gsub('A-*',"(O-",processed_input)
    }
    if(grepl('^P\\D\\([OP]\\d+',processed_input) ){      ########## ether lipids 20180219
      processed_input <- gsub('\\(O',"(O-",processed_input)
      processed_input <- gsub('\\(P',"(P-",processed_input)
    }
    
    if(grepl("[[:upper:]]{1}[[:lower:]]{2}.+(sterol|enone)", processed_input)){         ### in case of Lanosterol / Cholestenone
      UpperLower <- regmatches(processed_input, regexpr("[[:upper:]]{1}[[:lower:]]{2}", processed_input)) 
      processed_input <- gsub("[[:upper:]]{1}[[:lower:]]{2}", tolower(UpperLower),processed_input )
    }
    
    
    ## reorder FAs
    if (         
      (length(unlist(regmatches(processed_input, gregexpr("\\d+:\\d+",        processed_input)))) > 1) &                         ## when more than 1 FA)
      (length(unlist(regmatches(processed_input, gregexpr("[dtm](\\d+:\\d+)|[PO]-(\\d+:\\d+)", processed_input)))) < 1) 
      #  (!grepl("0:0", processed_input))                           ## unless it's containing a 0:0 FA (=lyso, should be on second position))
    ){  
      FAs <- t(as.data.frame(strsplit(unlist(
        regmatches(
          processed_input,
          gregexpr("\\d+:\\d+", processed_input)
        )
      ), ":")))
      FAs <- as.data.frame(sapply(as.data.frame(FAs), as.numeric))
      
      ## added to sort lysoPCs
      FAs[apply(FAs,1,function(i){paste(i, collapse = ":")}) == "0:0",1] <- 999  ## set 0:0 to 999:0 to order right
      
      FAs <- FAs[order(FAs$V1, FAs$V2), ]
      ## added to sort lysoPCs
      FAs[apply(FAs,1,function(i){paste(i, collapse = ":")}) == "999:0",1] <- 0  ## reset 
      
      FAs <- apply(FAs, 1, function(row) {
        paste(row[1], ":", row[2], sep = "")
      })
      
      parts <- unlist(strsplit(processed_input, "\\d+:\\d+"))
      text <- ""
      
      for (i in 1:length(parts)) {
        if (!is.na(parts[i])) {
          text <- paste(text, parts[i], sep = "")
        }
        
        if (!is.na(FAs[i])) {
          text <- paste(text, FAs[i], sep = "")
        }
      }
      text
    } else {
      processed_input
    }
    
  })
  
  
  
}       ## 20191106 updated in apps

simplifyLipidnames <- function(name){
  options(stringsAsFactors = FALSE)
  
  sapply(name, function(input) {
    processed_input <- input
    
    to_simplify <-
      unlist(regmatches(processed_input,    gregexpr("(iso-)*[dt]*\\d+:\\d+([/_]+\\d+:\\d+)+", processed_input)    ))
    to_simplify <- to_simplify[!grepl("^0:0|[_/]+0:0",to_simplify)]
    
    if (length(to_simplify) > 0) {
      simplified <-  sapply(to_simplify, function(part_i) {
        FAs <- unlist(regmatches(part_i,  gregexpr("\\d+:\\d+", part_i)))
        Cs <-  sum(sapply(strsplit(FAs, split = ":"), function(i) {  as.numeric(i[1])     }))
        DBs <- sum(sapply(strsplit(FAs, split = ":"), function(i) {  as.numeric(i[2])  }))
        paste(Cs, ":", DBs, sep = "")
      })
      
      nameLUT <-
        data.frame(original = paste('\\Q', to_simplify, '\\E', sep = ""),
                   simplified = simplified)
      
      for (i in 1:dim(nameLUT)[1]) {
        processed_input <-
          gsub(nameLUT$original[i], nameLUT$simplified[i], processed_input)
      }
      
      processed_input
    } else {processed_input}
  })
  
}   ### 20191001 

getFAs <- function(lipids){
  output <- sapply(lipids, function(lipid_i){
    FAs <- unlist(regmatches(lipid_i, gregexpr("(iso-)*[dt]*\\d+:\\d+", lipid_i)))
    FAs[grepl("^\\d+",FAs)] <- paste("C",FAs[grepl("^\\d+",FAs)], sep = "")
    FAs <- gsub("d18:0","d18:0 (dihydrosphingosine)",FAs)
    FAs <- gsub("d18:1","d18:1 (sphingosine)",FAs)
    
    FAs
  }, simplify = FALSE)
  if(length(output)<2){
    output <- unlist(output)
    names(output) <- NULL
    output
  } else {
    output
  }
}

getLipidHeadgroup <- function(lipids){      ### 20191030
  sapply(lipids, function(lipid_i){
    HG <- gsub("\\(|\\[.+\\]","",regmatches(lipid_i, regexpr("^.+\\(", lipid_i)) )
    if(length(HG)==0){HG <- ""}
    
    if(HG %in% c("PC","PE","PA","PG","PI","PS","CL") &
       any(getFAs(lipid_i) %in% "C0:0")){
      HG <- paste("L",HG,sep = "")
    }
    HG
  })
  
}

getLipidLinkage <- function(lipids){
  ifelse(grepl("P-\\d+", lipids), "P",
         ifelse(grepl("O-\\d+", lipids), "O", "D"))
}

predict_FAs <- function(headgroup = NULL, summedFA = NULL, composition_table = NULL){    ## v 20191101
  FA_n <-   ## number of FAs associated, based on headgroup
    switch(
      c(which(c(
        headgroup %in% c("PA", "PS", "PG", "PE", "PC", "PI", "BMP", "DG"),
        headgroup %in% c("TG","aPG","NAPE"),
        headgroup %in% c("CL")
      )),4)[1], 2, 3, 4,NA)
  
  FA_output <- composition_table$FAs[composition_table$nr_of_FAs == FA_n &   composition_table$FA_summed == summedFA][1]
  
  if(!is.na(FA_output)){
    return(paste("C",strsplit(x = FA_output, split = "\\|")[[1]],sep = ""))
  } else {
    return(NULL)
  }
  
}

sendEmail <- function(subject = "LION/web usage", mail_message = "empty", from = "system"){
  url <- "######"
  api_key <- "######"
  the_body <-
    list(
      from="####",
      to="aaa@xxx.nl",
      subject=subject,
      text=paste(mail_message,"\n\n=====================\nfrom: ",from,"\n=====================\nLION/web: www.lipidontology.com", sep= "")
    )
  req <- httr::POST(url, httr::authenticate("api", api_key),  encode = "form",  body = the_body)
  httr::stop_for_status(req)
  TRUE
}   ## mail function

associatedTerms <- function(lipid, ontologyObject, all = FALSE, reformat = FALSE) {    ## function: which LION-terms are associated with a lipid?
  lipidInTerms <- genesInTerm(ontologyObject)
  if (!reformat) {
    associations <- sapply(lipid, function(lipid_i) {
      terms <-
        ontologyObject@termName[names(ontologyObject@termName) %in% names(lipidInTerms)[sapply(lipidInTerms, function(term) {
          any(term == lipid_i)
        })]]   ## find terms
      
      if (!all) {
        terms <- terms[grepl("LION", names(terms))]
      }     ## is all == TRUE, than remove CAT and all  terms
      return(terms)
    })
    return(associations)
    
  } else if (reformat) {
    
    terms <- names(lipidInTerms)
    if (!all) {
      terms <- terms[grepl("LION", terms)]
    }     ## is all == TRUE, than remove CAT and all  terms
    
    lipid_term_table <- t(sapply(lipid, function(lipid_i) {
      terms_per_lipid <-
        names(ontologyObject@termName[names(ontologyObject@termName) %in% names(lipidInTerms)[sapply(lipidInTerms, function(term) {
          any(term == lipid_i)
        })]])
      return(ifelse(terms %in% terms_per_lipid, "x", ""))
    }))
    
    names <- ontologyObject@termName[match(terms, names(ontologyObject@termName))]
    
    lipid_term_table <- rbind(terms, names, lipid_term_table)
    lipid_term_table <- cbind(rownames(lipid_term_table), lipid_term_table)
    rownames(lipid_term_table) <- NULL
    colnames(lipid_term_table) <- NULL
    lipid_term_table <- as.data.frame(lipid_term_table)
    lipid_term_table[[1]][1:2] <- c("LION-term", "LION-name")
    return(lipid_term_table)
  }
  
} ## which LION-terms are associated with a lipid?

