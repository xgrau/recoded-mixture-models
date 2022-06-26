# input
aa_order = c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
model_list = c("C10","C20","C30","C40","C50","C60")

# recoding dictionaries
recoding_dictionary = list(
  "SR4" = list(
    "A" = "AGNPST",
    "C" = "CHWY",
    "G" = "DEKQR",
    "T" = "FILMV"
  ),
  "Dayhoff6" = list(
    "0" = "AGPST",
    "1" = "DENQ", 
    "2" = "HKR", 
    "3" = "ILMV",
    "4" = "FWY",
    "5" = "C"
  ),
  "Dayhoff9" = list(
    "0" = "DEHNQ",
    "1" = "ILMV",
    "2" = "FY",
    "3" = "AST",
    "4" = "KR",
    "5" = "G",
    "6" = "P",
    "7" = "C",
    "8" = "W"
  ),
  "SR6" = list(
    "0" = "APST",
    "1" = "DENG",
    "2" = "QKR",
    "3" = "MIVL",
    "4" = "WC",
    "5" = "FYH"
  ),
  "KGB6" = list(
    "0" = "AGPS",
    "1" = "DENQHKRT",
    "2" = "MIL",
    "3" = "W",
    "4" = "FY",
    "5" = "CV"
  )
)

# references:
# Dayhoff6: M.O. Dayhoff, R.M. Schwartz, B.C. Orcutt, A model of evolutionary change in proteins
# Dayhoff9: Alexandra M Hernandez, Joseph F Ryan,  Six-State Amino Acid Recoding is not an Effective Strategy to Offset Compositional Heterogeneity and Saturation in Phylogenetic Analyses
# SR6: On reduced amino acid alphabets for phylogenetic inference, Mol. Biol. Evol., 24 (2007), pp. 2139-2150
# KGB6: A new criterion and method for amino acid classification J. Theor. Biol., 228 (2004), pp. 97-106

for (model in model_list) {
  
  # read model
  mod_fn = sprintf(sprintf("data/models_iqtree_%s.txt", model))
  mod = read.table(mod_fn, fill = TRUE)
  
  # get nexus final line and drop from table
  nex_s = mod [ nrow(mod), ]
  nex_s = nex_s [ !is.na(nex_s) ]
  nex_s = paste(nex_s, collapse = " ")
  mod = mod [ 1:(nrow(mod) - 1), ]
  
  # remove nonnumeric columns
  mon = mod [ , 4:ncol(mod) ]
  mon[, ncol(mon)] = gsub(";","",mon[, ncol(mon)])
  mon = apply(mon, 2, as.numeric)
  
  # rownames are CXX categories, colnames are aa
  colnames(mon) = aa_order
  rownames(mon) = mod[,2]
  
  
  mon_t = t(mon)
  # mod_t$aa = rownames(mod_t)
  
  for (dic in names(recoding_dictionary)) {
    
    recoding_dictionary_v = unlist(lapply(1:length(recoding_dictionary[[dic]]), function(i) {
      
      new = names(recoding_dictionary[[dic]])[i]
      exp = unlist(strsplit(recoding_dictionary[[dic]][[new]], split = ""))
      new = rep(new, length(exp))
      names(new) = exp
      return(new)
      
    }  ))
    
    
    # vector of recoded categories
    recoded_categories = as.character(recoding_dictionary_v [ rownames(mon_t) ])
    
    # sum over recoded categories
    mon_t_r = aggregate(mon_t, by = list( recoded_categories = recoded_categories ), FUN=sum)
    mon_r = as.matrix(data.frame(t(mon_t_r)))
    colnames(mon_r) = mon_r[1,]
    mon_r = mon_r [2:nrow(mon_r),]
    mon_r = apply(mon_r, 2, as.numeric)
    rownames(mon_r) = rownames(mon)
    
    # normalise to sum 1
    mon_r_norm = mon_r / rowSums(mon_r)
    
    # format nexus
    nex_m = mon_r_norm
    rownames(nex_m) = sprintf("frequency %s%s_%s = ", model, dic, rownames(mon_r_norm))

    # write nexus for GTR
    nex_o_string = paste(sprintf("%s%s_%s", model, dic, rownames(mon_r_norm)), collapse = ",")
    nex_o = sprintf("model xm%s%s = GTR+G+FMIX{%s} ;", model, dic, nex_o_string)
    nex_fn = sprintf("recoded_models/xm%s%s.nex", model, dic)
    write("#nexus",        file = nex_fn)
    write("begin models ;", file = nex_fn, append = TRUE)
    for (i in 1:nrow(nex_m)) {
      cat_s = sprintf("%s %s ;", rownames(nex_m)[i], paste(nex_m[i,], collapse = " "))
      write(cat_s, file = nex_fn, append = TRUE)
    }
    write(nex_o, file = nex_fn, append = TRUE)
    write("end ;", file = nex_fn, append = TRUE)
    
    # write nexus for Poisson
    nex_o_string = paste(sprintf("%s%s_%s", model, dic, rownames(mon_r_norm)), collapse = ",")
    nex_o = sprintf("model xmPoi%s%s = Poisson+G+FMIX{%s} ;", model, dic, nex_o_string)
    nex_fn = sprintf("recoded_models/xmPoi%s%s.nex", model, dic)
    write("#nexus",        file = nex_fn)
    write("begin models ;", file = nex_fn, append = TRUE)
    for (i in 1:nrow(nex_m)) {
      cat_s = sprintf("%s %s ;", rownames(nex_m)[i], paste(nex_m[i,], collapse = " "))
      write(cat_s, file = nex_fn, append = TRUE)
    }
    write(nex_o, file = nex_fn, append = TRUE)
    write("end ;", file = nex_fn, append = TRUE)
    
  }
  
}

