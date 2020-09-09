library(crul)
library(rvest)
library(phylocomr)
library(dplyr)
library(phangorn)
library(stringr)

timetree_pair <- function(tax1, tax2){
  x <- HttpClient$new(url = "http://www.timetree.org")
  
  res_get0 <- x$get(paste0("ajax/names/", tax1, "/", tax2))
  res_get <- read_html(res_get0$parse(), encoding = "UTF-8")
  
  tax1_val <- res_get %>%
    html_node("#pairwise-resolve-taxon-a1 > option:nth-child(1)") %>%
    html_attr("value") 
  
  tax2_val <- res_get %>%
    html_node("#pairwise-resolve-taxon-b1 > option:nth-child(1)") %>%
    html_attr("value") 
  
  url <- paste0("http://www.timetree.org/ajax/pairwise/", tax1_val, "/", tax2_val)
  
  res <- read_html(url) %>% html_node("#pairwise-results") %>% html_text()
  if(is.na(res)){
    out = "No pairwise results found"
  } else {
    if(stringr::str_detect(res, "CI: [(]\\d")){
      out = stringr::str_replace(res, "^.*Median Time:([.0-9]+ MYA)[ ]?Estimated Time:([.0-9]+ MYA).*CI: .([0-9 -.]*)MYA.*derived.*from ([.0-9]+ [a-z]+).*$", 
                                 "Median Time: \\1, Estimated Time: \\2, From: \\4, CI: \\3")
    } else {
      out = stringr::str_replace(res, "^.*Median Time:([.0-9]+ MYA)[ ]?Estimated Time:([.0-9]+ MYA).*derived .*from ([.0-9]+ [a-z]+).*$", 
                                 "Median Time: \\1, Estimated Time: \\2, From: \\3,  - ")
    }
  }
  return(out)
}

# timetree_pair(tax1 = "Gagea_fragifera", tax2 = "Allium_strictum")

timetree_phylo <- function(phy){
  nt <- length(phy$tip.label)
  nn <- Nnode(phy)
  nodes = phy$node.label
  nodes[which(nodes == "")] = paste0("node", which(nodes == ""))
  l <- vector("list", nn)
  names(l) = nodes
  for (i in 1:nn){
    cat("Node ", nt+i, "\n")
    d <- unlist(phangorn::Descendants(phy, nt+i, type="tips"))
    taxa <- phy$tip.label[d][c(1,length(d))]
    div <- timetree_pair(taxa[1], taxa[2])
    if(div == "No pairwise results found"){
      l[[i]] = NA
    } else {
      l[[i]] = unlist(stringr::str_split(div, ","))
      names(l[[i]]) = c("median", "estimated", "n_studies", "CI")
    }
  }
  tmp_f = tempfile(fileext = ".rds")
  saveRDS(l, tmp_f)
  ll0 = data_frame(node = names(which(is.na(l))), median = NA, estimated = NA, n_studies = NA, CI = NA)
  ll1 = try(dplyr::bind_rows(!!! l[which(!is.na(l))], .id = "node"))
  if("try-error" %in% class(ll1)) {
    stop("Error when converting list to data.frame; the list saved as ", tmp_f)
  }
  ll = bind_rows(ll1, ll0)
  ll = mutate_at(ll, 2:4, ~as.numeric(stringr::str_extract(., "[0-9.]+")))
  ll = bind_cols(ll, 
                 stringr::str_replace_all(ll$CI, "CI: | ", "") %>% 
                   stringr::str_split("-") %>% do.call(rbind, .) %>% as_data_frame() %>% 
                   setNames(c("age.min", "age.max")) %>% 
                   mutate(age.min = as.numeric(age.min), age.max = as.numeric(age.max))
  )
  phy$node.label <- nodes
  unlink(tmp_f)
  return(list(info = ll, phy = phy))
}

itree_otl = timetree_phylo(insect_tree)

est_age <- itree_otl$info %>% filter(!is.na(median)) %>% 
  select(node, estimated) 

# if you have Phylocom in your computer:
system("cd data && ~/Documents/phylocom-4.2/src/phylocom bladj -f tree_otl_nodes.tre > tree_otl2.tre")
# if you use R package phylocomr:
phylocomr::ph_bladj(ages = est_age, phylo = itree_otl$phy)
tree_otl = read.tree("data/tree_otl2.tre")
tree_otl$root.edge = max(phytools::nodeHeights(tree_otl))
is.ultrametric(tree_otl)
plot(tree_otl,  show.tip.label = F)
plot(read.tree("data/tree_otl.tre"), show.tip.label = F)

# itree_otl$Nnode
# length(itree_otl$tip.label)

# url = "http://www.timetree.org/ajax/names/Acer%20rubrum/Acer%20saccharinum"
# tr = ape::read.tree("data/sp.nwk")
# tr$edge.length

# you may not need this...

phy = insect_tree
phy2 = timetree_phylo(insect_tree)

cali = phy2$info %>% select(node, estimated, age.min, age.max) %>% mutate(soft.bounds = NA)
wi = which(is.na(cali$age.min))
cali$age.min[wi] = cali$age.max[wi] = cali$estimated[wi]
cali = select(cali, -estimated)
cali = select(cali, -node) %>% tibble::rowid_to_column() %>% rename(node = rowid)
phy3 = compute.brlen(phy2$phy)
plot(phy3); tiplabels(); nodelabels()
phy3$edge.length
phy4 = ape::chronos(phy3, calibration = makeChronosCalib(phy2$phy), model = "relaxed")
plot(phy4); tiplabels(); nodelabels()
phy4$edge.length

agez <- select(phy2$info, node, estimated) %>% na.omit()
write.table(x = agez, file = "data/phylogeny/insect_ages", quote = F, row.names = F, col.names = F)
write.tree(phy2$phy, file = "data/phylogeny/insect_test.tre")

# if you use R package phylocomr:
phylocomr::ph_bladj(ages = agez, phylo = phy2$phy) 




phy2 = datelife::EstimateDates(input = phy$tip.label, output.format = "phylo.median")

