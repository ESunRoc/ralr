##### Required packages #####


#' robust additive log-ratio 
#' 
#' Computes the robust additive log-ratio (ralr) transformation for longitudinal/grouped compositional data. The returned data are free of a compositional effect and follow a MVN distribution.
#' 
#' @param dat Data frame; the count-level data to be normalized with columns (at least):
#' \itemize{
#' \item `pid`: patient/subject ID
#' \item `visit`: the visit or observation number
#' \item One column per feature/taxon
#' }
#' @param alpha Numeric; the confidence bandwidth parameter taking values in the interval \eqn{[0,1]}. Defaults to 0.05.
#' @param pseudo Numeric; a pseudocount that replaces zero feature values so that every ratio and log is well-defined. Defaults to 0.5.
#' @param min_prev Numeric; a minimum prevalence in the interval \eqn{[0,1)} required for a feature to enter the reference component. Features can 
#'                          be selected for the reference only if they are non-zero in at least `100*min_prev`\% of samples at \emph{both} visits.
#'                          Defaults to 0.1, with values closer to 0 allowing for rarer taxa to enter the reference.
#'                          
#' @returns A list with:
#' \itemize{
#' \item `ralr`: log-ratios of every post-baseline visit with respect to baseline and a CLR-transformed baseline value) normalized data at the subject and visit level
#' \item `refs`: the names of the features selected for each pair of visits, named "j:k" by internal visit
#' \item `ref_info`: per pair information on visits, sample sizes, candidate taxa, the selected angle, and the number of features in the reference set
#' \item `dropped`: the samples that were dropped for having no reads
#' } 
#' 
#' @importFrom tidyselect all_of ends_with where
#' @importFrom dplyr %>% across bind_rows desc distinct everything full_join group_by group_split pick rename_with vars
#' @importFrom utils combn
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom stats qnorm xtabs
#' @importFrom purrr reduce
#' @importFrom stringr str_detect str_replace
#' 
#' @export

ralr <- function(dat, alpha=0.05, pseudo=0.5, min_prev=0.1){
  # Based on initial coding of M. Sohn with modifications by E. Sundplyr
  ##### Feature columns, zero-read samples, and zeros
  feat.cols <- setdiff(names(dat), c("pid", "visit"))
  lib <- rowSums(dat[, feat.cols])
  dropped <- dat[lib == 0, c("pid", "visit")]
  dat <- dat[lib > 0, , drop = FALSE]
  observed <- as.matrix(dat[, feat.cols]) > 0                 # non-zero before the pseudo-count
  dat[, feat.cols][dat[, feat.cols] == 0] <- pseudo
  
  ##### Visits into 1, 2, ...
  dat$internal_visit <- as.numeric(as.factor(dat$visit))
  ##### PID and Visits Info
  vsts <- dat$internal_visit %>% unique() %>% sort()
  pid_visit_xtabs <- xtabs(~pid+internal_visit, data=dat)
  pid_one_visit <- names(which(rowSums(pid_visit_xtabs) == 1))
  ##### Remove pid with only one visit
  dat.no1v <- dat %>% dplyr::filter(pid %in% names(which(rowSums(pid_visit_xtabs) > 1)))
  
  pid_visit_df <- xtabs(~pid+internal_visit, data=dat.no1v) %>% as.data.frame()
  tmp.ref.vsts <- list()
  for(i in 1:length(vsts)){
    tmp.ref.vsts[[i]] <- with(pid_visit_df, pid[internal_visit==vsts[i] & Freq == 1]) %>% as.character()
  }
  ref.vsts <- list()
  ref.vsts[[1]] <- tmp.ref.vsts[[1]]
  ref.vsts.accum <- ref.vsts[[1]]
  for(i in 2:length(vsts)){
    ref.vsts[[i]] <- tmp.ref.vsts[[i]][!(tmp.ref.vsts[[i]] %in% ref.vsts.accum)]
    ref.vsts.accum <- c(ref.vsts.accum, ref.vsts[[i]])
  }
  
  pid_visit_mat_df <- xtabs(~pid+internal_visit, data=dat.no1v) %>% as.data.frame.matrix() %>%
    tibble::rownames_to_column(var="pid")
  
  ##### Normalization: Ratio wrt the earliest visit for pid
  n.visits <- length(unique(dat.no1v$internal_visit))
  comp.pairs <- combn(n.visits, 2) %>% t() %>% as.data.frame() %>%
    dplyr::mutate(V3 = V2-V1) %>% arrange(V3) %>% dplyr::select(-V3) %>%
    t() %>% as.matrix()
  n.pairs <- ncol(comp.pairs)
  
  tmp.ref.vsts <- NULL
  for(i in 1:n.pairs){
    tmp.ref.vsts <- cbind(tmp.ref.vsts, pid_visit_mat_df[,comp.pairs[,i]+1] %>% apply(1, prod))
  }
  colnames(tmp.ref.vsts) <- apply(comp.pairs, 2, function(x) paste(x[1], x[2], sep=":"))
  
  visit.pairs <- apply(comp.pairs, 2, function(x) vsts[x])
  tailing.ratio.id <- apply(visit.pairs, 2, function(x) paste0("_R", x[2], x[1]))
  
  ##### Get a set of non-DA features for each pair of visits
  ## All samples observed at a visit are used, including subjects seen at only one visit, so the
  ## two groups can differ in size. find_nda_pairwise takes the first n1 rows as the earlier visit
  ## and the next n2 rows as the later one, so the rows are ordered by visit first.
  ref.sets <- list(); ref.info <- list()
  for(i in 1:n.pairs){
    v1 <- visit.pairs[1,i]; v2 <- visit.pairs[2,i]
    in.pair <- which(dat$internal_visit %in% c(v1, v2))
    ord <- in.pair[order(dat$internal_visit[in.pair] != v1)]   # visit v1 rows first, then v2
    n1 <- sum(dat$internal_visit[ord] == v1); n2 <- sum(dat$internal_visit[ord] == v2)
    
    feat <- as.matrix(dat[ord, feat.cols])
    M <- feat / rowSums(feat)                                  # relative abundances
    obs <- observed[ord, , drop = FALSE]
    cand <- which(colMeans(obs[1:n1, , drop = FALSE]) >= min_prev &
                    colMeans(obs[(n1+1):(n1+n2), , drop = FALSE]) >= min_prev)
    if(length(cand) == 0) stop(sprintf("No features meet min_prev = %s at both visits %s and %s", min_prev, vsts[v1], vsts[v2]))
    
    rslt.tmp.ref <- find_nda_pairwise(M[, cand, drop = FALSE], n1 = n1, n2 = n2, alpha = alpha)
    ref.sets[[i]] <- feat.cols[cand][rslt.tmp.ref$id.nda]
    if(length(ref.sets[[i]]) == 0) stop(sprintf("Empty reference for visits %s and %s", vsts[v1], vsts[v2]))
    ref.info[[i]] <- data.frame(pair = paste(v1, v2, sep = ":"), visit1 = vsts[v1], visit2 = vsts[v2],
                                n1 = n1, n2 = n2, candidates = length(cand),
                                theta = rslt.tmp.ref$theta_selected, ref_size = length(ref.sets[[i]]))
  }
  names(ref.sets) <- apply(comp.pairs, 2, function(x) paste(x[1], x[2], sep = ":"))
  
  ##### Get each pair of visits using each set of non-DA 
  dat.vst <- dat.no1v %>% group_by(internal_visit) %>% group_split()
  paired.ratios <- list()
  for(i in 1:n.pairs){
    pw.ref.j <- dat.vst[[visit.pairs[1,i]]][,ref.sets[[i]]] %>% rowSums()
    pw.ref.k <- dat.vst[[visit.pairs[2,i]]][,ref.sets[[i]]] %>% rowSums()
    
    num.ratios <- dat.vst[[visit.pairs[1,i]]] %>%
      dplyr::select(-visit) %>%
      dplyr::mutate_at(vars(-pid, -internal_visit), ~./pw.ref.j) %>%
      dplyr::mutate(internal_visit = 1)
    
    den.ratios <- dat.vst[[visit.pairs[2,i]]] %>%
      dplyr::select(-visit) %>%
      dplyr::mutate_at(vars(-pid, -internal_visit), ~./pw.ref.k) %>%
      dplyr::mutate(internal_visit = 2)
    
    paired.ratios[[i]] <- bind_rows(num.ratios, den.ratios) %>% 
      pivot_wider(names_from = internal_visit, values_from = -c(pid, internal_visit)) %>%
      dplyr::mutate(across(ends_with("_2"), ~ .x / get(str_replace(cur_column(), "_2$", "_1")),
                           .names = "{col}_RT")) %>%
      dplyr::select(pid, ends_with("_RT")) %>%
      rename_with(~ str_replace(., "_2_RT", tailing.ratio.id[i]))
  }
  names(paired.ratios) <- apply(visit.pairs, 2, function(x) paste0("R", x[2], x[1]))
  
  missing_visits <- distinct(pid_visit_mat_df[,-1]) %>% arrange(desc(across(everything())))
  
  n.missing_visits <- nrow(missing_visits)
  full.ratios.pid <- list()
  for(i in 1:n.missing_visits){
    full.ratios.pid[[i]] <- pid_visit_mat_df$pid[apply(pid_visit_mat_df[,-1], 1, 
                                                       function(x) identical(as.numeric(x), 
                                                                             as.numeric(missing_visits[i,])))]
  }
  
  ##### Normalize each pair of visits using each set of non-DA 
  ratio.wrt.bl.l <- NULL
  for(i in 1:n.missing_visits){
    tmp.ratio.wrt.bl.l <- NULL
    tmp.ratio.wrt.bl.ls <- list()
    loc.bl <- which(missing_visits[i,]==1)
    target.ratios <- paste0("R", loc.bl[-1], loc.bl[-length(loc.bl)])
    n.target.ratios <- length(target.ratios)
    tmp.comp.pairs <- tmp.ratio.wrt.bl.ls[[1]] <- 
      paired.ratios[[target.ratios[1]]] %>% dplyr::filter(pid %in% full.ratios.pid[[i]])
    tailing.ratio.id.wrt.bl <- paste0("_R", loc.bl[2], loc.bl[1])
    if(n.target.ratios > 1){
      for(j in 2:n.target.ratios){
        tmp.paired.ratios <- paired.ratios[[target.ratios[j]]] %>% 
          dplyr::filter(pid %in% full.ratios.pid[[i]])
        tmp.comp.pairs <- tmp.ratio.wrt.bl.ls[[j]] <- 
          full_join(tmp.comp.pairs, tmp.paired.ratios, by="pid") %>%
          dplyr::mutate(across(ends_with(tailing.ratio.id.wrt.bl), 
                               ~ .x * get(str_replace(cur_column(),
                                                      paste0(tailing.ratio.id.wrt.bl, "$"),
                                                      paste0("_", target.ratios[j]))),
                               .names = "{col}_RT")) %>%
          dplyr::select(pid, ends_with("_RT")) %>%
          rename_with(~ str_replace(., paste0(tailing.ratio.id.wrt.bl, "_RT"), 
                                    paste0("_R", loc.bl[j+1], loc.bl[1])))
        tailing.ratio.id.wrt.bl <- paste0("_R", loc.bl[j+1], loc.bl[1])
      }
      tmp.ratio.wrt.bl.l <- tmp.ratio.wrt.bl.ls %>% reduce(full_join, by="pid") %>%
        pivot_longer(!pid,
                     cols_vary="slowest",
                     names_to=c(".value", "visit"),
                     names_pattern="(.+)_(.+)")
    } else{
      tmp.ratio.wrt.bl.l <- tmp.ratio.wrt.bl.ls[[1]] %>%
        pivot_longer(!pid,
                     cols_vary="slowest",
                     names_to=c(".value", "visit"),
                     names_pattern="(.+)_(.+)")
    }
    
    ratio.wrt.bl.l <- bind_rows(ratio.wrt.bl.l, tmp.ratio.wrt.bl.l)
  }
  
  ##### Only consider samples having baseline measures (visit 1)
  fratios <- ratio.wrt.bl.l %>% dplyr::filter(endsWith(visit, "1")) %>%
    dplyr::mutate(across(where(is.numeric), log),
                  visit = as.factor(visit))
  colnames(fratios)[str_detect(colnames(fratios), "^X\\d")] <-
    str_replace(colnames(fratios)[str_detect(colnames(fratios), "^X\\d")], "^X", "")
  
  ##### CLR transformed baseline measures
  fratios.bl.gm <- dat.vst[[1]] %>% dplyr::select(-pid, -visit, -internal_visit) %>%
    as.data.frame() %>%
    dplyr::mutate(log_gmean = rowMeans(log(pick(everything())))) %>%
    dplyr::mutate(across(-log_gmean, ~ log(.x) - log_gmean)) %>%
    dplyr::select(-log_gmean) %>%
    dplyr::mutate(pid=dat.vst[[1]]$pid, visit="R11") %>%
    dplyr::select(pid, visit, all_of(colnames(fratios)[-c(1:2)]))
  
  ##### Ratios from clr baseline to each post visit
  fratios <- bind_rows(fratios, fratios.bl.gm) %>%
    arrange(pid, visit)
  return(list(ralr=fratios, refs=ref.sets, ref_info=do.call(rbind, ref.info), dropped=dropped))
}
