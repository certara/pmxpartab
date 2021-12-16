#' @description  reads all relevant inputs (ie, parameter estimates, bootstrap, shrinkage, .) using directly model run or with the help of a YAML file. 
#' The parameter name must be provided in the model file. Other fields are optional and can be provided either in (1) a .mod file or a (2) YAML file.
#' Few fields needs to be inputted directly from user to complete the table:
#'   name (compulsory): name of parameter (eg, CL, Vd, nCL, nVd, corr(nCL, nVd), etc.)
#'   label:  description of parameter (eg, "Volume of distribution", "Apparent Clearance", etc.) 
#'   units: parameter unit (eg, L/h, 1/hr, etc)
#'   trans: parameter transformation (ie, %, exp, ilogit, sqrt, CV%)
#'   type: (ie, Typical Values, Between Subject Variability, Inter Occasion Variability, Residual Error, Covariates)
#' 
#' @param run_dir NONMEM model directory
#' @param run_prefix NONMEM run number prefix
#' @param runno NONMEM run number  
#' @param bootstrap flag for availability of bootstrap results. 
#' @param run_dir.boot bootstrap results directory
#' @param runno.boot bootstrap NONMEM run number (if changes from runno)
#' @param conf.level bootstrap results confidence interval. By default set to 95%.
#' @param min_succ filter bootstrap results on minimization successful. By default set to TRUE.
#' @param read.boot flag for reading directly the customized filtering bootstrap raw_results.csv
#' @param boot.obj customized filtering bootstrap raw_results.csv data frame passed by user 
#' @param yaml.file flag for using the yaml.file option. All the non-compulsory fields (ie, label, trans, units, type) will be read from yaml file. 
#' @param yaml.file.name name of the yaml file.
#' @return meta file information (df_m) and parameter information (prm)

parframe2setup <- function(run_dir, run_prefix, runno, bootstrap = NULL, run_dir.boot = NULL, runno.boot = NULL, conf.level = 0.95, min_suc = TRUE, read.boot = NULL, boot.obj = NULL, yaml.file = NULL, yaml.file.name = NULL) {  
  
  #Load xpose database
  xpdb   <- xpose_data(prefix = run_prefix, runno = runno, dir = run_dir)
  
  # set-up bootstrap info if available
  have.bootstrap <- !is.null(bootstrap)
  read.bootstrap <- !is.null(read.boot)
  have.boot.obj  <- !is.null(boot.obj) 										 
  if (have.bootstrap & !(read.bootstrap)) {
    boot_dir <- paste0(run_dir, run_dir.boot)  
    boot_res <- sprintf("/raw_results_%s%s.csv", run_prefix, ifelse(is.null(runno.boot),runno,runno.boot))
    boot   <- read.csv(paste0(boot_dir, boot_res), header = TRUE, check.names = FALSE)
  } else if (have.bootstrap & read.bootstrap & have.boot.obj) {
    boot <- boot.obj
  }
  
  # extract par info using xpose 
  yaml.file <- is.null(yaml.file)
  if (yaml.file){
    
    prm <- get_prm(xpdb, transform = FALSE) %>%
      rowwise() %>%
      rename(name.xpose = name) %>%
      mutate(meta = ifelse((label!=""),map(label,parse_parameter_description),list(list(name = '')))) 
    
    vec=prm %>% dplyr::select(m,n,diagonal,meta)
    off.diag = which(prm$diagonal == FALSE)
    for (i in 1:length(off.diag)) {
      
      prm$meta[off.diag[i]] <- list(list(name = paste0(vec$meta[[which(vec$m==prm$m[off.diag[i]] & vec$diagonal == TRUE)[1]]]$name,
                                                       ',',
                                                       vec$meta[[which(vec$n==prm$n[off.diag[i]] & vec$diagonal == TRUE)[1]]]$name),
                                         label = paste0('Correlation ',vec$meta[[which(vec$m==prm$m[off.diag[i]] & vec$diagonal == TRUE)[1]]]$label,
                                                        ',',
                                                        vec$meta[[which(vec$n==prm$n[off.diag[i]] & vec$diagonal == TRUE)[1]]]$label),
                                         type = 'IIV'))
      
    }
    
    prm$name = unlist(lapply(prm$meta, function(x) x[c('name')]))
    
    df_m = do.call(bind_rows, list(parmaters=prm$meta)) %>% as_tibble() %>% mutate(type = factor(type, levels = c('Structural','CovariateEffect','IIV', 'IOV','RUV'))) %>% arrange(type) %>% mutate(type = as.character(type))
    
  } else { #yaml file
    
    prm <- get_prm(xpdb, transform = FALSE) %>%
      rowwise() %>%
      rename(name.xpose = name) %>%
      mutate(meta = ifelse((label!=""),map(label,parse_parameter_description),list(list(name = '')))) 
    
    vec=prm %>% dplyr::select(m,n,diagonal,meta)
    off.diag = which(prm$diagonal == FALSE & prm$label=='')
    for (i in 1:length(off.diag)) {
      
      prm$meta[off.diag[i]] <- list(list(name = paste0(vec$meta[[which(vec$m==prm$m[off.diag[i]] & vec$diagonal == TRUE)[1]]]$name,
                                                       ',',
                                                       vec$meta[[which(vec$n==prm$n[off.diag[i]] & vec$diagonal == TRUE)[1]]]$name),
                                         label = paste0('Correlation ',vec$meta[[which(vec$m==prm$m[off.diag[i]] & vec$diagonal == TRUE)[1]]]$label,
                                                        ',',
                                                        vec$meta[[which(vec$n==prm$n[off.diag[i]] & vec$diagonal == TRUE)[1]]]$label),
                                         type = 'IIV'))
      
    }
    
    prm$name = unlist(lapply(prm$meta, function(x) x[c('name')])) 
    
    meta <- read_yaml(file = yaml.file.name)
    list(meta)
    
    meta$parameters
    do.call(bind_rows, meta$parameters) %>% as_tibble() -> tmp
    
    if(is.null(tmp$trans)) {tmp$trans = NA}
    if(is.null(tmp$units)) {tmp$units = NA}
    if(is.null(tmp$label)) {tmp$label = NA}
    if(is.null(tmp$type)) {tmp$type = NA}  										   
    df_m = tmp %>% select(name, label, units, trans, type)
    
  }

  # merge shrinkage information

  if (nrow(prm %>% filter(type %in% c("ome","sig") & !fixed))>0) {
    etashk = data.frame(type = "ome",
                        shk = strsplit(xpdb$summary %>% filter(label=="etashk") %>% pull(value), split=", ")[[1]]) %>%
      mutate(m = as.numeric(gsub(".*\\[(\\d+)\\].*", "\\1", shk)),
             n = m, # added to ensure that shrinkage only joins diagonal omega
             shk = as.numeric(gsub("\\[.*", "", shk)))

    epsshk = data.frame(type = "sig",
                        shk = strsplit(xpdb$summary %>% filter(label=="epsshk") %>% pull(value), split=", ")[[1]]) %>%
      mutate(m = as.numeric(gsub(".*\\[(\\d+)\\].*", "\\1", shk)),
             n = m, # added to ensure that shrinkage only joins diagonal omega
             shk = as.numeric(gsub("\\[.*", "", shk)))
    shk.tmp = rbind(etashk,epsshk)
    prm = prm %>%
      left_join(shk.tmp, by=c("type", "m", "n")) %>% rename(shrinkage = shk)
    if (prm$value[str_detect(prm$name.xpose,'SIGMA')]==1){ # error coded with fixed effects
      
      prm$shrinkage[str_detect(prm$label,'ERR') | str_detect(prm$label,'PROP') | str_detect(prm$label,'ADD')]=prm$shrinkage[str_detect(prm$name.xpose,'SIGMA')]
      
    }

  }
  
  # merge bootstrap information 
  if (have.bootstrap) {
    
    if (min_suc){
      boot = boot %>% filter(minimization_successful!=0) # default remove the unsuccessful runs    
    } 
    boot = boot[,(which(names(boot)=='ofv')+1):(which(regexpr("^se", names(boot), perl=T)==1)[1]-1)]  # select the columns of interest
    boot = boot %>%
      gather(key = "label_boot") %>% # long format
      group_by(label_boot) %>%
      mutate(bootstrap.median = median(value,na.rm = TRUE),
             bootstrap.uci    = quantile(value,probs=c((1+conf.level)/2), na.rm = TRUE),
             bootstrap.lci    = quantile(value,probs=c((1-conf.level)/2), na.rm = TRUE)) %>%
      ungroup() %>% 
      distinct(label_boot, .keep_all = TRUE) %>%
      dplyr::select(-value)
    
    prm = cbind(prm,boot)     
    
  }
  
  list(prm,df_m)
  
}

#' @description Get all relevant parameter information in a data.frame, input the desired transformation for parameters 
								 
#' and finally merge parameter estimates with their corresponding meta information  
#' 
#' @param out paramter data frame from parframe2setup
#' @param meta meta information data frame from parframe2setup
#' @param bootstrap flag for availability of bootstrap results
#' @param conf.level confidence interval of parameter estimates using normal distribution assumptions (in case no bootstrap rasults are available). By default set to 95%. 
parframe <- function(out, meta, bootstrap = NULL, conf.level = 0.95) {
  z <- meta
  
  z$fixed     <- as.logical(NA)
  z$value     <- as.numeric(NA)
  z$se        <- as.numeric(NA)
  z$rse       <- as.numeric(NA)
  z$lci       <- as.numeric(NA)
  z$uci       <- as.numeric(NA)
  z$pval      <- as.numeric(NA)
  z$shrinkage <- as.numeric(NA)
  
  have.bootstrap <- !is.null(bootstrap)
  if (have.bootstrap) {
    z$boot.median <- as.numeric(NA)
    z$boot.lci    <- as.numeric(NA)
    z$boot.uci    <- as.numeric(NA)
  }
  
  `%||%` <- function(x, y) { if (is.null(x) || is.na(x)) y else x }
  for (i in 1:nrow(z)) {
    name <- z$name[i] %||% NA
    trans <- z$trans[i] %||% NA
    
    if (have.bootstrap) {
      boot.median <- NA
      boot.lci <- NA
      boot.uci <- NA
    }
    
    # Check parameter type
    j <- which(out$name == name)
    if (length(j) == 1) {
      value <- as.numeric(out$value[j])
      se <- as.numeric(out$se[j])
      fixed <- as.logical(out$fixed[j])
      shrinkage <- if ("shrinkage" %in% names(out)) as.numeric(out$shrinkage[j]) else NA
      
      # !! The boostratp output may not be uniquely identified by a "name", but rather raw nonmem nm_names

      if (have.bootstrap) {
        boot.median <- out$bootstrap.median[j]
        boot.lci <- out$bootstrap.lci[j]
        boot.uci <- out$bootstrap.uci[j]
      }
    } else { #no parameters match
      
      value <- NA
      se <- NA
      fixed <- FALSE
      shrinkage <- NA
      
      #!! The boostratp output may not be uniquely identified by a "name", but rather raw nonmem nm_names
      if (have.bootstrap) {
        boot.median <- NA
        boot.lci <- NA
        boot.uci <- NA
      }
    }
    
    if (fixed || is.null(se) || is.na(se)) {
      se <- NA
      rse <- NA
      ci <- c(NA, NA)
      pval <- NA
    } else {
      rse <- 100*se/abs(value)
      ci <- value + c(-1,1) * qnorm((1+conf.level)/2) * se
      pval <- 2*(1 - pnorm(abs(value/se)))
    }
    
    # Check transformation
    if (!is.na(trans) && trans == "%") {
      value <- 100*value
      if (!fixed) {
        se  <- 100*se
        ci  <- 100*ci
      }
      if (have.bootstrap) {
        boot.median <- 100*boot.median
        boot.lci <- 100*boot.lci
        boot.uci <- 100*boot.uci
      }
    } else if (!is.na(trans) && trans == "exp") {
      value <- exp(value)
      if (!fixed) {
        rse <- 100*se
        se  <- (rse/100)*value
        ci  <- exp(ci)
      }
      if (have.bootstrap) {
        boot.median <- exp(boot.median)
        boot.lci <- exp(boot.lci)
        boot.uci <- exp(boot.uci)
      }
    } else if (!is.na(trans) && trans == "ilogit") {
      ilogit <- function(x) { 1 / (1 + exp(-x)) }
      value <- ilogit(value)
      if (!fixed) {
        rse <- 100*se*(1 - value)
        se  <- (rse/100)*value
        ci  <- ilogit(ci)
      }
      if (have.bootstrap) {
        boot.median <- ilogit(boot.median)
        boot.lci <- ilogit(boot.lci)
        boot.uci <- ilogit(boot.uci)
      }
    } else if (!is.na(trans) && trans == "sqrt") {
      g <- function(x) { sqrt(x) }
      dg <- function(x) { 1/(2*sqrt(x)) }
      x <- value
      value <- g(x)
      if (!fixed) {
        se  <- se*dg(x)
        rse <- 100*se/abs(value)
        ci  <- g(ci)
      }
      if (have.bootstrap) {
        boot.median <- g(boot.median)
        boot.lci <- g(boot.lci)
        boot.uci <- g(boot.uci)
      }
    } else if (!is.na(trans) && trans == "CV%") {
      g <- function(x) { 100 * sqrt(exp(x) - 1) }

      x <- value
      value <- g(x)
      if (!fixed) {
        dg <- function(x) { 100 * exp(x)/(2*sqrt(exp(x)-1)) }															 
        se  <- se*dg(x)
        rse <- 100*se/abs(value)
        ci  <- g(ci)
      }
      if (have.bootstrap) {
        boot.median <- g(boot.median)
        boot.lci <- g(boot.lci)
        boot.uci <- g(boot.uci)
      }
    } else if (!is.na(trans) && trans == "CV2%") {
      g <- function(x) { 100*sqrt(exp(x^2) - 1) }

      x <- value
      value <- g(x)
      if (!fixed) {
        dg <- function(x) { 100*0.5*(1/sqrt(exp(x^2) - 1))*exp(x^2)*2*x }																		 
        se  <- se*dg(x)
        rse <- 100*se/abs(value)
        ci  <- g(ci)
      }
      if (have.bootstrap) {
        boot.median <- g(boot.median)
        boot.lci <- g(boot.lci)
        boot.uci <- g(boot.uci)
      }
    }
    
    z$fixed[i] <- fixed
    z$value[i] <- value
    z$se[i]    <- se
    z$rse[i]   <- rse
    z$lci[i]   <- ci[1]
    z$uci[i]   <- ci[2]
    z$pval[i]  <- pval
    z$shrinkage[i] <- shrinkage
    
    if (have.bootstrap) {
      z$boot.median[i] <- boot.median
      z$boot.lci[i] <- boot.lci
      z$boot.uci[i] <- boot.uci
    }
  }
  z <- subset(z, !is.na(value) & name!='')
  as.data.frame(z)
}

#' @description Uses the output from parframe to generate a parameter estimates table in HTML 
#' 
#' @param parframe object from parframe
#' @param meta meta information data frame from parframe2setup
#' @param columns which information to display in the table
#' @param sections flag for having separated sections displayed
#' @param section.labels which sections will be displayed
#' @param show.fixed.to.zero flag to show or not parameter fixed to zero. By default will be not displayed.
#' @param na how NA will be displayed
#' @param digits how many digits will be displayed

pmxpartab <- function(
  parframe, meta, # added meta as argument
  
  columns=c(value="Estimate", rse="RSE%", ci95="95%CI", shrinkage="Shrinkage"), # changed est to value
  
  sections = TRUE,
  section.labels = c(
    Structural      = "Typical Values",
    CovariateEffect = "Covariate Effects",
    RUV             = "Residual Error",
    IIV             = "Between Subject Variability",
    IOV             = "Inter-Occasion Variability"),
  
  show.fixed.to.zero=F,
  na="n/a",
  digits=3) {
  
  if (isFALSE(show.fixed.to.zero)) {
    parframe <- subset(parframe, !(fixed & value==0)) #changed est to value
  }
  
  ncolumns <- length(columns) + 1
  
  thead <- paste0('<tr>\n<th>Parameter</th>\n',
                  paste0(paste0('<th>', columns, '</th>'), collapse="\n"), '\n</tr>')
  
  
  tbody <- ""
  for (i in 1:nrow(parframe)) {
    if (isTRUE(sections)) {
      newsection <- (!is.null(parframe$type) && !is.na(parframe$type[i]) && (i == 1 || parframe$type[i] != parframe$type[i-1]))
      if (newsection) {
        type <- parframe$type[i]
        if (type %in% names(meta$labels)) {
          label <- meta$labels[[type]]
        } else if (type %in% names(section.labels)) {
          label <- section.labels[[type]]
        } else {
          label <- type
        }
        
        tbody <- paste0(tbody, parameter.estimate.table.section(label, ncolumns=ncolumns), '\n')
      }
    }
    args <- c(parframe[i,], list(na=na, digits=digits, indent=sections))
    tbody <- paste0(tbody, do.call(parameter.estimate.table.row, args), '\n')
  }
  
  table <- paste0('<table>\n<thead>\n', thead, '\n</thead>\n<tbody>\n', tbody, '\n</tbody>\n</table>\n')
  structure(table, class=c("pmxpartab", "html", "character"), html=TRUE)
}

# Internal function to help format numbers
p <- function(x, digits=3, flag="", round.integers=FALSE){
  if (!is.numeric(x)) {
    return(x)
  }
  prefix <- ifelse(flag=="+" & x > 0, "+", "")
  paste0(prefix, table1::signif_pad(x, digits=digits, round.integers=round.integers))
}

parameter.estimate.table.section <- function(label, ncolumns) {
  paste0(c('<tr>',
           paste0(sprintf('<td class="paramsectionheading">%s</td>', c(label, rep("", ncolumns-1))), collapse='\n'),
           '</tr>'), collapse='\n')
}

parse_parameter_description <- function(string) {
  # Returns a structured object representing a description of a parameter
  # (in this example just a list with some attributes; only the name is mandatory)
  parameter_description <- function(name, label=NULL, units=NULL, trans=NULL, type=NULL) {
    list(name=name, label=label, units=units, trans=trans, type=type)
  }
  
  x <- str2lang(paste0("parameter_description(", string, ")"))
  x[[2]] <- as.character(x[[2]]) # Interpret the first element (name) as a string even if not quoted
  eval(x)
}
parameter.estimate.table.row <- function(
  name,
  label          = NULL,
  units          = NULL,
  type           = c("Structural", "CovariateEffect", "IIV", "IOV", "RUV", "Unspecified"),
  trans          = c("identity", "%", "exp", "ilogit", "CV%", "SD (CV%)"),
  expression     = NULL,
  relatedTo      = NULL,
  superscript    = NULL,
  fixed          = NULL,
  value          = NULL, # changed from est to value
  se             = NULL,
  rse            = NULL,
  lci95          = NULL,
  uci95          = NULL,
  boot.median    = NULL,
  boot.lci       = NULL,
  boot.uci       = NULL,
  shrinkage      = NULL,
  na             = "n/a",
  digits         = 3,
  indent         = TRUE,
  have.bootstrap = !is.null(boot.median),
  columns=c(value="Estimate", rse="RSE%", ci95="95%CI", shrinkage="Shrinkage"), # changed est to value
  ...) {
  
  # Check for superscript
  if (is.null(superscript) || is.na(superscript)) {
    superscript <- ""
  } else {
    superscript <- paste0("<sup>", superscript, "</sup>")
  }
  
  # Check for label
  if (is.null(label) || is.na(label)) {
    label <- name
  }
  
  # Check for units
  if (!is.null(units) && !is.na(units)) {
    label <- sprintf("%s (%s)", label, units)
  }
  
  # changed est to value --- 
  if (!is.null(trans) && !is.na(trans) && trans == "SD (CV%)") {
    g <- function(x) { 100*sqrt(exp(x^2) - 1) }
    x <- value
    value <- sprintf("%s (%s%%)", p(x, digits), p(g(x), digits))
  } else {
    value <- p(value, digits)
  }
  
  value <- paste0(value, superscript)
  
  if (fixed) {
    value <- sprintf('%s Fixed', value)
  }
  # up to here
  
  if (is.na(se)) {
    se <- na
    rse <- na
    ci95 <- na
  } else {
    rse <- p(rse, digits)
    ci95 <- sprintf('%s &ndash; %s', p(lci95, digits), p(uci95, digits))
  }
  
  if (have.bootstrap) {
    if (is.na(boot.median)) {
      boot.median <- na
      boot.ci95 <- na
    } else {
      boot.median <- p(boot.median, digits)
      boot.ci95 <- sprintf('%s &ndash; %s', p(boot.lci, digits), p(boot.uci, digits))
    }
  } else {
    boot.ci95 <- NULL
  }
  
  if (!is.null(shrinkage)) {
    if (is.na(shrinkage)) {
      shrinkage <- ""
    } else {
      shrinkage <- sprintf("%s%%", p(shrinkage, digits))
    }
  }
  
  all <- c(value=value, rse=rse, ci95=ifelse(have.bootstrap,boot.ci95,ci95), shrinkage=shrinkage) #changed est to value
  paste0(c('<tr>',
           sprintf('<td class="%s">%s</td>', ifelse(isTRUE(indent), "paramlabelindent", "paramlabelnoindent"), label),
           paste0(sprintf('<td>%s</td>', all[names(columns)]), collapse='\n'),
           # paste0(sprintf('<td>%s</td>', all[names(all)]), collapse='\n'),
           '</tr>'), collapse='\n')
}

