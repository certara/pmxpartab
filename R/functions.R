#' Setup model info function (@author: certara-alargajolli to check)
#' 
#' 

parframe2setup <- function(run_dir, run_prefix, runno, bootstrap = NULL, read.bootstrap = NULL, boot.obj = NULL, run_dir.boot = NULL, runno.boot = NULL, conf.level = 0.95, min_suc = TRUE, yaml.file = NULL) {  #modified
  
  #Load xpose database
  xpdb   <- xpose_data(prefix = run_prefix, runno = runno, dir = run_dir)
  
  # set-up bootstrap info if available
  have.bootstrap <- !is.null(bootstrap)
  read.bootstrap <- !is.null(read.bootstrap)
  if (have.bootstrap & !(read.bootstrap)) {
    boot_dir <- paste0(run_dir, run_dir.boot)  
    boot_res <- sprintf("/raw_results_%s%s.csv", run_prefix, ifelse(is.null(runno.boot),runno,runno.boot))
    boot   <- read.csv(paste0(boot_dir, boot_res), header = TRUE, check.names = FALSE)
  } else if (have.bootstrap & read.bootstrap) {
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
                                         type = 'IIV'))
      
    }
    
    prm$name = unlist(lapply(prm$meta, function(x) x[c('name')]))
    
    df_m = do.call(bind_rows, list(parmaters=prm$meta)) %>% as_tibble()
    
  } else { #yaml file
    
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
                                         type = 'IIV'))
      
    }
    
    prm$name = unlist(lapply(prm$meta, function(x) x[c('name')])) 
    
    meta <- read_yaml(file = "meta.yaml")
    list(meta)
    
    meta$parameters
    do.call(bind_rows, meta$parameters) %>% as_tibble() -> tmp
    
    df_m = tmp %>% select(name, label, units, trans, type)
    
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

#'
#' Get parameters in a data.frame
#'
#' @export
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
      # @certara-alargajolli, please confirm from here ---
      if (have.bootstrap && name %in% names(out$bootstrap$median)) {
        boot.median <- out$bootstrap$median[[name]]
        boot.lci <- out$bootstrap$lci[[name]][1]
        boot.uci <- out$bootstrap$uci[[name]][2]
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
    # @certara-alargajolli, please check up to here ---
    
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
      dg <- function(x) { 100 * exp(x)/(2*sqrt(exp(x)-1)) }
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
    } else if (!is.na(trans) && trans == "CV2%") {
      g <- function(x) { 100*sqrt(exp(x^2) - 1) }
      dg <- function(x) { 100*0.5*(1/sqrt(exp(x^2) - 1))*exp(x^2)*2*x }
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
  z <- subset(z, !is.na(value))
  as.data.frame(z)
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
  
  all <- c(value=value, rse=rse, ci95=ci95, boot.median=boot.median, boot.ci95=boot.ci95, shrinkage=shrinkage) #changed est to value
  paste0(c('<tr>',
           sprintf('<td class="%s">%s</td>', ifelse(isTRUE(indent), "paramlabelindent", "paramlabelnoindent"), label),
           # paste0(sprintf('<td>%s</td>', all[names(columns)]), collapse='\n'),
           paste0(sprintf('<td>%s</td>', all[names(all)]), collapse='\n'),
           '</tr>'), collapse='\n')
}

#' Generate a parameter estimates table in HTML
#' 
#' @examples
#' outputs <- list(
#'     th = list(CL = 0.482334, VC = 0.0592686),
#'     om = list(nCL = 0.315414, nVC = 0.536025),
#'     sg = list(ERRP = 0.0508497),
#'     se = list(
#'         th = list(CL = 0.0138646, VC = 0.00555121),
#'         om = list(nCL = 0.0188891, nVC = 0.0900352),
#'         sg = list(ERRP = 0.00182851)),
#'     fixed = list(
#'         th = list(CL = FALSE, VC = FALSE),
#'         om = list(nCL = FALSE, nVC = FALSE),
#'         sg = list(ERRP = FALSE)),
#'     shrinkage = list(nCL = 9.54556, nVC = 47.8771))
#' 
#' meta <- list(
#'     parameters = list(
#'         list(name="CL",   label="Clearance",    units="L/h", type="Structural"),
#'         list(name="VC",   label="Volume",       units="L",   trans="exp", type="Structural"),
#'         list(name="nCL",  label="On Clearance",              trans="SD (CV%)", type="IIV"),
#'         list(name="nVC",  label="On Volume",                 type="IIV"),
#'         list(name="ERRP", label="Proportional Error",        units="%", trans="%", type="RUV")))
#' 
#' parframe(outputs, meta)
#' 
#' pmxpartab(parframe(outputs, meta),
#'     columns=c(est="Estimate", rse="RSE%", ci95="95%CI", shrinkage="Shrinkage"))
#' @export
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
  
  thead <- paste0('<tr>\n<th rowspan="2">Parameter</th>\n',
                  paste0(paste0('<th rowspan="2">', columns, '</th>'), collapse="\n"), '\n</tr>\n')
  
  
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


#' Print \code{pmxpartab} object.
#' @param x An object returned by \code{\link{pmxpartab}}.
#' @param ... Further arguments passed on to other \code{print} methods.
#' @return Returns \code{x} invisibly.
#' @details In an interactive context, the rendered table will be displayed in
#' a web browser. Otherwise, the HTML code will be printed as text.
#' @export
print.pmxpartab <- function(x, ...) {
  if (interactive()) {
    z <- htmltools::HTML(x)
    default.style <- htmltools::htmlDependency("pmxpartab", "1.0",
                                               src=system.file(package="pmxpartab", "pmxpartab_defaults_1.0"),
                                               stylesheet="pmxpartab_defaults.css")
    z <- htmltools::div(class="Rpmxpartab", default.style, z)
    z <- htmltools::browsable(z)
    print(z, ...) # Calls htmltools:::print.html(z, ...)
  } else {
    cat(x)
  }
  invisible(x)
}


#' Method for printing in a \code{knitr} context.
#' @param x An object returned by \code{\link{pmxpartab}}.
#' @param ... Further arguments passed on to \code{knitr::knit_print}.
#' @importFrom knitr knit_print
#' @export
knit_print.pmxpartab <- function(x, ...) {
  knit_to_html <-
    !is.null(knitr::opts_knit$get("rmarkdown.pandoc.to")) &&
    grepl("^html", knitr::opts_knit$get("rmarkdown.pandoc.to"))
  
  if (knit_to_html) {
    z <- htmltools::HTML(x)
    default.style <- htmltools::htmlDependency("pmxpartab", "1.0",
                                               src=system.file(package="pmxpartab", "pmxpartab_defaults_1.0"),
                                               stylesheet="pmxpartab_defaults.css")
    z <- htmltools::div(class="Rpmxpartab", default.style, z)
    knitr::knit_print(z, ...)
  } else {
    knitr::knit_print(as.character(x), ...)
  }
}

#' Parse a parameter description string
#'
#' @examples
#' # Example 1: all elements present
#' x <- "CL, label='Clearance', units='L/h', trans=exp, type='Structural'"
#' parse_parameter_description(x)
#' 
#' # Example 2: Some elements missing (trans), will take default value (NULL)
#' x <- "CL, label='Clearance', units='L/h', type='Structural'"
#' parse_parameter_description(x)
#' 
#' # Example 3: Only the name is given
#' x <- "CL"
#' parse_parameter_description(x)
#' 
#' # Example 4: positional arguments
#' x <- "CL, 'Clearance', 'L/h', type='Structural'"
#' parse_parameter_description(x)
#' @export
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

