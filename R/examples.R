# Working examples

library(tidyverse)
library(knitr)
library(xpose)


# Convert "list of list" to data.frame, similar to xpose output

ex_outputs <- list(
  th = list(CL = 0.482334, VC = 0.0592686, `Age on CL` = -0.234, `Age on VC` = 0.234),
  om = list(nCL = 0.315414, nVC = 0.536025),
  sg = list(ERRP = 0.0508497),
  se = list(
    th = list(CL = 0.0138646, VC = 0.00555121, `Age on CL` = -.0234, `Age on VC` = 0.0234),
    om = list(nCL = 0.0188891, nVC = 0.0900352),
    sg = list(ERRP = 0.00182851)),
  fixed = list(
    th = list(CL = FALSE, VC = FALSE,`Age on CL` = FALSE, `Age on VC` = FALSE),
    om = list(nCL = FALSE, nVC = FALSE),
    sg = list(ERRP = FALSE)),
  shrinkage = list(nCL = 9.54556, nVC = 47.8771))

df_o = data.frame(value = as.numeric(cbind(with(ex_outputs, c(th, om, sg)))),
                  se = as.numeric(cbind(with(ex_outputs$se, c(th, om, sg)))),
                  fixed = as.numeric(cbind(with(ex_outputs$fixed, c(th, om, sg)))))
df_o$name = names(with(ex_outputs, c(th, om, sg)))

df_s = data.frame(shrinkage = as.numeric(cbind(ex_outputs$shrinkage)))
df_s$name = names(ex_outputs$shrinkage)

df_o = df_o %>%
  left_join(df_s, by="name") %>%
  as_tibble()

df_o = bind_rows(df_o,
                 df_o %>%
                   filter(name=="nCL") %>%
                   mutate(name="nCL2",
                          value=sqrt(value)))


#Convert "meta" to data.frame /// if this is called something else than meta then it has to be changed in all places in the function 
ex_meta <- list(parameters = list(
  list(name="CL",   label="Clearance",    units="L/h", type="Structural"),
  list(name="VC",   label="Volume",       units="L",   trans="exp", type="Structural"),
  list(name="nCL",  label="On Clearance",              trans=NA, type="IIV"),
  list(name="nCL",  label="On Clearance",              trans="sqrt", type="IIV"),
  list(name="nCL",  label="On Clearance",              trans="CV%", type="IIV"),
  list(name="nCL2", label="On Clearance",              trans="CV2%", type="IIV"),
  list(name="nVC",  label="On Volume",                 type="IIV"),
  list(name="ERRP", label="Proportional Error",        units="-", trans="", type="RUV"),
  list(name="Age on CL", label="Age effect on CL",        units="%", trans="%", type="CovariateEffects"),
  list(name="Age on VC", label="Age effect on VC", type="CovariateEffects")
)
)

df_m = do.call(bind_rows, ex_meta) %>% as_tibble()


parframe(out=df_o, meta=df_m)
pmxpartab(parframe(out=df_o, meta=df_m), meta=df_m,
           columns=c(value="Estimate", rse="RSE%", ci95="95%CI", shrinkage="Shrinkage"))

# use xpose output with parframe

