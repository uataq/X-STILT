# subroutine to get vertical transport error, DW, 08/31/2018

# no random error covariance on PBL heights used for now
# but one can assign/prescribe values to '*zierr' as below

# 'nhrs.zisf': # of hours for zi scaling, always positive
# 'const.zisf': constant PBL scaling factors for all receptors, during 'nhrs.zisf'

get.zierr <- function(run_ver_err, nhrs.zisf = NULL, const.zisf = NULL,
                      sigzierr = NA, TLzierr = NA, horcorzierr = NA){

  ### Besides horizontal wind error, do we want to account for PBL height?
  if (run_ver_err) {    # add vertical trans error via ziscale *****
    zicontroltf <- 1                         # 0: NO scaling; 1: with scaling

    # create ziscale as list, with list dimension of # of receptors
    ziscale  <- list(rep(const.zisf, nhrs.zisf))  # if const ziscale for all recp
    cat(paste('+++ PBL scaling of', const.zisf, 'when generating trajec +++\n'))

  } else {
    cat('NO PBL scaling when generating trajec...\n')
    zicontroltf <- 0
    const.zisf  <- 1      # 1 means no scaling
    ziscale     <- const.zisf
  } # end if run_ver_err

  err.stat <- list(run_ver_err = run_ver_err, zicontroltf = zicontroltf, 
                   ziscale = ziscale, sigzierr = sigzierr, TLzierr = TLzierr, 
                   horcorzierr = horcorzierr)
  return(err.stat)
}
