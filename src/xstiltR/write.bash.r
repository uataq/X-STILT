# subroutine to write bash script for running configuration file
# to run this .sh file, "sbatch *.sh"
# by DW, 05/09/2018
# clean codes, DW, 06/12/2018

write.bash <- function(namelist, job.time = "24:00:00", num.nodes, num.cores,
                       account = "lin-kp", partition = "lin-kp",
                       email = "dien.wu@utah.edu",
                       email.type = "FAIL,BEGIN,END"){

   if (namelist$rmTF) {   # if delete previous one
     if (length(list.files(namelist$workdir, ".sh")) > 0) {
        file.remove(list.files(namelist$workdir, ".sh"))
     }
   }

  namelist$BASH <- paste0("XSTILT_", namelist$timestr, ".sh")
  bash.file <- file.path(namelist$workdir, namelist$BASH) # store at workdir
  run.file <- file.path(namelist$workdir, namelist$RUN)
  write("#!/bin/bash", bash.file)

  sbatch.setup <- paste0("\n#SBATCH --time=", job.time,
                        "\n#SBATCH --nodes=", num.nodes,
                        "\n#SBATCH --ntasks=", num.cores,
                        "\n#SBATCH --account=", account,
                        "\n#SBATCH --partition=", partition,
                        "\n#SBATCH --mail-user=", email,
                        "\n#SBATCH --mail-type=", email.type,"\n")
  write(sbatch.setup, file = bash.file, append = T)
  write("rm -f multiprog.conf\n", bash.file, append = T)

  # loop over each job and output in configuration file
  job.loop <- paste0("JOBNR=0\nSLURM_NTASKS=", num.cores,
                     "\nwhile [[ $JOBNR -lt $SLURM_NTASKS ]]; do\n",
                     "  echo $JOBNR ", run.file,
                     " $JOBNR >> multiprog.conf\n",
                     "  let JOBNR=JOBNR+1\ndone\n")
  write(job.loop, bash.file, append = T)

  # submit job
  write("srun --multi-prog multiprog.conf\n", bash.file, append = T)
  cat(paste("write.bash(): Done with bash file...\n"))

  # turn on excutable permission
  system(paste0("chmod u+x ", bash.file))

  return(bash.file)
}
