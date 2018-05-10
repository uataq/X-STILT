# subroutine to write bash script for running configuration file
# to run this .sh file, "sbatch *.sh"
# by DW, 05/09/2018

write.bash <- function(namelist, job.time="24:00:00", num.nodes, num.cores,
                       account="lin-kp", partition="lin-kp",
                       email="dien.wu@utah.edu",
                       email.type=c("FAIL,BEGIN,END")){

   # delete previous one
   if(length(list.files(namelist$workdir, ".sh"))>0){
    file.remove(list.files(namelist$workdir, ".sh"))
   }

  namelist$BASH <- paste("XSTILT_", namelist$timestr, ".sh", sep="")
  bash.file <- file.path(namelist$workdir, namelist$BASH) # store at workdir
  run.file <- file.path(namelist$workdir, namelist$RUN)

  write("#!/bin/bash", bash.file)

  sbatch.setup <- paste("\n#SBATCH --time=", job.time,
                        "\n#SBATCH --nodes=", num.nodes,
                        "\n#SBATCH --ntasks=", num.cores,
                        "\n#SBATCH --account=", account,
                        "\n#SBATCH --partition=", partition,
                        "\n#SBATCH --mail-user=", email,
                        "\n#SBATCH --mail-type=", email.type,"\n", sep="")

  write(sbatch.setup, file=bash.file, append=T)

  write("rm -f multiprog.conf\n", bash.file, append=T)

  # loop over each job and output in configuration file
  job.loop <- paste("JOBNR=0\nSLURM_NTASKS=", num.cores,
                    "\nwhile [[ $JOBNR -lt $SLURM_NTASKS ]]; do\n",
                    "  echo $JOBNR ", run.file,
                    " $JOBNR >> multiprog.conf\n",
                    "  let JOBNR=JOBNR+1\ndone\n", sep="")
  write(job.loop, bash.file, append=T)

  # submit job
  write("srun --multi-prog multiprog.conf\n", bash.file, append=T)
  cat(paste("write.bash(): Done with bash file...\n"))

  return(bash.file)
}
