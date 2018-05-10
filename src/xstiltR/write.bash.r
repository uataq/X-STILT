# subroutine to write bash script for running configuration file
# to run this .sh file, "sbatch *.sh"
# by DW, 05/09/2018

write.bash <- function(job.time="24:00:00", num.nodes, num.cores,
                       account="lin-kp", partition="lin-kp",
                       email="dien.wu@utah.edu", email.type=c("FAIL,BEGIN,END"),
                       workdir, filenm="XSTILT.sh", run.filenm="run_XSTILT"){

  write("#!/bin/bash", filenm)

  sbatch.setup <- paste("\n#SBATCH --time=", job.time,
                        "\n#SBATCH --nodes=", num.nodes,
                        "\n#SBATCH --ntasks=", num.cores,
                        "\n#SBATCH --account=", account,
                        "\n#SBATCH --partition=", partition,
                        "\n#SBATCH --mail-user=", email,
                        "\n#SBATCH --mail-type=", email.type,"\n", sep="")

  write(sbatch.setup, file=filenm, append=T)

  write("rm -f multiprog.conf\n", filenm, append=T)

  # loop over each job
  job.loop <- paste("JOBNR=0\nSLURM_NTASKS=", num.cores,
                    "\nwhile [[ $JOBNR -lt $SLURM_NTASKS ]]; do\n",
                    "  echo $JOBNR ", file.path(workdir, run.filenm),
                    " $JOBNR >> multiprog.conf\n",
                    "  let JOBNR=JOBNR+1\ndone\n", sep="")

  write(job.loop, filenm, append=T)

  write("srun --multi-prog multiprog.conf\n", filenm, append=T)
  cat(paste("write.bash(): Done with bash file...\n"))
}
