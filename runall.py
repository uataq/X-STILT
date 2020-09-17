import os
import time
import glob
import datetime
import subprocess

cmd_template = ["Rscript", "--vanilla", "run_xstilt_covid.r"]

t = 60 * 5  # sec
pid = None
idx = 0
jobs = [[1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11], [12]]

print("Generating jobs[{}]: {}".format(len(jobs), jobs))
#exit()

start_time = None
end_time = None

def submit():
    global idx
    global pid
    global start_time

    cmd = cmd_template + list(map(str, jobs[idx]))
    job = subprocess.run(cmd, stdout=subprocess.PIPE)
    job_stdout = str(job.stdout)
    print(job_stdout)
    
    # search for PID after "submitted batch job" from running the main script
    pid_line_idx = job_stdout.find("Submitted batch job")
    print(pid_line_idx)

    if pid_line_idx != -1:
        pid = job_stdout[pid_line_idx:].split(" ")[-1][:-3]
        start_time = datetime.datetime.now()
        print("Caught {} at {} ".format(pid, start_time))
    else:
        print("Submit job failed!")
        print(" ".join(cmd))
    idx = idx + 1

submit()

while True:
    if pid == None:
        print("No job is currently running.")
        submit()
        continue

    current_jobs = str(subprocess.run(["squeue", "-A", "pow"], stdout=subprocess.PIPE).stdout)

    if current_jobs.find(pid) == -1:
        end_time = datetime.datetime.now()
        print("{} had finished at {}".format(pid, end_time))
        print("{} time diff: {}".format(pid, end_time - start_time))
        if idx >= len(jobs):
            print("all jobs had finished")
            break
        else:
            submit()
    else:
        print("{} is still running".format(pid))

    time.sleep(t)

