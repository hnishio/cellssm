
## Test environments
* macOS 11.6.4, R 4.1.1
* Windows 10, R 4.2.2
* Ubuntu 18.04, R 4.2.1

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE:

❯ checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    ‘01_dist_vis’ ‘03_ssm_KFAS’ ‘04_nomodel’ ‘05_lm_dist_beta’
    ‘06_lm_dist_start’ ‘07_lm_signal’ ‘11_dist_vis’ ‘13_ssm_KFAS’
    ‘14_nomodel’ ‘16_lm_dist_start’ ‘17_lm_signal’

These directories are created within functions or example codes.
For ease of use, the main functions are designed to save output files 
in the created directories. Thus, this is necessary for this package.

## Downstream dependencies
There is no downstream dependencies for this package.
