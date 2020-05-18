# EMS GWAS 
This repository contains the singularity environment to run the EMS GWAS pipeline.

Since compiling all the R packages takes a while (~45 mins), the singularity images are split into two parts.

Build part 1:
```
sudo singularity build EMS_GWA_PART1.sif singularity.def
```

Then build part 2:
```
sudo singularity build EMS_GWA_PART2.sif singularity_part2.def
```
