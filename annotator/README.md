# Description
This folder contains the nextflow script to invoke a docker to generate annotation.
Therefore there is no need to set up complete software/packages dependencies.
This script does not demand a lot of resource, so it is ok to run on laptops (except for very large amount of ChIP-seq peaks).

# Requirements
1. singularity or docker
2. nextflow

# Install
singularity: This is usually pre-installed on HPC. If not, then run https://docs.sylabs.io/guides/3.0/user-guide/installation.html.

java:  Install  Java 11 (or later, up to 20) as this is need by nextflow. Usually this is preinstall. Check installation by 

```
java --version
```

nextflow: Once java is installed,   run 
```
wget -qO- https://get.nextflow.io | bash

```

# Run

bash run.sh # set up nextflow.config before running


# Suggestion
1. Run on HPC for quicker performance. You can set multiple CPUs.
2. Adjust the nextflow.config to set up correct paths, e.g. path to hg19.fa
3. Use singularity over Docker, as later is usually not supported by HPC.


  
