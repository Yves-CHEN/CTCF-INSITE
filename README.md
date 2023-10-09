# CTCF-INSITE
CTCF-INSITE is a very useful tool to predict the persistently-bound CTCF-BS (P-CTCF-BS) out of all the hundreds of thousands of sites, by learning the features of experimentally-determined P-CTCF-BS using RandomForest/logistic model.

# How to run
Step1: generate an annotation file
```
cd annotator 
bash run.sh # set the parameter from nextflow.config before hiting the run
```

Step2: run prediction
Go to https://when.shinyapps.io/ctcf-insight/
Upload the annotation file to 

Cite:
