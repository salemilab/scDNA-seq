#<span style="color:gray">Tapestri pipeline sample pre-processing</span>
--

<b>[Step 1]</b> Assess end quality of fastq R2 reads using FastQC tool and determine whether to trim. Execute the following in terminal or run fastqc.sh shell file to run several in parallel.

```
(a) ml fastqc/0.11.7

(b) fastqc <name_of_fastq_file>
```
If ends require trimming (e.g. 17 bases), use the cutadapt parameters below and execute or run TrimR2array.sh to run several in parallel.

```
(a) ml cutadapt/3.4

(b) cutadapt -u -17 -o DNApooledTapestri_S2_L001_R2_001.fastq.gz DNApooledTapestri_S2_L001_R2_001_untrimmed.fastq.gz
```

<b>[Step 2]</b> While FastQC is running, modify the tapestri pipeline parameters for completeness and uniformity (recall this is for virus and not human derived cellular sequences).


<u>To change uniformity</u>:

Change the cell finder command in file \<installer>/bin/src/code\_GATK\_p1.sh - line 172 add additional parameter as --minor_threshold (default is 0.2) to 0.05.

<u>To change completeness</u>:

Edit the function get\_cell_count in the file \<installer>/lib/python3.6/site-packages/missionbio/core/cell\_utils.py. On line 84, change 80 to 5.

<b>[Step 3]</b> Modify config.yaml file and run the tapestri pipeline.