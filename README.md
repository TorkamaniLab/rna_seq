RNA Seq Pipeline
================

*A pipeline for basic RNA Seq analysis.*


About this pipeline
-------------------

This pipeline performs the introductory steps in any RNA Seq analysis:

- Trimming and aligning the reads for a sample.
- Counting the reads that align to various genomic features.

**Important**: This pipeline does not attempt any further analysis like DESeq
or GSEA.

This pipeline is designed to run on a cluster environment, though some slight
modification of the run script will cause it to run locally (see Running
Locally). 


Using this pipeline
-------------------

This pipeline is built using metapipe. Before running, please run the setup 
script to fetch all of the required tools. If you don't want to use the setup 
sciript or you have different tools you'd like to use, you'll need to edit 
the metapipe file directly.

```bash
$ sh setup.sh path_to_bowtie_index path_to_gtf_file
```

Once the environment is set up, it's time to run the pipeline on your given
samples.

To run the pipeline on a single sample, do the following:

```bash
$ sh run_one.sh <sample_id> <read_1_filename> <read_2_filename> 
```

To run a variety of samples, create a `samples.txt` file containing a line for
each sample.

```
<sample_id> <read_1_filename> <read_2_filename>
```

Then run the following script.

```bash 
$ sh run_all.sh samples.txt
```


Project Organization
--------------------

This pipeline relies on a variety of external tools. Please make sure the paths
are up to date for these tools.

Shortcut execution scripts are located on the top level and are the preferred
way to execute the pipeline.


Running Locally
---------------

In order to run this pipeline locally, edit the pipeline.mp file and uncomment
the paths to the local binaries. Then comment out or remove the cluster binary paths.

**Note**: The Trimmomatic settings has a path inbuilt which depends on the version 
you're using. You might have to change that if you're running it locally.
