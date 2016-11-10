# DKFZBiasFilter

## Description

The DKFZ bias filter flags SNVs that appear to be biased based on the variant read support. We differentiate between bias on the PCR template strand and bias on the forward/reverse strand.  To avoid massive over filtering variants are only considered for filtering if the whole context field (96 triplet) appears to be biased. A detailed description of the method can be found [here] (https://wiki.oicr.on.ca/display/PANCANCER/DKFZ+SNV+bias+filtering?preview=/66951525/66952026/description_dkfz_bias_filter.pdf "description of DKFZBiasFilter").

The bias filter is part of the DKFZ workflow. The version of the filter prsented here is a stand alone version to post filter somatic SNV calls from other workflows.

## Usage

First you need a running docker installation.

Clone the git repository and build the docker.

```
$ git clone https://github.com/eilslabs/DKFZBiasFilter.git
$ cd DKFZBiasFilter/
$ docker build -t "DKFZBiasFilter" .
```

## Running the Docker

To run the docker you need a somatic vcf file, a tumor bam and the corresponding bam index file (bai) as well as the reference genome your sample was aligned with.

Then just run the Docker:
```
$ docker run \
    -v /path/to/tumor.bam:/home/pcawg/tumor.bam \
    -v /path/to/tumor.bam.bai:/home/pcawg/tumor.bam.bai \
    -v /path/to/reference.fa:/home/pcawg/hs37d5.fa \
    -v /path/to/somatic.vcf:/home/pcawg/input.vcf \
    -v /path/to/results_directory/:/home/pcawg/results \
    DKFZBiasFilter
```

You should find the results in the `/path/to/results_directory/` directory.

Please be aware that the somatic VCF file should only contain filtered "PASS" SNVs.
