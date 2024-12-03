# dig-ldsc-methods

Houses the code for two scripts and two methods for use in generating s-LDSC results

## Scripts
### snpmap

To transform generic GWAS summary stats datasets into a sumstats input for exclusive use of the s-LDSC method a map is needed.
To start run the bootstrap script (`src/scripts/sumstats_inputs/sumstats_inputs.bootstrap.sh`)

Using the G1000 datasets for AFR, AMR, EAS, EUR, and SAS ancestries the script (`src/scripts/sumstats_snpmap/make_sumstats_snpmap.py`) 
will generate a map from varId (chromosome:position:reference:alternate) to rsID for both standard and flipped (alternate:reference)
varIDs and for both hg19 and hg38 human genome builds.

The resulting datasets should be placed either locally or in s3 under a directory `<bucket/director>/bin/snpmap/`.

These inputs are also available upon request.

### s-LDSC inputs

For fast s-LDSC regressions the baseline model and tissue (and eventually user) annotations need to be transformed into a set of 
overlap numpy arrays, ld numpy arrays, parameter names, and count files. To start run the bootstrap script 
(`src/scripts/sldsc_inputs/sldsc_inputs.bootstrap.sh`)

Running both `src/scripts/sldsc_inputs/make_overlap.py` and `src/scripts/sldsc_inputs/make_sldsc_inputs.py` takes around
1.5 hours per ancestry on a personal computer. Zip each individual ancestry results (e.g. for EUR navigate to the resulting
`src/scripts/sldsc_inputs/inputs/EUR` directury and run `zip -r ../sldsc_inputs.EUR.zip ./*`).

The resulting zip files should be placed either locally or in s3 under a directory `<bucket/director>/bin/sldsc_inputs/`

These inputs are also available upon request.

## Methods (sumstats, sldsc)

Once the necessary inputs are generated or downloaded the two available methods (run in series, sumstats -> sldsc) are 
run either remotely using AWS Batch or locally.

### local

Navigate to `src` where the `main.py` file is located. The command to run is:

`S3_BUCKET=<local or remote directory> INPUT_PATH=<local path for input files> python main.py --username=<username> --dataset=<dataset> --method=<sumstats|sldsc>`

S3_BUCKET is where the raw data and output data will be stored. The INPUT_PATH is where the input data is housed 
(either locally or it will be downloaded from the same directory as defined in S3_BUCKET, it will only be downloaded once
if you are running multiple analyses in series)

To run locally place your raw GWAS file in a directory e.g.
`<some directory>/userdata/<username>/genetic/<dataset>/raw/`

You will need to define a metadata file to be placed with the GWAS data like:
```
{
    "file": "<NAME OF FILE>",
    "ancestry": "<AFR|AMR|EAS|EUR|SAS>",
    "separator": "<\t|,>",
    "genome_build": "<GRCh37|GRCh38>",
    "col_map": {
        "chromosome": "<chromosome field as str>",
        "position": "<position_field as int>",
        "reference": "<reference allele field>",
        "alt": "<alternate allele field>",
        "pValue": "<pValue field as float>",
        "beta": "<beta field as float>",
        "n": "<n field as float>"
    }
}
```

If the GWAS has oddsRatio instead of beta that can be defined instead:

```
{
    "file": "<NAME OF FILE>",
    "ancestry": "<AFR|AMR|EAS|EUR|SAS>",
    "separator": "<\t|,>",
    "genome_build": "<GRCh37|GRCh38>",
    "col_map": {
        "chromosome": "<chromosome field as str>",
        "position": "<position_field as int>",
        "reference": "<reference allele field>",
        "alt": "<alternate allele field>",
        "pValue": "<pValue field as float>",
        "oddsRatio": "<odds ratio field as float>",
        "n": "<n field as float>"
    }
}
```

Similarly if a sample size field (`n`) is unavailable you can set an effective sample size for all rows

```
{
    "file": "<NAME OF FILE>",
    "ancestry": "<AFR|AMR|EAS|EUR|SAS>",
    "separator": "<\t|,>",
    "genome_build": "<GRCh37|GRCh38>",
    "col_map": {
        "chromosome": "<chromosome field as str>",
        "position": "<position_field as int>",
        "reference": "<reference allele field>",
        "alt": "<alternate allele field>",
        "pValue": "<pValue field as float>",
        "beta": "<beta field as float>"
    },
    "effective_n": <float>
}
```
### batch

This run is identical to local except all of the require data is required to be houses in an available s3 bucket.

TODO: How to define an appropriate IAM role, Batch Job Definition, Batch queue

To ultimately run a remote job AWS cli can be used e.g.
`aws batch submit-job --job-name <name> --job-queue <queue> --job-definition <definition> --parameters username=<username>,dataset=<dataset>,method=<sumstats|sldsc>    `

The resulting output will be stored in the s3 bucket defined in Job Definition.

## Server + GUI

see dig-job-server for code and information about using the job server and associated GUI to set up a remote system
for uploaded and running LDSC analyses in AWS.

## Tests

Requires pytest. To run:

`python -m pytest`
