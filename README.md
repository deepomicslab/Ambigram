# Ambigram
[![DOI](https://zenodo.org/badge/442341093.svg)](https://zenodo.org/badge/latestdoi/442341093)

Ambigram is a graph-based algorithm to reconstruct the BFB local haplotype during the evolution process. Ambigram currently supports general WGS sequencing (PE), 10X linked-reads (10x), PacBio SMRT (PB) and Oxford Nanopore (ONT).

### Prerequisites
- [Cbc](https://github.com/coin-or/Cbc) solves integer linear programming.
- [htslib](https://github.com/samtools/htslib) is used for accessing common file formats, such as SAM, CRAM and VCF, used for high-throughput sequencing data.
- [SVAS](https://github.com/paprikachan/SVAS) integrates SEG and SV files into a .lh file, which is the standard input file of Ambigram.
- [SpecHap](https://github.com/deepomicslab/SpecHap) extracts barcodes from 10x data.
- [hpvpipe](https://github.com/panguangze/hpvpipe) generate JUNCS files by processing PB or ONT data.
- [Python3](https://www.python.org/downloads/)
- [samtools](https://github.com/samtools/samtools) deals with .bam files. 
- [pysam](https://pysam.readthedocs.io/en/latest/) counts the coverage of genomic positions by reads in region.

## Installation

1. Clone the repository and enter the directory:

```
git clone https://github.com/deepomicslab/Ambigram
cd ./Ambigram/
```
2. Create a conda environment with all dependencies and enter the environment:
```
conda env create --prefix=./Ambigram_env -f environment.yml
conda activate ./Ambigram_env
```
3. Create a build directory and compile Amibgram under it (use **sudo**, if required):

```
mkdir build
cd build
cmake ..
make && make install
```
4. To see the option Ambigram supports, run the following command:

```
Ambigram --op bfb --help
```

## Using Ambigram

### Data preprocessing

Ambigram supports various data types including Illumina pair-end (PE) reads, Oxford Nanopore (ONT) long reads, Pacific Biosciences (PB) long reads, 10x Genomics linked-reads with varying tumor purity and sequencing depth. The procedures and methods to process various data sets are available on the other GitHub page (https://github.com/deepomicslab/Ambigram_Paper). 

To begin with, we need two files: (1) SV file (e.g. test_sv.txt) that consists of start and end chromosome names, breakpoint positions, and strands as well as copy number or depth of each SV junction; (2) BAM file that contains read information of all SVs in the SV file. 

* **test_sv.txt**

```
chrom_5p	bkpos_5p	strand_5p	chrom_3p	bkpos_3p	strand_3p	avg_cn
chr7	55282001	-	chr7	55282001	+	2
chr7	55283001	-	chr7	55283001	+	1
chr7	55285000	+	chr7	55285000	-	1
```

Then we need to divide the corresponding genomic regions into segments based on the breakpoints of SVs and read depth. The segment information is stored in a SEG file (e.g. test_seg.txt) that indicates genomic regions including chromosome names, start and end coordinates, and copy number or depth. We provide a script for users to generate the SEG file by inputting a SV file and the corresponding BAM file. To generate the SEG file, run the following command:

```
python Ambigram/script/bfb_scripts.py generate_seg --sv_file=[path to SV file] --bam_file=[path to bam file] --wgs_depth=[whole genome average depth (default: 100)] --tumor_purity=[sample tumor purity (default: 1)]
```

* **test_seg.txt**

```
chr7:55281001-55282000	1
chr7:55282001-55283000	5
chr7:55283001-55284000	8
chr7:55284001-55285000	8
chr7:55285001-55286000	6
chr7:55286001-55287000	4
```

Finally, we need to integrate the SV and segment information by converting the SV and SEG files into an LH file (e.g. test.lh), which is the standard input file of our algorithm. To convert the SEG and SV files into an LH file, run the following command:

```
python Ambigram/script/bfb_scripts.py generate_lh --sv_fn=[path to SV file] --seg_fn=[path to SEG file] --sample=[output sample name]
```

* **test.lh**

```
SAMPLE test
AVG_CHR_SEG_DP 30
AVG_WHOLE_HOST_DP 30
AVG_JUNC_DP 30
PURITY 1
AVG_TUMOR_PLOIDY 2
PLOIDY 2m1
VIRUS_START 7
SOURCE 1
SINK 6
SEG H:1:chr7:55281001:55282000 60.0 2.0
SEG H:2:chr7:55282001:55283000 180.0 6.0
SEG H:3:chr7:55283001:55284000 240.0 8.0
SEG H:4:chr7:55284001:55285000 240.0 8.0
SEG H:5:chr7:55285001:55286000 120.0 4.0
SEG H:6:chr7:55286001:55287000 120.0 4.0
JUNC H:2:- H:2:+ 60.0 2.0 U B
JUNC H:3:- H:3:+ 30.0 1.0 U B
JUNC H:4:+ H:4:- 60.0 2.0 U B
JUNC H:6:- H:6:+ 60.0 2.0 U B
```

The input file has three parts: header, segment, and junction. The header spans from the first line to the line starting with "VIRUS_START", which contain some default properties, e.g., "SAMPLE" refers to the sample name, and users do not need to change them.

The segment part spans from all the lines starting with "SOURCE", "SINK", and "SEG". They contains all the segment information, where "SOURCE" and "SINK" indicate the first and last segment indices, respectively. For each "SEG" line, the information, following "H:" and split by colons, represents segment index, chromosome name, start breakpoint, and end breakpoint. The two trailing numbers are segment depth and copy number. If the copy number is -1, Ambigram will automatically calculate the copy number by the segment depth.

The junction part consists of all the lines starting with "JUNC", containing all the junction information. For each "JUNC" line, the second and third columns represent two segments connected by the junction. The sign '+' indicates the forward strand of the segment, while '-' refers to the reverse strand. The following two columns are junction depth and copy number. If the copy number is -1, Ambigram will automatically calculate the copy number by the junction depth. The last two columns are defualt value, indicating the junction is not inferred by the program. 

### Run Ambigram
Ambigram requires only one input file (**.lh**) that contains CN profile and SV information of local regions invovled with complex BFB. To decipher BFB paths, run the following command:

``` 
Ambigram --op bfb --in_lh [path to your .lh file] --lp_prefix [sample name]
```
The output will be shown in the console as following. It is a sequence of segments, i.e., a haplotype cased by BFB. Each number is the segment index which is followed by a sign that indicates the strand. It is noted that the vertical bar '|' represents a fold-back inversion.
```
1+2+3+4+5+6+|6-5-4-3-2-|2+3+4+|4-3-|3+4+|4-3-2-|2+3+4+5+6+|6-5-4-3-2-1-
```
### Extra options
To decipher BFB events with translocation, we just need to add one line (shown as following) to the end of the LH file. Ambigram supports 4 complex cases involved with BFB and translocation. 

1. Insertion before BFB (TRX-BFB) (e.g., virus segments are inserted into chr1) - [test data](https://github.com/deepomicslab/Ambigram_paper/blob/master/simulation/test6.lh)

``` 
PROP I1:chr8:virus:chr8 M:chr8
```
* Expected output (double vertical bars "||" represent translocation)
```
1+2+3+||6+||4+|4-||6-||3-2-|2+3+||6+||4+|4-||6-||3-2-
```

2. Insertion after BFB (BFB-TRX) (e.g., segments of chr6 and chr13 are inserted into chr2, and insertion starts at Segment 3 on chr2) - [test data](https://github.com/deepomicslab/Ambigram_paper/blob/master/simulation/test3.lh)

``` 
PROP I2:chr2:chr6:chr13 M:chr2 S:3
```
* Expected output
```
1+2+3+||5+6+7+|7-6-||8+9+||4-3-2-|2+3+4+|4-3-
```

3. Concatenation before BFB (TRX-BFB) (e.g., segments of chr7 are linked to 260T-HBV_C3-RC by translocation) - [test data](https://github.com/deepomicslab/Ambigram_paper/blob/master/oncovirus/HBV/260T_chr1.lh)

``` 
PROP C1:chr1:260T-HBV_C3-RC
```
* Expected output
```
TRX-BFB mode: BFB path in the first stage:
8+||2+3+4+5+|7-6-5-4-3-2-||8-|8+||2+3+4+5+|7-6-5-4-3-2-|2+3+4+5+
TRX-BFB mode: BFB path in the second stage:
10-||4+5+|7-6-5-4-3-2-||8-|8+||2+3+4+5+|7-6-5-4-3-2-|2+3+4+5+
```

4. Concatenation after BFB (BFB-TRX) (e.g., segments of chr2 are concatenated with segments of chr6) - [test data](https://github.com/deepomicslab/Ambigram_paper/blob/master/simulation/test4.lh)

``` 
PROP C2:chr2:chr6
```
* Expected output
```
1+2+3+4+|4-3-2-|2+3+||6+7+|7-6-|6+7+|7-6-5-
```
Besides, we can input a JUNCS file (e.g. test.juncs) that contains extra information from linked reads (10x), long reads (PB and ONT), and optical mapping, which may help Ambigram resolve more accurate BFB paths. A JUNCS file that comprises groups of segments, which are possible fragments on BFB paths. 

* **test.juncs**

```
6+ 6- 5- 4- 3- 2- 2+
2- 2+ 3+ 4+ 5+ 6+ 6-
6+ 6- 5- 4- 3-
```
Firstly, we can use BarcodeExtractor in SpecHap to get barcodes from 10x data and run the script (process_barcode.py) to generate the JUNCS file. To generate a JUNCS file for 10x data ([SpecHap](https://github.com/deepomicslab/SpecHap) should be installed):

``` 
BarcodeExtract [path to 10x .bam file] [output path of .bed file]
python Ambigram/script/process_barcode.py --seg_file=[path to SEG file] --bed_file=[path to .bed file] --sample_name=[prefix of the output .juncs file]
```

Besides, we can use another script (hpvpipe/main.py process_tgs) to directly get the JUNCS file from PB or ONT data. To generate a JUNCS file for PB or ONT data ([hpvpipe](https://github.com/panguangze/hpvpipe) should be installed):

```
samtools bam2fq [path to PB/ONT .bam file] | seqtk seq -A > [output path of .fasta file]
python hpvpipe/main.py process_tgs -r [path to the reference genome] -l [path to .lh file] -t [path to .fasta file] -o [output path of .juncs file]
```

Moreover, we can use SegAligner in AmpliconReconstructor to align optical map contigs and run the script (bfb_scripts.py OM2juncs) to generate the JUNCS file. To generate a JUNCS file for OM data ([SegAligner](https://github.com/jluebeck/AmpliconReconstructor#segaligner) should be installed):

```
SegAligner [reference.cmap] [BioNano_contig.cmap]
python Ambigram/script/bfb_scripts.py OM2juncs -i [path to SA_output_score_thresholds.txt] -p [prefix of output .juncs file]
```

With a **JUNCS** file, we can decipher BFB paths by running the following command: ([test data with JUNCS files](https://github.com/deepomicslab/Ambigram_paper/tree/master/simulation/tumor%20purity))

``` 
Ambigram --op bfb --in_lh [path to your .lh file] --juncdb [path to your .juncs file] --junc_info true --lp_prefix [sample name]
```
<!-- 
If you have a very complicated sample, e.g. a sample with high and various copy numbers, try **Sequential Mode** that will resolve a BFB path with length-decreasing components (without nested loops):

``` 
Ambigram --op bfb --in_lh [path to your .lh file] --juncdb [path to your .juncs file] --seq_mode true --lp_prefix [sample name]
``` -->

We also provide an option to solve BFB paths for **single-cell** data, which can reconstruct BFB paths with similar components. You can input **multiple LH files** that represents different subclones with distinct SV and CNV profiles: ([test data of single-cell mode](https://github.com/deepomicslab/Ambigram_paper/tree/master/single_cell/simulation))

``` 
Ambigram --op sc_bfb --in_lh [paths to your .lh files (seperated by ,), e.g., test1.lh,test2.lh,test3.lh] --lp_prefix [sample name]
```
## Author

Ambigram is developed by DeepOmics lab under the supervision of Dr. Li Shuaicheng, City University of Hong Kong, Hong Kong, China.
Should you have any queries, please feel free to contact us by <chaohuili3-c@my.cityu.edu.hk>.

## License
This project is licensed under the MIT License - see the [LICENSE.txt](./LICENSE.txt) file for details
