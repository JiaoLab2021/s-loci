# s-loci

[![GitHub last commit](https://img.shields.io/github/last-commit/JiaoLab2021/s-loci.svg?label=Last%20commit&logo=github&style=flat)](https://github.com/JiaoLab2021/s-loci/releases)

## Introduction
Software for S-locus Genotyping in Rutaceae Samples Based on Second-Generation Sequencing Data

## Requirements

Please note the following requirements before building and running the software:

* `Linux` operating system
* cmake version `3.12` or higher
* C++ compiler that supports `C++17` or higher, and the `zlib` library installed (we recommend using GCC version `"7.3.0"` or newer) for building `s-loci`
* `SOAPdenovo-63mer` and `minimap2`

## Installation

### Building on Linux

Use the following script to build the software:

1. First, obtain the source code.

```shell
git clone https://github.com/JiaoLab2021/s-loci.git
cd s-loci
```

2. Next, compile the software and add the current directory to your system's `PATH` environment variable.

```shell
g++ s-loci.cpp -o s-loci -lpthread -lz -lhts -O2 -std=c++17
echo 'export PATH="$PATH:'$(pwd)'"' >> ~/.bashrc
source ~/.bashrc
```

## Usage

### Running

For convenience, let's assume the following file names for the input:

* `s-locus.fa`
* `s-rnase.dedup.cds`

**1. Filter:**

```shell
s-loci kmerfilter -i s-locus.fa -r sample.1.fq.gz -r sample.2.fq.gz -m 0.85 --prefix 0.85
```

**2. Performing Genotyping:**

```shell
s-loci genotype -t 1 -f s-rnase.dedup.cds -r 0.85.sample.1.fq.gz -r 0.85.sample.2.fq.gz --prefix sample
```

**3. Assembly:**

```shell
s-loci assembly -i 0.85.sample.1.fq.gz -i 0.85.sample.2.fq.gz --avg_ins 500 -f s-rnase.dedup.cds --prefix sample --SOAPdenovo-63mer /path/SOAPdenovo-63mer --minimap2 /path/minimap2  -t 1
```