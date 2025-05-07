# MicroDNA
An algorithm to detect and quantify microDNAs from alignment data.

## microDNA
### Usage
```
usage: testSuffix.py [-h] [--bam_file BAM_FILE] [--fasta_file FASTA_FILE]

optional arguments:
  -h, --help                    show this help message and exit
  --bam_file BAM_FILE           BAM file to read
  --fasta_file FASTA_FILE       FASTA (.fna) file to read
```
### Installs/Downloads
First, ensure pysam is installed:
```
$ pip install pysam
```

Then, make sure the BAM and FASTA folders from the shared Google Drive are uploaded into the /data/ folder


### Example (using given data files)
First, ensure pysam is installed:
```
$ pip install pysam
```
To output results:
```
$ python microDNA.py
```
```
POSITION                  LENGTH       SCORE       
956749-956934             185 bp       87.91
2066285-2066645           360 bp       84.32
2403270-2403652           382 bp       68.17
2403298-2403643           345 bp       85.32
11458983-11459179         196 bp       62.99
17998305-17998648         343 bp       95.25
18062885-18063265         380 bp       66.88
21698974-21699333         359 bp       93.31
26053913-26054114         201 bp       75.80
27478693-27479037         344 bp       82.95
31775822-31776185         363 bp       70.98
35928489-35928760         271 bp       86.55
42105303-42105483         180 bp       88.56
53792372-53792739         367 bp       80.94
116247456-116247806       350 bp       86.98
121485089-121485434       345 bp       68.34
144058864-144059048       184 bp       78.50
161940384-161940765       381 bp       74.36
167461877-167462071       194 bp       86.83
181943749-181944106       357 bp       76.56
201093850-201094038       188 bp       89.35
206521552-206521736       184 bp       78.20
230653575-230653895       320 bp       65.57
230653581-230653895       314 bp       73.78
232062795-232063012       217 bp       62.91
236595941-236596163       222 bp       84.10

Number of circles: 26
```
