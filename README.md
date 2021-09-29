## clipperhelper
Script to create files required to run [clipper](https://github.com/YeoLab/clipper) on custom references.

### pre-requisites
clipperhelper only needs pandas and argparse python libraries. Details provided in the `requirements.txt` file, if you wish to create a conda environment.

### usage
```
% python clipperhelper/util/make_custom_species_files.py --help
usage: make_custom_species_files.py [-h] --gtf GTF --species SPECIES

Create clipper input files for custom species

optional arguments:
  -h, --help         show this help message and exit
  --gtf GTF          input GTF file
  --species SPECIES  species name to be used as output file prefix
```

### example
```
% python clipperhelper/util/make_custom_species_files.py --gtf hg19+virus.gtf --species hybrid
```
The above command will create the following files:

  * `hybrid_genes.bed`
  * `hybrid_exons.bed`
  * `hybrid.AS.STRUCTURE.COMPILED.gff`

All of these are required to run clipper on a custom genome with customized annotations (GTF).

#### comments etc.
Feel free to contact me at vishal.koparde@nih.gov