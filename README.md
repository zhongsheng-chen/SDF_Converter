# A SDF Conversion Tool
This utility is intended to handle SDF files from [Mass of North America (MoNA)](https://mona.fiehnlab.ucdavis.edu)

## Package requirement
+ [Tensorflow 1.14.0](https://www.tensorflow.org/install/pip)
+ [RDKit 2018.09.3](https://www.rdkit.org/docs/Install.html)
+ [Openbabel 3.0.0](http://openbabel.org/wiki/Category:Installation)


## Dataset description
The dataset in SDF format used to test this code is download from 
[Vaniya-Fiehn Natural Products Library of MoNA](https://mona.fiehnlab.ucdavis.edu/downloads). It is failed to 
load as molecules from its blocks, each of which starts with molecule's title line and ends at the four dollar signs
(`$$$$`). 

## Motivation
For some reasons unknown or unveiled yet by [MoNA](https://mona.fiehnlab.ucdavis.edu), all the SDF datasets
provided are not work well as them are loaded. It seems to I've gotten the reasons why they can not works well after 
aligning them to a standard SDF file. Following SDF format specifications, comparisons find that those dataset files 
are missing some necessary lines and `M  END` ahead of atom's coordination and bond's connections. Please see SDF 
file specifications from a overview on [Chemical Table File](https://en.wikipedia.org/wiki/Chemical_table_file). 
In order to convert those bad SDF datasets to their good counterparts, both required lines and `M  END` have been
append to the raw files properly.

## Usage
For converting a bad SDF dataset, specify its path, a file which save failed blocks and a directory where the 
converted SDF file should be saved.
```
$ python convert_sdf_utils.py \

            --path_to_bad_sdf=/sdf/like/file/path \

            --failed_block_file_name=/save/failed/block/to/file \

            --output_dir=/save/path/to/converted/sdf \

            --alsologtostderr
```
