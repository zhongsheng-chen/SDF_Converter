# A SDF Conversion Tool
This tool is intended to handle SDF files from [Mass of North America (MoNA)]

## Package requirement
+ [Tensorflow 1.14.0]
+ [RDKit 2018.09.3]
+ [Openbabel 3.0.0]

## Dataset description
The dataset in SDF format used to test this tool was download from 
[Vaniya-Fiehn Natural Products Library of MoNA](https://mona.fiehnlab.ucdavis.edu/downloads). 
It was failed to retrieve molecules from its blocks, 
each of which starts with molecule's title line and ends at the four dollar signs (`$$$$`). 

## Motivation
For some reasons unknown or unveiled yet by [MoNA](https://mona.fiehnlab.ucdavis.edu), all the SDFs
provided are not work well when loaded. It seems to me that I've gotten the reasons why they can not works well after 
aligning them to a standard SDF file. Following SDF format specifications, comparisons show that those dataset files 
are missing some necessary lines and `M  END` ahead of atom's coordination and bond's connections. Please see SDF 
file specifications from a overview on [Chemical Table File](https://en.wikipedia.org/wiki/Chemical_table_file). 
In order to convert those bad SDFs to their good counterparts, both lines and `M  END` required have been
append to the raw files properly.

## Usage
For converting a bad SDF dataset, followed by specifying the path to it, a
converted SDF file would be written, 
as do a file for keeping that failed blocks if `failed_block_file_name` is designated.
```
$ python convert_sdf_utils.py \

            --path_to_bad_sdf=/sdf/like/file/path \

            --failed_block_file_name=/save/failed/block/to/file \

            --output_dir=/save/path/to/converted/sdf \

            --alsologtostderr
```

When loading molecules from the converted SDF, it is worth mentioning that you can reset 
the global constant variable `MAX_ATOMS` in `mass_spec_constants.py` to a proper value to 
passe out any molecule whose number of atoms is below `MAX_ATOMS`. 

For example, a maximum number of atoms of the converted SDF for `MoNA-export-HMDB.sdf` 
is `92`, so I assigned here `MAX_ATOMS` to `1000` so that make ensure all the molecules stored
in that converted SDF can be fully loaded.
```python
# Assign MAX_ATOMS to 1000 in mass_spec_constants.py
MAX_ATOMS = 1000
MAX_ATOM_ID = 1000
```
After that settings, you can load molecules from the converted SDF with
``` python
import parse_sdf_utils

def main():
    converted_sdf_name = 'path/to/converted/SDF'
    mol_list = parse_sdf_utils.get_sdf_to_mol(converted_sdf_name)
```

## Acknowledgement
Some modules are imported from google brain team's efforts on [deep-molecular-massspec],
which give a easy way to parse molecules from SDFs. 


[//]: <> (Refences)
[Mass of North America (MoNA)]: https://mona.fiehnlab.ucdavis.edu
[Tensorflow 1.14.0]: https://www.tensorflow.org/install/pip
[RDKit 2018.09.3]: https://www.rdkit.org/docs/Install.html
[Openbabel 3.0.0]: http://openbabel.org/wiki/Category:Installation
[deep-molecular-massspec]: https://github.com/brain-research/deep-molecular-massspec


