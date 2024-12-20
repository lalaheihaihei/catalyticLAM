# Scripts

This section contains some useful scripts to generate cluster structures and convert the format of files.

- `cif2pos.py`: Convert CIF file to POSCAR.
- `get-cluster.py`: Generate the structures of metal clusters in xyz format.
- `json2cif.py`: Convert JSON file to CIF file.
- `xyz2pos.py`: Convert XYZ file to POSCAR.
- `sim_model.py`: For deleting the unnecessary keys in checkpoint files (oc22).
- `cal_nframes.py`: Calculate the number of frames in a dataset with dp (deepmd-kit) format.
- `make_test.py`: Make a dataset test with lmdb format.

## Usage

1. cif2pos.py

```
python cif2pos.py --input_dir CIF --output_dir POSCAR   # Convert the cif files in the CIF folder to corresponding POSCAR files in the POSCAR folder
python cif2pos.py --input_dir CIF
# --input_dir represents the path to the input directory containing CIF files
# --output_dir represents the path to the output directory for POSCAR files, default is 'POSCAR'
```

2. get-cluster.py

```
python get-cluster.py   # generate all cluster structures in xyz format
python get-cluster.py --type all --element all   # generate all cluster structures in xyz format
python get-cluster.py --type all --element Cu   # generate all Cu-related clusters
python get-cluster.py --type s_fcc --element all   # generate all supported metal clusters
python get-cluster.py --type fcc --element Au   # generate Cu-related clusters in FCC lattice
python get-cluster.py --type hcp --element Zn   # generate Zn-related clusters in HCP lattice
# --type represents which type of structure to generate, including 's_fcc', 'fcc', 'bcc', 'hcp', and 'all', default is 'all', where 's_fcc' reprents the supported cluster
# --type represents which element or metal to use, default is 'all'
```

3. json2cif.py

```
python json2cif.py --file input.json --outdir CIF   # Convert the json file to corresponding cif files in the CIF folder
python json2cif.py --file input.json
# --file represents the path to the JSON file
# --outdir represents the path to the output directory for CIF files, default is 'CIF'
```

4. xyz2pos.py

```
python xyz2pos.py --input_dir XYZ --output_dir POSCAR --padding 5.0
python xyz2pos.py --input_dir XYZ
# --input_dir represents the path to the input directory containing XYZ files
# --output_dir represents the path to the output directory for POSCAR files, default is 'POSCAR'
# --padding reprensts the distance of the atoms from the box boundary in three directions, default is '5.0'
```

5. sim_model.py

```
python sim_model.py
# You should revise the path of checkpoint files and specify the keys you want to delete by yourself.
```

6. cal_nframes.py

```
python cal_nframes.py --dataset "dataset path"
return the number of frames in a dataset
# --dataset represents the path to the dataset
```

7. make_test.py

```
python make_test.py --checkpoint "checkpoint path" --dataset "dataset path" 
return thr rmse and mae error of a dataset
# --checkpoint represents the path to the checkpoint file
# --dataset represents the path to the dataset file
```
