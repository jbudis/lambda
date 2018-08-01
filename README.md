# Lambda method

Method for reducing number of uninformative calls in non-invasive prenatal testing.

Detailed description of statistical models is in [preprint manuscript](https://www.biorxiv.org/content/early/2018/06/22/353862).

## Requirements
Scripts were tested on following configuration
```
Python3 3.5.2
R 3.2.3
```

Python requires libraries
```
numpy (1.11.3)
pandas (0.19.2)
rpy2 (2.8.6)
```

## Usage

Train parameters
```bash
python3 train.py --help
usage: train.py [-h] train_dir param_file

Trains parameters for FL a NCV models.Parameters are then utilized in the
evaluate.py script to infer predictions for aneuploidy of diagnosed fetus

positional arguments:
  train_dir   Directory with TSV count files for sample with healthy euploid
              fetus. All files with ".tsv" suffix in the directory would be
              used for training. Count file contains the corrected number or
              reads for each assumed fragment length (50-220, organised in
              columns) and autosomes (chr1..chr22, organised in rows)
  param_file  Output YAML file with trained parameters

```

Evaluate sample
```bash
python3 evaluate.py --help
usage: evaluate.py [-h] counts param_file

Calculates standard and FL-based score for diagnosed sampleand predicts
diagnosis, i.e. positive, false positive, uninformative, negative

positional arguments:
  counts      TSV file with the corrected number or reads for each assumed
              fragment length (50-220, organised in columns) and autosomes
              (chr1..chr22, organised in rows)
  param_file  YAML file with trained parameters (output of the train.py
              script)

```

# Example

Set of test files is located at the examples/ directory.

First train parameters 
```bash
python3 train.py example/train example/my_params.yaml
```

Output file example/my_params.yaml should correspond to the pre-computed parameters in example/parameters.yaml.

Now use trained parameters to predict trisomy of chromosome 13, 18 and 21.

Healthy samples should have resulting scores lower than 2.5.
```bash
python3 evaluate.py example/test/neg_101.tsv example/my_params.yaml
```

Trisomic samples should have scores of affected chromosomes higher than 4. 
```bash
python3 evaluate.py example/test/t21_001.tsv example/my_params.yaml
```

# Citation
If you find method useful, please cite 

Budis J, Gazdarica J, Radvanszky J, Szucs G, Kucharik M, Strieskova L, Gazdaricova I, Harsanyova M, Duris F, Minarik G, Sekelska M. Innovative method for reducing uninformative calls in non-invasive prenatal testing. arXiv preprint arXiv:1806.08552. 2018 Jun 22.
