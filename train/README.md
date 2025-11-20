# Training input files, configurations and checkpoints by DPA1 and Gemnet-OC
## Table of Contents

- [Table of Contents](#table-of-contents)
   - [0. CLAM dataset (version: 2024Q2)](#0-clam-dataset-version-2024q2)
   - [1. Finetune model](#1-finetune-model)
     - [1.1 Finetune by DPA2 in DeePMD-kit](#11-finetune-by-dpa-in-deepmd-kit)
     - [1.2 Finetune by Gemnet-OC](#12-finetune-by-gemnet-oc)
   - [2. Checkpoints and configurations](#2-checkpoints-and-configurations)
   - [3. Test Results](#3-test-results)

## 0. CLAM dataset (version: 2024Q2)

All bulk, slab, and molecules dataset are obtained by AIMD simulations and labeled per 10 steps. The initial structures are discussed in [README.md](../generation/README.md)

|                       |  2D   |  metal   |  QM9       |
| --------------------- | ------- | ------- | ---------------- |
|  Training dataset     |  183,215  |  411,125  |  596, 160         |
|  Validation dataset   |  10,180   |  22,835   |  33,120          |
|  Test dataset         |  10,178   |  22,839   |  33,120          |

## 1. finetune model

### 1.1 finetune by DPA in DeePMD-kit

Example of finetune based on pretrained `DPA2_medium_28_10M_beta4.pt`, note that `--model-branch` is needed to specify the task head.

```sh
dp --pt train finetune.json --finetune DPA2_medium_28_10M_beta4.pt --model-branch Domains_OC2M > finetune.out
```

### 1.2 finetune by Gemnet-oc

Example of finetune based on OC pretrained `gnoc_oc22_oc20_all_s2ef.pt`, and if you want to finetune the model from your self-made checkpoint, the `epoch`, `steps` information need to be removed by the [script](../scripts/sim_model.py).

```sh
python main.py --mode train --config-yml finetune.yml --checkpoint gnoc_oc22_oc20_all_s2ef.pt --print-every 1000 >> out
```

## 2. Checkpoints and configurations

The pretrained checkpoints by Gemnet-OC can be found [here](https://fair-chem.github.io/core/model_checkpoints.html).

The pretrained checkpoints by DPA2 can be found [here](https://www.aissquare.com/models?page=1&type=models).

## 3. Test Results

Performance of the trained model on the test dataset by DPA2 and Gemnet-OC:

| model      | metrics | energy (meV/atom) |||| force (meV/Å) ||||
|------------|---------|:-----------------:|:-------:|:----:|:----:|:-------------:|:------:|:----:|:----:|
|            |         | OC22              | metal   | QM9  | 2D   | OC22          | metal  | QM9  | 2D   |
| GNOC-CLAM  | MAE     | 24.5              | 1.8     | 2.9  | 3.2  | 45.5          | 28.5   | 21.5 | 27.1 |
| DPA2-CLAM  | MAE     | 68.1              | 5.8     | 7.8  | 5.5  | 277.0         | 72.2   | 78.0 | 84.6 |
