# CONETT
Combinatorial optimization method for identification of conserved evolutionary trajectories in tumors

### System Requirements
- make (version 3.81 or higher)
- g++ (GCC version 4.1.2 or higher)
- Gurobi Optimizer

### Compiling `CONETT`
In the `Makefile`, set `GUROBIROOT` to the path of your root Gurobi folder.

Simply run `make` command in the root CONETT folder. It will create the executables.

### Running `CONETT`
**Usage:**
```sh
./conett -p [input tumor graph] -g [gene sets to be used as seeds] -t [minimum seed recurrence] -e [subgraph tumor conservation rate] -i [permutation p-value iterations] -a [edge weight cutoff] -f [outputFolder]
```

| Parameters | Description |
| ------ | ------ |
| `-p` | Tumor graph file which shows temporal order of alterations |
| `-g` | (optional) File of gene (sets) that should be used as seeds |
| `-t` | Minimum recurrence frequency of the seeds |
| `-e` | Minimum number of tumors that each node on the trajectory has to follow the seed in |
| `-i` | (optional) The number of iterations for experimental estimation of the p-value of the size of the largest subgraph |
| `-a` | Threshold for the edge confidence in the ILP formulation |
| `-f` | Output folder name; all output files are placed into this folder |

#### Example

Command to run:
```sh
./conett -p data/TRACERxDAGs_inac_singleRoot.txt -g SeedPathways.txt -t 10 -e 10 -f TRACERx_ccRCC_20191030 -a 0.85
```
Tumor graph file is a transitive precedence graph, describing each **directed** edge with the name of the patient/sample and two pairs of node names and the corresponding alteration types. For example:
```sh
T001  GL  - gene1 SNV
T001  GL  - gene2 CNGAIN
T001  gene1 SNV gene2 CNGAIN
T002  GL  - gene1 CNLOSS
...
```
