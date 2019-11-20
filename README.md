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
| `-t` | (used if -g is not set) Minimum recurrence frequency of the seeds |
| `-e` | Minimum number of tumors that each node on the trajectory has to follow the seed in |
| `-i` | (optional) The number of iterations for experimental estimation of the p-value of the size of the largest subgraph |
| `-a` | Threshold for the edge confidence in the ILP formulation |
| `-f` | Output folder name; all output files are placed into this folder |
