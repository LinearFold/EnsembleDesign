# EnsembleDesign: Messenger RNA Design Minimizing Ensemble Free Energy via Probabilistic Lattice Parsing
This repository  hosts the source code for the research paper titled *EnsembleDesign: Messenger RNA Design Minimizing Ensemble Free Energy via Probabilistic Lattice Parsing* by [Ning Dai](https://ndai.ai/), [Tianshuo Zhou](https://eecs.oregonstate.edu/~zhoutian), [Wei Yu Tang](https://eecs.oregonstate.edu/~tangwe), [David H. Mathews](https://rna.urmc.rochester.edu/), and [Liang Huang](https://eecs.oregonstate.edu/~huanlian), to appear in the *Proceedings of ISMB 2025* (*Bioinformatics*, 2025).

If you use this repository in your research, please cite the above paper as well as [the LinearDesign paper](https://www.nature.com/articles/s41586-023-06127-z).


## Dependencies
- Clang 11.0.0 (or above) or GCC 4.8.5 (or above)
- Python 3

## To Compile
Run the following commands:

```
chmod +x setup.sh
./setup.sh
```

## Usage

The mRNA Design program is designed to process protein sequences provided in a FASTA file and apply our mRNA Design algorithm to each sequence. The results will be stored in `output_dir`, with each sequence having its own subfolder containing the raw outputs from all executions. The optimal solution for each sequence, determined across all runs, will be displayed to the user.

### Command Syntax

Use the following command format to run the mRNA Design program:

```
python EnsembleDesign.py [--fasta <path>] [--output_dir <path>] [--beam_size <int>] [--lr <float>] [--epsilon <float>] [--num_iters <int>] [--num_runs <int>] [--num_threads <int>]
```

- `--fasta <path>`: Specifies the path to the input protein FASTA file. The default is `./examples.fasta`.
- `--output_dir <path>`: Sets the directory for saving output files. The default is `./outputs`.
- `--beam_size <int>`: Determines the beam size for beam pruning. The default is `200`. A smaller beam size can speed up the process but may result in search errors.
- `--lr <float>`: Sets the learning rate for projected gradient descent. The default is `0.03`.
- `--epsilon <float>`: Defines the epsilon parameter for soft-MFE initialization. The default is `0.5`.
- `--num_iters <int>`: Specifies the number of optimization steps. The default is `30`.
- `--num_runs <int>`: Indicates the number of execution runs per sample file. The default is `20`. More runs increase the likelihood of finding optimal solutions but require more computational resources.
- `--num_threads <int>`: Sets the number of threads in the thread pool. The default is `8`. Adjust this based on your CPU's core count to optimize parallel processing without overloading your system.

### Example Command and Expected Output

To run the program with a specific set of parameters, you can use a command similar to the following (note that we are using smaller parameters to make it run faster):

```
python EnsembleDesign.py --fasta data/examples.fasta --output_dir ./outputs --beam_size 100 --lr 0.03 --epsilon 0.5 --num_iters 3 --num_runs 2 --num_threads 4
```

The expected output is:

```
>seq1|Ensemble Free Energy: -18.0 kcal/mol
AUGAAUACCUAUCAUAUUACUCUGCCGUGGCCGCCGAGUAAUAAUAGGUAU
>seq2|Ensemble Free Energy: -43.71 kcal/mol
AUGAAUGAGUACCAGUUUGUGCUCCCCUAUCCGCCGUCAUUGAAUACUUAUUGGCGGCGGAGGGGGAGCCAGUACUACAUU
>seq3|Ensemble Free Energy: -70.29 kcal/mol
AUGGCGACCGUUCUCCUGGCGCUUCUGGUUUAUUUGGGUGCGCUGGUGGAUGCGUACCCAAUUAAGCCAGAAGCGCCAGGAGAGGACGCCUUCUUAGGG
```

Note that following LinearDesign, we use `-V -d0 -b0` options in LinearPartition to evaluate the ensemble free energy, which uses the Vienna energy model with no dangling ends and no beam search (although the mRNA design code does use beam search).

## Other Files

Alongside the main application, this repository includes additional files that are useful for testing and understanding the capabilities of the mRNA Design tool:

- `data/examples.fasta`: This file contains 3 example protein sequences. It is designed to offer a quick and straightforward way to test the functionality of the software with pre-defined input. This file serves as the default input for the tool.

- `data/uniprot.fasta`: Contains 20 protein sequences from the UniProt database used in our experiments.

- `data/covid_spike.fasta`: Contains the SARS-CoV-2 spike protein sequence used in our experiments, taken from the LinearDesign paper (Zhang et al., Nature 2023).

- `coding_wheel.txt`: RNA codon table (i.e., the genetic code).

- `codon_usage_freq_table_human.csv`: A human codon usage frequency table from the LinearDesign paper, which is used by LinearDesign code to compute the codon adapation index (CAI). It is originally from [GenScript](https://www.genscript.com/tools/codon-frequency-table) (choose human instead of yeast).

- `tools/LinearDesign`: [LinearDesign codebase](https://github.com/LinearDesignSoftware/LinearDesign) (baseline, for min MFE design).

- `tools/LinearPartition`: [LinearPartition codebase](https://github.com/LinearFold/LinearPartition) (to calculate ensemble free energy). Note that we use `-V -d0 -b0` options to evaluate energies, following LinearDesign.