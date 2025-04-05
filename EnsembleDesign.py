import os
import re
import sys
import subprocess
import threading
import argparse
import math
import shutil
from queue import Queue
from multiprocessing import Value, Lock

prog_path = "./bin/EnsembleDesign"


def read_fasta(file_path):
    records = []
    with open(file_path, "r") as file:
        seq_id = None
        seq_lines = []
        for line in file:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_id is not None:
                    records.append((seq_id, "".join(seq_lines)))
                seq_id = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)

        if seq_id is not None:
            records.append((seq_id, "".join(seq_lines)))
    return records

def get_mfe_solutoin(seq):
    full_command = f"cd ./tools/LinearDesign && echo {seq} | ./lineardesign"
    result = subprocess.run(full_command, shell=True, capture_output=True, text=True)
    mfe_solution = result.stdout.strip()
    return mfe_solution.split("\n")[-3].strip().split(" ")[-1]

def eval_partition(seq):
    full_command = f"cd ./tools/LinearPartition && echo {seq} | ./linearpartition -V -d0 -b0"
    result = subprocess.run(full_command, shell=True, capture_output=True, text=True)
    partition_result = result.stderr.strip()
    return float(partition_result.split(' ')[-2].strip())

def process_run_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    seq = re.search(r'Final mRNA sequence: (.*)', lines[-1]).group(1)
    efe = eval_partition(seq)
    return seq, efe

def execute_mrna_design(protein, output_dir, run_number, args, progress_tracker, lock):

    mfe_solutoin = get_mfe_solutoin(protein)
    command = f"echo {protein} | {prog_path} {args.beam_size} {args.num_iters} {args.lr} {args.epsilon} {run_number} {mfe_solutoin}"

    output_path = os.path.join(output_dir, f"run_{run_number}.txt")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    with open(output_path, "w") as output_file:
        subprocess.run(command, shell=True, stdout=output_file, stderr=subprocess.STDOUT)

    with lock:
        progress_tracker.value += 1
        sys.stderr.write(f"Completed: {progress_tracker.value} runs\n")
        sys.stderr.flush()

def worker(task_queue, progress_tracker, lock):
    while True:
        task = task_queue.get()
        if task is None:
            break
        execute_mrna_design(*task, progress_tracker=progress_tracker, lock=lock)
        task_queue.task_done()

def run_mrna_design(args):
    records = read_fasta(args.fasta)
    num_runs = args.num_runs
    num_threads = args.num_threads

    task_queue = Queue()
    progress_tracker = Value('i', 0)
    lock = Lock()

    # Initialize progress tracker
    total_runs = len(records) * num_runs
    sys.stderr.write(f"{len(records)} input sequences.\n{num_runs} runs for each sequence.\nTotal {total_runs} runs started.\n")

    # Create and start threads
    threads = []
    for _ in range(num_threads):
        t = threading.Thread(target=worker, args=(task_queue, progress_tracker, lock))
        t.start()
        threads.append(t)

    # Clean log folders
    for seq_id, _ in records:
        log_dir = os.path.join(args.output_dir, seq_id)
        if os.path.exists(log_dir):
            shutil.rmtree(log_dir)
            
    # Enqueue tasks
    for run_number in range(1, num_runs + 1):
        for seq_id, protein in records:
            log_dir = os.path.join(args.output_dir, seq_id)
            task_queue.put((protein, log_dir, run_number, args))

    # Block until all tasks are done
    task_queue.join()

    # Stop workers
    for _ in range(num_threads):
        task_queue.put(None)
    for t in threads:
        t.join()

    sys.stderr.write(f"All {total_runs} runs completed.\n")

    width = int(math.log10(num_runs)) + 1;
    for seq_id, protein in records:
        log_dir = os.path.join(args.output_dir, seq_id)
        results = []
        for run_file in sorted(os.listdir(log_dir)):
            if not run_file.endswith(".txt"):
                continue
            run_id = int(run_file[:-4].split("_")[-1])
            seq, efe = process_run_file(os.path.join(log_dir, run_file))
            results.append((run_id, seq, efe))

        logs = []
        best_seq, best_efe = None, 0
        for run_id, seq, efe in sorted(results, key=lambda x: x[0]):
            logs.append(f"[{run_id:>{width}}] {seq} | Ensemble Free Energy: {efe} kcal/mol")

            if efe < best_efe:
                best_seq, best_efe = seq, efe
        logs.append(f"[Best] {best_seq} | Ensemble Free Energy: {best_efe} kcal/mol")

        print(f">{seq_id}|Ensemble Free Energy: {best_efe} kcal/mol\n{best_seq}")

        result_file = os.path.join(log_dir, "results.txt")
        with open(result_file, 'w') as f:
            f.write("\n".join(logs))

def main():
    parser = argparse.ArgumentParser(description='Run EnsembleDesign on protein fasta file.')
    parser.add_argument('--fasta', type=str, default='./examples.fasta', help='Path to the input protein fasta file (default: ./examples.fasta)')
    parser.add_argument('--output_dir', type=str, default='./outputs', help='Directory to save output files (default: ./outputs)')
    parser.add_argument('--beam_size', type=int, default=200, help='Beam size for beam pruning (default: 200)')
    parser.add_argument('--lr', type=float, default=0.03, help='Learning rate for projected gradient decent (default: 0.03)')
    parser.add_argument('--epsilon', type=float, default=0.5, help='The epsilon paramter of soft-MFE initialization (default: 0.5)')
    parser.add_argument('--num_iters', type=int, default=30, help='Number of optimization steps (default: 30)')
    parser.add_argument('--num_runs', type=int, default=20, help='Number of runs per sample file (default: 20)')
    parser.add_argument('--num_threads', type=int, default=8, help='Number of threads in the thread pool (default: 16)')

    args = parser.parse_args()
    run_mrna_design(args)

if __name__ == "__main__":
    main()