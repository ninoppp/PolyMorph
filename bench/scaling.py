import pandas as pd
import matplotlib.pyplot as plt

ensemble = pd.read_csv('benchmark_ensemble.csv', sep='\s+')
solver128 = pd.read_csv('benchmark_solver_128.csv', sep='\s+')
solver256 = pd.read_csv('benchmark_solver_256.csv', sep='\s+')
full_sim = pd.read_csv('benchmark_full_simulation.csv', sep='\s+')
ensemble_box = pd.read_csv('benchmark_ensemble_box.csv', sep='\s+')

P = 0.93 # reference: Polyhoop paper
def amdahl(t):
    return 1.0 / (1.0 - P + P / t)

# ensemble.step() runtime vs total vertices
ensemble_1_thread = ensemble[ensemble['num_threads'] == 1]
ensemble_2_thread = ensemble[ensemble['num_threads'] == 2]
ensemble_box_1_thread = ensemble_box[ensemble_box['num_threads'] == 1]
ensemble_box_2_thread = ensemble_box[ensemble_box['num_threads'] == 2]
plt.figure(figsize=(10, 5))
plt.scatter(ensemble_1_thread['num_vertices'], ensemble_1_thread['time'], label='1 Thread')
plt.scatter(ensemble_2_thread['num_vertices'], ensemble_2_thread['time'], label='2 Threads')
plt.scatter(ensemble_box_1_thread['num_vertices'], ensemble_box_1_thread['time'], label='1 Thread, box')
plt.scatter(ensemble_box_2_thread['num_vertices'], ensemble_box_2_thread['time'], label='2 Threads, box')
plt.title('Frame Time vs Total Vertices (single cell)')
plt.xlabel('Total Vertices')
plt.ylabel('Time (seconds)')
plt.xscale('log')
plt.yscale('log')
plt.grid(True)
plt.legend()
plt.savefig("ensemble_runtime.pdf")

# ensemble.step() scaling
threads = ensemble['num_threads'].unique()
ensemble['time_per_vertex'] = ensemble['time'] / ensemble['num_vertices']
ensemble_box['time_per_vertex'] = ensemble_box['time'] / ensemble_box['num_vertices']
ensemble_box_filtered = ensemble_box[ensemble_box['num_vertices'] > 50000]
ensemble_box_filtered['time_per_vertex'] = ensemble_box_filtered['time'] / ensemble_box_filtered['num_vertices']
avg_ratio_by_thread = ensemble.groupby('num_threads')['time_per_vertex'].mean().reset_index()
avg_ratio_by_thread_box = ensemble_box.groupby('num_threads')['time_per_vertex'].mean().reset_index()
avg_ratio_by_thread_box_filtered = ensemble_box_filtered.groupby('num_threads')['time_per_vertex'].mean().reset_index()
speedup = avg_ratio_by_thread['time_per_vertex'][0] / avg_ratio_by_thread['time_per_vertex']
speedup_box = avg_ratio_by_thread_box['time_per_vertex'][0] / avg_ratio_by_thread_box['time_per_vertex']
speedup_box_filtered = avg_ratio_by_thread_box_filtered['time_per_vertex'][0] / avg_ratio_by_thread_box_filtered['time_per_vertex']
speedup_reference = amdahl(avg_ratio_by_thread['num_threads'])
plt.figure(figsize=(10, 5))
plt.plot(threads, speedup, marker='o', label='Single cell, Ns=1e3')
plt.plot(threads, speedup_box, marker='o', label='Single cell + box, Ns=3e3')
plt.plot(threads, speedup_box_filtered, marker='o', label='Single cell + box, only 50k+ vertices')
plt.plot(threads, speedup_reference, label='Amdahl\'s Law, P=93%')
plt.title('Average Frame Time per Vertex vs Number of Threads')
plt.xlabel('Number of Threads')
plt.ylabel('Speedup S = T(1) / T(p)')
plt.xscale('log', base=2)
plt.yscale('log')
plt.grid(True)
plt.legend()
plt.savefig("ensemble_scaling.pdf")

# solver.step() scaling
avg_runtime_by_thread_128 = solver128.groupby('num_threads')['time'].mean().reset_index()
avg_runtime_by_thread_256 = solver256.groupby('num_threads')['time'].mean().reset_index()
plt.figure(figsize=(10, 5))
plt.plot(avg_runtime_by_thread_128['num_threads'], avg_runtime_by_thread_128['time'], marker='o', label='128x128 grid')
plt.plot(avg_runtime_by_thread_256['num_threads'], avg_runtime_by_thread_256['time'], marker='o', label='256x256 grid')
plt.title('Solver Average Frame Time vs Number of Threads. 256x256 grid')
plt.xlabel('Number of Threads')
plt.ylabel('Average Time per Vertex (seconds)')
plt.xscale('log', base=2)
plt.yscale('log')
plt.grid(True)
plt.legend()
plt.savefig("solver_scaling.pdf")

