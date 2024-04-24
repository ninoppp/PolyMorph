import pandas as pd
import matplotlib.pyplot as plt

ensemble = pd.read_csv('benchmark_ensemble.csv', sep='\s+')
solver128 = pd.read_csv('benchmark_solver_128.csv', sep='\s+')
solver256 = pd.read_csv('benchmark_solver_256.csv', sep='\s+')
full_sim = pd.read_csv('benchmark_full_simulation.csv', sep='\s+')

# ensemble.step() runtime vs total vertices
ensemble_1_thread = ensemble[ensemble['num_threads'] == 1]
ensemble_2_thread = ensemble[ensemble['num_threads'] == 2]
plt.figure(figsize=(10, 5))
plt.scatter(ensemble_1_thread['total_vertices'], ensemble_1_thread['time'], label='1 Thread')
plt.scatter(ensemble_2_thread['total_vertices'], ensemble_2_thread['time'], label='2 Threads')
plt.title('Frame Time vs Total Vertices')
plt.xlabel('Total Vertices')
plt.ylabel('Time (seconds)')
plt.grid(True)
plt.legend()
plt.savefig("ensemble_runtime.pdf")

# ensemble.step() scaling
ensemble['time_per_vertex'] = ensemble['time'] / ensemble['total_vertices']
avg_ratio_by_thread = ensemble.groupby('num_threads')['time_per_vertex'].mean().reset_index()
plt.figure(figsize=(10, 5))
plt.plot(avg_ratio_by_thread['num_threads'], avg_ratio_by_thread['time_per_vertex'], marker='o')
plt.title('Average Frame Time per Vertex vs Number of Threads. Starting with 1 cell')
plt.xlabel('Number of Threads')
plt.ylabel('Time (seconds)')
plt.xscale('log', base=2)
plt.yscale('log')
plt.grid(True)
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