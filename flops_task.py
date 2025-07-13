import numpy as np
import matplotlib.pyplot as plt
import itertools

# Número de threads utilizados en los experimentos
threads = [1, 2, 4, 8, 16, 32, 64, 128]

# Tiempos de ejecución para la implementación con 'tasks'
task_times = {
    "N=1024": [342.322, 171.763, 87.5689, 48.4029, 27.4464, 18.3951, 13.1698, 10.9407],
    "N=512": [44.7552, 22.8137, 11.9034, 6.17732, 3.87999, 2.67343, 1.80716, 1.65459],
    "N=256": [7.16534, 3.52782, 2.27189, 1.22987, 0.704408, 0.449432, 0.349796, 0.302823],
    "N=128": [1.44816, 0.731304, 0.376581, 0.194846, 0.12436, 0.0749535, 0.0673379, 0.0733135],
    "N=64": [0.203752, 0.111712, 0.0603501, 0.0344061, 0.0219477, 0.016642, 0.0153573, 0.0215591],
    "N=32": [0.0373969, 0.021276, 0.0115092, 0.00765573, 0.0050895, 0.00420818, 0.00805823, 0.0107166],
}

# FLOPs totales calculados para cada tamaño de N
flops_data = {
    32: 14849032,
    64: 118252136,
    128: 943810888,
    256: 7541676584,
    512: 60298644232,
    1024: 482249623496,
}

# Configuración del gráfico
plt.figure(figsize=(12, 7))
colors = itertools.cycle(plt.cm.tab10.colors)

# Iterar sobre cada tamaño de problema N para graficar sus resultados
for N_key, color in zip(task_times.keys(), colors):
    N_val = int(N_key.split('=')[1])
    total_flops = flops_data[N_val]

    # Calcular GFLOP/s para la versión 'tasks'
    gflops_tasks = [(total_flops / t) / 1e9 if t >
                    1e-6 else 0 for t in task_times[N_key]]

    # Graficar los resultados
    plt.plot(threads, gflops_tasks, 'x--',
             label=f'{N_key}', color=color, markersize=7)

# Configuración de los ejes y etiquetas
plt.xscale('log', base=2)
plt.xlabel('Número de Threads (p)')
plt.ylabel('Rendimiento (GFLOP/s)')
plt.title("Rendimiento (GFLOP/s) para Implementación con 'Tasks'")
plt.legend(title='Tamaño del Problema',
           bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True, which="both", linestyle='--', alpha=0.7)
plt.xticks(threads, labels=[str(p) for p in threads])
plt.tight_layout()
plt.show()
