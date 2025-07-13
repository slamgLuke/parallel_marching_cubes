import numpy as np
import matplotlib.pyplot as plt
import itertools

# Número de threads utilizados en los experimentos
threads = [1, 2, 4, 8, 16, 32, 64, 128]

# Tiempos de ejecución para la implementación con bucles 'for'
for_times = {
    "N=1024": [192.849, 140.965, 97.1488, 70.6403, 62.0982, 60.8236, 64.5615, 62.2575],
    "N=512": [25.23, 19.0662, 12.6847, 10.0845, 8.44556, 7.5536, 7.98846, 7.74935],
    "N=256": [3.8631, 2.63748, 1.96705, 1.25413, 1.09923, 1.06948, 1.0931, 1.07937],
    "N=128": [0.8688, 0.688, 0.512, 0.35, 0.2967, 0.301, 0.219, 0.1456],
    "N=64": [0.1422, 0.0989832, 0.0700724, 0.0468122, 0.0399077, 0.0400399, 0.0425398, 0.0540],
    "N=32": [0.0292, 0.015, 0.01099, 0.0819, 0.0689, 0.073, 0.01152, 0.02198],
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
for N_key, color in zip(for_times.keys(), colors):
    N_val = int(N_key.split('=')[1])
    total_flops = flops_data[N_val]

    # Calcular GFLOP/s para la versión 'for'
    gflops_for = [(total_flops / t) / 1e9 if t >
                  1e-6 else 0 for t in for_times[N_key]]

    # Graficar los resultados
    plt.plot(threads, gflops_for, 'o-',
             label=f'{N_key}', color=color, markersize=7)

# Configuración de los ejes y etiquetas
plt.xscale('log', base=2)
plt.xlabel('Número de Threads (p)')
plt.ylabel('Rendimiento (GFLOP/s)')
plt.title("Rendimiento (GFLOP/s) para Implementación con 'For'")
plt.legend(title='Tamaño del Problema',
           bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True, which="both", linestyle='--', alpha=0.7)
plt.xticks(threads, labels=[str(p) for p in threads])
plt.tight_layout()
plt.show()
