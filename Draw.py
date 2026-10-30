import pandas as pd
import matplotlib.pyplot as plt
import os

def plot_benchmark_results():
    files = {
        'Binet': 'binet_results.csv',
        'Strassen': 'strassen_results.csv',
        'Laderman': 'laderman_results.csv'
    }

    for name, filepath in files.items():
        if not os.path.exists(filepath):
            print(f"Błąd: Nie znaleziono pliku {filepath}")
            return

    df_binet = pd.read_csv(files['Binet'])
    df_strassen = pd.read_csv(files['Strassen'])
    df_laderman = pd.read_csv(files['Laderman'])


    plt.figure(figsize=(10, 6))
    plt.plot(df_binet['Size'], df_binet['AvgTime_ms'], label='Binet', color='blue', linewidth=2)
    plt.plot(df_strassen['Size'], df_strassen['AvgTime_ms'], label='Strassen', color='orange', linewidth=2)
    plt.plot(df_laderman['Size'], df_laderman['AvgTime_ms'], label='Laderman', color='green', linewidth=2)

    plt.title('Czas wykonania algorytmów mnożenia macierzy', fontsize=14, fontweight='bold')
    plt.xlabel('Rozmiar macierzy (N)', fontsize=12)
    plt.ylabel('Średni czas (ms)', fontsize=12)
    plt.legend(fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()

    plt.savefig('img/wykres_czas.png', dpi=300)
    plt.show()


    plt.figure(figsize=(10, 6))
    plt.plot(df_binet['Size'], df_binet['Flops'], label='Binet', color='blue', linewidth=2)
    plt.plot(df_strassen['Size'], df_strassen['Flops'], label='Strassen', color='orange', linewidth=2)
    plt.plot(df_laderman['Size'], df_laderman['Flops'], label='Laderman', color='green', linewidth=2)

    plt.title('Złożoność obliczeniowa (FLOPs)', fontsize=14, fontweight='bold')
    plt.xlabel('Rozmiar macierzy (N)', fontsize=12)
    plt.ylabel('Ilość operacji (FLOPs)', fontsize=12)
    plt.legend(fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.7)

    # Skala logarytmiczna dla osi Y jest świetna do pokazywania różnic rzędu wielkości
    # Opcjonalnie odkomentuj poniższą linię, jeśli krzywe zlewają się na dole:
    # plt.yscale('log')

    plt.tight_layout()

    # Zapis do pliku i wyświetlenie
    plt.savefig('img/wykres_flops.png', dpi=300)
    plt.show()

if __name__ == "__main__":
    plot_benchmark_results()