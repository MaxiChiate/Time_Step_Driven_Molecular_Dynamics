#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
import re

def load_analysis_file(filename: str):
    """Carga archivo de análisis con columnas: time, energy, rhm, tStar"""
    data = np.genfromtxt(filename, delimiter=",", names=["time", "energy", "rhm", "tStar"])
    return data["time"], data["rhm"]

def analyze_directory(path: Path):
    """Promedia los archivos *_analysis.csv de una carpeta N*"""
    files = sorted(path.glob("*_analysis.csv"))
    if not files:
        return None, None

    rhm_all = []
    times_ref = None
    for f in files:
        t, rhm = load_analysis_file(f)
        if times_ref is None:
            times_ref = t
            rhm_all.append(rhm)
        else:
            tmin, tmax = max(t[0], times_ref[0]), min(t[-1], times_ref[-1])
            mask_ref = (times_ref >= tmin) & (times_ref <= tmax)
            mask_t = (t >= tmin) & (t <= tmax)
            rhm_interp = np.interp(times_ref[mask_ref], t[mask_t], rhm[mask_t])
            times_ref = times_ref[mask_ref]
            rhm_all = [rhm_i[:len(times_ref)] for rhm_i in rhm_all]
            rhm_all.append(rhm_interp)

    all_rhm = np.vstack(rhm_all)
    mean_rhm = np.mean(all_rhm, axis=0)
    return times_ref, mean_rhm

def slope_stationary(times, mean_rhm, frac=0.7):
    """Calcula la pendiente lineal en el último 30% del tiempo total."""
    start = int(len(times) * frac)
    coeffs = np.polyfit(times[start:], mean_rhm[start:], 1)
    return coeffs[0]

def main():
    ap = argparse.ArgumentParser(description="Calcular pendiente estacionaria de rhm(t) vs N")
    ap.add_argument("--base", required=True, help="Directorio base con subcarpetas N####")
    ap.add_argument("--save", action="store_true", help="Guardar CSV y gráfico en images/")
    args = ap.parse_args()

    base = Path(args.base)
    out_dir = Path("images")
    if args.save:
        out_dir.mkdir(exist_ok=True)

    Ns, slopes = [], []

    for d in sorted(base.glob("N*")):
        m = re.match(r"N(\d+)", d.name)
        if not m:
            continue
        N = int(m.group(1))

        times, mean_rhm = analyze_directory(d)
        if times is None:
            continue

        slope = slope_stationary(times, mean_rhm, frac=0.7)
        Ns.append(N)
        slopes.append(slope)
        print(f"N={N:4d}  pendiente estacionaria = {slope:.4e}")

    if not Ns:
        print("No se encontraron resultados válidos.")
        return

    Ns, slopes = np.array(Ns), np.array(slopes)
    plt.figure(figsize=(7, 5))
    plt.plot(Ns, slopes, "o-", color="C0")
    plt.xlabel("Número de partículas")
    plt.ylabel("Pendiente del rhm")
    plt.grid(True, alpha=0.4)
    plt.tight_layout()

    if args.save:
        plt.savefig(out_dir / "slope_vs_N.png", dpi=150)
        with open(out_dir / "slopes.csv", "w") as f:
            f.write("N,slope\n")
            for n, s in zip(Ns, slopes):
                f.write(f"{n},{s}\n")
        print("[OK] Guardados: images/slope_vs_N.png y images/slopes.csv")
    else:
        plt.show()

if __name__ == "__main__":
    main()
