#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
import re

#### USAGE:
# python3 rhm_plot_shadow.py --base MolecularDynamic/outputs/single/gear_order_5 --save

def load_analysis_file(filename: str):
    """Carga archivo de análisis con columnas: time, energy, rhm, tStar"""
    try:
        data = np.genfromtxt(filename, delimiter=",", names=["time", "energy", "rhm", "tStar"])
        return data["time"], data["rhm"]
    except Exception:
        return None, None


def analyze_directory(path: Path):
    """Promedia los archivos *_analysis.csv de una carpeta N*"""
    files = sorted(path.glob("*_analysis.csv"))
    if not files:
        print(f"[WARN] {path}: no hay archivos *_analysis.csv")
        return None, None, None

    rhm_all = []
    times_ref = None

    for f in files:
        t, rhm = load_analysis_file(f)
        if t is None or rhm is None or len(rhm) == 0:
            continue

        if times_ref is None:
            times_ref = t
            rhm_all.append(rhm)
        else:
            # Interpolar para alinear tiempos
            tmin = max(t[0], times_ref[0])
            tmax = min(t[-1], times_ref[-1])
            mask_ref = (times_ref >= tmin) & (times_ref <= tmax)
            mask_t = (t >= tmin) & (t <= tmax)
            if np.sum(mask_ref) == 0:
                continue
            rhm_interp = np.interp(times_ref[mask_ref], t[mask_t], rhm[mask_t])
            times_ref = times_ref[mask_ref]
            rhm_all = [rhm_i[:len(times_ref)] for rhm_i in rhm_all]
            rhm_all.append(rhm_interp)

    if len(rhm_all) == 0:
        print(f"[WARN] {path}: no se pudo leer ningún archivo válido")
        return None, None, None

    rhm_all = np.vstack(rhm_all)
    mean_rhm = np.mean(rhm_all, axis=0)
    std_rhm = np.std(rhm_all, axis=0)
    return times_ref, mean_rhm, std_rhm


def main():
    ap = argparse.ArgumentParser(description="Graficar rhm(t) con banda de desvío estándar")
    ap.add_argument("--base", required=True, help="Directorio base con subcarpetas N####")
    ap.add_argument("--save", action="store_true", help="Guardar gráficos en images/")
    args = ap.parse_args()

    base = Path(args.base)
    out_dir = Path("images")
    if args.save:
        out_dir.mkdir(exist_ok=True)

    for d in sorted(base.glob("N*")):
        m = re.match(r"N(\d+)", d.name)
        if not m:
            continue
        N = int(m.group(1))

        times, mean_rhm, std_rhm = analyze_directory(d)
        if times is None:
            continue

        plt.figure(figsize=(7, 5))
        color = "C0"
        plt.plot(times, mean_rhm, color=color, lw=1, label="Datos")
        plt.fill_between(times, mean_rhm - std_rhm, mean_rhm + std_rhm,
                         color=color, alpha=0.25, label="Desvío estándar")

        plt.xlabel("Tiempo")
        plt.ylabel("Radio de media masa")
        plt.grid(True, alpha=0.4)
        plt.legend()
        plt.tight_layout()

        if args.save:
            plt.savefig(out_dir / f"rhm_shadow_N{N}.png", dpi=150)
            plt.close()
            print(f"[OK] Guardado: images/rhm_shadow_N{N}.png")
        else:
            plt.show()


if __name__ == "__main__":
    main()
