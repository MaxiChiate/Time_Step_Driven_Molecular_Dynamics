#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from pathlib import Path
import argparse

def load_frames(filename: str):
    times, frames = [], []
    with open(filename, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    i = 0
    while i < len(lines):
        if "," not in lines[i]:
            t = float(lines[i]); i += 1
            rows = []
            while i < len(lines) and "," in lines[i]:
                parts = lines[i].split(",")
                pid = int(parts[0])
                x, y, z = map(float, parts[1:4])
                vx, vy, vz = map(float, parts[4:7])
                m = float(parts[7])
                rows.append([pid, x, y, z, vx, vy, vz, m])
                i += 1
            if rows:
                frames.append(np.array(rows, dtype=float))
                times.append(t)
        else:
            i += 1
    return np.array(times), frames

def animate_particles(filename, out=None, fps=30, interval=50, mode="single"):
    times, frames = load_frames(filename)

    # Rango de ejes
    all_x = [p[1] for fr in frames for p in fr]
    all_y = [p[2] for fr in frames for p in fr]
    all_z = [p[3] for fr in frames for p in fr]
    margin = 0.1 * max(
        max(all_x) - min(all_x),
        max(all_y) - min(all_y),
        max(all_z) - min(all_z),
        )
    xmin, xmax = min(all_x) - margin, max(all_x) + margin
    ymin, ymax = min(all_y) - margin, max(all_y) + margin
    zmin, zmax = min(all_z) - margin, max(all_z) + margin

    # Figura 3D
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d")
    ax.set_facecolor("black")
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_zlim(zmin, zmax)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

    scat = ax.scatter([], [], [], s=5)
    time_text = ax.text2D(0.02, 0.95, "", transform=ax.transAxes, color="green")

    def init():
        scat._offsets3d = ([], [], [])
        time_text.set_text("")
        return scat, time_text

    def update(i):
        data = frames[i]
        xs = data[:, 1]; ys = data[:, 2]; zs = data[:, 3]
        ids = data[:, 0].astype(int)

        if mode == "clusters":
            # N se puede inferir como la mitad de la cantidad de partículas
            N = len(ids) // 2
            colors = ["pink" if pid < N else "blue" for pid in ids]
        else:
            colors = ["white"] * len(ids)

        scat._offsets3d = (xs, ys, zs)
        scat.set_color(colors)
        time_text.set_text(f"t = {times[i]:.2f}")
        return scat, time_text

    ani = animation.FuncAnimation(
        fig, update, frames=len(frames),
        init_func=init, interval=interval, blit=False
    )

    if out:
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        ani.save(out, writer="ffmpeg", fps=fps)
        print(f"Animación guardada en {out}")
    else:
        plt.show()

def main():
    ap = argparse.ArgumentParser(description="Visualizador 3D interactivo de simulaciones SDS4.")
    ap.add_argument("file", help="Archivo de simulación (.csv)")
    ap.add_argument("--out", help="Archivo de salida .mp4 (opcional)")
    ap.add_argument("--fps", type=int, default=30, help="FPS del video (default=30)")
    ap.add_argument("--interval", type=int, default=50, help="Intervalo entre frames en ms (default=50)")
    ap.add_argument("--mode", choices=["single", "clusters"], default="single",
                    help="Modo de visualización: single (blanco) o clusters (rosa/azul)")
    args = ap.parse_args()

    animate_particles(args.file, out=args.out, fps=args.fps, interval=args.interval, mode=args.mode)

if __name__ == "__main__":
    main()
