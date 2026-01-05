#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FFMpegWriter
from pathlib import Path
import argparse
import re


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


def infer_output_path(csv_path: Path) -> Path:
    parts = [p.lower() for p in csv_path.parts]
    mode = "single" if any("single" in p for p in parts) else ("clusters" if any("cluster" in p for p in parts) else "mode")
    scheme_candidates = ["verlet", "beeman", "gear", "euler", "leapfrog", "velocity-verlet"]
    scheme = next((s for s in scheme_candidates if any(s in p for p in parts)), "sim")

    name = csv_path.name
    mN  = re.search(r"N_?(\d+)", name, re.I)
    mdt = re.search(r"dt([0-9]*\.?[0-9]+(?:e-?\d+)?)", name, re.I)
    mt  = re.search(r"_t([0-9]*\.?[0-9]+(?:e-?\d+)?)_", name, re.I)
    mver= re.search(r"_([0-9]{3,})\.csv$", name, re.I)

    N   = mN.group(1)  if mN  else "X"
    dt  = mdt.group(1) if mdt else "X"
    t   = mt.group(1)  if mt  else "X"
    ver = mver.group(1)if mver else "0000"

    outdir = Path("animations"); outdir.mkdir(exist_ok=True)
    return outdir / f"{scheme}_{mode}_N{N}_dt{dt}_t{t}_{ver}.mp4"


def animate_particles(filename, out=None, fps=30, interval=50, mode="single", export_mp4=False):
    times, frames = load_frames(filename)

    # Rango de ejes
    all_x = [p[1] for fr in frames for p in fr]
    all_y = [p[2] for fr in frames for p in fr]
    all_z = [p[3] for fr in frames for p in fr]
    span = max(
        (max(all_x) - min(all_x)) if all_x else 1.0,
        (max(all_y) - min(all_y)) if all_y else 1.0,
        (max(all_z) - min(all_z)) if all_z else 1.0,
    )
    margin = 0.10 * span
    xmin, xmax = min(all_x) - margin, max(all_x) + margin
    ymin, ymax = min(all_y) - margin, max(all_y) + margin
    zmin, zmax = min(all_z) - margin, max(all_z) + margin

    # Figura HD 1080p
    fig = plt.figure(figsize=(12.8, 7.2), dpi=150)
    fig.patch.set_facecolor("black")
    ax = fig.add_subplot(111, projection="3d")
    ax.set_facecolor("black")

    # Límites
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_zlim(zmin, zmax)

    # Ejes visibles (sin ticks)
    ax.set_xlabel("X", color="white", labelpad=8)
    ax.set_ylabel("Y", color="white", labelpad=8)
    ax.set_zlabel("Z", color="white", labelpad=8)
    ax.tick_params(colors="white", labelsize=0, length=0)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.fill = False
    ax.xaxis.line.set_color("white")
    ax.yaxis.line.set_color("white")
    ax.zaxis.line.set_color("white")
    ax.grid(False)
    plt.subplots_adjust(0, 0, 1, 1)

    # Scatter inicial (sin depthshade)
    scat = ax.scatter([], [], [], s=10, c="white", marker=".", edgecolors="none", depthshade=False)
    time_text = ax.text2D(0.02, 0.95, "", transform=ax.transAxes, color="#7CFC00")

    def init():
        scat._offsets3d = ([], [], [])
        time_text.set_text("")
        return scat, time_text

    def update(i):
        data = frames[i]
        xs, ys, zs = data[:, 1], data[:, 2], data[:, 3]
        ids = data[:, 0].astype(int)

        # --- Zoom fijo fuerte sobre el inicio ---
        # centrado en el centro inicial y con recorte ajustado
        cx, cy, cz = np.mean(frames[0][:, 1]), np.mean(frames[0][:, 2]), np.mean(frames[0][:, 3])
        zoom_factor = 0.1  # cuanto menor, más zoom (0.3 ≈ fuerte, 1.0 = normal)

        dx = (xmax - xmin) * zoom_factor / 2
        dy = (ymax - ymin) * zoom_factor / 2
        dz = (zmax - zmin) * zoom_factor / 2

        ax.set_xlim(cx - dx, cx + dx)
        ax.set_ylim(cy - dy, cy + dy)
        ax.set_zlim(cz - dz, cz + dz)

        # --- Profundidad y color ---
        z_norm = (zs - zmin) / (zmax - zmin + 1e-9)
        sizes = 3 + 25 * (1 - z_norm)**1.5
        brightness = 0.3 + 0.7 * (1 - z_norm)**1.2

        if mode == "clusters":
            Nhalf = len(ids) // 2
            base_colors = np.array([[1.0, 0.31, 0.64] if pid < Nhalf else [0.25, 0.65, 1.0] for pid in ids])
        else:
            base_colors = np.ones((len(ids), 3))

        colors = np.clip(base_colors * brightness[:, None], 0, 1)

        # --- Actualizar scatter ---
        scat._offsets3d = (xs, ys, zs)
        scat.set_color(colors)
        scat.set_sizes(sizes)

        # --- Texto tiempo ---
        time_text.set_text(f"t = {times[i]:.2f}")
        return scat, time_text





    ani = animation.FuncAnimation(
        fig, update, frames=len(frames),
        init_func=init, interval=interval, blit=False
    )

    if export_mp4 or out:
        out_path = Path(out) if out else infer_output_path(Path(filename))
        out_path.parent.mkdir(parents=True, exist_ok=True)

        writer = FFMpegWriter(
            fps=fps,
            bitrate=16000,
            codec="libx264",
            extra_args=["-pix_fmt", "yuv420p"]
        )
        ani.save(
            out_path,
            writer=writer,
            savefig_kwargs={"facecolor": "black", "bbox_inches": "tight", "pad_inches": 0}
        )
        print(f"Animación guardada en {out_path}")
    else:
        plt.show()


def main():
    ap = argparse.ArgumentParser(description="Visualizador 3D de simulaciones SDS4.")
    ap.add_argument("file", help="Archivo de simulación (.csv)")
    ap.add_argument("--out", help="Archivo de salida .mp4 (opcional)")
    ap.add_argument("--mp4", action="store_true", help="Exportar animación a .mp4 con nombre deducido del .csv")
    ap.add_argument("--fps", type=int, default=30, help="FPS del video (default=30)")
    ap.add_argument("--interval", type=int, default=50, help="Intervalo entre frames en ms (default=50)")
    ap.add_argument("--mode", choices=["single", "clusters"], default="single",
                    help="Modo de visualización: single (blanco) o clusters (rosa/azul)")
    args = ap.parse_args()

    export_mp4 = bool(args.mp4 or args.out)
    animate_particles(
        args.file,
        out=args.out,
        fps=args.fps,
        interval=args.interval,
        mode=args.mode,
        export_mp4=export_mp4
    )


if __name__ == "__main__":
    main()
