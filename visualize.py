#!/usr/bin/env python3
"""
1D Electromagnetic Wave Visualizer
Plays the generated simulation frames as an animated graph.
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


def load_data(filepath: str = "output.csv"):
    """Load simulation data from CSV (tabular or matrix format)."""
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {filepath}")

    with open(path, "r", encoding="utf-8") as f:
        first_line = f.readline().strip().lower()

    # Tabular CSV with headers: timestep,cell,Ex,Hy
    if "cell" in first_line or "ex" in first_line:
        raw = np.loadtxt(path, delimiter=",", skiprows=1)
        num_steps = int(raw[:, 0].max()) + 1
        num_cells = int(raw[:, 1].max()) + 1
        ex = raw[:, 2].reshape(num_steps, num_cells)
        hy = raw[:, 3].reshape(num_steps, num_cells)
    else:
        # Matrix format (each row = 1 timestep, each col = 1 cell)
        ex = np.loadtxt(path, delimiter=",")
        if ex.ndim == 1:
            ex = ex.reshape(1, -1)
        hy = np.zeros_like(ex)

    return ex, hy


def play_simulation(ex: np.ndarray, hy: np.ndarray, interval: int = 20):
    """Plays the frames as an animated 2-panel graph (Ex and Hy)."""
    num_steps, num_cells = ex.shape
    cells = np.arange(num_cells)

    fig, (ax_e, ax_h) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)

    e_max = max(float(np.max(np.abs(ex))), 1e-6) * 1.15
    h_max = max(float(np.max(np.abs(hy))), 1e-9) * 1.15

    # Electric Field subplot
    ax_e.set_xlim(0, num_cells - 1)
    ax_e.set_ylim(-e_max, e_max)
    ax_e.set_ylabel("Electric Field $E_x$ (V/m)", fontweight="bold")
    ax_e.grid(True, linestyle="--", alpha=0.6)
    (line_e,) = ax_e.plot(cells, ex[0], color="#007acc", lw=1.8, label="$E_x$")
    ax_e.legend(loc="upper right")

    # Magnetic Field subplot
    ax_h.set_xlim(0, num_cells - 1)
    ax_h.set_ylim(-h_max, h_max)
    ax_h.set_xlabel("Grid Cell Index", fontweight="bold")
    ax_h.set_ylabel("Magnetic Field $H_y$ (A/m)", fontweight="bold")
    ax_h.grid(True, linestyle="--", alpha=0.6)
    (line_h,) = ax_h.plot(cells, hy[0], color="#d9534f", lw=1.8, label="$H_y$")
    ax_h.legend(loc="upper right")

    title = fig.suptitle(f"Timestep: 0 / {num_steps}", fontsize=12, fontweight="bold")

    is_paused = False

    def update(frame):
        line_e.set_ydata(ex[frame])
        line_h.set_ydata(hy[frame])
        title.set_text(f"Timestep: {frame} / {num_steps}")
        return line_e, line_h, title

    anim = FuncAnimation(fig, update, frames=num_steps, interval=interval)

    def on_key(event):
        nonlocal is_paused
        if event.key == " ":
            if is_paused:
                anim.resume()
            else:
                anim.pause()
            is_paused = not is_paused

    fig.canvas.mpl_connect("key_press_event", on_key)
    plt.tight_layout()
    plt.show()


def main():
    filepath = sys.argv[1] if len(sys.argv) > 1 else "output.csv"
    try:
        ex, hy = load_data(filepath)
    except Exception as err:
        print(f"Error loading '{filepath}': {err}", file=sys.stderr)
        sys.exit(1)

    print(f"Playing '{filepath}' ({ex.shape[0]} timesteps, {ex.shape[1]} cells)...")
    print("Press [Spacebar] in the window to Play/Pause.")
    play_simulation(ex, hy)


if __name__ == "__main__":
    main()
