#!/usr/bin/env python3
from __future__ import annotations

"""
1D Electromagnetic Wave Simulation Visualizer
============================================
Visualizes 1D FDTD electromagnetic wave simulation results:
- X-axis: Spatial Grid Cells (Cell Index: 0, 1, 2, ..., N-1)
- Y-axis: Magnitude / Amplitude of the Electric (Ex) and Magnetic (Hy) Fields

Supports:
- Live/interactive wave propagation animation (with Play/Pause & Step controls)
- Single-plot overlay view (--overlay) with dual y-axes for Ex and Hy
- Absolute magnitude mode (--abs) to plot |E| and |H|
- Spatiotemporal heatmaps / waterfall diagrams (Cells vs Timesteps)
- Multi-step snapshot comparison plots across cells
- Exporting to MP4, GIF, and PNG
- Auto-loading CSV, NumPy (.npz/.npy), and JSON formats
- Built-in FDTD demo mode (--demo) for instant testing

Dependencies:
    pip install numpy matplotlib

Usage Examples:
    # 1. Run built-in demo simulation (X: Cells, Y: Field Magnitude)
    python3 visualize.py --demo

    # 2. Plot absolute magnitude |E| and |H|
    python3 visualize.py --demo --abs

    # 3. Overlay both E and H fields on a single plot
    python3 visualize.py --demo --overlay

    # 4. Spatiotemporal heatmap (X: Cells, Y: Timesteps)
    python3 visualize.py --demo --mode heatmap

    # 5. Progression snapshots across cells
    python3 visualize.py --demo --mode snapshots

    # 6. Load results from simulation file (e.g. CSV or NPZ)
    python3 visualize.py --input sim_output.csv
"""

import argparse
import json
import os
import sys
from pathlib import Path
from dataclasses import dataclass
from typing import Optional, Tuple

# Physical Constants
C0 = 299_792_458.0  # Speed of light in vacuum (m/s)
MU0 = 4e-7 * 3.141592653589793  # Vacuum permeability (H/m)
EPS0 = 1.0 / (MU0 * C0 * C0)  # Vacuum permittivity (F/m)

np = None
plt = None
animation = None


def ensure_dependencies():
    """Validates and imports required visualization libraries."""
    global np, plt, animation
    missing = []
    try:
        import numpy as _np
        np = _np
    except ImportError:
        missing.append("numpy")

    try:
        import matplotlib.pyplot as _plt
        import matplotlib.animation as _animation
        plt = _plt
        animation = _animation
    except ImportError:
        missing.append("matplotlib")

    if missing:
        print(f"\n[Error] Missing required packages: {', '.join(missing)}")
        print("Please install them using:")
        print("    pip install -r requirements.txt\n  or:\n    pip install " + " ".join(missing) + "\n")
        sys.exit(1)


@dataclass
class SimulationData:
    """Container for 1D EM simulation results."""
    ex: np.ndarray  # Shape: (num_timesteps, num_cells)
    hy: Optional[np.ndarray] = None  # Shape: (num_timesteps, num_cells)
    x: Optional[np.ndarray] = None  # Shape: (num_cells,) in meters
    t: Optional[np.ndarray] = None  # Shape: (num_timesteps,) in seconds
    dx: float = 1.0
    dt: float = 1.0
    frequency: Optional[float] = None
    title: str = "1D EM Wave Simulation"

    def __post_init__(self):
        num_timesteps, num_cells = self.ex.shape
        if self.x is None:
            self.x = np.arange(num_cells) * self.dx
        if self.t is None:
            self.t = np.arange(num_timesteps) * self.dt
        if self.hy is None:
            self.hy = np.zeros_like(self.ex)


def generate_demo_simulation(
    num_cells: int = 300,
    num_steps: int = 500,
    frequency: float = 2.4e9,
    source_type: str = "gaussian"
) -> SimulationData:
    """
    Generates a 1D FDTD simulation of an electromagnetic wave pulse.
    Demonstrates pulse injection, free-space propagation, and boundary interaction.
    """
    wavelength = C0 / frequency
    dx = wavelength / 20.0  # 20 grid cells per wavelength
    courant_factor = 0.99   # Courant stability limit for 1D FDTD is S <= 1.0
    dt = (courant_factor * dx) / C0

    # Initialize field arrays
    ex = np.zeros(num_cells, dtype=np.float64)
    hy = np.zeros(num_cells, dtype=np.float64)

    # History buffers to store data over time
    ex_history = np.zeros((num_steps, num_cells), dtype=np.float64)
    hy_history = np.zeros((num_steps, num_cells), dtype=np.float64)

    # Source location (center of the domain)
    source_pos = num_cells // 4
    pulse_width = 25 * dt
    t0 = 40 * dt

    # Boundary condition buffers (Simple Mur ABC)
    ex_left_prev = 0.0
    ex_right_prev = 0.0

    print(f"Running Demo 1D FDTD Simulation:")
    print(f"  - Grid cells : {num_cells} cells")
    print(f"  - Timesteps  : {num_steps} steps")
    print(f"  - dx         : {dx * 1e3:.3f} mm/cell")
    print(f"  - dt         : {dt * 1e12:.3f} ps/step")
    print(f"  - Frequency  : {frequency / 1e9:.2f} GHz")

    for n in range(num_steps):
        # 1. Update Magnetic Field Hy
        hy[:-1] += (dt / (MU0 * dx)) * (ex[1:] - ex[:-1])

        # 2. Update Electric Field Ex
        ex[1:] += (dt / (EPS0 * dx)) * (hy[1:] - hy[:-1])

        # 3. Source Injection
        current_time = n * dt
        if source_type == "gaussian":
            pulse = np.exp(-((current_time - t0) / pulse_width) ** 2)
            ex[source_pos] += pulse
        elif source_type == "sine":
            if current_time >= 0:
                ex[source_pos] += np.sin(2.0 * np.pi * frequency * current_time)

        # 4. Simple First-Order Mur ABC at domain boundaries
        c_factor = (C0 * dt - dx) / (C0 * dt + dx)
        ex[0] = ex_left_prev + c_factor * (ex[1] - ex[0])
        ex_left_prev = ex[1]

        ex[-1] = ex_right_prev + c_factor * (ex[-2] - ex[-1])
        ex_right_prev = ex[-2]

        # Record state
        ex_history[n, :] = ex
        hy_history[n, :] = hy

    x_coords = np.arange(num_cells) * dx
    t_coords = np.arange(num_steps) * dt

    return SimulationData(
        ex=ex_history,
        hy=hy_history,
        x=x_coords,
        t=t_coords,
        dx=dx,
        dt=dt,
        frequency=frequency,
        title="1D FDTD Wave Simulation"
    )


def load_simulation_data(filepath: str, hy_filepath: Optional[str] = None) -> SimulationData:
    """
    Loads simulation results from .npz, .npy, .csv, or .json files.
    Supports both Grid Matrix CSVs and Tabular (column-based) CSVs.
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"Simulation file not found: {filepath}")

    ext = path.suffix.lower()

    if ext == ".npz":
        data = np.load(path)
        ex = data["ex"] if "ex" in data else data["Ex"]
        hy = data.get("hy", data.get("Hy", None))
        x = data.get("x", None)
        t = data.get("t", None)
        dx = float(data.get("dx", 1.0))
        dt = float(data.get("dt", 1.0))
        freq = float(data.get("frequency", data.get("freq", 0.0))) or None
        return SimulationData(ex=ex, hy=hy, x=x, t=t, dx=dx, dt=dt, frequency=freq, title=path.stem)

    elif ext == ".npy":
        arr = np.load(path)
        if arr.ndim == 1:
            arr = arr.reshape(1, -1)
        return SimulationData(ex=arr, title=path.stem)

    elif ext == ".csv":
        # Check if CSV is Tabular (has named headers) or Matrix (pure numbers / comments)
        with open(filepath, "r", encoding="utf-8") as f:
            first_line = f.readline().strip()

        # Check for metadata comments (e.g. # dx=0.001 dt=1e-12)
        dx = 1.0
        dt = 1.0
        if first_line.startswith("#"):
            for part in first_line.lstrip("#").split():
                if "=" in part:
                    k, v = part.split("=", 1)
                    if k.strip().lower() == "dx":
                        dx = float(v)
                    elif k.strip().lower() == "dt":
                        dt = float(v)

        # Check for tabular header columns like "timestep,cell,ex,hy"
        header_tokens = [t.strip().lower() for t in first_line.split(",")]
        is_tabular = any(h in header_tokens for h in ["cell", "timestep", "step", "ex", "hy", "mag", "magnitude"])

        if is_tabular:
            # Parse tabular columns
            data_arr = np.genfromtxt(filepath, delimiter=",", names=True, dtype=float)
            names = [n.lower() for n in data_arr.dtype.names]

            step_col = next((n for n in data_arr.dtype.names if n.lower() in ["timestep", "step", "t_step", "t"]), None)
            cell_col = next((n for n in data_arr.dtype.names if n.lower() in ["cell", "cell_index", "x_cell", "x"]), None)
            ex_col = next((n for n in data_arr.dtype.names if n.lower() in ["ex", "e", "efield", "mag", "magnitude"]), None)
            hy_col = next((n for n in data_arr.dtype.names if n.lower() in ["hy", "h", "hfield"]), None)

            if step_col and cell_col and ex_col:
                steps = data_arr[step_col].astype(int)
                cells = data_arr[cell_col].astype(int)
                max_step = steps.max() + 1
                max_cell = cells.max() + 1

                ex_mat = np.zeros((max_step, max_cell), dtype=float)
                ex_mat[steps, cells] = data_arr[ex_col]

                hy_mat = None
                if hy_col:
                    hy_mat = np.zeros((max_step, max_cell), dtype=float)
                    hy_mat[steps, cells] = data_arr[hy_col]

                return SimulationData(ex=ex_mat, hy=hy_mat, dx=dx, dt=dt, title=path.stem)

        # Standard Matrix Format (each row = 1 timestep, each column = 1 grid cell)
        matrix = np.loadtxt(filepath, delimiter=",", comments="#")
        if matrix.ndim == 1:
            matrix = matrix.reshape(1, -1)

        hy_mat = None
        if hy_filepath:
            hy_path = Path(hy_filepath)
            if hy_path.exists():
                hy_mat = np.loadtxt(hy_filepath, delimiter=",", comments="#")
                if hy_mat.ndim == 1:
                    hy_mat = hy_mat.reshape(1, -1)

        return SimulationData(ex=matrix, hy=hy_mat, dx=dx, dt=dt, title=path.stem)

    elif ext == ".json":
        with open(path, "r", encoding="utf-8") as f:
            raw = json.load(f)
        ex = np.array(raw["ex" if "ex" in raw else "Ex"], dtype=np.float64)
        if ex.ndim == 1:
            ex = ex.reshape(1, -1)
        hy = np.array(raw["hy" if "hy" in raw else "Hy"], dtype=np.float64) if ("hy" in raw or "Hy" in raw) else None
        if hy is not None and hy.ndim == 1:
            hy = hy.reshape(1, -1)
        dx = float(raw.get("dx", 1.0))
        dt = float(raw.get("dt", 1.0))
        freq = raw.get("frequency", None)
        return SimulationData(ex=ex, hy=hy, dx=dx, dt=dt, frequency=freq, title=path.stem)

    else:
        raise ValueError(f"Unsupported file format '{ext}'. Expected .npz, .npy, .csv, or .json")


def plot_animation(
    data: SimulationData,
    save_path: Optional[str] = None,
    fps: int = 30,
    interval: int = 25,
    overlay: bool = False,
    abs_magnitude: bool = False,
    use_distance: bool = False
):
    """
    Renders an interactive animated plot showing Ex and Hy field magnitudes along the cell grid.
    - X-axis: Cell indices (0, 1, ..., N-1)
    - Y-axis: Magnitude of E and H fields
    """
    num_steps, num_cells = data.ex.shape
    cells_axis = np.arange(num_cells)

    if use_distance:
        x_axis = data.x * 1e3 if np.max(data.x) < 1.0 else data.x
        x_label = "Spatial Position (mm)" if np.max(data.x) < 1.0 else "Spatial Position (m)"
    else:
        x_axis = cells_axis
        x_label = "Grid Cell (Index)"

    # Field data (absolute magnitude or signed amplitude)
    ex_data = np.abs(data.ex) if abs_magnitude else data.ex
    hy_data = np.abs(data.hy) if abs_magnitude else data.hy

    e_label_name = r"$|E_x|$" if abs_magnitude else r"$E_x$"
    h_label_name = r"$|H_y|$" if abs_magnitude else r"$H_y$"

    # Set up styling
    plt.style.use("seaborn-v0_8-darkgrid" if "seaborn-v0_8-darkgrid" in plt.style.available else "default")

    if overlay:
        # Single plot with dual y-axes for E and H field magnitudes
        fig, ax_ex = plt.subplots(figsize=(10, 6))
        ax_hy = ax_ex.twinx()

        ex_max = max(np.max(np.abs(ex_data)), 1e-6) * 1.15
        hy_max = max(np.max(np.abs(hy_data)), 1e-9) * 1.15

        if abs_magnitude:
            ax_ex.set_ylim(0, ex_max)
            ax_hy.set_ylim(0, hy_max)
        else:
            ax_ex.set_ylim(-ex_max, ex_max)
            ax_hy.set_ylim(-hy_max, hy_max)

        ax_ex.set_xlim(x_axis[0], x_axis[-1])
        ax_ex.set_xlabel(x_label, fontsize=12, fontweight="bold")
        ax_ex.set_ylabel(f"Electric Field Magnitude {e_label_name} (V/m)", color="#007acc", fontsize=11, fontweight="bold")
        ax_hy.set_ylabel(f"Magnetic Field Magnitude {h_label_name} (A/m)", color="#d9534f", fontsize=11, fontweight="bold")

        ax_ex.tick_params(axis='y', labelcolor="#007acc")
        ax_hy.tick_params(axis='y', labelcolor="#d9534f")
        ax_ex.grid(True, linestyle="--", alpha=0.6)

        (line_ex,) = ax_ex.plot([], [], color="#007acc", lw=2, label=f"Electric Field ({e_label_name})")
        (line_hy,) = ax_hy.plot([], [], color="#d9534f", lw=2, linestyle="--", label=f"Magnetic Field ({h_label_name})")

        # Combine legends
        lines = [line_ex, line_hy]
        labels = [l.get_label() for l in lines]
        ax_ex.legend(lines, labels, loc="upper right", frameon=True)

    else:
        # Stacked subplots: Top = E-Field Magnitude, Bottom = H-Field Magnitude
        fig, (ax_ex, ax_hy) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

        ex_max = max(np.max(np.abs(ex_data)), 1e-6) * 1.15
        hy_max = max(np.max(np.abs(hy_data)), 1e-9) * 1.15

        if abs_magnitude:
            ax_ex.set_ylim(0, ex_max)
            ax_hy.set_ylim(0, hy_max)
        else:
            ax_ex.set_ylim(-ex_max, ex_max)
            ax_hy.set_ylim(-hy_max, hy_max)

        ax_ex.set_xlim(x_axis[0], x_axis[-1])
        ax_ex.set_ylabel(f"Electric Field Magnitude {e_label_name} (V/m)", fontsize=11, fontweight="bold")
        ax_ex.grid(True, linestyle="--", alpha=0.6)
        (line_ex,) = ax_ex.plot([], [], color="#007acc", lw=2, label=f"Electric Field ({e_label_name})")
        ax_ex.legend(loc="upper right", frameon=True)

        ax_hy.set_xlim(x_axis[0], x_axis[-1])
        ax_hy.set_xlabel(x_label, fontsize=12, fontweight="bold")
        ax_hy.set_ylabel(f"Magnetic Field Magnitude {h_label_name} (A/m)", fontsize=11, fontweight="bold")
        ax_hy.grid(True, linestyle="--", alpha=0.6)
        (line_hy,) = ax_hy.plot([], [], color="#d9534f", lw=2, label=f"Magnetic Field ({h_label_name})")
        ax_hy.legend(loc="upper right", frameon=True)

    title_text = fig.suptitle(f"{data.title} | Timestep: 0 / {num_steps}", fontsize=13, fontweight="bold")

    # Animation state
    is_paused = [False]
    current_frame = [0]

    def init():
        line_ex.set_data([], [])
        line_hy.set_data([], [])
        return line_ex, line_hy, title_text

    def update(frame):
        current_frame[0] = frame
        t_val = data.t[frame]
        t_str = f"{t_val * 1e12:.2f} ps" if t_val < 1e-6 else f"{t_val * 1e9:.2f} ns"

        line_ex.set_data(x_axis, ex_data[frame, :])
        line_hy.set_data(x_axis, hy_data[frame, :])

        title_text.set_text(f"{data.title} | Timestep {frame}/{num_steps} (t = {t_str})")
        return line_ex, line_hy, title_text

    save_ext = Path(save_path).suffix.lower() if save_path else None
    if save_ext in [".png", ".jpg", ".jpeg", ".pdf", ".svg"]:
        update(min(num_steps // 2, num_steps - 1))
        plt.savefig(save_path, dpi=200)
        plt.close(fig)
        print(f"Snapshot frame successfully saved to {save_path}")
        return

    anim = animation.FuncAnimation(
        fig,
        update,
        init_func=init,
        frames=num_steps,
        interval=interval,
        blit=False,
        repeat=True
    )

    # Interactive keyboard shortcuts
    def on_key(event):
        if event.key == " ":
            if is_paused[0]:
                anim.resume()
                is_paused[0] = False
            else:
                anim.pause()
                is_paused[0] = True
        elif event.key == "right":
            if is_paused[0]:
                current_frame[0] = (current_frame[0] + 1) % num_steps
                update(current_frame[0])
                fig.canvas.draw_idle()
        elif event.key == "left":
            if is_paused[0]:
                current_frame[0] = (current_frame[0] - 1 + num_steps) % num_steps
                update(current_frame[0])
                fig.canvas.draw_idle()

    fig.canvas.mpl_connect("key_press_event", on_key)
    plt.tight_layout()

    if save_path:
        print(f"Saving animation to '{save_path}' (FPS: {fps})...")
        if save_ext == ".gif":
            anim.save(save_path, writer="pillow", fps=fps)
        elif save_ext in [".mp4", ".avi", ".mov"]:
            anim.save(save_path, writer="ffmpeg", fps=fps)
        else:
            anim.save(save_path, fps=fps)
        plt.close(fig)
        print(f"Animation successfully saved to {save_path}")
    else:
        print("\nInteractive controls:")
        print("  - [Spacebar] : Play / Pause")
        print("  - [Left / Right Arrow] : Step frame (when paused)")
        plt.show()


def plot_heatmap(
    data: SimulationData,
    save_path: Optional[str] = None,
    field: str = "ex",
    abs_magnitude: bool = False
):
    """
    Renders a 2D Spatiotemporal heatmap:
    - X-axis: Grid Cells (Cell Index: 0 to N-1)
    - Y-axis: Timestep index (0 to T-1)
    - Color : Field Magnitude
    """
    num_steps, num_cells = data.ex.shape
    raw_field = data.ex if field.lower() == "ex" else data.hy
    field_arr = np.abs(raw_field) if abs_magnitude else raw_field

    field_name = r"$|E_x|$" if abs_magnitude else r"$E_x$"
    unit_str = "V/m" if field.lower() == "ex" else "A/m"

    fig, ax = plt.subplots(figsize=(10, 7))

    if abs_magnitude:
        cmap = "inferno"
        vmin = 0.0
        vmax = np.max(field_arr)
    else:
        cmap = "RdBu_r"
        vbound = max(np.max(np.abs(field_arr)), 1e-6)
        vmin = -vbound
        vmax = vbound

    im = ax.pcolormesh(
        np.arange(num_cells),
        np.arange(num_steps),
        field_arr,
        shading="auto",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax
    )

    cbar = plt.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label(f"Field Magnitude {field_name} ({unit_str})", fontsize=11, fontweight="bold")

    ax.set_title(f"Spatiotemporal Propagation Map: {data.title}", fontsize=13, fontweight="bold", pad=12)
    ax.set_xlabel("Grid Cell (Index)", fontsize=12, fontweight="bold")
    ax.set_ylabel("Timestep Index", fontsize=12, fontweight="bold")

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=200)
        print(f"Heatmap plot saved to {save_path}")
    else:
        plt.show()


def plot_snapshots(
    data: SimulationData,
    save_path: Optional[str] = None,
    num_snapshots: int = 5,
    abs_magnitude: bool = False
):
    """
    Plots multi-step spatial snapshots:
    - X-axis: Grid Cells (Cell Index: 0 to N-1)
    - Y-axis: Field Magnitude
    """
    num_steps, num_cells = data.ex.shape
    step_indices = np.linspace(0, num_steps - 1, num_snapshots, dtype=int)
    cells_axis = np.arange(num_cells)

    ex_data = np.abs(data.ex) if abs_magnitude else data.ex
    hy_data = np.abs(data.hy) if abs_magnitude else data.hy

    e_label_name = r"$|E_x|$" if abs_magnitude else r"$E_x$"
    h_label_name = r"$|H_y|$" if abs_magnitude else r"$H_y$"

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, num_snapshots))

    for idx, color in zip(step_indices, colors):
        t_val = data.t[idx]
        t_str = f"{t_val * 1e12:.1f} ps" if t_val < 1e-6 else f"{t_val * 1e9:.2f} ns"
        pct = int((idx / (num_steps - 1)) * 100)

        ax1.plot(cells_axis, ex_data[idx, :], color=color, lw=1.8, label=f"Step {idx} ({pct}% | {t_str})")
        ax2.plot(cells_axis, hy_data[idx, :], color=color, lw=1.8, label=f"Step {idx} ({pct}% | {t_str})")

    ax1.set_title(f"Simulation Progression Snapshots: {data.title}", fontsize=13, fontweight="bold")
    ax1.set_ylabel(f"Electric Field Magnitude {e_label_name} (V/m)", fontsize=11, fontweight="bold")
    ax1.grid(True, linestyle="--", alpha=0.6)
    ax1.legend(loc="upper right", fontsize=9)

    ax2.set_xlabel("Grid Cell (Index)", fontsize=12, fontweight="bold")
    ax2.set_ylabel(f"Magnetic Field Magnitude {h_label_name} (A/m)", fontsize=11, fontweight="bold")
    ax2.grid(True, linestyle="--", alpha=0.6)
    ax2.legend(loc="upper right", fontsize=9)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=200)
        print(f"Snapshots plot saved to {save_path}")
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="1D Electromagnetic Wave Simulation Visualizer (X: Cells, Y: Field Magnitude)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 visualize.py --demo
  python3 visualize.py --demo --abs
  python3 visualize.py --demo --overlay
  python3 visualize.py --demo --mode heatmap
  python3 visualize.py --demo --mode snapshots --save snapshots.png
  python3 visualize.py --input sim_data.csv --mode animation
        """
    )

    parser.add_argument(
        "-i", "--input",
        type=str,
        default=None,
        help="Path to simulation data file (.npz, .npy, .csv, .json)"
    )
    parser.add_argument(
        "--hy-input",
        type=str,
        default=None,
        help="Optional path to separate Hy field CSV file (if Ex and Hy are in separate files)"
    )
    parser.add_argument(
        "-m", "--mode",
        choices=["animation", "heatmap", "snapshots"],
        default="animation",
        help="Visualization mode: 'animation' (default), 'heatmap' (Cells vs Timesteps), 'snapshots'"
    )
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Run built-in 1D FDTD simulation demo"
    )
    parser.add_argument(
        "--abs",
        action="store_true",
        dest="abs_magnitude",
        help="Plot absolute field magnitude (|E| and |H|) on the y-axis"
    )
    parser.add_argument(
        "--overlay",
        action="store_true",
        help="Overlay both E and H field magnitudes on a single plot with dual y-axes"
    )
    parser.add_argument(
        "--distance",
        action="store_true",
        dest="use_distance",
        help="Use physical distance (meters/mm) instead of grid cell index for x-axis"
    )
    parser.add_argument(
        "--save",
        type=str,
        default=None,
        help="Save plot/animation to file (.gif, .mp4, .png, .pdf)"
    )
    parser.add_argument(
        "--fps",
        type=int,
        default=30,
        help="Frames per second for saved animation (default: 30)"
    )
    parser.add_argument(
        "--interval",
        type=int,
        default=25,
        help="Animation frame delay in milliseconds (default: 25)"
    )
    parser.add_argument(
        "--field",
        choices=["ex", "hy"],
        default="ex",
        help="Field component to display in heatmap mode ('ex' or 'hy')"
    )
    parser.add_argument(
        "--source",
        choices=["gaussian", "sine"],
        default="gaussian",
        help="Source type for demo mode ('gaussian' or 'sine')"
    )

    args = parser.parse_args()

    # Ensure numpy and matplotlib are available before proceeding
    ensure_dependencies()

    # Load data or generate demo
    if args.demo or args.input is None:
        if args.input is None and not args.demo:
            print("No input file provided. Launching built-in 1D FDTD simulation demo...\n")
        data = generate_demo_simulation(source_type=args.source)
    else:
        print(f"Loading simulation results from '{args.input}'...")
        data = load_simulation_data(args.input, hy_filepath=args.hy_input)
        print(f"Loaded data with {data.ex.shape[0]} timesteps and {data.ex.shape[1]} spatial cells.")

    # Execute selected visualization mode
    if args.mode == "animation":
        plot_animation(
            data,
            save_path=args.save,
            fps=args.fps,
            interval=args.interval,
            overlay=args.overlay,
            abs_magnitude=args.abs_magnitude,
            use_distance=args.use_distance
        )
    elif args.mode == "heatmap":
        plot_heatmap(
            data,
            save_path=args.save,
            field=args.field,
            abs_magnitude=args.abs_magnitude
        )
    elif args.mode == "snapshots":
        plot_snapshots(
            data,
            save_path=args.save,
            abs_magnitude=args.abs_magnitude
        )


if __name__ == "__main__":
    main()
