#!/usr/bin/env python3
"""
Electromagnetic Wave Visualizer (2D & 1D)
Plays generated FDTD simulation frames as animated 2D heatmaps (or 1D line graphs).
"""

import argparse
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


def load_data(
    filepath: str = "output.csv",
    shape: tuple[int, int] | None = None,
    force_1d: bool = False,
):
    """
    Load simulation data from CSV (tabular or matrix format).
    Returns (e_field, h_field, is_2d, (rows, cols), timesteps).
    For 2D: fields have shape (num_steps, rows, cols)
    For 1D: fields have shape (num_steps, num_cells)
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {filepath}")

    with open(path, "r", encoding="utf-8") as f:
        first_line = f.readline().strip()

    headers = (
        [h.strip().lower() for h in first_line.split(",")]
        if ("," in first_line and not first_line.replace("-", "").split(",")[0].isdigit())
        else []
    )

    is_tabular = any(
        k in headers
        for k in ("cell", "ex", "ey", "ez", "hx", "hy", "hz", "row", "col", "x", "y")
    )

    rows, cols = None, None

    if is_tabular:
        raw = np.loadtxt(path, delimiter=",", skiprows=1)
        if raw.ndim == 1:
            raw = raw.reshape(1, -1)

        # Check for explicit 2D coordinate columns (e.g. timestep, row, col, ...)
        if ("row" in headers and "col" in headers) or ("y" in headers and "x" in headers):
            r_idx = headers.index("row") if "row" in headers else headers.index("y")
            c_idx = headers.index("col") if "col" in headers else headers.index("x")
            rows = int(raw[:, r_idx].max()) + 1
            cols = int(raw[:, c_idx].max()) + 1
            num_cells = rows * cols
            e_col = max(r_idx, c_idx) + 1
            h_col = e_col + 1 if raw.shape[1] > e_col + 1 else None
        elif "cell" in headers:
            cell_idx = headers.index("cell")
            num_cells = int(raw[:, cell_idx].max()) + 1
            e_col = cell_idx + 1
            h_col = cell_idx + 2 if raw.shape[1] > cell_idx + 2 else None
        else:
            num_cells = int(raw[:, 1].max()) + 1
            e_col = 2
            h_col = 3 if raw.shape[1] > 3 else None

        num_steps = raw.shape[0] // num_cells
        if num_steps == 0:
            raise ValueError(
                f"File contains {raw.shape[0]} rows, which is less than 1 frame ({num_cells} cells)."
            )

        valid_len = num_steps * num_cells
        raw = raw[:valid_len]
        timesteps = raw[::num_cells, 0].astype(int)

        e_raw = raw[:, e_col]
        h_raw = (
            raw[:, h_col]
            if (h_col is not None and h_col < raw.shape[1])
            else np.zeros_like(e_raw)
        )
    else:
        # Matrix format (each row = 1 timestep, each col = 1 cell)
        raw = np.loadtxt(path, delimiter=",")
        if raw.ndim == 1:
            raw = raw.reshape(1, -1)
        num_steps, num_cells = raw.shape
        timesteps = np.arange(num_steps)
        e_raw = raw.flatten()
        h_raw = np.zeros_like(e_raw)

    is_2d = False
    if not force_1d:
        if shape is not None:
            rows, cols = shape
            if rows * cols != num_cells:
                raise ValueError(
                    f"Specified shape ({rows}, {cols}) = {rows * cols} cells does not match data ({num_cells} cells)."
                )
            is_2d = rows > 1 and cols > 1
        elif rows is not None and cols is not None:
            is_2d = rows > 1 and cols > 1
        else:
            # Check if num_cells is a perfect square
            side = int(round(np.sqrt(num_cells)))
            if side * side == num_cells and side > 1:
                rows, cols = side, side
                is_2d = True

    if is_2d and rows is not None and cols is not None:
        e_field = e_raw.reshape(num_steps, rows, cols)
        h_field = h_raw.reshape(num_steps, rows, cols)
        grid_dim = (rows, cols)
    else:
        e_field = e_raw.reshape(num_steps, num_cells)
        h_field = h_raw.reshape(num_steps, num_cells)
        grid_dim = (1, num_cells)
        is_2d = False

    return e_field, h_field, is_2d, grid_dim, timesteps


def play_simulation_2d(
    ez: np.ndarray,
    h: np.ndarray,
    timesteps: np.ndarray | None = None,
    interval: int = 20,
    save_path: str | None = None,
    fps: int = 30,
    field_mode: str = "both",
    cmap_e: str = "RdBu_r",
    cmap_h: str = "PuOr_r",
):
    """Plays the frames as animated 2D heatmaps (Ez and H)."""
    num_steps, ny, nx = ez.shape
    if timesteps is None:
        timesteps = np.arange(num_steps)

    e_max = float(np.max(np.abs(ez)))
    e_max = 1.0 if (e_max == 0.0 or np.isnan(e_max)) else e_max * 1.15

    h_max = float(np.max(np.abs(h)))
    h_max = 1.0 if (h_max == 0.0 or np.isnan(h_max)) else h_max * 1.15

    if field_mode == "e":
        fig, ax_e = plt.subplots(figsize=(8, 6.5))
        ax_h = None
    elif field_mode == "h":
        fig, ax_h = plt.subplots(figsize=(8, 6.5))
        ax_e = None
    else:
        fig, (ax_e, ax_h) = plt.subplots(1, 2, figsize=(13, 6))

    im_e, im_h = None, None

    if ax_e is not None:
        im_e = ax_e.imshow(
            ez[0],
            cmap=cmap_e,
            vmin=-e_max,
            vmax=e_max,
            origin="lower",
            aspect="equal",
            extent=[0, nx, 0, ny],
        )
        ax_e.set_title("Electric Field $E_z$ (V/m)", fontweight="bold")
        ax_e.set_xlabel("X (cells)", fontweight="bold")
        ax_e.set_ylabel("Y (cells)", fontweight="bold")
        cbar_e = fig.colorbar(im_e, ax=ax_e, shrink=0.8, pad=0.03)
        cbar_e.set_label("V/m", rotation=270, labelpad=15)

    if ax_h is not None:
        im_h = ax_h.imshow(
            h[0],
            cmap=cmap_h,
            vmin=-h_max,
            vmax=h_max,
            origin="lower",
            aspect="equal",
            extent=[0, nx, 0, ny],
        )
        ax_h.set_title("Magnetic Field $H$ (A/m)", fontweight="bold")
        ax_h.set_xlabel("X (cells)", fontweight="bold")
        ax_h.set_ylabel("Y (cells)", fontweight="bold")
        cbar_h = fig.colorbar(im_h, ax=ax_h, shrink=0.8, pad=0.03)
        cbar_h.set_label("A/m", rotation=270, labelpad=15)

    t0 = timesteps[0] if len(timesteps) > 0 else 0
    title = fig.suptitle(
        f"Timestep: {t0}  (Frame 1 / {num_steps})  [Grid: {nx}×{ny}]",
        fontsize=12,
        fontweight="bold",
    )

    is_paused = False
    current_frame = 0

    def update(frame):
        nonlocal current_frame
        current_frame = frame
        artists = [title]
        if im_e is not None:
            im_e.set_data(ez[frame])
            artists.append(im_e)
        if im_h is not None:
            im_h.set_data(h[frame])
            artists.append(im_h)
        status = " [PAUSED]" if is_paused else ""
        t_val = timesteps[frame] if frame < len(timesteps) else frame
        title.set_text(
            f"Timestep: {t_val}  (Frame {frame + 1} / {num_steps})  [Grid: {nx}×{ny}]{status}"
        )
        return tuple(artists)

    anim = FuncAnimation(fig, update, frames=num_steps, interval=interval, blit=False)

    def on_key(event):
        nonlocal is_paused, current_frame
        if event.key == " ":
            if is_paused:
                anim.resume()
            else:
                anim.pause()
            is_paused = not is_paused
            update(current_frame)
            fig.canvas.draw_idle()
        elif event.key == "right":
            current_frame = min(current_frame + 1, num_steps - 1)
            update(current_frame)
            fig.canvas.draw_idle()
        elif event.key == "left":
            current_frame = max(current_frame - 1, 0)
            update(current_frame)
            fig.canvas.draw_idle()
        elif event.key == "up":
            current_frame = min(current_frame + 10, num_steps - 1)
            update(current_frame)
            fig.canvas.draw_idle()
        elif event.key == "down":
            current_frame = max(current_frame - 10, 0)
            update(current_frame)
            fig.canvas.draw_idle()
        elif event.key == "r":
            current_frame = 0
            update(0)
            fig.canvas.draw_idle()

    fig.canvas.mpl_connect("key_press_event", on_key)

    if save_path:
        print(f"Saving animation to '{save_path}' (fps={fps})...")
        anim.save(save_path, fps=fps)
        print("Save complete.")
    else:
        plt.tight_layout()
        plt.show()


def play_simulation_1d(
    ex: np.ndarray,
    hy: np.ndarray,
    timesteps: np.ndarray | None = None,
    interval: int = 20,
    save_path: str | None = None,
    fps: int = 30,
):
    """Plays the frames as an animated 2-panel 1D line graph (Ex and Hy)."""
    num_steps, num_cells = ex.shape
    cells = np.arange(num_cells)
    if timesteps is None:
        timesteps = np.arange(num_steps)

    fig, (ax_e, ax_h) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)

    e_max = max(float(np.max(np.abs(ex))), 1e-6) * 1.15
    h_max = max(float(np.max(np.abs(hy))), 1e-9) * 1.15

    # Electric Field subplot
    ax_e.set_xlim(0, num_cells - 1)
    ax_e.set_ylim(-e_max, e_max)
    ax_e.set_ylabel("Electric Field $E$ (V/m)", fontweight="bold")
    ax_e.grid(True, linestyle="--", alpha=0.6)
    (line_e,) = ax_e.plot(cells, ex[0], color="#007acc", lw=1.8, label="$E$")
    ax_e.legend(loc="upper right")

    # Magnetic Field subplot
    ax_h.set_xlim(0, num_cells - 1)
    ax_h.set_ylim(-h_max, h_max)
    ax_h.set_xlabel("Grid Cell Index", fontweight="bold")
    ax_h.set_ylabel("Magnetic Field $H$ (A/m)", fontweight="bold")
    ax_h.grid(True, linestyle="--", alpha=0.6)
    (line_h,) = ax_h.plot(cells, hy[0], color="#d9534f", lw=1.8, label="$H$")
    ax_h.legend(loc="upper right")

    t0 = timesteps[0] if len(timesteps) > 0 else 0
    title = fig.suptitle(
        f"Timestep: {t0}  (Frame 1 / {num_steps})", fontsize=12, fontweight="bold"
    )

    is_paused = False
    current_frame = 0

    def update(frame):
        nonlocal current_frame
        current_frame = frame
        line_e.set_ydata(ex[frame])
        line_h.set_ydata(hy[frame])
        status = " [PAUSED]" if is_paused else ""
        t_val = timesteps[frame] if frame < len(timesteps) else frame
        title.set_text(
            f"Timestep: {t_val}  (Frame {frame + 1} / {num_steps}){status}"
        )
        return line_e, line_h, title

    anim = FuncAnimation(fig, update, frames=num_steps, interval=interval)

    def on_key(event):
        nonlocal is_paused, current_frame
        if event.key == " ":
            if is_paused:
                anim.resume()
            else:
                anim.pause()
            is_paused = not is_paused
            update(current_frame)
            fig.canvas.draw_idle()
        elif event.key == "right":
            current_frame = min(current_frame + 1, num_steps - 1)
            update(current_frame)
            fig.canvas.draw_idle()
        elif event.key == "left":
            current_frame = max(current_frame - 1, 0)
            update(current_frame)
            fig.canvas.draw_idle()
        elif event.key == "r":
            current_frame = 0
            update(0)
            fig.canvas.draw_idle()

    fig.canvas.mpl_connect("key_press_event", on_key)

    if save_path:
        print(f"Saving animation to '{save_path}' (fps={fps})...")
        anim.save(save_path, fps=fps)
        print("Save complete.")
    else:
        plt.tight_layout()
        plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Electromagnetic Wave Visualizer (2D Heatmap & 1D Line Graph)"
    )
    parser.add_argument(
        "filepath",
        nargs="?",
        default="output.csv",
        help="Path to CSV simulation output file (default: output.csv)",
    )
    parser.add_argument(
        "--shape",
        "-s",
        nargs=2,
        type=int,
        metavar=("ROWS", "COLS"),
        help="Explicit 2D grid shape (Ny Nx / rows cols). Auto-detected if square.",
    )
    parser.add_argument(
        "--field",
        "-f",
        choices=["both", "e", "h"],
        default="both",
        help="Field panel to display: both, e, or h (default: both)",
    )
    parser.add_argument(
        "--interval",
        "-i",
        type=int,
        default=20,
        help="Delay between frames in milliseconds (default: 20)",
    )
    parser.add_argument(
        "--fps",
        type=int,
        default=30,
        help="Frames per second when saving animation (default: 30)",
    )
    parser.add_argument(
        "--save",
        type=str,
        default=None,
        help="Save animation to file (e.g. wave.mp4 or wave.gif) instead of showing GUI",
    )
    parser.add_argument(
        "--cmap-e",
        type=str,
        default="RdBu_r",
        help="Colormap for Electric Field (default: RdBu_r)",
    )
    parser.add_argument(
        "--cmap-h",
        type=str,
        default="PuOr_r",
        help="Colormap for Magnetic Field (default: PuOr_r)",
    )
    parser.add_argument(
        "--1d",
        dest="force_1d",
        action="store_true",
        help="Force 1D line plot visualization",
    )

    args = parser.parse_args()

    shape_tuple = tuple(args.shape) if args.shape else None

    try:
        e_field, h_field, is_2d, grid_dim, timesteps = load_data(
            args.filepath, shape=shape_tuple, force_1d=args.force_1d
        )
    except Exception as err:
        print(f"Error loading '{args.filepath}': {err}", file=sys.stderr)
        sys.exit(1)

    num_steps = e_field.shape[0]

    if is_2d:
        rows, cols = grid_dim
        print(
            f"Playing 2D simulation '{args.filepath}' ({num_steps} steps, grid: {cols}×{rows})..."
        )
        print("Controls: [Space] Play/Pause | [← / →] Step ±1 | [↑ / ↓] Step ±10 | [r] Reset")
        play_simulation_2d(
            e_field,
            h_field,
            timesteps=timesteps,
            interval=args.interval,
            save_path=args.save,
            fps=args.fps,
            field_mode=args.field,
            cmap_e=args.cmap_e,
            cmap_h=args.cmap_h,
        )
    else:
        num_cells = grid_dim[1]
        print(
            f"Playing 1D simulation '{args.filepath}' ({num_steps} steps, {num_cells} cells)..."
        )
        print("Controls: [Space] Play/Pause | [← / →] Step ±1 | [r] Reset")
        play_simulation_1d(
            e_field,
            h_field,
            timesteps=timesteps,
            interval=args.interval,
            save_path=args.save,
            fps=args.fps,
        )


if __name__ == "__main__":
    main()

