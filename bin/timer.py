from collections import namedtuple
from contextlib import contextmanager
import os
import threading
import time

# Never reset: one process is one step, so this holds that step's whole run.
_timings = []

# vram is None for a step that never touched a GPU, which is most of them.
Measurement = namedtuple("Measurement", "label seconds rss vram")

# Resolved on first use and kept: importing cupy costs a second or two and every sample
# asks. False means unresolved, None means this process has no GPU.
_cuda = False

_PAGE_SIZE = os.sysconf("SC_PAGE_SIZE")

# Fine enough to catch a short allocation, coarse enough to be free next to a step that
# runs for minutes.
_SAMPLE_SECONDS = 0.05


def _rss():
    """Bytes resident right now. Read from /proc rather than resource.getrusage, whose
    ru_maxrss is a process-wide high-water mark with no way to reset it -- every step after
    a heavy one would report the heavy one's number."""
    with open("/proc/self/statm") as f:
        return int(f.read().split()[1]) * _PAGE_SIZE


def _vram():
    """Bytes in use on the GPU right now, or None in a process without one.

    Asked of the driver, not of an allocator: rapids-singlecell allocates through RMM,
    which leaves cupy's default memory pool empty however much of the card is in use. This
    counts every allocator plus the CUDA context -- what nvidia-smi shows for the device.
    """
    global _cuda
    if _cuda is False:
        try:
            import cupy

            cupy.cuda.runtime.getDeviceCount()   # raises without a driver or a card
            _cuda = cupy.cuda.runtime
        except Exception:
            _cuda = None

    if _cuda is None:
        return None

    free, total = _cuda.memGetInfo()
    return total - free


@contextmanager
def _peaks():
    """Sample memory while a block runs, yielding the largest readings it reaches."""
    peak = {"rss": _rss(), "vram": _vram()}
    done = threading.Event()

    def sample():
        peak["rss"] = max(peak["rss"], _rss())
        vram = _vram()
        if vram is not None:
            peak["vram"] = max(peak["vram"], vram)

    def poll():
        while not done.wait(_SAMPLE_SECONDS):
            sample()

    # Daemon, so a step that dies mid-block cannot leave the process hanging on it.
    watcher = threading.Thread(target=poll, daemon=True)
    watcher.start()
    try:
        yield peak
    finally:
        done.set()
        watcher.join()
        sample()   # a block shorter than the interval would otherwise never be sampled


def _gb(size):
    return f"{size / 1024 ** 3:.1f} GB"


@contextmanager
def timer(label):
    """Context manager that times a block and prints elapsed time and peak memory inline.

    The memory is the block's own peak, so a column names the step that needs it.
    """
    start = time.perf_counter()
    with _peaks() as peak:
        yield
    elapsed = time.perf_counter() - start
    measurement = Measurement(label, elapsed, peak["rss"], peak["vram"])
    _timings.append(measurement)

    minutes, seconds = divmod(elapsed, 60)
    if minutes > 0:
        line = f"[{label}] {int(minutes)}m {seconds:.1f}s"
    else:
        line = f"[{label}] {seconds:.2f}s"
    line += f" | rss {_gb(measurement.rss)}"
    if measurement.vram is not None:
        line += f" | vram {_gb(measurement.vram)}"
    print(line)


def write_timings_tsv(path):
    """Write recorded measurements to a TSV with columns: step, seconds, rss_gb, vram_gb.

    Appends a final 'Total' row: seconds sum, memory takes the largest seen. vram_gb is
    empty for a step that never touched a GPU.
    """
    with open(path, "w") as f:
        f.write("step\tseconds\trss_gb\tvram_gb\n")
        total = 0.0
        for label, elapsed, rss, vram in _timings:
            vram_gb = "" if vram is None else f"{vram / 1024 ** 3:.3f}"
            f.write(f"{label}\t{elapsed:.4f}\t{rss / 1024 ** 3:.3f}\t{vram_gb}\n")
            total += elapsed

        peak_rss = max(measurement.rss for measurement in _timings)
        vrams = [measurement.vram for measurement in _timings if measurement.vram is not None]
        peak_vram = f"{max(vrams) / 1024 ** 3:.3f}" if vrams else ""
        f.write(f"Total\t{total:.4f}\t{peak_rss / 1024 ** 3:.3f}\t{peak_vram}\n")


def timing_summary(path=None):
    """Print a formatted table of all recorded measurements and a total.

    If path is provided, also write them to a TSV file.
    """
    if not _timings:
        print("No timings recorded.")
        return

    show_vram = any(measurement.vram is not None for measurement in _timings)
    col = max(max(len(measurement.label) for measurement in _timings), 4)

    header = f"\n{'Step':<{col}}  {'Time':>10}  {'RSS':>9}"
    if show_vram:
        header += f"  {'VRAM':>9}"
    print(header)
    print("-" * (len(header) - 1))

    total = 0
    for label, elapsed, rss, vram in _timings:
        minutes, seconds = divmod(elapsed, 60)
        time_str = f"{int(minutes)}m {seconds:.1f}s" if minutes > 0 else f"{seconds:.2f}s"
        row = f"{label:<{col}}  {time_str:>10}  {_gb(rss):>9}"
        if show_vram:
            row += f"  {'' if vram is None else _gb(vram):>9}"
        print(row)
        total += elapsed

    print("-" * (len(header) - 1))
    minutes, seconds = divmod(total, 60)
    total_str = f"{int(minutes)}m {seconds:.1f}s" if minutes > 0 else f"{seconds:.2f}s"
    # Memory is the largest step's rather than a sum: the peaks do not coincide.
    peak_rss = max(measurement.rss for measurement in _timings)
    vrams = [measurement.vram for measurement in _timings if measurement.vram is not None]
    row = f"{'Total':<{col}}  {total_str:>10}  {_gb(peak_rss):>9}"
    if show_vram:
        row += f"  {_gb(max(vrams)) if vrams else '':>9}"
    print(row)

    if path is not None:
        write_timings_tsv(path)
