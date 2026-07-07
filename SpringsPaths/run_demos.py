import subprocess
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

MAX_PARALLEL = 8

root_path = Path(__file__).parent
demos_path = root_path / "demos" / "demos_realSubstrat"
assert demos_path.is_dir()

results_path = root_path / "demos_results"
results_path.mkdir(exist_ok=True)

bin_file = root_path / "bin" / "springPaths"
assert bin_file.is_file()

commands = []
for mol_dir in demos_path.iterdir():
    if not mol_dir.is_dir():
        continue
    results_mol_dir = results_path / mol_dir.name
    results_mol_dir.mkdir(parents=True, exist_ok=True)
    for file_name in mol_dir.iterdir():
        if "_moc" in file_name.name:
            number = file_name.name.split("_moc")[1].split(".")[0]
            output_path = results_mol_dir / f"moc_{number}.mol2"
            cmd = [
                str(bin_file.relative_to(root_path)),
                "-i",
                str(mol_dir.relative_to(root_path)),
                "-n",
                number,
                "-o",
                str(output_path.relative_to(root_path)),
                "-g",
                "1",
                "-c",
                "1",
            ]
            commands.append(cmd)

total = len(commands)
started = 0
completed = 0
lock = threading.Lock()


def run_cmd(cmd):
    global started, completed

    with lock:
        started += 1
        idx = started
        print(f"[{idx}/{total} started] {' '.join(cmd)}")

    result = subprocess.run(cmd, cwd=root_path)

    with lock:
        completed += 1
        status = (
            "OK" if result.returncode == 0 else f"FAILED (code {result.returncode})"
        )
        print(
            f"[{completed}/{total} done | {total - completed} left] {status}: {' '.join(cmd)}"
        )

    return cmd, result.returncode


print(f"Running {total} commands with max {MAX_PARALLEL} in parallel\n")

with ThreadPoolExecutor(max_workers=MAX_PARALLEL) as executor:
    for cmd, returncode in executor.map(run_cmd, commands):
        pass  # errors already printed inside run_cmd
