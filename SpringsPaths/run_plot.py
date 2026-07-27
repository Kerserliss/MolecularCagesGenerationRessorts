import matplotlib.pyplot as plt
import subprocess
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np

def variable_computing_raw(directory: Path):
    """
    Same traversal as variable_computing, but keeps every individual
    file's value instead of averaging per molecule.
    Returns: list_dict where each dict maps molecule_name -> list of raw values.
    """
    dictonary_mol_angle = {}
    dictonary_mol_dist = {}
    dictonary_mol_time = {}
    dictonary_mol_collision = {}
    list_dict = [dictonary_mol_dist, dictonary_mol_angle, dictonary_mol_time, dictonary_mol_collision]

    for mol_dir in directory.iterdir():
        if not mol_dir.is_dir():
            continue
        for i in range(len(list_dict)):
            list_dict[i][mol_dir.name] = []

        for file_name in mol_dir.iterdir():
            if mol_dir.name == "BOQJUM" and file_name.name == "moc_2.mol2":
                continue
            if "moc_" in file_name.name:
                with open(file_name, "r") as f:
                    f.readline()  # skip first line
                    for i in range(len(list_dict)):
                        value = float(f.readline().strip().split()[2])
                        list_dict[i][mol_dir.name].append(value)

    return list_dict

def variable_computing(directory : Path):
    """Get all RMSD for an directory and all mol2 file and return the average."""
    dictonary_mol_angle = {}
    dictonary_mol_dist = {}
    dictonary_mol_time = {}
    dictonary_mol_collision = {}
    list_dict = [dictonary_mol_dist,dictonary_mol_angle,dictonary_mol_time,dictonary_mol_collision]

    for mol_dir in directory.iterdir():
        compteur = 0
        if not mol_dir.is_dir():
            continue
        for i in range(len(list_dict)):
            list_dict[i][mol_dir.name] = 0

        for file_name in mol_dir.iterdir():
            if mol_dir.name == "BOQJUM" and file_name.name == "moc_2.mol2":
                continue
            if "moc_" in file_name.name :
                with open(file_name,"r") as f :
                    # We skip the first line
                    f.readline()
                    # And we start filling our variable.
                    for i in range(len(list_dict)):
                        list_dict[i][mol_dir.name] += float(f.readline().strip().split()[2])
                    compteur += 1
        if compteur == 0:
            raise ValueError(f"No matching mol2 files found in {mol_dir}")
        for i in range(len(list_dict)):
            list_dict[i][mol_dir.name] /= compteur

    return list_dict
root_path = Path(__file__).parent
results_path_euclidien = root_path / "demos_results_euclidien"
results_path_astar = root_path / "demos_results_astar"
list_euclidien = variable_computing(results_path_euclidien)
list_astar = variable_computing(results_path_astar)

list_euclidien_raw = variable_computing_raw(results_path_euclidien)
list_astar_raw = variable_computing_raw(results_path_astar)

metric_labels = ["RMSD_dist", "RMSD_angle", "Time", "Collision"]

def print_results(methods):
    """
    methods: dict mapping method name -> list_dict
             e.g. {"Euclidian": list_euclidien, "Astar": list_astar}
    """
    mol_names = set()
    for list_dict in methods.values():
        for d in list_dict:
            mol_names.update(d.keys())

    for mol_name in sorted(mol_names):
        print(mol_name)
        for method_name, list_dict in methods.items():
            print(f"\t{method_name}")
            for label, metric_dict in zip(metric_labels, list_dict):
                value = metric_dict.get(mol_name, "N/A")
                print(f"\t\t{label}: {value}")

def plot_results(methods):
    """
    methods: dict mapping method name -> list_dict
             e.g. {"Euclidian": list_euclidien, "Astar": list_astar}
             where each list_dict = [dict_dist, dict_angle, dict_time, dict_collision]
             and each of those dicts is keyed by molecule name.
    """
    mol_names = sorted(set().union(*[d.keys() for list_dict in methods.values() for d in list_dict]))
    method_names = list(methods.keys())

    fig, axes = plt.subplots(
        len(metric_labels), 1,
        figsize=(max(8, len(mol_names) * 0.6), 4 * len(metric_labels))
    )

    x = np.arange(len(mol_names))
    bar_width = 0.8 / len(method_names)

    for ax, label, metric_index in zip(axes, metric_labels, range(len(metric_labels))):
        for m_idx, method_name in enumerate(method_names):
            list_dict = methods[method_name]
            metric_dict = list_dict[metric_index]
            values = [metric_dict.get(mol, 0) for mol in mol_names]
            offset = (m_idx - len(method_names) / 2) * bar_width + bar_width / 2
            ax.bar(x + offset, values, bar_width, label=method_name)

        ax.set_title(label)
        ax.set_xticks(x)
        ax.set_xticklabels(mol_names, rotation=45, ha="right")
        ax.legend()
        ax.set_ylabel(label)

    plt.tight_layout()
    plt.savefig("results_plot10.png", dpi=150, bbox_inches="tight")

def plot_boxplot_per_molecule(methods):
    """
    methods: dict mapping method name -> list_dict (from variable_computing_raw)
             each dict inside list_dict maps molecule_name -> list of raw values

    For each metric, draws one subplot with grouped boxplots:
    one x-position per molecule, one box per method at that position.
    """
    mol_names = sorted(set().union(*[d.keys() for list_dict in methods.values() for d in list_dict]))
    method_names = list(methods.keys())

    fig, axes = plt.subplots(
        len(metric_labels), 1,
        figsize=(max(8, len(mol_names) * 0.8), 4 * len(metric_labels))
    )

    x = np.arange(len(mol_names))
    box_width = 0.8 / len(method_names)
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    for ax, label, metric_index in zip(axes, metric_labels, range(len(metric_labels))):
        for m_idx, method_name in enumerate(method_names):
            metric_dict = methods[method_name][metric_index]
            data = [metric_dict.get(mol, []) for mol in mol_names]

            offset = (m_idx - len(method_names) / 2) * box_width + box_width / 2
            positions = x + offset

            bp = ax.boxplot(data, positions=positions, widths=box_width * 0.9, patch_artist=True)
            for patch in bp['boxes']:
                patch.set_facecolor(colors[m_idx % len(colors)])
                patch.set_alpha(0.6)

        ax.set_title(label)
        ax.set_xticks(x)
        ax.set_xticklabels(mol_names, rotation=45, ha="right")
        ax.set_ylabel(label)

        handles = [plt.Rectangle((0, 0), 1, 1, facecolor=colors[i % len(colors)], alpha=0.6)
                   for i in range(len(method_names))]
        ax.legend(handles, method_names)

    plt.tight_layout()
    plt.savefig("boxplot_per_molecule10.png", dpi=150, bbox_inches="tight")

def plot_percentage_change(reference, comparison, reference_name="Euclidian", comparison_name="Astar"):
    """
    reference:  list_dict used as the baseline (e.g. list_euclidien)
    comparison: list_dict to compare against the reference (e.g. list_astar)
    Both are [dict_dist, dict_angle, dict_time, dict_collision], keyed by molecule name.

    For each molecule and each metric, computes:
        pct_change = (comparison - reference) / reference * 100
    Negative = comparison is lower than reference; positive = higher.
    """
    mol_names = sorted(set(reference[0].keys()) & set(comparison[0].keys()))

    fig, axes = plt.subplots(
        len(metric_labels), 1,
        figsize=(max(8, len(mol_names) * 0.6), 4 * len(metric_labels))
    )

    x = np.arange(len(mol_names))

    for ax, label, metric_index in zip(axes, metric_labels, range(len(metric_labels))):
        ref_dict = reference[metric_index]
        comp_dict = comparison[metric_index]

        pct_changes = []
        for mol in mol_names:
            ref_val = ref_dict.get(mol)
            comp_val = comp_dict.get(mol)
            if ref_val in (None, 0):
                pct_changes.append(0)
            else:
                pct_changes.append((comp_val - ref_val) / ref_val * 100)

        colors = ["tab:red" if v > 0 else "tab:green" for v in pct_changes]
        ax.bar(x, pct_changes, color=colors)
        ax.axhline(0, color="black", linewidth=0.8)
        ax.set_title(f"{label}: % change ({comparison_name} vs {reference_name})")
        ax.set_xticks(x)
        ax.set_xticklabels(mol_names, rotation=45, ha="right")
        ax.set_ylabel("% change")

    plt.tight_layout()
    plt.savefig("pct_change10.png", dpi=150, bbox_inches="tight")


#plot_percentage_change(list_euclidien, list_astar)

#plot_boxplot_per_molecule({"Euclidian": list_euclidien_raw, "Astar": list_astar_raw})



#plot_results({"Euclidian": list_euclidien, "Astar": list_astar})



results_path_euclidien10 = root_path / "Demos_temp10" / "demos_results_euclidien"
results_path_astar10 = root_path / "Demos_temp10" /"demos_results_astar"

list_euclidien_10 = variable_computing(results_path_euclidien10)
list_astar_10 = variable_computing(results_path_astar10)

plot_results({"Euclidian": list_euclidien_10, "Astar": list_astar_10})
plot_percentage_change(list_euclidien_10, list_astar_10)
