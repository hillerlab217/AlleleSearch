"""
This script:

- Loads a structure file (.pdb, .cif, or .pse)
- Highlights variant positions from a CSV file (based on species percentage)
- Optionally labels those positions on a specified chain
- Optionally shows highlighted residues as spheres with adjustable scale
- Renders a high-resolution transparent image
- Saves both the image and a PyMOL session (.pse)

Usage:
pymol -cq highlight_variants.py -- structure.pdb variants.csv threshold chain output.png [--adjust_pos N] [--highlight_color COLOR] [--label] [--sphere_scale F] [--no_spheres]
"""

import sys
import csv
import os
import re
from pymol import cmd

def parse_variant_positions(variant_str):
    if not variant_str or variant_str.strip() == "":
        return []

    raw_entries = re.split(r'[;,]', variant_str)
    variants = []
    for entry in raw_entries:
        try:
            entry = entry.strip()
            if not entry:
                continue
            pos_str, change = entry.split(':')
            ref, alt = change.split('→')
            pos = int(pos_str.strip())
            label = f"{pos_str}:{ref}->{alt}"
            variants.append((pos, label))
        except ValueError:
            continue  # Skip malformed entries

    return variants

def save_variant_locations_csv(variant_dict, output_image_path):
    """
    variant_dict is expected to be {position(int): label(str)}
    """
    csv_path = os.path.splitext(output_image_path)[0] + "_positions.csv"
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Residue_Position", "Label"])
        for pos in sorted(variant_dict.keys()):
            writer.writerow([pos, variant_dict[pos]])
    print(f"✔ Variant positions saved to: {csv_path}")

def highlight_variants(
    pdb_path,
    csv_path,
    species_pct_threshold,
    chain,
    adjust_pos=0,
    highlight_color="red",
    label=False,
    sphere_scale=1.0,
    spheres=True,
):
    # Load structure file or session
    pdb_name = os.path.basename(pdb_path)
    obj_name = os.path.splitext(pdb_name)[0]

    if pdb_path.endswith(".pse"):
        cmd.load(pdb_path)
        print(f"✔ Loaded session file: {pdb_path}")
    else:
        cmd.load(pdb_path, obj_name)
        print(f"✔ Loaded structure: {pdb_path}")

    # Collect variant positions
    # NOTE: store as {position: label} so CSV matches headers
    variant_dict = {}
    with open(csv_path, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            species_pct = float(row['species_percentage'])
            if species_pct >= species_pct_threshold:
                variants = parse_variant_positions(row['variant_positions(position:reference→changed)'])
                for pos, label_text in variants:
                    adjusted_pos = pos - adjust_pos
                    if adjusted_pos <= 0:
                        continue
                    variant_dict[adjusted_pos] = label_text

    # Optional label settings
    if label:
        cmd.set("label_connector", "on")
        cmd.set("label_color", "black")
        cmd.set("label_size", 16)
        cmd.set("label_position", [24, 4, -2])

    # Highlight & show spheres
    for pos in sorted(variant_dict.keys()):
        selection = f"chain {chain} and resi {pos}"
        # color
        cmd.color(highlight_color, selection)
        if spheres:
            cmd.show("spheres", selection)
            cmd.set("sphere_scale", float(sphere_scale), selection)
        else:
            # ensure no spheres even if session had them on
            cmd.hide("spheres", selection)

        # optional label (on CA to avoid clutter)
        if label:
            cmd.label(f"({selection} and name CA)", f'"{variant_dict[pos]}"')

    return variant_dict

def render_and_save_output(image_path):
    output_dir = os.path.dirname(image_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    cmd.set("ray_opaque_background", 0)
    cmd.set("ray_trace_fog", 0)
    cmd.set("antialias", 2)
    # cmd.orient()

    cmd.ray(4000, 3000)
    cmd.png(image_path, dpi=1200)
    print(f"Image saved to: {image_path}")

    session_path = os.path.splitext(image_path)[0] + ".pse"
    cmd.save(session_path)
    print(f"Session saved to: {session_path}")

# Main entry point for PyMOL
if __name__ == '__main__':
    if len(sys.argv) < 6:
        print("Usage: pymol -cq highlight_variants.py -- structure.pdb variants.csv threshold chain output.png [--adjust_pos N] [--highlight_color COLOR] [--label] [--sphere_scale F]")
        sys.exit(1)

    structure_file = sys.argv[1]
    csv_file = sys.argv[2]
    threshold = float(sys.argv[3])
    chain = sys.argv[4]
    output_image = sys.argv[5]

    # Optional args
    adjust_pos = 0
    highlight_color = "red"
    label = False
    spheres = True
    sphere_scale = 1.0

    i = 6
    while i < len(sys.argv):
        if sys.argv[i] == "--adjust_pos":
            try:
                adjust_pos = int(sys.argv[i + 1])
                i += 2
            except (IndexError, ValueError):
                print("Error: --adjust_pos must be followed by an integer")
                sys.exit(1)
        elif sys.argv[i] == "--highlight_color":
            try:
                highlight_color = sys.argv[i + 1]
                i += 2
            except IndexError:
                print("Error: --highlight_color must be followed by a valid PyMOL color name")
                sys.exit(1)
        elif sys.argv[i] == "--label":
            label = True
            i += 1
        elif sys.argv[i] == "--sphere_scale":
            try:
                sphere_scale = float(sys.argv[i + 1])
                i += 2
            except (IndexError, ValueError):
                print("Error: --sphere_scale must be followed by a float (e.g., 1.0)")
                sys.exit(1)
        elif sys.argv[i] == "--no_spheres":
            spheres = False
            i += 1
        else:
            print(f"Unknown argument: {sys.argv[i]}")
            sys.exit(1)

    print(" Running analysis...")
    variant_dict = highlight_variants(
        structure_file,
        csv_file,
        threshold,
        chain,
        adjust_pos=adjust_pos,
        highlight_color=highlight_color,
        label=label,
        sphere_scale=sphere_scale,
        spheres=spheres,
    )
    save_variant_locations_csv(variant_dict, output_image)
    print(" Rendering and saving...")
    render_and_save_output(output_image)
    print(" Done.")