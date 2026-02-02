# highlight_variant_pos_protein_structure.py takes 
# a input .pse file, an input csv file, a filtering percentage, a chain name, 
# an output file name end with .png
# an optional highlight color, an optional sphere_scale
# it read in the .pse file and ,csv file and label the vaiable position based on information of .csv file
# it output a .png file, a modified .pse file, and an csv file in output/protein_structure

pymol -cq -r code/highlight_variant_pos_protein_structure.py \
    -- data/protein_structure/PrsA.pse \
    output/allele_info/PrsA_D39V/self_reference/allele_info_no_strain_name.csv \
    1.0 \
    A \
    output/protein_structure/PrsA.png \
    --highlight_color yellow \
    --sphere_scale 1.0

pymol -cq -r code/highlight_variant_pos_protein_structure.py \
    -- data/protein_structure/alia_colored.pse \
    output/allele_info/AliA_D39V/self_reference/allele_info_no_strain_name.csv \
    1.0 \
    A \
    output/protein_structure/AliA.png \
    --highlight_color yellow \
    --sphere_scale 1.0

pymol -cq -r code/highlight_variant_pos_protein_structure.py \
    -- data/protein_structure/alib_colored.pse \
    output/allele_info/AliB_D39V/self_reference/allele_info_no_strain_name.csv \
    1.0 \
    A \
    output/protein_structure/AliB.png \
    --highlight_color yellow \
    --sphere_scale 1.0

pymol -cq -r code/highlight_variant_pos_protein_structure.py \
    -- data/protein_structure/amia_colored.pse \
    output/allele_info/AmiA_D39V/self_reference/allele_info_no_strain_name.csv \
    1.0 \
    A \
    output/protein_structure/AmiA.png \
    --highlight_color yellow \
    --sphere_scale 1.0


pymol -cq -r code/highlight_variant_pos_protein_structure.py \
    -- data/protein_structure/SlrA_colored.pse \
    output/allele_info/SlrA_D39V/self_reference/allele_info_no_strain_name.csv \
    1.0 \
    A \
    output/protein_structure/SlrA.png \
    --highlight_color yellow \
    --sphere_scale 1.0


pymol -cq -r code/highlight_variant_pos_protein_structure.py \
    -- data/protein_structure/SlrA_colored.pse \
    output/allele_info/SlrA_D39V/self_reference/allele_info_no_strain_name.csv \
    1.0 \
    A \
    output/protein_structure/SlrA.png \
    --highlight_color yellow \
    --sphere_scale 1.0

pymol -cq -r code/highlight_variant_pos_protein_structure.py \
    -- data/protein_structure/PrsA.pse \
    output/allele_info/PrsA_D39V/self_reference/allele_info_no_strain_name.csv \
    1.0 \
    A \
    output/protein_structure/PrsA.png \
    --highlight_color yellow \
    --sphere_scale 1.0

pymol -cq -r code/highlight_variant_pos_protein_structure.py \
  -- data/protein_structure/SlrA_colored.pse \
  output/allele_info/SlrA_D39V/self_reference/allele_info_no_strain_name.csv \
  1.0 \
  A \
  output/protein_structure/SlrA_old.png \
  --highlight_color yellow \
  --no_spheres

pymol -cq -r code/highlight_variant_pos_protein_structure.py \
  -- data/protein_structure/PrsA.pse \
  output/allele_info/PrsA_D39V/self_reference/allele_info_no_strain_name.csv \
  1.0 \
  A \
  output/protein_structure/PrsA_old.png \
  --highlight_color yellow \
  --no_spheres
