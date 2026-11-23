reinitialize
bg_color white
set ray_opaque_background, on
set antialias, 2
set ray_shadows, 0
set label_color, black
set label_size, 22

load /Users/jxs794/Library/CloudStorage/Box-Box/PROTAC_Degrader_DB_Project_2025/Kif20a/Manuscript Pieces/KIF20A_Final_Analysis/MolecularDynamics/UMF607-Z/Analysis_Results/Figure5_Replacement/UMF607_250ns_snapshot.pdb, complex
hide everything

select protein, polymer.protein
select ligand, resn UNK
select pocket, byres (protein within 5.0 of ligand)
select key_residues, resi 152+337+338+339+386+400+401+402+469+470
remove elem H

show cartoon, protein
color gray80, protein
set cartoon_transparency, 0.42, protein

show surface, pocket
util.cbag pocket
color purple, (pocket and elem C)
set transparency, 0.35, pocket
show sticks, pocket
set stick_radius, 0.10, pocket

show sticks, ligand
util.cbag ligand
color purple, (ligand and elem C)
show spheres, ligand
set sphere_scale, 0.22, ligand
set stick_radius, 0.22, ligand

zoom pocket, 4.5
ray 1800,1800
png /Users/jxs794/Library/CloudStorage/Box-Box/PROTAC_Degrader_DB_Project_2025/Kif20a/Manuscript Pieces/KIF20A_Final_Analysis/MolecularDynamics/UMF607-Z/Analysis_Results/Figure5_Replacement/Fig5A_UMF607_KIF20A_pocket.png, dpi=300

hide cartoon, protein
hide sticks, pocket
show surface, pocket
set transparency, 0.60, pocket
show sticks, key_residues
color tv_orange, key_residues
show sticks, ligand
color cyan, ligand
distance arg386_contact, ligand, (resn ARG and resi 386), 4.0, mode=2
distance ile400_contact, ligand, (resn ILE and resi 400), 4.0, mode=2
distance glu402_contact, ligand, (resn GLU and resi 402), 4.0, mode=2
distance phe470_contact, ligand, (resn PHE and resi 470), 4.0, mode=2
color red, arg386_contact
color yellow, ile400_contact
color magenta, glu402_contact
color orange, phe470_contact
set dash_radius, 0.10
set dash_gap, 0.30
set dash_length, 0.20
hide labels
ray 2400,1800
png /Users/jxs794/Library/CloudStorage/Box-Box/PROTAC_Degrader_DB_Project_2025/Kif20a/Manuscript Pieces/KIF20A_Final_Analysis/MolecularDynamics/UMF607-Z/Analysis_Results/Figure5_Replacement/Fig5B_UMF607_key_contacts.png, dpi=300
quit
