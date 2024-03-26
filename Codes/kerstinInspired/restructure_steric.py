import sys
from pymol import cmd
import numpy as np
import scipy
import matplotlib
import Bio.PDB
import pickle

# global var
# locations save/read
# pickle save for dict
# remove water fun
# parse pdb using biopython

# is stecric clash
#   radii calc

# Work flow
#     parse pdb
#     declare start and end
#         chain dependant DNMT Atom calculation
#     superposition while Loop 

## Where did we get these values from??
vdw_radii = {"C" : 1.70,
             "O" : 1.52,
             "N" : 1.55,
             "H" : 1.20,
             "P" : 1.80,
             "S" : 1.80
            }

## who decided this is the suger backbone
sugar_backbone_dna_atoms = ["P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'"]

main_path = 'D:\\Work\\DNA-Methylation-patterns\\'
steric_path = main_path + 'Results\\Steric_Clash\\'

pdb_nuclesome = '1KX5'
pdb_dnmt = '3PTA'

def pickle_dump(data,filei):
    filo = open(filei,"wb")
    pickle.dump(data,filo,-1)
    filo.close()

def remove_water_from_pdb(pdb_file, pdb_name):
    cmd.load(pdb_file, pdb_name)
    cmd.remove('resn hoh')
    cmd.save(pdb_file.replace('.pdb', '_water_rem.pdb'))
    cmd.quit()

def parse_PDB(pdb_file, pdb_id):
    pdb_parser = Bio.PDB.PDBParser(QUIET = True)
    pdb_structure = pdb_parser.get_structure(pdb_id, pdb_file)
                        
    return pdb_structure

def is_steric_clash(atom1,atom2,distance):
    r1 = vdw_radii[atom1]
    r2 = vdw_radii[atom2]
    d_radius = r1+r2
    
    if distance < d_radius:
        return True
    else:
        return False
    
def save_fun():
    if I_start != 1:
        cmd.delete("all")

    cmd.load(dnmt_moved_file, output_superpos_name)
    cmd.load(nuc_file, pdb_nuclesome)
    cmd.show_as("cartoon")
    cmd.set("cartoon_fancy_helices","1")

                
    cmd.color("white","all")
    cmd.color("palecyan","chain I")#palecyan
    cmd.color("palecyan","chain J")
        
    #color dnmt 3pta
    cmd.color("bluewhite","chain X")
    cmd.color("pink","chain X and resi 755-880")
    cmd.color("lightorange","chain X and resi 972-1100")
    cmd.color("lightblue","chain X and resi 1139-1599")
    cmd.color("firebrick","chain X and resi 647-692")
    cmd.color("slate","chain X and resi 693-754")
    
    cmd.color("white","chain Y")
    cmd.color("white","chain Z")

    #color aligned residues nucleosome
    for resi in res_nbrs_I_color:
        if resi < 0:
            cmd.color("marine","chain I and resi \\"+str(resi))
        else:
            cmd.color("marine","chain I and resi "+str(resi))
    
    for resi in res_nbrs_J_color:
        if resi < 0:
            cmd.color("marine","chain J and resi \\"+str(resi))#green
        else:
            cmd.color("marine","chain J and resi "+str(resi))#green
    
    
    methyl_C_pos = res_nbrs_I_color[0]
    # print methyl_C_pos
    if methyl_C_pos < 0:
        cmd.color("cyan","chain I and resi \\"+str(methyl_C_pos))
    else:
        cmd.color("cyan","chain I and resi "+str(methyl_C_pos))   
                
    #dnmt
    for resi in range_dnmt1:
        #unmethylted C
        if resi == 4:
            cmd.color("cyan","chain Y and resi "+str(resi))
        elif resi == 5 or resi == 14 or resi == 15:
            cmd.color("yellow","chain Y and resi "+str(resi))
        else:
            cmd.color("orange","chain Y and resi "+str(resi))

    for resi in range_dnmt2:
        if resi == 24 or resi == 25 or resi == 34 or resi == 35:
            cmd.color("yellow","chain Z and resi "+str(resi))
        else:
            cmd.color("orange","chain Z and resi "+str(resi))
            
    cmd.bg_color("black")
    cmd.set('''seq_view''','''1''',quiet=1)
    
    
    cmd.save(steric_path+"superposition\\superpos_dnmt_nuc_pse\\superposition_dnmt_nucleosome_"+str(I_start)+".pse")
    fname = steric_path+"superposition\\superpos_dnmt_nuc_pse\\png\\superposition_dnmt_nucleosome_"+str(I_start)+"_1.png"
    cmd.png(fname, width=800, height=700, dpi=300,ray=1,quiet=1)

    
if __name__ == '__main__':
    dnmt_file = steric_path + 'pdb_files\\' + pdb_dnmt + '.pdb'
    nuc_file = steric_path + 'pdb_files\\' + pdb_nuclesome + '.pdb'
    
    #remove water
    # remove_water_from_pdb(nuc_file, pdb_nuclesome)
    # remove_water_from_pdb(dnmt_file, pdb_dnmt)

    # parse structure
    # nuc and dnmt are the id of the structure object
    nucleosome_structure = parse_PDB(nuc_file, 'nuc')
    dnmt_structure = parse_PDB(dnmt_file, 'dnmt')

    ## https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ
    ## Go Through above link for structure object documentation

    ## for now we work with only model id = 0
    assert len(nucleosome_structure.get_list()) == 1
    assert len(dnmt_structure.get_list()) == 1

    ## model object - child of the structure object
    model_nuc = nucleosome_structure[0]
    model_dnmt = dnmt_structure[0]

    ## IN DNMT1
    ## Chain 1 = DNMT1
    ## Chain 2 = DNA Strand 1
    ## Chain 3 = DNA Strand 3

    len_DNA_part = 12

    ## from the looks of it we are comparing 12 residues in the DNA -- > Why not 19??
    
    # DNA strand 1 (B,Y): numbers 1-19
    # DNA strand 2 (C,Z): numbers 20-38
    B_start = 4
    B_end = B_start + len_DNA_part
    C_start = 36 - len_DNA_part
    C_end = 36

    dnmt_atoms = []
    range_dnmt1 = range(B_start, B_end)
    range_dnmt2 = range(C_start, C_end)

    ## append DNMT DNA atoms to list if in backbone and in middle of DNA (4-15, and 24 to 35)
    ## Why are not using 1-19 of chain B and 20-38 of chain C
    nbr_atoms_chain_Y = 0
    for chain in model_dnmt:

        ## rename for coloring since nucleosome structure has same chain ids
        rename = {'A': 'X', 'B': 'Y', 'C': 'Z'}
        chain.id = rename[chain.id]

        for residue in chain:
            res_id = residue.id
            
            if res_id[0] == ' ':
                ## We can combine both these range values
                if res_id[1] in range_dnmt1 or res_id[1] in range_dnmt2:
                    for atom in residue:
                        if atom.id in sugar_backbone_dna_atoms:
                            dnmt_atoms.append(atom)
                            if chain.id == 'Y':
                                nbr_atoms_chain_Y += 1

    print(len(dnmt_atoms))
    print(nbr_atoms_chain_Y)
    # print(dnmt_atoms)

    # posiition 1 and 136 missing 3 atoms
    dnmt_atoms_pos1 = []
    for a_pos in range(0,len(dnmt_atoms)):
        ## This is removing the 1st 3 indices -- But why
        if a_pos not in [0,1,2]:
            dnmt_atoms_pos1.append(dnmt_atoms[a_pos])

    # print(len(dnmt_atoms_pos1))
    # print(dnmt_atoms_pos1)
    
    ## 136 next chain begins - 12 residues * 11 atoms == 132
    ## nbr_atoms_chain_Y = 132 -- what is 136???

    dnmt_atoms_pos136 = []
    for a_pos in range(0,len(dnmt_atoms)):
        ## Removing 3 atoms at pos 132,133,134 -- again why??
        if a_pos not in [nbr_atoms_chain_Y,nbr_atoms_chain_Y+1,nbr_atoms_chain_Y+2]:
            dnmt_atoms_pos136.append(dnmt_atoms[a_pos])

    ## Nucleosome
    ## What are these values and why
    
    I_start = 1 #2#12
    I_end = I_start + len_DNA_part
    J_start = 148 - len_DNA_part
    J_end = 148

    ## DNA length 147 in pdb -- why 148??
    ## There are two strands we always read for 5'-3' -- is this correct even in pdb??
    ## 5'-3' read first strand to get its complementatry take the respective 3'-5' for the second stand
    ## sub DNA 1st strand say [0, 12] second stand [147, 135]
    ## Why are we working with 12 residues

    clash_dict = dict()
    while I_end <= 13: #####148 for superpos of every position, 14 for first two positions
        clash_dict[I_start] = dict()
        clash_dict[I_start]["steric_clash_list"] = dict()

        #chain I: DNA strand 1
        #chain J: DNA strand 2
        #chain A: nucleosome complex
        nucleosome_atoms = []
        range_53 = range(I_start, I_end)
        range_35 = range(J_start, J_end)

        # do not count them in steric clash, atoms in DNA of nucleosme that are in supoerposition with DNA od DNMT
        ## what is _all_superpos range??
        ## from debugging this is bringing it back to the 19 original DNA residues and in DNMT
        range_53_all_superpos = range(I_start-3, I_end+4)
        range_35_all_superpos = range(J_start-4, J_end+3)
        nuc_atoms_all_superpos = set()

        res_nbrs_I_color = []
        res_nbrs_J_color = []
        for chain in model_nuc:
            pos_res_1based = 1
            ## I and J are DNA chains /// these can be soft coded
            ## residues here are 934 for I and 955 for J
            ## in pdb its 147 reidues for DNA ???
            ## What do the numbers 934 and 955 mean??

            for residue in chain:
                if chain.id == 'I':
                    if pos_res_1based in range_53:
                        res_id = residue.id
                        res_nbrs_I_color.append(res_id[1])
                        if res_id[0] == " ":
                            for atom in residue:
                                if atom.id in sugar_backbone_dna_atoms:
                                    nucleosome_atoms.append(atom)
                    ## 1614 nucleosome atoms in nucleosome_atoms

                    if pos_res_1based in range_53_all_superpos:
                        res_id = residue.id
                        if res_id[0] == " ":
                            for atom in residue:
                                ## why not checking sugar_backbone_atoms
                                nuc_atoms_all_superpos.add(atom)

                if chain.id == "J":
                    if pos_res_1based in range_35:
                        res_id = residue.id
                        res_nbrs_J_color.append(res_id[1])
                        if res_id[0] == " ":
                            for atom in residue:
                                if atom.id in sugar_backbone_dna_atoms:
                                    nucleosome_atoms.append(atom)
                    
                    if pos_res_1based in range_35_all_superpos:
                        res_id = residue.id
                        if res_id[0] == " ":
                            for atom in residue:
                                nuc_atoms_all_superpos.add(atom)

                pos_res_1based += 1

        sup = Bio.PDB.Superimposer()
        #if nuc position 1, nuc atoms <Atom P>, <Atom OP1>, <Atom OP2>, are missing -> delete Atom P>, <Atom OP1>, <Atom OP2>, (pos 0-2) from dnmt atoms
        if I_start == 1:
            sup.set_atoms(nucleosome_atoms, dnmt_atoms_pos1)
        elif I_start == 136:
            sup.set_atoms(nucleosome_atoms, dnmt_atoms_pos136)
        else:
            sup.set_atoms(nucleosome_atoms, dnmt_atoms)

        sup.apply(model_dnmt.get_atoms())

        clash_dict[I_start]["rmsd"] = sup.rms

        io = Bio.PDB.PDBIO()
        io.set_structure(dnmt_structure)    #write PDB file
        output_superpos_name = "dnmt_superpos_coords_" + str(I_start)
        dnmt_moved_file = steric_path + "superposition\\dnmt_coords\\" + output_superpos_name + ".pdb"
        io.save(dnmt_moved_file)

        save_fun()

        dnmt_moved_structure = parse_PDB(dnmt_moved_file,"dnmt_moved")
        atom_list_dnmt = Bio.PDB.Selection.unfold_entities(dnmt_moved_structure, 'A')
        atom_list_nuc = Bio.PDB.Selection.unfold_entities(nucleosome_structure, 'A') 

        ns = Bio.PDB.NeighborSearch(atom_list_nuc)
        i = 0
        dnmt_atoms_consider_clash = 0

        for atom_dnmt in atom_list_dnmt:
            residue_dnmt = atom_dnmt.get_parent()
            chain_dnmt = residue_dnmt.get_parent()

            if i%1000==0:
                print(str(i)+"/"+str(len(atom_list_dnmt)))
            i += 1

            if atom_dnmt.element != "ZN" and atom_dnmt.element != "MN" and atom_dnmt.element != "CL" and chain_dnmt.id == "X":
                dnmt_atoms_consider_clash += 1
                center_coords = atom_dnmt.get_coord()
                neighbors = ns.search(center_coords, 5.0) # 5.0 for distance in angstrom

                for neighbor_atom in neighbors:
                    residue_neighbor = neighbor_atom.get_parent()
                    if residue_neighbor.id[0] != "W" and neighbor_atom.element != "ZN" and neighbor_atom.element != "MN" and neighbor_atom.element != "CL" and neighbor_atom not in nuc_atoms_all_superpos:
                        distance = atom_dnmt-neighbor_atom
                        steric_clash = is_steric_clash(atom_dnmt.element,neighbor_atom.element,distance)
                        if steric_clash:
                            chain_neighbor = residue_neighbor.get_parent()
                            dnmt_res_str = chain_dnmt.id+"_"+residue_dnmt.get_resname()+str(residue_dnmt.id[1])
                            dnmt_str = dnmt_res_str+"_"+atom_dnmt.id
                            nuc_str = chain_neighbor.id+"_"+residue_neighbor.get_resname()+str(residue_neighbor.id[1])+"_"+neighbor_atom.id
                            
                            if dnmt_res_str not in clash_dict[I_start]["steric_clash_list"].keys():
                                clash_dict[I_start]["steric_clash_list"][dnmt_res_str] = [[dnmt_str,nuc_str,distance]]
                            else:
                                clash_dict[I_start]["steric_clash_list"][dnmt_res_str].append([dnmt_str,nuc_str,distance])


        ## I from 0 to 147/148 and J from 147/148 to 0
        I_start += 1
        I_end += 1
        J_start -= 1
        J_end -= 1 













