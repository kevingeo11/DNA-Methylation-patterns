import sys
from pymol import cmd
import numpy as np
import scipy
import matplotlib
import Bio.PDB
import pickle

print(sys.version)
# print(pymol.get_version_message())
print(np.__version__)
print(scipy.__version__)
print(matplotlib.__version__)
# print(Bio.PDB.__version__)

vdw_radii = {"C":1.70,
            "O":1.52,
            "N":1.55,
            "H":1.20,
            "P":1.80,
            "S":1.80
            }

sugar_backbone_dna_atoms = ["P",
                            "OP1",
                            "OP2",
                            "O5'",
                            "C5'",
                            "C4'",
                            "O4'",
                            "C3'",
                            "O3'",
                            "C2'",
                            "C1'"
                            ]

main_path = '/home/kevin/DNA-Methylation-patterns/'
data_path = main_path + 'downstream/chr1_pstrand/'
steric_path = main_path + 'downstream/steric_clash/'

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

def make_superimpose(steric_path, pdb_nuclesome, pdb_dnmt):
    dnmt_file = steric_path + 'pdb_files/' + pdb_dnmt + '.pdb'
    nuc_file = steric_path + 'pdb_files/' + pdb_nuclesome + '.pdb'
    
    #remove water
    # remove_water_from_pdb(nuc_file, pdb_nuclesome)
    # remove_water_from_pdb(dnmt_file, pdb_dnmt)

    # parse structure
    nucleosome_structure = parse_PDB(nuc_file, 'nuc')
    dnmt_structure = parse_PDB(dnmt_file, 'dnmt')

    model_nuc = nucleosome_structure[0]
    model_dnmt = dnmt_structure[0]

    ################ ATOMS ######################
    len_DNA_part = 12
    
    # DNA strand 1 (B,Y): numbers 1-19
    # DNA strand 2 (C,Z): numbers 20-38
    B_start = 4
    B_end = B_start + len_DNA_part
    C_start = 36 - len_DNA_part  #21 38
    C_end = 36
    
    #################### DNMT ATOMS ##########################
    # chain A = DNMT1
    # chain B = DNA strand 1
    # chain C = DNA strand 2
    
    dnmt_atoms = []
    range_dnmt1 = range(B_start, B_end)
    range_dnmt2 = range(C_start, C_end)

    #append DNMT DNA atoms to list if in backbone and in middle of DNA (4-15, and 24 to 35)
    nbr_atoms_chain_Y = 0
    for chain in model_dnmt:
        chain_id = chain.id
        #rename for coloring since nucleosome structure has same chain ids
        if chain_id == "A":
            chain.id = "X"
            chain_id = "X"
        if chain_id == "B":
            chain.id = "Y"
            chain_id = "Y"
        if chain_id == "C":
            chain.id = "Z"
            chain_id = "Z"
        for residue in chain:
            res_id = residue.id
            res_info = res_id[0]
            res_pos = res_id[1]
            if res_info == " ":
                #DNA strand 1
                if chain_id == "Y":
                    if res_pos in range_dnmt1:
                        res_name = residue.resname #DA
                        for atom in residue:
                            atom_id = atom.id
                            if atom_id in sugar_backbone_dna_atoms:
                                dnmt_atoms.append(atom)
                                nbr_atoms_chain_Y += 1
                
                #DNA strand 2    
                if chain_id == "Z":
                    if res_pos in range_dnmt2:
                        res_name = residue.resname #DA
                        for atom in residue:
                            atom_id = atom.id
                            if atom_id in sugar_backbone_dna_atoms:
                                dnmt_atoms.append(atom)
    
    # posiition 1 and 136 missing 3 atoms
    dnmt_atoms_pos1 = []
    for a_pos in range(0,len(dnmt_atoms)):
        if a_pos not in [0,1,2]:
            dnmt_atoms_pos1.append(dnmt_atoms[a_pos])
    
    dnmt_atoms_pos136 = []
    for a_pos in range(0,len(dnmt_atoms)):
        if a_pos not in [nbr_atoms_chain_Y,nbr_atoms_chain_Y+1,nbr_atoms_chain_Y+2]:
            dnmt_atoms_pos136.append(dnmt_atoms[a_pos])
    
    # print(dnmt_atoms)
    # print()
    # print(dnmt_atoms_pos1)
    # print()
    # print(dnmt_atoms_pos136)

    '''
    dnmt_atoms [<Atom P>, <Atom OP1>, <Atom OP2>, <Atom O5'>, <Atom C5'>, <Atom C4'>, <Atom O4'>, <Atom C3'>, <Atom O3'>, <Atom C2'>, <Atom C1'>, <Atom P>, <Atom OP1>, <Atom OP2>, <Atom O5'>, <Atom C5'>, <Atom C4'>, <Atom O4'>, <Atom C3'>, <Atom O3'>]
    nucleosome_atoms [<Atom P>, <Atom OP1>, <Atom OP2>, <Atom O5'>, <Atom C5'>, <Atom C4'>, <Atom O4'>, <Atom C3'>, <Atom O3'>, <Atom C2'>, <Atom C1'>, <Atom P>, <Atom OP1>, <Atom OP2>, <Atom O5'>, <Atom C5'>, <Atom C4'>, <Atom O4'>, <Atom C3'>, <Atom O3'>]
    
    dnmt_atoms_pos1 [<Atom O5'>, <Atom C5'>, <Atom C4'>, <Atom O4'>, <Atom C3'>, <Atom O3'>, <Atom C2'>, <Atom C1'>, <Atom P>, <Atom OP1>, <Atom OP2>, <Atom O5'>, <Atom C5'>, <Atom C4'>, <Atom O4'>, <Atom C3'>, <Atom O3'>, <Atom C2'>, <Atom C1'>, <Atom P>]    
    nucleosome_atoms [<Atom O5'>, <Atom C5'>, <Atom C4'>, <Atom O4'>, <Atom C3'>, <Atom O3'>, <Atom C2'>, <Atom C1'>, <Atom P>, <Atom OP1>, <Atom OP2>, <Atom O5'>, <Atom C5'>, <Atom C4'>, <Atom O4'>, <Atom C3'>, <Atom O3'>, <Atom C2'>, <Atom C1'>, <Atom P>]
    '''

    #################### NUCLEOSOME ##########################
    
    I_start = 1#2#12
    I_end = I_start+len_DNA_part
    J_start = 148-len_DNA_part
    J_end = 148


    ######## START WHILE LOOP ITERATION DNA NUCLEOSOME ##########
    clash_dict = dict()
    ### ?? What is 148??
    while I_end <= 148:#####148 for superpos of every position, 14 for first two positions
        print("pos nucleosome", I_start)
        clash_dict[I_start] = dict()
        clash_dict[I_start]["steric_clash_list"] = dict()
        
        ################ NUCLEOSOME ATOMS ######################
        #chain I: DNA strand 1
        #chain J: DNA strand 2
        #chain A: nucleosome complex
        nucleosome_atoms = []
        range_53 = range(I_start,I_end)
        range_35 = range(J_start,J_end)
        
        #do not count them in steric clash, atoms in DNA of nucleosme that are in supoerposition with DNA od DNMT
        range_53_all_superpos = range(I_start-3,I_end+4)
        range_35_all_superpos = range(J_start-4,J_end+3)
        nuc_atoms_all_superpos = set()
        
        # print(range_53_all_superpos)
        # print(range_35_all_superpos)
        # print(range_53)
        # print(range_35)

        res_nbrs_I_color = []
        res_nbrs_J_color = []
        for chain in model_nuc:
            chain_id = chain.id
    
            pos_res_1based = 1
            for residue in chain:
                if chain_id == "I":
                    if pos_res_1based in range_53:
                        res_id = residue.id
                        res_nbrs_I_color.append(res_id[1])
                        
                        res_info = res_id[0]
                        res_pos = res_id[1]
                        if res_info == " ":
                            #if chain_id == "B" and res_pos in range():
                            res_name = residue.resname #DA
                            #print res_name
                            for atom in residue:
                                atom_id = atom.id
                                if atom_id in sugar_backbone_dna_atoms:
                                    nucleosome_atoms.append(atom)
                                    
                    if pos_res_1based in range_53_all_superpos:
                        res_id = residue.id
                        res_info = res_id[0]
                        res_pos = res_id[1]
                        if res_info == " ":
                            res_name = residue.resname #DA
                            for atom in residue:
                                atom_id = atom.id
                                #if atom_id in sugar_backbone_dna_atoms:
                                nuc_atoms_all_superpos.add(atom)
                                
                if chain_id == "J":
                    if pos_res_1based in range_35:
                        #print pos_res_1based
                        res_id = residue.id
                        res_nbrs_J_color.append(res_id[1])
                        res_info = res_id[0]
                        res_pos = res_id[1]
                        if res_info == " ":
                            #if chain_id == "B" and res_pos in range():
                            res_name = residue.resname #DA
                            #print res_name
                            for atom in residue:
                                atom_id = atom.id
                                if atom_id in sugar_backbone_dna_atoms:
                                    nucleosome_atoms.append(atom)
                    
                    if pos_res_1based in range_35_all_superpos:
                        res_id = residue.id
                        res_info = res_id[0]
                        res_pos = res_id[1]
                        if res_info == " ":
                            res_name = residue.resname #DA
                            for atom in residue:
                                atom_id = atom.id
                                #if atom_id in sugar_backbone_dna_atoms:
                                nuc_atoms_all_superpos.add(atom)
                pos_res_1based += 1  
    
        
        #print(len(nucleosome_atoms))
        #print(len(nuc_atoms_all_superpos))
        #print(nuc_atoms_all_superpos)
        #if I_start == 136 or I_start == 1 or I_start == 2:
        #    print("dnmt_atoms",dnmt_atoms#[-20:-1])
        #    print("dnmt_atoms_pos1",dnmt_atoms_pos1#[-20:-1])
        #    print("dnmt_atoms_pos136",dnmt_atoms_pos136#[-20:-1])
        #    print("nucleosome_atoms",nucleosome_atoms#[-20:-1])
        
        
        #print(dnmt_atoms)

        ######### SUPERIMPOSITION ###########
        
        sup = Bio.PDB.Superimposer()
        #if nuc position 1, nuc atoms <Atom P>, <Atom OP1>, <Atom OP2>, are missing -> delete Atom P>, <Atom OP1>, <Atom OP2>, (pos 0-2) from dnmt atoms
        if I_start == 1:
            sup.set_atoms(nucleosome_atoms, dnmt_atoms_pos1)
        elif I_start == 136:
            sup.set_atoms(nucleosome_atoms, dnmt_atoms_pos136)
        else:
            sup.set_atoms(nucleosome_atoms, dnmt_atoms)
        #print(sup.rotran)
        
        #all_dnmt_atoms = dnmt_structure.get_atoms()
        #all_nucleosome_atoms = nucleosome_structure.get_atoms()
        
        sup.apply(model_dnmt.get_atoms())
        #sup.apply(dnmt_atoms)
        
        clash_dict[I_start]["rmsd"] = sup.rms
         
        io = Bio.PDB.PDBIO()
        io.set_structure(dnmt_structure)    #write PDB file
        output_superpos_name = "dnmt_superpos_coords_" + str(I_start)
        dnmt_moved_file = steric_path + "superposition/dnmt_coords/" + output_superpos_name + ".pdb"
        io.save(dnmt_moved_file)

        ############ save both structures in one file ############
        if True:
            if I_start != 1:    #first round
                cmd.delete("all")
    
            cmd.load(dnmt_moved_file, output_superpos_name)
            cmd.load(nuc_file, pdb_nuclesome)
            #pymol.cmd.hide("all")
            cmd.show_as("cartoon")
            cmd.set("cartoon_fancy_helices","1")
            #cmd.set("cartoon_ring_mode","3")
            #cmd.set("cartoon_ring_finder","1")
            #cmd.set("cartoon_ladder_mode","1")
            #cmd.set("cartoon_nucleic_acid_mode","4")
            #cmd.set("cartoon_ring_transparency","0.5")
                        
            cmd.color("white","all")
            '''
            #color nucleosome 1kx5
            cmd.color("lightpink","chain A")#deepsalmon
            cmd.color("lightpink","chain E")
            
            cmd.color("wheat","chain B")#yelloworange
            cmd.color("wheat","chain F")
            
            cmd.color("palegreen","chain C")#splitpea
            cmd.color("palegreen","chain G")
            
            cmd.color("lightblue","chain D")#slate
            cmd.color("lightblue","chain H")
            '''
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
            '''
            cmd.color("tv_yellow","chain Y and resi 4")
            cmd.color("tv_yellow","chain Y and resi 5")
            cmd.color("tv_yellow","chain Y and resi 14")
            cmd.color("tv_yellow","chain Y and resi 15")
            
            cmd.color("tv_yellow","chain Z and resi 24")
            cmd.color("tv_yellow","chain Z and resi 25")
            cmd.color("tv_yellow","chain Z and resi 34")
            cmd.color("tv_yellow","chain Z and resi 35")
            
            
            
            cmd.color("lightblue","all")
            cmd.color("palecyan","resn DA")
            cmd.color("palecyan","resn DC")
            cmd.color("palecyan","resn DG")
            cmd.color("palecyan","resn DT")
            
            cmd.color("orange","chain X") #lightorange
            cmd.color("lightpink","chain Y")
            cmd.color("lightpink","chain Z")
            
            #cmd.color("palegreen","resn DA")
            #cmd.color("paleyellow","resn DC")
            #cmd.color("lightpink","resn DG")
            #cmd.color("lightblue","resn DT")
        '''
            
            
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
            
            
            cmd.save(steric_path+"superposition/superpos_dnmt_nuc_pse/superposition_dnmt_nucleosome_"+str(I_start)+".pse")
            fname = steric_path+"superposition/superpos_dnmt_nuc_pse/png/superposition_dnmt_nucleosome_"+str(I_start)+".png"
            cmd.png(fname, width=800, height=700, dpi=300,ray=1,quiet=1)

        ######## get new coordinates #############
        dnmt_moved_structure = parse_PDB(dnmt_moved_file,"dnmt_moved")
        
        #print len(list(dnmt_structure.get_atoms()))
        #print len(list(nucleosome_structure.get_atoms()))
        #print len(list(dnmt_moved_structure.get_atoms()))
        #print dnmt_moved_pdb_file_dict["Z"]["ARG1259"]
        #print nucleosome_pdb_file_dict["A"]["ALA31"]
        
        ########################## CALCULATE STERIC CLASH ###############################
        atom_list_dnmt = Bio.PDB.Selection.unfold_entities(dnmt_moved_structure, 'A') # A for atoms
        atom_list_nuc = Bio.PDB.Selection.unfold_entities(nucleosome_structure, 'A') 
        #print len(atom_list_dnmt)
        #print len(atom_list_nuc)
        
        ns = Bio.PDB.NeighborSearch(atom_list_nuc)
        i = 0
        dnmt_atoms_consider_clash = 0
        for atom_dnmt in atom_list_dnmt:#[0:1000]:
            residue_dnmt = atom_dnmt.get_parent()
            chain_dnmt = residue_dnmt.get_parent()
            if i%1000==0:
                print(str(i)+"/"+str(len(atom_list_dnmt)))
            i += 1
            #chain X: use only protein atoms, not DNA atoms
            if  atom_dnmt.element != "ZN" and atom_dnmt.element != "MN" and atom_dnmt.element != "CL" and chain_dnmt.id == "X":
                dnmt_atoms_consider_clash += 1
                #print atom_dnmt
                #print atom_dnmt.element
                center_coords = atom_dnmt.get_coord()
                neighbors = ns.search(center_coords, 5.0) # 5.0 for distance in angstrom
                
                #if len(neighbors) == 0:
                #    residue_dnmt = atom_dnmt.get_parent()
                #    chain_dnmt = residue_dnmt.get_parent()
                #    print chain_dnmt.id+"_"+residue_dnmt.get_resname()+str(residue_dnmt.id[1])+"_"+atom_dnmt.id
                for neighbor_atom in neighbors:
                    residue_neighbor = neighbor_atom.get_parent()
                    #if is_aa(residue_neighbor):
                    if residue_neighbor.id[0] != "W" and neighbor_atom.element != "ZN" and neighbor_atom.element != "MN" and neighbor_atom.element != "CL" and neighbor_atom not in nuc_atoms_all_superpos:
                    #if True:
                        distance = atom_dnmt-neighbor_atom
                        steric_clash = is_steric_clash(atom_dnmt.element,neighbor_atom.element,distance)
                        if steric_clash:
                            chain_neighbor = residue_neighbor.get_parent()
                            dnmt_res_str = chain_dnmt.id+"_"+residue_dnmt.get_resname()+str(residue_dnmt.id[1])
                            dnmt_str = dnmt_res_str+"_"+atom_dnmt.id
                            nuc_str = chain_neighbor.id+"_"+residue_neighbor.get_resname()+str(residue_neighbor.id[1])+"_"+neighbor_atom.id
                            #if I_start == 5:
                            #    print dnmt_str,nuc_str
                            if dnmt_res_str not in clash_dict[I_start]["steric_clash_list"].keys():
                                clash_dict[I_start]["steric_clash_list"][dnmt_res_str] = [[dnmt_str,nuc_str,distance]]
                            else:
                                clash_dict[I_start]["steric_clash_list"][dnmt_res_str].append([dnmt_str,nuc_str,distance])
                
        I_start += 1
        I_end += 1
        J_start -= 1
        J_end -= 1 

    ######## END WHILE LOOP ITERATION DNA NUCLEOSOME ##########
    
    #print clash_dict
    info_nbr_dict = dict()
    info_nbr_dict["model_dnmt_nbr_residues"] = len(list(model_dnmt.get_residues()))
    info_nbr_dict["model_dnmt_nbr_atoms"] = len(list(model_dnmt.get_atoms()))
    info_nbr_dict["dnmt_atoms_consider_clash"] = dnmt_atoms_consider_clash

    # print("dnmt residues",info_nbr_dict["model_dnmt_nbr_residues"])
    # print("dnmt atoms" ,info_nbr_dict["model_dnmt_nbr_atoms"])
    # print("dnmt considered clash atoms ", info_nbr_dict["dnmt_atoms_consider_clash"])
    # print("dump")
    if True:
        print('file writing')
        pickle_dump(clash_dict, steric_path + "clash_dict")
        pickle_dump(info_nbr_dict, steric_path + "info_nbr_dict")


make_superimpose(steric_path, pdb_nuclesome, pdb_dnmt)

