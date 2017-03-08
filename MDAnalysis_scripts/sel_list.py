# NOTATION:
# sel.append([file_descripter,atom_descriptor,selection1,selection2,donor_atoms,accepter_atoms])

sel = []
#sel.append(['Protein_Protein','protein','protein',('O*','N*','S*'),('O*','N*','S*')])
sel.append(['Protein_Nucleic','protein','nucleic or resname A5 C5 G5 U5 A3 C3 G3 U3',('O*','N*','S*'),('O*','N*','S*')])
#sel.append(['Nucleic_Nucleic','nucleic or resname A5 C5 G5 U5 A3 C3 G3 U3','nucleic or resname A5 C5 G5 U5 A3 C3 G3 U3 and not resname atp adp',('O*','N*','S*'),('O*','N*','S*')])

#res_list = [452,453,454,455,456,457,458]	# nucleic residues

