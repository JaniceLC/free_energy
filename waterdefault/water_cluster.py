'''
    example MolecularCluster subclass
'''

from molecular_cluster import MolecularCluster, Molecule

class WaterCluster(MolecularCluster):
    def define_molecule(self):
        # define atom types, internal permutations, bonds
        types=["O","H","H"]
        permlist=[[1,2]]
        bonds=[[0,1],[0,2]]
        return Molecule(types,permlist,bonds)

def main():
    nmol=3
    system=WaterCluster(nmol)
    print system.get_masses()
    print system.get_permlist()
    #print system.get_potential()
    print system.get_potential(eps=0.1550, sig =3.1536 )

if __name__=="__main__":
    main()