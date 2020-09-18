import numpy as np
import math
import os
from src.mulliken import *
from src.anIres import *
import pandas as pd
from ast import literal_eval

class Molecule:
    def __init__(self,mysystem,target_mol):
        self.atoms=mysystem.atoms_per_molecule[target_mol]
        self.species=[i for j,i in enumerate(mysystem.species) if mysystem.mol[j]-1==target_mol]
        self.x=[i for j,i in enumerate(mysystem.x) if mysystem.mol[j]-1==target_mol]
        self.y=[i for j,i in enumerate(mysystem.y) if mysystem.mol[j]-1==target_mol]
        self.z=[i for j,i in enumerate(mysystem.z) if mysystem.mol[j]-1==target_mol]
        self.id_native=mysystem.clusters[target_mol]
        self.B1_native=[i for i in mysystem.B1 if i in self.id_native]
        self.B2_native=[i for i in mysystem.B2 if i in self.id_native]
        self.id_map_dict={i:j+1 for j,i in enumerate(self.id_native)}
        self.B1=[self.id_map_dict[i] for i in self.B1_native]
        self.B2=[self.id_map_dict[i] for i in self.B2_native]
        self.connectivity=resolve_connectivity(self.atoms,self.B1,self.B2)
        self.calculate_geom_center()
        self.dist_unit='Ang'
        self.unique_species=list(set(self.species))
        self.AOM_dict={}
        self.orb_compl_dict={}
    def calculate_geom_center(self):
        """Geometric center calculation"""
        self.X,self.Y,self.Z=[np.array(self.x).mean(),np.array(self.y).mean(),np.array(self.z).mean()]
    def recenter(self,repeat=1):
        """Calculate geometric center; then resolve the radius R of a sphere that encompasses the whole molecule.
        Finally, center the molecule at (0,0,0) and shift by (1+2*repeat)*R.
        This way, the molecule is isolated in a CUBIC supercell with edge length 2*(1+2*repeat)*R."""
        self.calculate_geom_center()
        R=max([np.linalg.norm(np.array([self.x[i],self.y[i],self.z[i]])-np.array([self.X,self.Y,self.Z])) for i in range(self.atoms)])
        self.x=[i-self.X+(1+2*repeat)*R for i in self.x]
        self.y=[i-self.Y+(1+2*repeat)*R for i in self.y]
        self.z=[i-self.Z+(1+2*repeat)*R for i in self.z]
        self.calculate_geom_center()
        self.supercell=[2*(1+2*repeat)*R for i in range(3)]
    def angtobohr(self):
        """Convert Angstrom to Bohr"""
        AngToBohr=1.8897259886
        if self.dist_unit=='Ang':
            self.x=[i*AngToBohr for i in self.x]
            self.y=[i*AngToBohr for i in self.y]
            self.z=[i*AngToBohr for i in self.z]
            self.dist_unit='Bohr'
    def bohrtoang(self):
        """Convert Bohr to Angstrom"""
        AngToBohr=1.8897259886
        if self.dist_unit=='Bohr':
            self.x=[i/AngToBohr for i in self.x]
            self.y=[i/AngToBohr for i in self.y]
            self.z=[i/AngToBohr for i in self.z]
            self.dist_unit='Ang'
    def resolve_pvecs(self,STO_matrix=None):
        """Calculate normal vectors. If the STO matrix is given as input, calculate pi projection coeffs"""
        # convert ang to bohr
        self.angtobohr()
        # initialize normal vector components list
        self.px,self.py,self.pz=[np.zeros(self.atoms) for i in range(3)]
        aoinum=0
        atomlist=[]
        for j,i in enumerate(self.species):
            if i!='H':
                aoinum+=1
                atomlist.append(j)
        for i in range(len(atomlist)):
            neighbourlist=[[0 for j in range(4)] for i in range(4)]
            j=0
            for k in range(self.atoms):
                if atomlist[i]!=k:
                    r=np.linalg.norm(np.array([self.x[atomlist[i]],self.y[atomlist[i]],self.z[atomlist[i]]])-np.array([self.x[k],self.y[k],self.z[k]]))
                    if r<3.5:
                        neighbourlist[0][j]=k+1
                        neighbourlist[1][j]=self.x[k]
                        neighbourlist[2][j]=self.y[k]
                        neighbourlist[3][j]=self.z[k]
                        j+=1
            if j==2:
                neighbourlist[0][2]=atomlist[i]+1
                neighbourlist[1][2]=self.x[atomlist[i]]
                neighbourlist[2][2]=self.y[atomlist[i]]
                neighbourlist[3][2]=self.z[atomlist[i]]
            if j==1:
                # use as a second column entry the only 1-2 neighbor
                neighbourlist[0][1]=atomlist[i]+1
                neighbourlist[1][1]=self.x[atomlist[i]]
                neighbourlist[2][1]=self.y[atomlist[i]]
                neighbourlist[3][1]=self.z[atomlist[i]]
                for k in range(self.atoms):
                    if neighbourlist[0][0]!=k+1 and neighbourlist[0][1]!=k+1:
                        r=np.linalg.norm(np.array([self.x[neighbourlist[0][0]-1],self.y[neighbourlist[0][0]-1],self.z[neighbourlist[0][0]-1]])-np.array([self.x[k],self.y[k],self.z[k]]))
                        if r<3.5:
#                             print(f'r={r/1.8897259886} {neighbourlist[0][0]}-{k+1}')
                            neighbourlist[0][2]=k+1
                            neighbourlist[1][2]=self.x[k]
                            neighbourlist[2][2]=self.y[k]
                            neighbourlist[3][2]=self.z[k]
                            break
            veca_x,veca_y,veca_z=[neighbourlist[i][0]-neighbourlist[i][2] for i in [1,2,3]]
            vecb_x,vecb_y,vecb_z=[neighbourlist[i][1]-neighbourlist[i][2] for i in [1,2,3]]
            u=np.cross([veca_x,veca_y,veca_z],[vecb_x,vecb_y,vecb_z])
            u/=np.linalg.norm(u)
            self.px[atomlist[i]],self.py[atomlist[i]],self.pz[atomlist[i]]=u
            
#             print(f'[{atomlist[i]+1}] j={j}')
#             for c,k in enumerate(neighbourlist):
#                 if c>0:
#                     for l in k:
#                         print(f'{l/1.8897259886:>.3}\t',end='')
#                     print()
#                 else:
#                     for l in k:
#                         print(f'{l}\t',end='')
#                     print()
#             print('--------------------------------')
        if STO_matrix is not None:
            return [self.px[j]*i[1]+self.py[j]*i[2]+self.pz[j]*i[3] for j,i in enumerate(STO_matrix)]
    def initialize_STOs(self,STO_dict):
        self.angtobohr()
        self.STOs,self.STO_id_array,self.STO_type_array,self.STO_mu_array=initialize_STOs(self.species,self.x,self.y,self.z,STO_dict)
        self.Smatrix=calculate_overlap_S_matrix(self.x,self.y,self.z,
                                   self.STOs,self.STO_id_array,self.STO_type_array,self.STO_mu_array)
    def get_cp2k_info(self,filename,MO,driver):
        self.MO=MO
        self.MOcoeffs=get_cp2k_MO(filename,MO)
        self.basis_dict=read_basis(driver['cp2k_basis_file'],driver['basis'],self.unique_species)
        self.pcoeff,self.palpha,self.pqn,self.bfnPerAtom,self.GTO_depth=read_CP2K_GTOs(self.species,self.basis_dict)
    def project(self):
        self.STO_matrix,self.orb_compl=STO_GTO_projection(self.x,self.y,self.z,
                                        self.STO_id_array,self.STO_type_array,self.STO_mu_array,
                                        self.pcoeff,self.palpha,self.pqn,self.bfnPerAtom,self.GTO_depth,self.MOcoeffs,self.Smatrix)
        self.resolve_pvecs()
        self.AOM_pi_coeffs=[self.px[c]*i[1]+self.py[c]*i[2]+self.pz[c]*i[3] for c,i in enumerate(self.STO_matrix)]
        self.AOM_dict[self.MO]=self.AOM_pi_coeffs
        self.orb_compl_dict[self.MO]=self.orb_compl
    def cube_preview(self,STO_dict,filename):
        create_cube_file(self.species,self.x,self.y,self.z,self.STO_matrix,STO_dict,filename)
    def read_AOM_from_file(self,driver,mol_type,prefix='AOM_COEFF_molecule_type_',suffix='.dat',out_dir='out_dir'):
        filename=f'{driver[out_dir]}/{prefix}{mol_type}{suffix}'
        fp=open(filename,mode='r')
        aom_from_file=fp.read()
        fp.close
        self.AOM_dict=literal_eval(aom_from_file)

class System:
    """This is the class describing an atomistic system
    Basic attributes:
    - number of atoms/bonds/molecules
    - atomic info: species, xyz coordinates, and molecular ids
    - bonding info: B1 and B2 arrays
    - optional: periodic supercell as a/b/c vectors"""
    def __init__(self,filename):
        # upper values from DREIDING (Mayo1990)
        self.bond_radii={'H':0.330,'C':0.770,'N':0.702,'O':0.660,'F':0.611,'S':1.040}
        # reading xyz file
        if filename.split('.')[-1]=='dat' or filename.split('.')[-1]=='xyz':
            [self.species,self.x,self.y,self.z]=read_coords(filename)
            #[self.B1,self.B2,self.connectivity]=resolve_bonds(self.x,self.y,self.z)
            head,mylist,neighbors=initialize_linked_list_cells(self.x,self.y,self.z,cutoff=self.species_cutoff(self.bond_radii))
            [self.B1,self.B2,self.connectivity]=resolve_bonds_on_cells(head,mylist,neighbors,self.species,self.x,self.y,self.z,self.bond_radii)
            self.clusters=resolve_molecules(self.connectivity)
            self.atoms=len(self.species)
            self.bonds=len(self.B1)
            self.molecules=len(self.clusters)
            self.mol=[0 for i in range(self.atoms)]
            self.atoms_per_molecule=[0 for i in range(self.molecules)]
            for counter,i in enumerate(self.clusters):
                for j in i:
                    self.mol[j-1]=counter+1
                self.atoms_per_molecule[counter]=len(i)
        self.system_to_molecules()
        self.resolve_unique_molecule_types()
        self.name=filename.split('.'+filename.split('.')[-1])[0]
    # cutoff selector
    def species_cutoff(self,bond_radii,delta=0.1):
        unique_species=list(set(self.species))
        res=[]
        for i in range(len(unique_species)):
            for j in range(i,len(unique_species)):
                res.append(bond_radii[unique_species[i]]+bond_radii[unique_species[j]]+delta)
        return max(res)
    # breakdown the system to individual molecules        
    def system_to_molecules(self):
        self.molecule=[]
        for i in range(self.molecules):
            self.molecule.append(Molecule(self,i))
    # resolve unique molecular types based on the connectivity matrix and the species array
    def create_molecular_type_identifier(self,mol):
        return (tuple(self.molecule[mol].species),tuple(map(tuple,self.molecule[mol].connectivity)))
    def resolve_unique_molecule_types(self):
        unique=()
        for i in (self.create_molecular_type_identifier(i) for i in range(self.molecules)):
            if i not in unique:
                unique+=(i,)
        self.unique_dict={i:j+1 for j,i in enumerate(unique)}
        self.mol_type=[self.unique_dict[self.create_molecular_type_identifier(i)] for i in range(self.molecules)]
        self.mol_type_dict={i:j+1 for j,i in enumerate(self.mol_type)}
    def isolate_unique_mol_types(self):
        if os.path.exists(f'{self.name}')==False:
            os.mkdir(f'{self.name}')
        for key,value in self.mol_type_dict.items():
            with open(f'{self.name}/molecule_type_{key}.xyz',mode='w') as fp:
                print(f'{self.molecule[value-1].atoms}\n',file=fp)
                for i in range(self.molecule[value-1].atoms):
                    print(f'{self.molecule[value-1].species[i]}\t{self.molecule[value-1].x[i]}\t{self.molecule[value-1].y[i]}\t{self.molecule[value-1].z[i]}',file=fp)
    def prep_cp2k_single(self,repeat=1,template='sp_pbe_template.inp',basis='DZVP-GTH',potential='GTH-PBE'):
        for i in self.molecule:
            i.bohrtoang()
            i.recenter(repeat)
        self.isolate_unique_mol_types()
        if os.path.exists(f'{self.name}')==False:
            os.mkdir(f'{self.name}')
        fp=open(f'templates/{template}',mode='r')
        cp2k_input_ref=fp.readlines()
        fp.close()
        for key,value in self.mol_type_dict.items():
            filename=f'{self.name}/molecule_type_{key}_subsys.include'
            with open(filename,mode='w') as fp:
                print('&SUBSYS',file=fp)
                print(' &CELL',file=fp)
                print(f'  ABC {self.molecule[value-1].supercell[0]:.6f} {self.molecule[value-1].supercell[1]:.6f} {self.molecule[value-1].supercell[2]:.6f}',file=fp)
                print('  ALPHA_BETA_GAMMA 90.0 90.0 90.0',file=fp)
                print('  PERIODIC none',file=fp)
                print(' &END CELL',file=fp)
                unique_species=list(set(self.molecule[value-1].species))
                for i in unique_species:
                    print(f' &KIND {i}',file=fp)
                    print(f'  BASIS_SET {basis}',file=fp)
                    print(f'  POTENTIAL {potential}',file=fp)
                    print(' &END KIND',file=fp)
                print(' &COORD',file=fp)
                for i in range(self.molecule[value-1].atoms):
                    print(f'  {self.molecule[value-1].species[i]} {self.molecule[value-1].x[i]:.6f} {self.molecule[value-1].y[i]:.6f} {self.molecule[value-1].z[i]:.6f}',file=fp)
                print(' &END COORD',file=fp)
                print('&END SUBSYS',file=fp)
            print(f'Generated {filename}')
            cp2k_input=cp2k_input_ref.copy()
            for j,i in enumerate(cp2k_input):
                if i.find('_GLOBAL_PROJECT_NAME_')!=-1:
                    cp2k_input[j]=i.replace('_GLOBAL_PROJECT_NAME_',f'molecule_type_{key}')
                if i.find('_MGRID_CUTOFF_')!=-1:
                    cp2k_input[j]=i.replace('_MGRID_CUTOFF_',f'{450.0}')
                if i.find('_MGRID_REL_CUTOFF_')!=-1:
                    cp2k_input[j]=i.replace('_MGRID_REL_CUTOFF_',f'{75.0}')
                if i.find('_SCF_MAX_SCF_')!=-1:
                    cp2k_input[j]=i.replace('_SCF_MAX_SCF_',f'{200}')
                if i.find('_MO_CUBES_NHOMO_')!=-1:
                    cp2k_input[j]=i.replace('_MO_CUBES_NHOMO_',f'{5}')
                if i.find('_MO_CUBES_NLUMO_')!=-1:
                    cp2k_input[j]=i.replace('_MO_CUBES_NLUMO_',f'{5}')
                if i.find('_MO_CUBES_STRIDE_')!=-1:
                    cp2k_input[j]=i.replace('_MO_CUBES_STRIDE_',f'{5}')
            filename=f'{self.name}/molecule_type_{key}.inp'
            with open(filename,mode='w') as fp:
                for i in cp2k_input:
                    print(i,end='',file=fp)
            print(f'Generated {filename}')
    def get_unique_molecular_types(self):
        return self.unique_dict        
    def apply_AOM_types(self):
        for j,i in enumerate(self.mol_type):
            if self.molecule[j].AOM_dict=={}:
                self.molecule[j].AOM_dict=self.molecule[self.mol_type_dict[i]-1].AOM_dict.copy()
    def store_AOM_types(self,driver):
        for key,value in self.mol_type_dict.items():
            filename=f'{driver["out_dir"]}/AOM_COEFF_molecule_type_{key}.dat'
            with open(filename,mode='w') as fp:
                print(self.molecule[value-1].AOM_dict,file=fp)
    def resolve_mol_type_driver(self,driver,out_dir="out_dir",prefix="molecule_type_",suffix=".out"):
        self.mol_type_driver={}
        for key,value in self.mol_type_dict.items():
            filename=f'{driver[out_dir]}/{prefix}{key}{suffix}'
            fp=open(filename,mode='r')
            cp2k_out=fp.readlines()
            fp.close()
            cp2k_out=[i for i in cp2k_out if i!='\n']
            electrons=[int(i.strip().split(':')[-1]) for i in cp2k_out if i.find('Number of electrons')!=-1][0]
            self.mol_type_driver[key]={
                'id':value-1,
                'cp2k_out':filename,
                'MO':[electrons,electrons+1]
            }
        return self.mol_type_driver
    def apply_STO_projection(self,STO_dict,driver,cube_preview=False):
        for key,value in self.mol_type_driver.items():
            for MO in value['MO']:
                self.molecule[value['id']].initialize_STOs(STO_dict)
                self.molecule[value['id']].get_cp2k_info(value['cp2k_out'],MO,driver)
                self.molecule[value['id']].project()
                if cube_preview==True:
                    self.molecule[value['id']].cube_preview(STO_dict,f'{value["cp2k_out"].split(".out")[0]}_MO_{MO}.cube')
    def resolve_pvecs(self):
        for i in range(self.molecules):
            self.molecule[i].resolve_pvecs()
    def calculate_AOM_overlap(self,AOM_overlap_dict,MO1,MO2,mol1=1,mol2=2):
        self.resolve_pvecs()
        mol1-=1
        mol2-=1
        STOs_f1,STO_id_array_f1,STO_type_array_f1,STO_mu_array_f1=initialize_STOs(self.molecule[mol1].species,
                                                                          self.molecule[mol1].x,self.molecule[mol1].y,self.molecule[mol1].z,
                                                                          AOM_overlap_dict)
        STO_matrix_f1=resolve_STO_matrix(len(self.molecule[mol1].x),STOs_f1,STO_id_array_f1,STO_type_array_f1,
                           self.molecule[mol1].px,
                           self.molecule[mol1].py,
                           self.molecule[mol1].pz,
                        self.molecule[mol1].AOM_dict[MO1])
        STOs_f2,STO_id_array_f2,STO_type_array_f2,STO_mu_array_f2=initialize_STOs(self.molecule[mol2].species,
                                                                          self.molecule[mol2].x,self.molecule[mol2].y,self.molecule[mol2].z,
                                                                          AOM_overlap_dict)
        STO_matrix_f2=resolve_STO_matrix(len(self.molecule[mol2].x),STOs_f2,STO_id_array_f2,STO_type_array_f2,
                           self.molecule[mol2].px,
                           self.molecule[mol2].py,
                           self.molecule[mol2].pz,
                        self.molecule[mol2].AOM_dict[MO2])
        S_f1=AOM_overlap_calculation(0,STOs_f1,
                        0,STOs_f1,
                        self.molecule[mol1].x+self.molecule[mol2].x,
                        self.molecule[mol1].y+self.molecule[mol2].y,
                        self.molecule[mol1].z+self.molecule[mol2].z,
                        STO_id_array_f1+[i+len(self.molecule[mol1].x) for i in STO_id_array_f2],
                        STO_type_array_f1+STO_type_array_f2,
                        STO_mu_array_f1+STO_mu_array_f2,
                        STO_matrix_f1+STO_matrix_f2)   
        for ci,i in enumerate(STO_matrix_f1):
            for j in range(1,3+1):
                STO_matrix_f1[ci][j]/=math.sqrt(abs(S_f1))
        S_f1_norm=AOM_overlap_calculation(0,STOs_f1,
                        0,STOs_f1,
                        self.molecule[mol1].x+self.molecule[mol2].x,
                        self.molecule[mol1].y+self.molecule[mol2].y,
                        self.molecule[mol1].z+self.molecule[mol2].z,
                        STO_id_array_f1+[i+len(self.molecule[mol1].x) for i in STO_id_array_f2],
                        STO_type_array_f1+STO_type_array_f2,
                        STO_mu_array_f1+STO_mu_array_f2,
                        STO_matrix_f1+STO_matrix_f2)   
        S_f2=AOM_overlap_calculation(STOs_f1,STOs_f1+STOs_f2,
                        STOs_f1,STOs_f1+STOs_f2,
                        self.molecule[mol1].x+self.molecule[mol2].x,
                        self.molecule[mol1].y+self.molecule[mol2].y,
                        self.molecule[mol1].z+self.molecule[mol2].z,
                        STO_id_array_f1+[i+len(self.molecule[mol1].x) for i in STO_id_array_f2],
                        STO_type_array_f1+STO_type_array_f2,
                        STO_mu_array_f1+STO_mu_array_f2,
                        STO_matrix_f1+STO_matrix_f2)   
        for ci,i in enumerate(STO_matrix_f2):
            for j in range(1,3+1):
                STO_matrix_f2[ci][j]/=math.sqrt(abs(S_f2))
        S_f2_norm=AOM_overlap_calculation(STOs_f1,STOs_f1+STOs_f2,
                        STOs_f1,STOs_f1+STOs_f2,
                        self.molecule[mol1].x+self.molecule[mol2].x,
                        self.molecule[mol1].y+self.molecule[mol2].y,
                        self.molecule[mol1].z+self.molecule[mol2].z,
                        STO_id_array_f1+[i+len(self.molecule[mol1].x) for i in STO_id_array_f2],
                        STO_type_array_f1+STO_type_array_f2,
                        STO_mu_array_f1+STO_mu_array_f2,
                        STO_matrix_f1+STO_matrix_f2)
        Sab=AOM_overlap_calculation(0,STOs_f1,
                        STOs_f1,STOs_f1+STOs_f2,
                        self.molecule[mol1].x+self.molecule[mol2].x,
                        self.molecule[mol1].y+self.molecule[mol2].y,
                        self.molecule[mol1].z+self.molecule[mol2].z,
                        STO_id_array_f1+[i+len(self.molecule[mol1].x) for i in STO_id_array_f2],
                        STO_type_array_f1+STO_type_array_f2,
                        STO_mu_array_f1+STO_mu_array_f2,
                        STO_matrix_f1+STO_matrix_f2)
        return S_f1,S_f1_norm,S_f2,S_f2_norm,Sab
        
def read_coords(filename):
    """This is a generic coordinates reader
    The type of file is deduced from the suffix:
    - .dat or .xyz are read as xyz files
    - .mol2 are read as mol2 files ** currently NOT implemented
    For the mol2 files, the function must look for the @<TRIPOS>ATOM and 
    @<TRIPOS>BOND sections, read the number of atoms/bonds/molecules from
    the third line, and - if present - read the @<TRIPOS>CRYSIN section"""
    suffix=filename.split('.')[-1]
    if suffix=='dat' or suffix=='xyz':
        with open(filename,mode='r') as fp:
            data=fp.readlines()
        atoms=int(data[0])
        species=[data[i+2].split()[0] for i in range(atoms)]
        [x,y,z]=[[float(data[i+2].split()[j]) for i in range(atoms)] for j in [1,2,3]]
        [x,y,z]=[np.array(x),np.array(y),np.array(z)]
        return species,x,y,z
def resolve_connectivity(atoms,B1,B2):
    bonds=len(B1)
    connectivity=[[] for i in range(atoms)]
    for i in range(bonds):
        connectivity[B1[i]-1].append(B2[i])
        connectivity[B2[i]-1].append(B1[i])
    for i in range(atoms):
        connectivity[i].sort()
    return connectivity
def resolve_bonds(x,y,z,cutoff=2.0):
    atoms=len(x)
    B1=[]
    B2=[]
    for i in range(atoms-1):
        for j in range(i,atoms):
            if i!=j:
                r=np.linalg.norm(np.array([x[j],y[j],z[j]])-np.array([x[i],y[i],z[i]]))
                if r<cutoff:
                    #print(f'{i+1}-{j+1}: {r}')
                    B1.append(i+1)
                    B2.append(j+1)
    connectivity=resolve_connectivity(atoms,B1,B2)
    return B1,B2,connectivity    
def resolve_molecules(connectivity):
    atoms=len(connectivity)
    ignore=[False for i in range(atoms)]
    clusters=[]
    c_counter=-1
    while ignore != [True for i in range(atoms)]:
        # initialize the clusters list
        for counter,i in enumerate(ignore):
            if i==False:
                clusters.append(connectivity[counter].copy())
                ignore[counter]=True
                c_counter+=1
                break
        # find entries and store    
        for i in clusters[c_counter]:
            if ignore[i-1]==False:
                clusters[c_counter]+=[j for j in connectivity[i-1] if j not in clusters[c_counter]]
                ignore[i-1]=True
    molecules=len(clusters)
    for i in range(molecules):
        clusters[i].sort()
    return clusters
def initialize_linked_list_cells(x,y,z,cutoff,Mx=None,My=None,Mz=None):
    # this is a simple implementation of the linked-list cell distance calculation method
    
    #---- no PBCs ---------------------------------------------------
    
    # for non-periodic systems, we need to define a cell offset: deltar
    deltar=0.001
    xlo,xhi=[x.min()-deltar,x.max()+deltar]
    ylo,yhi=[y.min()-deltar,y.max()+deltar]
    zlo,zhi=[z.min()-deltar,z.max()+deltar]
    lx,ly,lz=[xhi-xlo,yhi-ylo,zhi-zlo]
    
    # resolve cells per cartesian direction
    if Mx is None and My is None and Mz is None:
        Mx,My,Mz=[math.floor(lx/cutoff),math.floor(ly/cutoff),math.floor(lz/cutoff)]
        if Mx<3: Mx=3
        if My<3: My=3
        if Mz<3: Mz=3
    lcx,lcy,lcz=[lx/Mx,ly/My,lz/Mz]
    
    # augment boundaries for small values
    if lcx<cutoff:
        xlo,xhi=[x.min()-(deltar+cutoff/2),x.max()+deltar+cutoff/2]
        lx=xhi-xlo
        Mx=math.floor(lx/cutoff)
        lcx=lx/Mx
    if lcy<cutoff:
        ylo,yhi=[y.min()-(deltar+cutoff/2),y.max()+deltar+cutoff/2]
        ly=yhi-ylo
        My=math.floor(ly/cutoff)
        lcy=ly/My
    if lcz<cutoff:
        zlo,zhi=[z.min()-(deltar+cutoff/2),z.max()+deltar+cutoff/2]
        lz=zhi-zlo
        Mz=math.floor(lz/cutoff)
        lcz=lz/Mz
    
    # build the neighbors
    M2=Mx*My
    M3=M2*Mz
    
    # introduce a basal x-y plane with cell indices, surrounded by ghost cells for PBCs
    ghostsx=Mx+2
    ghostsy=My+2
    cells=[[0 for j in range(ghostsx)] for i in range(ghostsy)]
    # populate the inner part - not the ghost cells
    l=0
    for i in range(My-1,-1,-1):
        for j in range(Mx):
            l+=1
            cells[i+1][j+1]=l
    
    # periodicity
    # for j in range(1,ghostsx-2+1):
    #     cells[0][j]=cells[ghostsy-2][j]
    #     cells[ghostsy-1][j]=cells[1][j]
    # for i in range(0,ghostsy-1+1):
    #     cells[i][ghostsx-1]=cells[i][1]
    #     cells[i][0]=cells[i][ghostsx-2]

    # no periodicity: designate ghost cells as 'na' (not available)
    for j in range(1,ghostsx-2+1):
        cells[0][j]='na'
        cells[ghostsy-1][j]='na'
    for i in range(0,ghostsy-1+1):
        cells[i][ghostsx-1]='na'
        cells[i][0]='na'
        
    # prepare the neighbors 2D matrix: will store the indices of neighboring cells per active cell
    # ! neighboring cell info is stored row-wise !
    neighbors=[[0 for j in range(M3)] for i in range(26)]
    
    # work on the basal plane
    l=0
    for i in range(My-1,-1,-1):
        for j in range(Mx):
            row=i+2
            column=j+2
            k=0
            for jj in range(3):
                for ii in range(3):
                    if ii!=1 or jj!=1:
                        neighbors[k][l]=cells[ii-2+row][jj-2+column]
                        k+=1
            l+=1
    
    # print for debugging purposes
#     for i in range(M3):
#         print(f'{i+1:>3} ',end='')
#     print()
#     for i in range(M3):
#         print(f'{"---":>3} ',end='')
#     print()
#     for i in neighbors:
#         for j in i:
#             print(f'{j:>3} ',end='')
#         print()
#     for i in range(M3):
#         print(f'{"===":>3} ',end='')
#     print()
    # end of print
    
    # populate all other planes up to top
    for j in range(M2,M3-1+1):
        for i in range(0,7+1):
            if neighbors[i][j-M2]!='na':
                neighbors[i][j]=neighbors[i][j-M2]+M2
            else:
                neighbors[i][j]='na'
    
    # print for debugging purposes
#     for i in range(M3):
#         print(f'{i+1:>3} ',end='')
#     print()
#     for i in range(M3):
#         print(f'{"---":>3} ',end='')
#     print()
#     for i in neighbors:
#         for j in i:
#             print(f'{j:>3} ',end='')
#         print()
#     for i in range(M3):
#         print(f'{"===":>3} ',end='')
#     print()
    # end of print
    
    for j in range(M3-1+1):
        for i in range(7+1):
            # upper peripheral
            if neighbors[i][j]!='na':
                if neighbors[i][j]<=M3-M2:
                    neighbors[i+8][j]=neighbors[i][j]+M2
                else:
                    # boundary
                    #neighbors[i+8][j]=neighbors[i][j]-(Mz-1)*M2
                    neighbors[i+8][j]='na'
            else:
                neighbors[i+8][j]='na'
            # lower peripheral
            if neighbors[i][j]!='na':                
                if neighbors[i][j]>M2:
                    neighbors[i+16][j]=neighbors[i][j]-M2
                else:
                    # boundary
                    #neighbors[i+16][j]=neighbors[i][j]+(Mz-1)*M2
                    neighbors[i+16][j]='na'
            else:
                neighbors[i+16][j]='na'
        # directly above
        if j+1+M2<=M3:
            neighbors[24][j]=j+1+M2
        else:
            # boundary
            #neighbors[24][j]=j+1-(Mz-1)*M2
            neighbors[24][j]='na'
        # directly below
        if j+1-M2>0:
            neighbors[25][j]=j+1-M2
        else:
            # boundary 
            #neighbors[25][j]=j+1+(Mz-1)*M2
            neighbors[25][j]='na'
    
    # print for debugging purposes
#     for i in range(M3):
#         print(f'{i+1:>3} ',end='')
#     print()
#     for i in range(M3):
#         print(f'{"---":>3} ',end='')
#     print()
#     for i in neighbors:
#         for j in i:
#             print(f'{j:>3} ',end='')
#         print()
#     for i in range(M3):
#         print(f'{"===":>3} ',end='')
#     print()
    # end of print
    
    # eliminate double-counting
    
    for j in range(M3):
        for i in range(26):
            if neighbors[i][j]!='na':
                if neighbors[i][j]<j+1:
                    neighbors[i][j]='na'
    
    # print for debugging purposes
#     for i in range(M3):
#         print(f'{i+1:>3} ',end='')
#     print()
#     for i in range(M3):
#         print(f'{"---":>3} ',end='')
#     print()
#     for i in neighbors:
#         for j in i:
#             print(f'{j:>3} ',end='')
#         print()
#     for i in range(M3):
#         print(f'{"===":>3} ',end='')
#     print()
    # end of print
    
    for i in range(26):
        for j in range(M3):
            if neighbors[i][j]!='na':
                neighbors[i][j]=neighbors[i][j]-1
    particles=len(x)
    head=[-1 for i in range(M3)]
    mylist=[-1 for i in range(particles)]
    for i in range(particles):
        cindex=math.floor((x[i]-xlo)/lcx)+math.floor((y[i]-ylo)/lcy)*Mx+math.floor((z[i]-zlo)/lcz)*Mx*My
        mylist[i]=head[cindex]
        head[cindex]=i
    return head,mylist,neighbors

def resolve_bonds_on_cells(head,mylist,neighbors,species,x,y,z,bond_radii,delta=0.1):
    
    B1=[]
    B2=[]
    
    particles=len(mylist)
    M3=len(head)
    for cellindex in range(M3):
        cellpartarray=[]
        cellpart=0
        headpointer=head[cellindex]
        if headpointer!=-1:
            cellpart+=1
            cellpartarray.append(headpointer)
            listpointer=mylist[headpointer]
            while listpointer!=-1:
                if listpointer!=-1:
                    cellpart+=1
                    cellpartarray.append(listpointer)
                    listpointer=mylist[listpointer]
        if cellpart!=0:
            cellpartarray.sort()
            for i in range(cellpart-1):
                for j in range(i,cellpart):
                    if i!=j:
                        locali=cellpartarray[i]
                        localj=cellpartarray[j]
                        r=np.linalg.norm(np.array([x[localj],y[localj],z[localj]])-np.array([x[locali],y[locali],z[locali]]))
                        cutoff=bond_radii[species[locali]]+bond_radii[species[localj]]+delta
                        if r<cutoff:
                            B1.append(locali+1)
                            B2.append(localj+1)

        neighborcellpartarray=[]
        neighborcellpart=0
        for neighborcellindex in range(26):
            if neighbors[neighborcellindex][cellindex]!='na':
                headpointer=head[neighbors[neighborcellindex][cellindex]]
                if headpointer!=-1:
                    neighborcellpart+=1
                    neighborcellpartarray.append(headpointer)
                    listpointer=mylist[headpointer]
                    while listpointer!=-1:
                        if listpointer!=-1:
                            neighborcellpart+=1
                            neighborcellpartarray.append(listpointer)
                            listpointer=mylist[listpointer]
        if neighborcellpart!=0:
            neighborcellpartarray.sort()
            for locali in cellpartarray:
                for localj in neighborcellpartarray:
                    r=np.linalg.norm(np.array([x[localj],y[localj],z[localj]])-np.array([x[locali],y[locali],z[locali]]))
                    cutoff=bond_radii[species[locali]]+bond_radii[species[localj]]+delta
                    if r<cutoff:
                        B1.append(locali+1)
                        B2.append(localj+1)

    connectivity=resolve_connectivity(particles,B1,B2)
    return B1,B2,connectivity

# initialize STOs
def initialize_STOs(species,x,y,z,STO_dict,debug=0):
    one_species=[]
    two_species=[]
    three_species=[]
    for key,value in STO_dict.items():
        for i in value:
            if i.find('1')!=-1:
                one_species.append(key)
                break
        for i in value:
            if i.find('2')!=-1:
                two_species.append(key)
                break
        for i in value:
            if i.find('3')!=-1:
                three_species.append(key)
                break
    atoms=len(species)
    STO_type_string={1:'1s',2:'2s',3:'2px',4:'2py',5:'2pz',6:'3s',7:'3px',8:'3py',9:'3pz'}
    STOs=sum([STO_dict[i]['STOs'] for i in species])
    STO_id_array=[]
    STO_type_array=[]
    STO_mu_array=[]
    j=-1
    for i in range(atoms):
        if species[i] in one_species:
            j+=1
            STO_id_array.append(i+1)
            STO_type_array.append(1)
            STO_mu_array.append(STO_dict[species[i]][STO_type_string[1][0:2]])
        if species[i] in two_species:
            j+=1
            STO_id_array.append(i+1)
            STO_type_array.append(2)
            STO_mu_array.append(STO_dict[species[i]][STO_type_string[2][0:2]])
            j+=1
            STO_id_array.append(i+1)
            STO_type_array.append(3)
            STO_mu_array.append(STO_dict[species[i]][STO_type_string[3][0:2]])
            j+=1
            STO_id_array.append(i+1)
            STO_type_array.append(4)
            STO_mu_array.append(STO_dict[species[i]][STO_type_string[4][0:2]])
            j+=1
            STO_id_array.append(i+1)
            STO_type_array.append(5)
            STO_mu_array.append(STO_dict[species[i]][STO_type_string[5][0:2]])
        if species[i] in three_species:
            j+=1
            STO_id_array.append(i+1)
            STO_type_array.append(6)
            STO_mu_array.append(STO_dict[species[i]][STO_type_string[6][0:2]])
            j+=1
            STO_id_array.append(i+1)
            STO_type_array.append(7)
            STO_mu_array.append(STO_dict[species[i]][STO_type_string[7][0:2]])
            j+=1
            STO_id_array.append(i+1)
            STO_type_array.append(8)
            STO_mu_array.append(STO_dict[species[i]][STO_type_string[8][0:2]])
            j+=1
            STO_id_array.append(i+1)
            STO_type_array.append(9)
            STO_mu_array.append(STO_dict[species[i]][STO_type_string[9][0:2]])
    if debug==1:
        print(f'Number of atoms: {atoms}')
        print('[Atomic id], species and coordinates (Ang|Bohr):')
        for i in range(atoms):
            print(f'[{i+1}] {species[i]} {x[i]/1.8897259886:>10.6f} {y[i]/1.8897259886:>10.6f} {z[i]/1.8897259886:>10.6f} | {x[i]:>10.6f} {y[i]:>10.6f} {z[i]:>10.6f}')
        print('--------------------------------------------------------------------------')
        print(f'Number of total STO basis functions: {STOs}\nAtomic decomposition:\n[id]\tSTO\torbital\tmu\ttype')
        for j,i in enumerate(STO_id_array):
            print(f'[{i}]\t{j+1}\t{STO_type_string[STO_type_array[j]]}\t{STO_mu_array[j]}\t{STO_type_array[j]}')
    return STOs,STO_id_array,STO_type_array,STO_mu_array
def calculate_overlap_S_matrix(x,y,z,STOs,STO_id_array,STO_type_array,STO_mu_array):
    Smatrix=np.identity(STOs)
    for i in range(STOs):
        for j in range(STOs):
            if STO_id_array[i]!=STO_id_array[j]:
                Smatrix[i][j]=overlap(x[STO_id_array[i]-1],y[STO_id_array[i]-1],z[STO_id_array[i]-1],x[STO_id_array[j]-1],y[STO_id_array[j]-1],z[STO_id_array[j]-1],STO_mu_array[i],STO_mu_array[j],STO_type_array[i],STO_type_array[j])
    return Smatrix
def get_cp2k_MO(filename,MO):
    # open cp2k output file and store MO info: alpha channel
    fp=open(filename,mode='r')
    cp2k_out=fp.readlines()
    fp.close()
    cp2k_out=[i for i in cp2k_out if i!='\n']
    cart_basis_functions=[int(i.strip().split(':')[-1]) for i in cp2k_out if i.find('Cartesian basis functions')!=-1][0]
    orbitals=[int(i.strip().split(':')[-1]) for i in cp2k_out if i.find('Number of molecular orbitals')!=-1]
    header=[j for j,i in enumerate(cp2k_out) if i.find('MO EIGENVALUES')!=-1]
    blocks=math.ceil(orbitals[0]/4)
    alphaMOs=cp2k_out[header[0]+1:header[0]+1+blocks*(cart_basis_functions+3)]
    alphaMOs=[i.strip().split() for i in alphaMOs]
    frames=[[] for i in range(blocks)]
    basis_info_block=[['MO'],['Energy'],['Occupation']]
    for i in range(blocks):
        start=3+i*(cart_basis_functions+3)
        stop=start+cart_basis_functions-1
        for j in range(3):
            frames[i].append(alphaMOs[start-3+j].copy())
        for j in range(start,stop+1):
            c=alphaMOs[j].copy()
            if i==0:
                basis_info_block.append(' '.join(c[0:4]))
            del c[0:4]
            frames[i].append(c)
    [i for i in frames[0]]
    alpha_records=[]
    for i in range(3,len(frames[0])):
        line=[]
        for j in range(blocks):
            line+=frames[j][i]
        alpha_records.append(line)
    alpha_dataframe=pd.DataFrame(alpha_records)
    MOcoeffs=[float(i) for i in alpha_dataframe[MO-1].to_list()]
    return MOcoeffs
def read_basis(filename,basis,unique_species,debug=0):
    fp=open(filename,mode='r')
    cp2k_basis_sets=fp.readlines()
    fp.close

    basis_dict={}

    for element in unique_species:
        basis_dict[element]={}
        if debug==1:
            print(f'*** Species: {element}\n*** Basis set breakdown:\n')
        counter=0
        for j,i in enumerate(cp2k_basis_sets):
            line=i.strip().split()
            if line[0]==element and line[2]==basis:
                current=j
                nset=int(cp2k_basis_sets[j+1])
                current+=1
                for k in range(nset):
                    # read set header
                    current+=1
                    if debug==1:
                        print(cp2k_basis_sets[current],end='')
                    lmin,lmax,nexp=[int(cp2k_basis_sets[current].split()[i]) for i in [1,2,3]]
                    nshell=[int(cp2k_basis_sets[current].split()[i+4]) for i in range(lmax-lmin+1)]
                    set_matrix=[]
                    for l in range(nexp):
                        current+=1
                        set_matrix.append([float(i) for i in cp2k_basis_sets[current].split()])
                    if debug==1:
                        for l in set_matrix:
                            for m in l:
                                print(f'{m}\t',end='')
                            print()
                        print(f'lmin={lmin}')
                        print(f'lmax={lmax}')
                    c_counter=0
                    for l in range(lmin,lmax+1):
                        if debug==1:
                            print(f'l={l} has nshell={nshell[l-lmin]} contractions')
                        for n in range(nshell[l-lmin]):
                            c_counter+=1
                            for m in range(-l,l+1):
                                counter+=1
                                if debug==1:
                                    print(f'gto_{element}_{counter} = r^{l} * (',end='')
                                coeff=[]
                                alpha=[]
                                for kk in range(nexp):
                                    if abs(set_matrix[kk][c_counter])>1.0e-9:
                                        if debug==1:
                                            print(f'{set_matrix[kk][c_counter]} * exp(-{set_matrix[kk][0]}*r^2) + ',end='')
                                    coeff.append(set_matrix[kk][c_counter])
                                    alpha.append(set_matrix[kk][0])
                                if debug==1:
                                    print(f') * Y({l},{m})')
                                basis_dict[element][counter]={'l':l,'m':m,'coeff':coeff,'alpha':alpha}
        if debug==1:
            print('\n--------------------------------------------------------------------------\n')
        # GTOs: spherical
        basis_dict[element]['GTOs']=counter
        # CGTOs: cartesian - up to l=2 supported right now...
        basis_dict[element]['CGTOs']=basis_dict[element]['GTOs']
        for i in range(counter):
            if basis_dict[element][i+1]['m']==-2:
                basis_dict[element]['CGTOs']+=1
    return basis_dict
def read_CP2K_GTOs(species,basis_dict):
    atoms=len(species)
    GTOs=sum([basis_dict[i]['GTOs'] for i in species])
    CGTOs=sum([basis_dict[i]['CGTOs'] for i in species])
    GTO_depth=[]
    pcoeff=[]
    palpha=[]
    pqn=[]
    bfnPerAtom=[]
    for i in species:
        bfnPerAtom.append(basis_dict[i]['CGTOs'])
        for j in range(basis_dict[i]['GTOs']):
            if basis_dict[i][j+1]['m']==0:
                # s type
                if basis_dict[i][j+1]['l']==0:
                    GTO_depth.append(len(basis_dict[i][j+1]['alpha']))
                    for l in basis_dict[i][j+1]['coeff']:
                        pcoeff.append(l)
                    for l in basis_dict[i][j+1]['alpha']:
                        palpha.append(l)
                    for l in basis_dict[i][j+1]['alpha']:
                        pqn.append([0,0,0])
                # p type
                if basis_dict[i][j+1]['l']==1:
                    for k in range(3):
                        GTO_depth.append(len(basis_dict[i][j+1]['alpha']))
                        for l in basis_dict[i][j+1]['coeff']:
                            pcoeff.append(l)
                        for l in basis_dict[i][j+1]['alpha']:
                            palpha.append(l)
                    for l in basis_dict[i][j+1]['alpha']:
                        pqn.append([1,0,0])
                    for l in basis_dict[i][j+1]['alpha']:
                        pqn.append([0,1,0])
                    for l in basis_dict[i][j+1]['alpha']:
                        pqn.append([0,0,1])
                # d type
                if basis_dict[i][j+1]['l']==2:
                    for k in range(6):
                        GTO_depth.append(len(basis_dict[i][j+1]['alpha']))
                        for l in basis_dict[i][j+1]['coeff']:
                            pcoeff.append(l)
                        for l in basis_dict[i][j+1]['alpha']:
                            palpha.append(l)
                    for l in basis_dict[i][j+1]['alpha']:
                        pqn.append([2,0,0])
                    for l in basis_dict[i][j+1]['alpha']:
                        pqn.append([1,1,0])
                    for l in basis_dict[i][j+1]['alpha']:
                        pqn.append([1,0,1])
                    for l in basis_dict[i][j+1]['alpha']:
                        pqn.append([0,2,0])
                    for l in basis_dict[i][j+1]['alpha']:
                        pqn.append([0,1,1])
                    for l in basis_dict[i][j+1]['alpha']:
                        pqn.append([0,0,2])
    return pcoeff,palpha,pqn,bfnPerAtom,GTO_depth
def STO_GTO_projection(x,y,z,STO_id_array,STO_type_array,STO_mu_array,pcoeff,palpha,pqn,bfnPerAtom,GTO_depth,MOcoeffs,Smatrix):
    projection_dict={
        1:{
            'oohatac':[0.0154199, 0.203216, 0.401002, 0.313305, 0.144647, 0.049609, 0.0139344, 0.00318156, 0.000622761, 0.0000886225],
            'oohataa':[0.0372356, 0.0848076, 0.191073, 0.458596, 1.18834, 3.34954, 10.514, 37.9276, 156.411, 1188.35],
            'Lx':0,
            'Ly':0,
            'Lz':0,
        },
        2:{
            'oohatac':[0.00578359, 0.225661, 0.527807, 0.265263, 0.0457229, -0.0488578, -0.0223909, -0.00600006, -0.00120626, -0.000163319],
            'oohataa':[0.0214081, 0.0497851, 0.0991734, 0.2010660, 0.263514, 1.3073, 3.98351, 14.1623, 60.6176, 432.099],
            'Lx':0,
            'Ly':0,
            'Lz':0,
        },
        3:{
            'oohatac':[0.162114, 0.490328, 0.362480, 0.118961, 0.0237179, 0.00306629],
            'oohataa':[0.0554971, 0.127302, 0.31, 0.834323, 2.56, 10.3],
            'Lx':1,
            'Ly':0,
            'Lz':0,
        },
        4:{
            'oohatac':[0.162114, 0.490328, 0.362480, 0.118961, 0.0237179, 0.00306629],
            'oohataa':[0.0554971, 0.127302, 0.31, 0.834323, 2.56, 10.3],
            'Lx':0,
            'Ly':1,
            'Lz':0,
        },
        5:{
            'oohatac':[0.162114, 0.490328, 0.362480, 0.118961, 0.0237179, 0.00306629],
            'oohataa':[0.0554971, 0.127302, 0.31, 0.834323, 2.56, 10.3],
            'Lx':0,
            'Ly':0,
            'Lz':1,
        },
        6:{
            'oohatac':[0.00621804, 0.434003, 0.678264, 0.00507343, -0.163049, -0.049145, -0.00574937, -0.000273434],
            'oohataa':[0.014179, 0.0408735, 0.0803231, 0.201715, 0.376225, 0.995934, 3.2689, 15.8],
            'Lx':0,
            'Ly':0,
            'Lz':0,
        },
        7:{
            'oohatac':[0.0196080, 0.324519, 0.521925, 0.207248, 0.00504818, -0.0139732, -0.00394290, -0.000523609],
            'oohataa':[0.0247242, 0.0510310, 0.100108, 0.207164, 0.424917, 1.24038, 3.81650, 15.8360],
            'Lx':1,
            'Ly':0,
            'Lz':0,
        },
        8:{
            'oohatac':[0.0196080, 0.324519, 0.521925, 0.207248, 0.00504818, -0.0139732, -0.00394290, -0.000523609],
            'oohataa':[0.0247242, 0.0510310, 0.100108, 0.207164, 0.424917, 1.24038, 3.81650, 15.8360],
            'Lx':0,
            'Ly':1,
            'Lz':0,
        },
        9:{
            'oohatac':[0.0196080, 0.324519, 0.521925, 0.207248, 0.00504818, -0.0139732, -0.00394290, -0.000523609],
            'oohataa':[0.0247242, 0.0510310, 0.100108, 0.207164, 0.424917, 1.24038, 3.81650, 15.8360],
            'Lx':0,
            'Ly':0,
            'Lz':1,
        },
    }
    V_array=[]
    for isto in range(len(STO_id_array)):
        S=0
        mm=-1
        for ll in range(len(projection_dict[STO_type_array[isto]]['oohatac'])):
            mm+=1
            k=-1
            m=-1
            for i in range(len(x)):
                for j in range(bfnPerAtom[i]):
                    k+=1
                    for l in range(GTO_depth[k]):
                        m+=1
                        integral=anIres(
                                x[STO_id_array[isto]-1],
                                y[STO_id_array[isto]-1],
                                z[STO_id_array[isto]-1],
                                STO_mu_array[STO_type_array[isto]-1]**2*projection_dict[STO_type_array[isto]]['oohataa'][mm],
                                projection_dict[STO_type_array[isto]]['Lx'],
                                projection_dict[STO_type_array[isto]]['Ly'],
                                projection_dict[STO_type_array[isto]]['Lz'],
                                x[i],
                                y[i],
                                z[i],
                                palpha[m],
                                pqn[m][0],
                                pqn[m][1],
                                pqn[m][2])
                        S+=projection_dict[STO_type_array[isto]]['oohatac'][mm]*MOcoeffs[k]*pcoeff[m]*integral
        V_array.append(S)
    STOs=len(STO_id_array)
    res=np.linalg.solve(Smatrix,V_array).round(decimals=12)
    orb_compl=sum([res[i]*res[j]*Smatrix[i][j] for i in range(STOs) for j in range(STOs)])
    res/=math.sqrt(abs(orb_compl))
    STO_matrix=[[0,0,0,0] for i in range(len(x))]
    for i in range(STOs):
        if STO_type_array[i]==1:
            STO_matrix[STO_id_array[i]-1][0]=res[i]
        else:
            STO_matrix[STO_id_array[i]-1][(STO_type_array[i]-2)%4]=res[i]
    return STO_matrix,orb_compl
def create_cube_file(species,x,y,z,STO_matrix,STO_dict,filename='test.cube',cube_grid=0.5,offset=5.0,print_thres=1.0e-20):
    def atomic_contrib(x,y,z,X,Y,Z,species,line,STO_matrix,STO_dict):
        c1s,c2s,c2px,c2py,c2pz,c3s,c3px,c3py,c3pz=[0 for i in range(9)]
        mu1s,mu2s,mu2p,mu3s,mu3p=[0 for i in range(5)]
        if '1s' in STO_dict[species]:
            c1s=STO_matrix[line][0]
            mu1s=STO_dict[species]['1s']
        if '2s' in STO_dict[species]:
            c2s,c2px,c2py,c2pz=[i for i in STO_matrix[line]]
            mu2s=STO_dict[species]['2s']
            mu2p=STO_dict[species]['2p']
        if '3s' in STO_dict[species]:
            c3s,c3px,c3py,c3pz=STO_matrix[line]
            mu3s=STO_dict[species]['3s']
            mu3p=STO_dict[species]['3p']
        R=math.sqrt((x-X)**2+(y-Y)**2+(z-Z)**2)
        res=(0.5641895835477563*c1s*mu1s**1.5)/math.exp(mu1s*R) \
        + (0.32573500793527993*c2s*mu2s**2.5*R)/math.exp(mu2s*R) \
        + (0.11894160774351806*c3s*mu3s**3.5*R*R)/math.exp(mu3s*R) \
        + (0.5641895835477563*mu2p**2.5*(c2px*(x - X) + c2py*(y - Y) + c2pz*(z - Z)))/math.exp(mu2p*R) \
        + (0.20601290774570113*mu3p**3.5*R*(c3px*(x - X) + c3py*(y - Y) + c3pz*(z - Z)))/math.exp(mu3p*R)
        return res
    atoms=len(species)
    xmin,ymin,zmin=[min(i)-offset for i in [x,y,z]]
    xmax,ymax,zmax=[max(i)+offset for i in [x,y,z]]
    resx,resy,resz=[int((xmax-xmin)/cube_grid+1),int((ymax-ymin)/cube_grid+1),int((zmax-zmin)/cube_grid+1)]
    with open(filename,mode='w') as fp:
        print('My CUBE\n',file=fp)
        print(f'{atoms}\t{xmin}\t{ymin}\t{zmin}',file=fp)
        print(f'{resx}\t{cube_grid}\t{0.0}\t{0.0}',file=fp)
        print(f'{resy}\t{0.0}\t{cube_grid}\t{0.0}',file=fp)
        print(f'{resz}\t{0.0}\t{0.0}\t{cube_grid}',file=fp)
        species_dict={'H':1,'C':6,'N':7,'O':8,'F':9,'S':16}
        for j,i in enumerate(species):
            print(f'{species_dict[i]}\t{float(species_dict[i]):>.2}\t{x[j]}\t{y[j]}\t{z[j]}',file=fp)
        for ix in range(resx):
            for iy in range(resy):
                for iz in range(resz):
                    X,Y,Z=[xmin+ix*cube_grid,ymin+iy*cube_grid,zmin+iz*cube_grid]
                    mysum=0
                    for i in range(atoms):
                        mysum+=atomic_contrib(X,Y,Z,x[i],y[i],z[i],species[i],i,STO_matrix,STO_dict)
                    if mysum**2>print_thres:
                        if mysum<0.0:
                            print(f'{mysum**2:>.2} ',end='',file=fp)
                        else:
                            print(f'{-mysum**2:>.2} ',end='',file=fp)
                    else:
                        print('0.0 ',end='',file=fp)
                    if iz % 6 == 5:
                        print(file=fp)
                print(file=fp)

def resolve_STO_matrix(atoms,STOs,STO_id_array,STO_type_array,px,py,pz,AOM_array):
    STO_matrix=[[0,0,0,0] for i in range(atoms)]
    for i in range(STOs):
        if STO_type_array[i]==1:
            STO_matrix[STO_id_array[i]-1][0]=0
        else:
            STO_matrix[STO_id_array[i]-1][0]=0
            STO_matrix[STO_id_array[i]-1][1]=px[STO_id_array[i]-1]*AOM_array[STO_id_array[i]-1]
            STO_matrix[STO_id_array[i]-1][2]=py[STO_id_array[i]-1]*AOM_array[STO_id_array[i]-1]
            STO_matrix[STO_id_array[i]-1][3]=pz[STO_id_array[i]-1]*AOM_array[STO_id_array[i]-1]
    return STO_matrix

def AOM_overlap_calculation(istart,istop,jstart,jstop,x,y,z,STO_id_array,STO_type_array,STO_mu_array,STO_matrix):
    S=0
    for i in range(istart,istop):
        locali=(STO_type_array[i]-2)%4
        for j in range(jstart,jstop):
            localj=(STO_type_array[j]-2)%4
            if locali>0 and localj>0:
                if STO_id_array[i]!=STO_id_array[j]:
                    S=S+STO_matrix[STO_id_array[i]-1][locali]*STO_matrix[STO_id_array[j]-1][localj]\
                        *overlap(x[STO_id_array[i]-1],y[STO_id_array[i]-1],z[STO_id_array[i]-1],
                                 x[STO_id_array[j]-1],y[STO_id_array[j]-1],z[STO_id_array[j]-1],
                                 STO_mu_array[i],STO_mu_array[j],
                                 STO_type_array[i],STO_type_array[j])
                else:
                    if STO_type_array[i]==STO_type_array[j]:
                        S=S+STO_matrix[STO_id_array[i]-1][locali]*STO_matrix[STO_id_array[j]-1][localj]
    return S


