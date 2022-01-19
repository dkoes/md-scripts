# https://github.com/openforcefield/openff-toolkit/blob/stable/examples/swap_amber_parameters/swap_existing_ligand_parameters.ipynb
# Run prepare amber as normal (antechamber will be used), then take the system apart, reparameterize the ligands and put it back together
import parmed
try:
    from openmm.app import PDBFile
except ImportError:
    from simtk.openmm.app import PDBFile
from openff.toolkit.topology import Molecule, Topology
from openff.toolkit.typing.engines.smirnoff import ForceField
try:
    from plumbum.cmd import obabel
except ImportError:
    raise ImportError('Check that obabel is on your path')
import pdb_util as util


def reparm(ligands, base): 
    print('**Running reparameterization of ligand(s) using open force fields\'s SMIRNOFF with openff 2.0.0**')
    # Load already parm'd system
    in_prmtop = base + '.prmtop'
    in_crd = base + '.inpcrd'

    # Create parmed strucuture
    orig_structure = parmed.amber.AmberParm(in_prmtop, in_crd)

    # Split orig_stucuture into unique structure instances e.g. protein, water, ligand, etc.
    pieces = orig_structure.split()
    for piece in pieces:
        # TODO: Figure out how to know which piece is which
        print(f"There are {len(piece[1])} instance(s) of {piece[0]}")

    # Generate an openff topology for the ligand
    # Openff Molecule does not support mol2 so conversion is needed
    ligs_w_sdf= []
    for ligand in ligands:
        obabel[ligand[0], '-O', util.get_base(ligand[0]) + '.sdf']()
        ligs_w_sdf.append((ligand[0], ligand[1], util.get_base(ligand[0]) + '.sdf'))

    # Keep track of ligands that were successfully reparmed so we know to skip them when putting the pieces back together
    reparmed_pieces = []
    complex_structure = parmed.Structure()
    force_field = ForceField("openff_unconstrained-2.0.0.offxml")

    for lig in ligs_w_sdf:
        # Set up openff topology
        ligand_off_molecule = Molecule(lig[2])
        ligand_pdbfile = PDBFile(lig[0])
        ligand_off_topology = Topology.from_openmm(
            ligand_pdbfile.topology,
            unique_molecules=[ligand_off_molecule],
        )
        # Parameterizing the ligand
        # Find ligand "piece", reparm, add to the new structure
        for piece in pieces:
            new_ligand_structure = None
            # TODO: Figure out how to know which piece is which
            if (ligand_off_molecule.n_atoms == len(piece[0].atoms)):
                if (ligand_off_molecule.n_bonds == len(piece[0].bonds)):
                    if ([atom.atomic_number for atom in ligand_off_molecule.atoms] == [atom.element for atom in piece[0].atoms]):
                        print('Found ligand piece', piece)
                        try:
                            # Since the method of matching the piece to ligand is imperfect, ligands that are isomers could mess things up.
                            # So try any piece that matches and see if we get an error
                            print('Reparameterizing ligand using SMIRNOFF')
                            ligand_system = force_field.create_openmm_system(ligand_off_topology)
                            new_ligand_structure = parmed.openmm.load_topology(
                                ligand_off_topology.to_openmm(),
                                ligand_system,
                                xyz=piece[0].positions,
                            )
                            # A quick check to make sure things were not messed up during param
                            if check_discrepencies(new_ligand_structure, piece):
                                # Add the newly parameterized ligand the complex structure
                                reparmed_pieces.append(piece)
                                new_ligand_structure *= len(piece[1])
                                complex_structure += parmed.amber.AmberParm.from_structure(new_ligand_structure)
                            break
                        except:
                            pass

    # Stick all the pieces back together 
    for piece in pieces:
        if (piece not in reparmed_pieces):
            curr_structure = parmed.Structure()
            curr_structure += piece[0]
            curr_structure *= len(piece[1])
            complex_structure += parmed.amber.AmberParm.from_structure(curr_structure)

    # print("Unique atom names:",sorted(list({atom.atom_type.name for atom in complex_structure})),)
    # print("Number of unique atom types:", len({atom.atom_type for atom in complex_structure}))
    # print("Number of unique epsilons:", len({atom.epsilon for atom in complex_structure}))
    # print("Number of unique sigmas:", len({atom.sigma for atom in complex_structure}))

    # # Copy over the original coordinates and box vectors
    complex_structure.coordinates = orig_structure.coordinates
    complex_structure.box_vectors = orig_structure.box_vectors

    # Save the newly parameterized system
    complex_structure.save(base + ".prmtop", overwrite=True)
    complex_structure.save(base + ".inpcrd", overwrite=True)


def check_discrepencies(new_ligand_structure, old_lig_struc):
    # Check how many atoms and which order elements are in the new ligand
    n_atoms_new = len(new_ligand_structure.atoms)
    elements_new = [atom.element for atom in new_ligand_structure.atoms]

    # Check how many atoms and which order elements are in the old ligand
    old_ligand_structure, n_copies = old_lig_struc
    n_atoms_old = len(old_ligand_structure.atoms)
    elements_old = [atom.element for atom in old_ligand_structure.atoms]

    # Print out error message if number of atoms doesn't match
    if n_atoms_new != n_atoms_old:
        print("Error: Number of atoms in input ligand doesn't match number extracted from prmtop file. Using antechamber parameters instead.")
        return False
    elif elements_new != elements_old:
        print("Error: Elements in input ligand don't match elements in the ligand from the prmtop file. Using antechamber parameters instead.")
        print(f"Old elements: {elements_old}")
        print(f"New elements: {elements_new}")
        return False
    else:
        print(f"There are {n_atoms_old} in the old ligand structure and {n_atoms_new} atoms in the new ligand structure")
        return True

