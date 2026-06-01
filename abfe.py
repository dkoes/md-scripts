#!/usr/bin/env python3

import argparse
import os
import sys
from rdkit.Chem import AllChem as Chem
import openfe
from openfe.protocols.openmm_utils.omm_settings import OpenFFPartialChargeSettings
from openfe.protocols.openmm_utils.charge_generation import assign_offmol_partial_charges
from openfe.protocols.openmm_afe import AbsoluteBindingProtocol
from openfe.protocols.openmm_md import PlainMDProtocol
from openff.units import unit

from gufe.protocols import execute_DAG
import pathlib
import gzip
import json
import gufe

def parse_args():
    parser = argparse.ArgumentParser(
        description="OpenFE Absolute Binding Free Energy"
    )

    parser.add_argument(
        "-l", "--ligand",
        required=True,
        help="Ligand file (e.g., ligand.sdf)"
    )

    parser.add_argument(
        "-r", "--receptor",
        required=True,
        help="Receptor file (e.g., receptor.pdb)"
    )

    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Output prefix (default: output)"
    )
    
    parser.add_argument("--plain", action='store_true',help='Execute a plain MD simulation')
    parser.add_argument("--prod_length", type=float,help='Length of production simulation (ns)')

    return parser.parse_args()


def check_file(path, label):
    if not os.path.exists(path):
        sys.exit(f"ERROR: {label} file not found: {path}")
    if not os.path.isfile(path):
        sys.exit(f"ERROR: {label} is not a file: {path}")


if __name__ == "__main__":
    args = parse_args()

    check_file(args.ligand, "Ligand")
    check_file(args.receptor, "Receptor")

    print(f"Reading {args.ligand}")
    suppl = Chem.SDMolSupplier(args.ligand, removeHs=False)
    mol = next(suppl)
    mol = Chem.AddHs(mol,addCoords=True)
    name = mol.GetProp('_Name').split()[0]
    ligand = openfe.SmallMoleculeComponent.from_rdkit(mol)

    print("Calculating partial charges")
    charge_settings = OpenFFPartialChargeSettings(partial_charge_method="am1bcc", off_toolkit_backend="ambertools")
    
    mol_with_charge = assign_offmol_partial_charges(
        offmol=ligand.to_openff(),
        overwrite=False,
        method=charge_settings.partial_charge_method,
        toolkit_backend=charge_settings.off_toolkit_backend,
        generate_n_conformers=charge_settings.number_of_conformers,
        nagl_model=charge_settings.nagl_model
    )

    ligand = openfe.SmallMoleculeComponent.from_openff(mol_with_charge)

    print("Setting up systems")
    solvent = openfe.SolventComponent()
    protein = openfe.ProteinComponent.from_pdb_file(args.receptor)

    # In state A the ligand is fully interacting in the solvent
    systemA = openfe.ChemicalSystem({
        'ligand': ligand,
        'protein': protein,
        'solvent': solvent,
    }, name=ligand.name)
    # In state B the ligand is fully decoupled in the solvent, therefore we are only defining the solvent here
    systemB = openfe.ChemicalSystem({
        'protein': protein,
        'solvent': solvent,
    })


    if args.plain:
        settings = PlainMDProtocol.default_settings()
        if args.prod_length:
            settings.simulation_settings.production_length = args.prod_length*unit.nanosecond
        protocol = PlainMDProtocol(settings=settings)
        dag = protocol.create(stateA=systemA, stateB=systemA, mapping=None)

    else: #free energy
        settings = AbsoluteBindingProtocol.default_settings()
        if args.prod_length:
            asettings.solvent_simulation_settings.production_length = args.prod_length*unit.nanosecond
            asettings.complex_simulation_settings.production_length = args.prod_length*unit.nanosecond
        protocol = AbsoluteBindingProtocol(settings=settings)
        dag = protocol.create(stateA=systemA, stateB=systemB, mapping=None)


    outdir = args.output
    if outdir is None:
        outdir = name

    path = pathlib.Path(outdir)
    path.mkdir(exist_ok=True)
    
    print("Running...")
    # Execute the DAG
    dag_results = execute_DAG(dag, scratch_basedir=path, shared_basedir=path, n_retries=3,keep_shared=True)

    protocol_results = protocol.gather([dag_results])
    print(f"ABFE dG: {protocol_results.get_estimate()}, err {protocol_results.get_uncertainty()}")

    # Save the results in a json file

    outdict = {
         "estimate": protocol_results.get_estimate(),
         "uncertainty": protocol_results.get_uncertainty(),
         "protocol_result": protocol_results.to_dict(),
         "unit_results": {
             unit.key: unit.to_keyed_dict()
             for unit in dag_results.protocol_unit_results
         }
     }
    
    with open(f"{outdir}/results.json",'wt') as stream:
        print(f"Saving {outdir}/results.json")
        json.dump(outdict, stream, cls=gufe.tokenization.JSON_HANDLER.encoder)
