# md-scripts

Scripts for preparing, running, and analyzing molecular dynamics simulations.

## `abfe.py`: absolute binding free energy and plain MD

`abfe.py` runs an OpenFE `AbsoluteBindingProtocol` for a protein–ligand complex,
or a `PlainMDProtocol` when `--plain` is supplied. ABFE includes both the complex
and solvent legs of the calculation. Protocol settings use OpenFE defaults
unless a production length is supplied.

### Dependencies and inputs

Run the script in a compatible OpenFE environment containing RDKit, GUFE,
OpenFF Toolkit and Units, OpenMM, and AmberTools. The script explicitly uses
AmberTools for AM1-BCC partial charges. This repository does not pin package versions;
use the [OpenFE installation instructions](https://docs.openfree.energy/en/stable/installation.html)
to set up an environment with compatible dependencies.

Provide:

- A prepared receptor PDB file.
- A ligand SDF containing the intended bound pose, with 3D coordinates in the
  same coordinate frame as the receptor. Only the first SDF record is used.
  Existing hydrogens are retained, and missing hydrogens are added with
  coordinates. A molecule title is optional.

Before setting up a simulation, the script validates the receptor topology
against OpenMM's Amber14 force field. It temporarily adds hydrogens to an
in-memory model for validation but never alters the supplied receptor. If the
topology is invalid, it stops and suggests fixing the structure with PDBFixer.
The script does not dock the ligand or otherwise prepare the receptor.
Protonation states, the binding pose, and suitability for the selected OpenFE
protocol must still be established before running it. Check the
[OpenFE ABFE guide](https://docs.openfree.energy/en/v1.12.0/guide/protocols/absolutebinding.html)
and the restrictions of your installed protocol when choosing inputs and
assessing convergence.

### Usage

Run ABFE with OpenFE's default production lengths:

```bash
python abfe.py --ligand ligand.sdf --receptor receptor.pdb --output abfe_run
```

Set the production length for both the solvent and complex ABFE simulations:

```bash
python abfe.py -l ligand.sdf -r receptor.pdb -o abfe_run --prod_length 10
```

Run plain MD of the solvated protein–ligand complex:

```bash
python abfe.py -l ligand.sdf -r receptor.pdb -o md_run --plain --prod_length 5
```

| Option | Meaning |
| --- | --- |
| `-l`, `--ligand` | Required ligand SDF path. |
| `-r`, `--receptor` | Required receptor PDB path. |
| `-o`, `--output` | Output directory. Defaults to the first word of the ligand's SDF title; for a blank or missing title, uses the ligand filename without its extension. |
| `--plain` | Run plain MD instead of ABFE. |
| `--prod_length` | Production length in nanoseconds. Supply a positive value; omit to use OpenFE defaults. Does not change equilibration lengths or the number of protocol repeats. |
| `-h`, `--help` | Show command-line help. Requires the Python dependencies to be installed. |

The output directory is created if needed; its parent directory must already
exist. Use a separate output directory for each run: an existing `results.json`
is overwritten when a run succeeds.

### Results and failures

On success, `<output>/results.json` contains:

- `estimate`: the binding free energy estimate for ABFE.
- `uncertainty`: the corresponding uncertainty.
- `protocol_result`: the gathered OpenFE result.
- `unit_results`: individual protocol unit results, keyed by their GUFE keys.

The JSON uses GUFE's encoder to preserve quantities and other supported types;
consumers should retain the associated units when interpreting estimates.
For plain MD, `estimate` and `uncertainty` are `null`; use the saved simulation
files and protocol results instead. Shared simulation files are retained in
subdirectories of the output directory.

Each failed protocol unit can be retried three times, for up to four attempts.
If the DAG still fails, the script prints the recorded errors and tracebacks,
exits with a nonzero status, and skips gathering and writing `results.json`.
Any `results.json` left from an earlier run is not removed on failure.
