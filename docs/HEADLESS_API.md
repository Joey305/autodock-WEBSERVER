# Headless API

AutoDock-Vina PrepServer includes a versioned JSON API for scripted workspace creation, receptor intake, docking-box generation, receptor preparation, ligand upload, and package generation.

The API is additive. Existing browser routes under `/api/...` and the `/build` workflow remain unchanged.

## Route Table

| Method | Route | Purpose |
| --- | --- | --- |
| GET | `/api/v1/health` | Check service and API version. |
| GET | `/api/v1/clients/headless_redock_bound_ligand.py` | Download the command-line redocking client script. |
| POST | `/api/v1/workspaces` | Create or reuse a workspace. |
| GET | `/api/v1/workspaces/<jobname>` | Read workspace state. |
| GET | `/api/v1/workspaces/<jobname>/summary` | Read workflow status and artifact summary. |
| POST | `/api/v1/workspaces/<jobname>/receptors/upload` | Upload receptor file, ZIP, or folder-style multipart files. |
| POST | `/api/v1/workspaces/<jobname>/receptors/fetch` | Fetch a receptor from RCSB by PDB ID. |
| GET | `/api/v1/workspaces/<jobname>/receptors` | List receptors. |
| GET | `/api/v1/workspaces/<jobname>/hetatms` | List grouped HETATM ligand/heterogen instances and centroids. |
| POST | `/api/v1/workspaces/<jobname>/centers/resolve` | Compute a center without saving it. |
| POST | `/api/v1/workspaces/<jobname>/centers/save` | Compute and save a center to `vina_centers.csv`. |
| GET | `/api/v1/workspaces/<jobname>/centers` | List saved centers. |
| POST | `/api/v1/workspaces/<jobname>/prep/start` | Run receptor preparation using the existing conversion path. |
| GET | `/api/v1/workspaces/<jobname>/prep/status` | Read receptor-prep status and log tail. |
| POST | `/api/v1/workspaces/<jobname>/ligands/upload` | Upload ligand file, ZIP, or folder-style multipart files. |
| POST | `/api/v1/workspaces/<jobname>/ligands/extract` | Extract a selected receptor HETATM instance into `Ligands/` as SDF. |
| GET | `/api/v1/workspaces/<jobname>/ligands` | List ligand files and upload metadata. |
| POST | `/api/v1/workspaces/<jobname>/build` | Build a portable or optional LSF package. |
| GET | `/api/v1/workspaces/<jobname>/artifacts` | List generated ZIP artifacts. |
| GET | `/api/v1/workspaces/<jobname>/download` | Download the latest or selected workspace artifact. |
| POST | `/api/v1/headless/package` | One-call headless receptor, center, bound-ligand extraction, prep, and package workflow. |

## Response Shape

Success:

```json
{"ok": true, "data": {}, "warnings": []}
```

Error:

```json
{"ok": false, "error": "machine_readable_error", "message": "Human-readable explanation.", "details": {}}
```

## Workspace Lifecycle

Create a workspace:

```bash
curl -s -X POST "$BASE/api/v1/workspaces" \
  -H "Content-Type: application/json" \
  -d '{"workspace_name":"api-test","reuse":true}'
```

If `reuse` is true and the workspace exists, the existing workspace is returned. Otherwise a new workspace is created.

## Receptor Input

Fetch from RCSB:

```json
{"pdb_id": "3EKY", "chains": ["A"]}
```

Upload receptors with multipart form data:

- `mode=single`, `file=@receptor.pdb`
- `mode=zip`, `file=@Receptors.zip`
- `mode=folder`, repeated `files=@...`

The browser accepts `.pdb`, `.cif`, `.mmcif`, `.ent`, and `.pdbqt`. Headless center resolution supports PDB, ENT, PDBQT-like fixed-column coordinate records, and common mmCIF atom-site loops.

## HETATM Discovery

List candidate bound ligands and heterogens after receptor intake:

```bash
curl -s "$BASE/api/v1/workspaces/api-example/hetatms?receptor=9g94.cif" \
  | python -m json.tool
```

Water is hidden by default. Add `include_water=true` when water positions are relevant.

## Center Selection Modes

Explicit XYZ:

```json
{"method": "xyz", "receptor": "3eky.pdb", "center": [12.3, 45.6, 7.89], "size": 20}
```

Residue centroid:

```json
{"method": "residue", "receptor": "3eky.pdb", "chain": "A", "resi": "82", "resname": "GLU", "size": 20}
```

HETATM or ligand-instance centroid:

```json
{"method": "hetatm", "receptor": "3eky.pdb", "het": "DR7", "chain": "A", "resi": "100", "size": 20}
```

Accepted HETATM aliases are `het`, `hetatm`, `ligand`, and `resname`.

Atom-specific center:

```json
{"method": "atom", "receptor": "3eky.pdb", "record": "HETATM", "resname": "DR7", "chain": "A", "resi": "100", "atom_name": "C12"}
```

Selection object:

```json
{"method": "selection", "receptor": "3eky.pdb", "selection": {"record": "HETATM", "resname": "DR7", "chain": "A", "resi": "100"}}
```

`centers/resolve` computes without persistence. `centers/save` computes and writes the canonical `PDB_ID,X,Y,Z,SIZE` row used by the browser workflow and package builder.

## Ambiguity Handling

The resolver does not silently choose the first ligand instance. If a HETATM selector matches multiple ligand instances, the API returns `ambiguous_selection` with candidates. Provide `chain` and `resi` to disambiguate.

Zero-match selectors return `no_atoms_matched`.

## Bound Ligand Extraction

After uploading or fetching a receptor, a bound HETATM ligand instance can be reused as a ligand input:

```json
{"receptor": "3eky.pdb", "resname": "DR7", "chain": "A", "resi": "100"}
```

`ligands/extract` writes an SDF into `Ligands/`, marks ligand input as uploaded, and returns `ambiguous_selection` if the selector matches multiple HETATM instances. Open Babel must be available on the server for SDF conversion.

## One-Call Redocking Package

For a local deployment or the hosted site:

```bash
BASE="${BASE_URL:-https://autodockvina.com}"
```

Create a portable redocking package for RCSB entry `9G94`, using bound ligand `A1D73` from chain `A`, residue `101` as both the ligand input and docking-box center:

```bash
curl -s -X POST "$BASE/api/v1/headless/package" \
  -H "Content-Type: application/json" \
  -d '{
    "workspace_name": "9g94-redock",
    "reuse": true,
    "receptor": {"pdb_id": "9G94"},
    "bound_ligand": {"resname": "A1D73", "chain": "A", "resi": "101"},
    "center": {"method": "same_as_bound_ligand", "size": 20},
    "prep": {"remove_hets": "all", "remove_chains": ["B", "C"], "altloc": "collapse"},
    "package": {"package_mode": "portable", "poses_conf": 64, "poses_vina": 9}
  }' | python -m json.tool
```

The response includes `artifact.download_url`. Download it with that returned path:

```bash
curl -L -o 9g94-redock.zip "$BASE<artifact.download_url>"
```

For an interactive command-line wrapper:

```bash
curl -fsSLo autodock-redock-bound-ligand.py \
  "$BASE/api/v1/clients/headless_redock_bound_ligand.py"

python autodock-redock-bound-ligand.py \
  --base-url "$BASE" \
  --pdb-id 9G94 \
  --ligand A1D73 \
  --download-dir "$HOME/Docking"
```

In bash or zsh, you can run it without saving a local copy:

```bash
python <(curl -fsSL "$BASE/api/v1/clients/headless_redock_bound_ligand.py") \
  --base-url "$BASE" \
  --pdb-id 9G94 \
  --ligand A1D73 \
  --download-dir "$HOME/Docking"
```

## Prep And Build

After centers are saved for every receptor:

```bash
curl -s -X POST "$BASE/api/v1/workspaces/$JOB/prep/start" \
  -H "Content-Type: application/json" \
  -d '{"remove_hets":["HOH"],"remove_chains":[],"altloc":"collapse"}'
```

Poll:

```bash
curl -s "$BASE/api/v1/workspaces/$JOB/prep/status"
```

Build a portable package:

```json
{"package_mode": "portable", "poses_conf": 64, "poses_vina": 9}
```

Build Joey's legacy-compatible Miami LSF package:

```json
{"package_mode": "joey_lsf", "poses_conf": 64, "poses_vina": 9}
```

Build a custom LSF package for another cluster:

```json
{
  "package_mode": "custom_lsf",
  "queue": "general",
  "project": "",
  "workers": 16,
  "mem_per_core": 2000,
  "confgen_walltime": "48:00",
  "vina_walltime": "96:00",
  "lsf_email": "cluster-user@example.org",
  "python_command": "$(command -v python3 || command -v python)",
  "conda_sh": "/shared/miniconda3/etc/profile.d/conda.sh",
  "conda_env": "vina_env",
  "vina_path": "/shared/miniconda3/envs/vina_env/bin/vina",
  "setup_commands": "module load miniconda3"
}
```

## Security Note

This API follows the app's existing public/no-login behavior when public mode is enabled. Do not expose write/build endpoints publicly without deployment controls if abuse is a concern. Do not upload confidential structures or ligand libraries to public instances unless the deployment is configured for that use. Future public deployments should consider authentication, rate limiting, quotas, cleanup, and monitoring.

## Limitations

- No production queue system is introduced; prep currently follows the app's lightweight synchronous conversion behavior.
- One-call JSON packaging supports fetched or existing workspace receptors. Local ligand files still use the staged multipart upload endpoint or the example CLI wrapper.
- Authentication, rate limiting, and quotas are not implemented by this task.
