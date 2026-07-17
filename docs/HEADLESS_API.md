# Headless API

AutoDock-Vina PrepServer includes a versioned JSON API for scripted workspace creation, receptor intake, docking-box generation, receptor preparation, ligand upload, and package generation.

The API is additive. Existing browser routes under `/api/...` and the `/build` workflow remain unchanged.

## Route Table

| Method | Route | Purpose |
| --- | --- | --- |
| GET | `/api/v1/health` | Check service and API version. |
| POST | `/api/v1/workspaces` | Create or reuse a workspace. |
| GET | `/api/v1/workspaces/<jobname>` | Read workspace state. |
| GET | `/api/v1/workspaces/<jobname>/summary` | Read workflow status and artifact summary. |
| POST | `/api/v1/workspaces/<jobname>/receptors/upload` | Upload receptor file, ZIP, or folder-style multipart files. |
| POST | `/api/v1/workspaces/<jobname>/receptors/fetch` | Fetch a receptor from RCSB by PDB ID. |
| GET | `/api/v1/workspaces/<jobname>/receptors` | List receptors. |
| POST | `/api/v1/workspaces/<jobname>/centers/resolve` | Compute a center without saving it. |
| POST | `/api/v1/workspaces/<jobname>/centers/save` | Compute and save a center to `vina_centers.csv`. |
| GET | `/api/v1/workspaces/<jobname>/centers` | List saved centers. |
| POST | `/api/v1/workspaces/<jobname>/prep/start` | Run receptor preparation using the existing conversion path. |
| GET | `/api/v1/workspaces/<jobname>/prep/status` | Read receptor-prep status and log tail. |
| POST | `/api/v1/workspaces/<jobname>/ligands/upload` | Upload ligand file, ZIP, or folder-style multipart files. |
| GET | `/api/v1/workspaces/<jobname>/ligands` | List ligand files and upload metadata. |
| POST | `/api/v1/workspaces/<jobname>/build` | Build a portable or optional LSF package. |
| GET | `/api/v1/workspaces/<jobname>/artifacts` | List generated ZIP artifacts. |
| GET | `/api/v1/workspaces/<jobname>/download` | Download the latest or selected workspace artifact. |
| POST | `/api/v1/headless/package` | Reserved. Returns a staged-workflow response in this release. |

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

The browser accepts `.pdb`, `.cif`, `.mmcif`, `.ent`, and `.pdbqt`. Headless center resolution currently supports PDB, ENT, and PDBQT-like fixed-column coordinate records. CIF/mmCIF upload remains accepted by the workflow, but selector-based headless center resolution for CIF/mmCIF is documented as a current limitation.

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
- CIF/mmCIF selector-based center resolution is not implemented in the headless resolver.
- `/api/v1/headless/package` is reserved and returns a staged-workflow response.
- Authentication, rate limiting, and quotas are not implemented by this task.
