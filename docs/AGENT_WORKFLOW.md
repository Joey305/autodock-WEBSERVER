# Agent Workflow

This guide describes how a command-line automation or future AI agent should use the headless API without relying on browser clicks.

## Recommended Staged Flow

1. Call `/api/v1/health`.
2. Create or reuse a workspace with `/api/v1/workspaces`.
3. Fetch or upload receptors.
4. List receptors and verify the expected file names.
5. Resolve centers with `/centers/resolve`.
6. Save centers with `/centers/save`.
7. Confirm `/centers` has one row for each expected receptor PDBQT name.
8. Start receptor prep.
9. Poll `/prep/status` or `/summary`.
10. Upload ligands.
11. Build a package.
12. List artifacts and download the ZIP.

## Validation Before Build

Before calling `/build`, verify:

- `receptors_total` is greater than zero.
- `have_rows_for_all` is true.
- `ligands_uploaded` is true.
- receptor prep has generated the expected PDBQT outputs if your workflow requires prepared receptors.

## Polling Strategy

Use moderate polling, for example every 2 to 5 seconds for local development. Stop polling after a deployment-appropriate timeout and surface the last log tail from `/prep/status`.

## Failure Handling

Treat `ok=false` as a machine-readable failure. Use `error` for branching and `message` for human display.

Important errors:

- `workspace_missing`
- `receptor_not_found`
- `unsupported_structure_format`
- `invalid_xyz`
- `no_atoms_matched`
- `ambiguous_selection`
- `centers_incomplete`
- `ligands_missing`

If a ligand selector returns `ambiguous_selection`, inspect `details.candidates` and retry with `chain` and `resi`.

## Safe Boundaries

The API can write files, run receptor conversion, and create downloadable packages. Agents should avoid uploading confidential data to public deployments, should avoid repeated large uploads, and should not assume authentication, rate limiting, quotas, or cleanup are present unless a deployment explicitly adds them.

## One-Call Package Endpoint

`/api/v1/headless/package` is reserved but intentionally returns a staged-workflow response in this release. A future one-call endpoint should still expose intermediate validation and artifact status rather than hiding failures behind a single opaque job.
