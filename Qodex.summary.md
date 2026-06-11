# Qodex.summary

## Task
Add a headless programmatic API for AutoDock-Vina PrepServer.

## Original Goal
The user wants a background/scriptable/API system so users, command-line workflows, Python scripts, and future AI agents can run the docking-preparation workflow without using the visual website. A key requirement is programmatic docking-box generation from residue/chain/residue ID or HETATM ligand instances such as DR7, chain A, residue 100.

## Assumptions
- The current dirty working tree is intentional user work and must be preserved.
- The active browser workflow is `templates/build.html` and the legacy browser endpoints in `app.py` are the compatibility surface.
- New endpoints should be additive under `/api/v1` and should not replace existing `/api/...` routes.
- Headless center resolution can safely support PDB, ENT, and PDBQT fixed-column coordinate records in this pass.
- CIF/mmCIF uploads remain accepted by the app, but selector-based headless center resolution for CIF/mmCIF should be documented as a limitation rather than implemented without a parser dependency.
- Full one-call package orchestration is riskier than a staged workflow because prep/build failures should remain visible, so `/api/v1/headless/package` is reserved and returns a documented staged-workflow response.

## Files Inspected
- `app.py` - inspected Flask pages, state helpers, existing `/api/...` endpoints, workspace handling, receptor upload/fetch, center saving, prep, ligand upload, build, and download behavior.
- `templates/build.html` - inspected the active browser workflow and JavaScript API usage.
- `templates/documentation.html` - inspected public documentation hub before adding the API mention.
- `README.md` - inspected existing setup, workflow, package-mode, and security notes.
- `packager.py` - inspected workspace, upload, fetch, center CSV, package assembly, and ZIP helpers.
- `tests/test_public_access.py` and `tests/test_packaging_modes.py` - inspected existing Flask and packaging test patterns.
- `tests/` - inspected available test suite layout.
- `docs/` - checked and confirmed it did not exist before this task.

## Files Changed
- `app.py` - added `center_resolver` imports and a versioned `/api/v1` JSON API layer for health, workspaces, receptors, centers, prep, ligands, build, artifacts, download, and a reserved headless package endpoint.
- `README.md` - added a concise Headless API section and feature bullet.
- `templates/documentation.html` - added a restrained documentation-card mention of the headless API.
- `Qodex.summary.md` - updated with this task summary.

## Files Created
- `center_resolver.py` - reusable center-resolution utility with PDB/ENT/PDBQT parsing, centroid calculation, ambiguity handling, and structured errors.
- `docs/HEADLESS_API.md` - route table, workflow, center modes, error shape, limitations, and security notes.
- `docs/API_EXAMPLES.md` - curl and Python examples.
- `docs/AGENT_WORKFLOW.md` - staged CLI/agent workflow, polling strategy, validation, and failure handling.
- `docs/examples/headless_center_by_residue.py`
- `docs/examples/headless_center_by_hetatm.py`
- `docs/examples/headless_build_from_xyz.py`
- `docs/examples/headless_create_package.py`
- `tests/test_center_resolver.py`
- `tests/test_headless_api_contract.py`

## Implementation Summary
The app now exposes an additive `/api/v1` API for creating or reusing workspaces, uploading or fetching receptors, resolving and saving docking centers, preparing receptors, uploading ligands, building packages, listing artifacts, and downloading generated ZIPs. Center generation can be done headlessly from explicit XYZ coordinates, residue selectors, HETATM/ligand instance selectors, atom-level selectors, or a simple selection object. Saved centers write to the same canonical center CSV used by the browser workflow and package builder.

## Key Decisions
- The `/api/v1` routes use a consistent JSON shape with `ok`, `data`, and `warnings` on success and `ok`, `error`, `message`, and `details` on failure.
- Center resolution lives in `center_resolver.py` so it can be tested without Flask.
- HETATM aliases `het`, `hetatm`, `ligand`, and `resname` are accepted.
- Ligand-instance ambiguity is explicit; the resolver returns `ambiguous_selection` with candidates instead of choosing the first match.
- The browser workflow was preserved by leaving existing `/api/...` endpoints and `templates/build.html` JavaScript unchanged.
- The staged `/api/v1` workflow is implemented; full one-call orchestration is deferred and documented through the reserved `/api/v1/headless/package` response.
- Agent/CLI usage is documented as a staged workflow with polling and validation rather than as a complete autonomous platform.

## Commands Run
- `git status --short` - confirmed the repository already had uncommitted changes before this task.
- `sed`, `find`, and `rg` inspections - reviewed routes, templates, API usage, docs, tests, and helper code.
- `python -m flask --app app routes` - attempted before and after implementation; blocked by missing `flask_login` in the ambient Python environment.
- `python -m py_compile app.py` - passed.
- `python -m py_compile center_resolver.py` - passed.
- `python -m unittest tests.test_center_resolver -v` - passed 7 resolver tests.
- `python -m unittest tests.test_headless_api_contract -v` - skipped because Flask app dependencies are not installed.
- `python -m unittest discover -s tests -v` - ran the broader suite; new resolver tests passed, API contract tests skipped, and unrelated existing app/PyMOL-builder tests failed.
- Source route searches with `rg` - confirmed existing browser endpoints and new `/api/v1` route declarations are present.
- Unsafe claim search - no new unsafe overclaims were introduced; existing safe “not a claim” text matched in informational templates.
- `git diff --check` - passed.

## Validation Results
- Syntax checks passed for `app.py` and `center_resolver.py`.
- New center resolver unit tests passed.
- New Flask API contract tests are present but skipped locally because `flask_login` is missing.
- Flask route listing could not run locally because `flask_login` is missing.
- Full test discovery had unrelated failures:
  - app-import tests fail due to missing `flask_login`;
  - existing PyMOL-builder tests fail due to pre-existing function/CLI expectation mismatches.
- Source checks confirmed the old browser endpoints remain declared in `app.py`.
- Source checks confirmed all new `/api/v1` route declarations are present.
- `git diff --check` passed.

## Known Issues
- Headless selector-based center resolution supports PDB, ENT, and PDBQT-like fixed-column records, not CIF/mmCIF.
- Receptor preparation still relies on local `obabel` availability.
- The API uses the app's existing lightweight synchronous prep behavior and does not add a production queue.
- `/api/v1/headless/package` is reserved and currently returns a staged-workflow response.
- Authentication, rate limiting, quotas, cleanup, and abuse controls were intentionally not implemented.
- Local dynamic Flask validation could not run until dependencies such as `flask_login` are installed in the active environment.

## Manual Verification
1. Open `/build` and confirm the browser workflow still works.
2. Call `/api/v1/health`.
3. Create a workspace through `/api/v1/workspaces`.
4. Fetch or upload a receptor through `/api/v1`.
5. Resolve a center by residue.
6. Resolve a center by HETATM ligand instance.
7. Save a center.
8. Start prep and poll status.
9. Upload ligands.
10. Build a portable package.
11. Confirm generated download URL works.
12. Review `docs/HEADLESS_API.md`, `docs/API_EXAMPLES.md`, and `docs/AGENT_WORKFLOW.md`.

## Suggested Next Prompt
Add an OpenAPI schema and a tiny Python client package for the `/api/v1` staged workflow, including fixture-backed API smoke tests.
