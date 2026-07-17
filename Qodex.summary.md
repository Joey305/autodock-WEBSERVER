# Qodex Summary

## Objective

Implement a third packaging mode for non-Miami LSF/HPC users without changing Joey's existing University of Miami workflow.

Date inspected: July 17, 2026.

## Current Behavior Before Changes

### Package mode flow

- `portable`
  - Normalized in [app.py](/Users/jxs794/Documents/autodock-WEBSERVER/app.py).
  - `assemble_job_tree(...)` stages the shared runtime files and excludes the LSF helpers.
  - `build_portable_runners(...)` writes the local shell runners.

- Legacy `lsf`
  - Previously normalized in [app.py](/Users/jxs794/Documents/autodock-WEBSERVER/app.py) as the only LSF-like mode.
  - `assemble_job_tree(...)` staged `1B_confgen_batch.py`, `3B_ServerDocks.py`, `4B_LSFbatch.py`, `lsf_templates.py`, and related HPC assets.
  - The web/UI and v1 API both called `build_confgen_lsfs(...)` and `build_vina_lsfs(...)`.

### Web form submission flow

- The build UI in [templates/build.html](/Users/jxs794/Documents/autodock-WEBSERVER/templates/build.html) offered:
  - `portable`
  - `lsf`
- The UI submitted:
  - `package_mode`
  - `include_lsf`
  - `queue`
  - `project`
  - `walltime`
  - `workers`
  - `mem_per_core`
  - `env_line`
  - `vina_path`
- There was no dedicated LSF email field.
- The environment override field was a single-line `<input>`.

### Headless/API flow

- `/api/build` and `/api/v1/workspaces/<jobname>/build` both normalized mode in [app.py](/Users/jxs794/Documents/autodock-WEBSERVER/app.py).
- Both routes built LSF files only when the normalized mode equaled `lsf`.
- Both routes always set `email = _public_email()` for generated LSF files instead of reading a custom LSF notification email.

### ZIP assembly flow

- [packager.py](/Users/jxs794/Documents/autodock-WEBSERVER/packager.py):
  - `assemble_job_tree(...)` creates `job/`
  - copies receptors, ligands, and canonicalized `vina_centers.csv`
  - stages runtime assets from `RUNTIME_ROOT_FILES`
  - stages `RUNTIME_LSF_FILES` only when `package_mode == "lsf"`
  - zips the finished tree with `zip_job_tree(...)`

### Generated `.lsf` files

- [lsf_templates.py](/Users/jxs794/Documents/autodock-WEBSERVER/lsf_templates.py) previously hard-coded:
  - `DEFAULT_EMAIL = "jxs794@miami.edu"`
  - `DEFAULT_QUEUE = "general"`
  - `DEFAULT_PROJECT = "brd"`
  - Joey's Conda activation path
  - Joey's `vina_env`
  - Joey's default `VINA_EXE`
  - `#BSUB -P`, `#BSUB -B`, `#BSUB -N`, and `#BSUB -u` in every header
  - mixed Python selection behavior across scripts

### Included helper scripts that could regenerate `.lsf`

- [1B_confgen_batch.py](/Users/jxs794/Documents/autodock-WEBSERVER/1B_confgen_batch.py)
  - hard-coded Joey defaults in CLI args and template
  - default env activation: `/nethome/jxs794/.../conda.sh && conda activate vina_env`

- [3B_ServerDocks.py](/Users/jxs794/Documents/autodock-WEBSERVER/3B_ServerDocks.py)
  - hard-coded Joey queue/project/email defaults
  - hard-coded Joey Conda path and `VINA_EXE`
  - generated docking `.lsf` files with `python 3_Complete_batch_docking.py`

- [4B_LSFbatch.py](/Users/jxs794/Documents/autodock-WEBSERVER/4B_LSFbatch.py)
  - hard-coded Joey prompt defaults for queue/project/email
  - hard-coded Joey Conda activation

## Joey-Specific Values Found

### Direct literal values

- `jxs794@miami.edu`
  - [lsf_templates.py](/Users/jxs794/Documents/autodock-WEBSERVER/lsf_templates.py)
  - [1B_confgen_batch.py](/Users/jxs794/Documents/autodock-WEBSERVER/1B_confgen_batch.py)
  - [3B_ServerDocks.py](/Users/jxs794/Documents/autodock-WEBSERVER/3B_ServerDocks.py)
  - [4B_LSFbatch.py](/Users/jxs794/Documents/autodock-WEBSERVER/4B_LSFbatch.py)

- `/nethome/jxs794/miniconda3/etc/profile.d/conda.sh`
  - [lsf_templates.py](/Users/jxs794/Documents/autodock-WEBSERVER/lsf_templates.py)
  - [1B_confgen_batch.py](/Users/jxs794/Documents/autodock-WEBSERVER/1B_confgen_batch.py)
  - [3B_ServerDocks.py](/Users/jxs794/Documents/autodock-WEBSERVER/3B_ServerDocks.py)
  - [4B_LSFbatch.py](/Users/jxs794/Documents/autodock-WEBSERVER/4B_LSFbatch.py)

- `vina_env`
  - [lsf_templates.py](/Users/jxs794/Documents/autodock-WEBSERVER/lsf_templates.py)
  - [1B_confgen_batch.py](/Users/jxs794/Documents/autodock-WEBSERVER/1B_confgen_batch.py)
  - [3B_ServerDocks.py](/Users/jxs794/Documents/autodock-WEBSERVER/3B_ServerDocks.py)
  - [4B_LSFbatch.py](/Users/jxs794/Documents/autodock-WEBSERVER/4B_LSFbatch.py)

- `$HOME/miniconda3/envs/vina_env/bin/vina`
  - [lsf_templates.py](/Users/jxs794/Documents/autodock-WEBSERVER/lsf_templates.py)
  - [3B_ServerDocks.py](/Users/jxs794/Documents/autodock-WEBSERVER/3B_ServerDocks.py)

- `general`
  - [lsf_templates.py](/Users/jxs794/Documents/autodock-WEBSERVER/lsf_templates.py)
  - [1B_confgen_batch.py](/Users/jxs794/Documents/autodock-WEBSERVER/1B_confgen_batch.py)
  - [3B_ServerDocks.py](/Users/jxs794/Documents/autodock-WEBSERVER/3B_ServerDocks.py)
  - [4B_LSFbatch.py](/Users/jxs794/Documents/autodock-WEBSERVER/4B_LSFbatch.py)
  - [templates/build.html](/Users/jxs794/Documents/autodock-WEBSERVER/templates/build.html)
  - [app.py](/Users/jxs794/Documents/autodock-WEBSERVER/app.py)

- `brd`
  - [lsf_templates.py](/Users/jxs794/Documents/autodock-WEBSERVER/lsf_templates.py)
  - [1B_confgen_batch.py](/Users/jxs794/Documents/autodock-WEBSERVER/1B_confgen_batch.py)
  - [3B_ServerDocks.py](/Users/jxs794/Documents/autodock-WEBSERVER/3B_ServerDocks.py)
  - [4B_LSFbatch.py](/Users/jxs794/Documents/autodock-WEBSERVER/4B_LSFbatch.py)
  - [templates/build.html](/Users/jxs794/Documents/autodock-WEBSERVER/templates/build.html)
  - [app.py](/Users/jxs794/Documents/autodock-WEBSERVER/app.py)

### Derived or structural Joey-specific behavior

- Notification directives were always emitted whenever LSF files were generated.
- `#BSUB -P` was always emitted, even for clusters that do not use an account/project directive.
- Python invocation was inconsistent:
  - shared template chose `python3 || python`
  - `3B_ServerDocks.py` and `4B_LSFbatch.py` used `python`
- `3B_ServerDocks.py`, `1B_confgen_batch.py`, and `4B_LSFbatch.py` could regenerate LSF files using Joey defaults even after a package had been customized earlier.

## Changes Implemented

### Centralized mode and profile model

- Added [hpc_profiles.py](/Users/jxs794/Documents/autodock-WEBSERVER/hpc_profiles.py)
  - canonical modes: `portable`, `joey_lsf`, `custom_lsf`
  - legacy alias: `lsf -> joey_lsf`
  - immutable `JOEY_LSF_PROFILE`
  - custom profile builder and packaged profile persistence
  - shared LSF header and setup rendering

### Route and API changes

- [app.py](/Users/jxs794/Documents/autodock-WEBSERVER/app.py)
  - now normalizes package modes through the shared profile module
  - preserves backward compatibility for incoming `lsf`
  - builds Joey packages from the built-in Joey profile
  - builds custom packages from user-provided LSF/HPC fields

### Package assembly changes

- [packager.py](/Users/jxs794/Documents/autodock-WEBSERVER/packager.py)
  - stages LSF assets for `joey_lsf` and `custom_lsf`
  - includes `hpc_profiles.py` in LSF packages

### LSF generation changes

- [lsf_templates.py](/Users/jxs794/Documents/autodock-WEBSERVER/lsf_templates.py)
  - uses the shared profile model
  - writes `hpc_profile.json` into generated LSF packages
  - emits `#BSUB -P` only when a project/account is configured
  - emits email directives only when enabled
  - uses one shared Python command field consistently

### Packaged helper script changes

- [1B_confgen_batch.py](/Users/jxs794/Documents/autodock-WEBSERVER/1B_confgen_batch.py)
- [3B_ServerDocks.py](/Users/jxs794/Documents/autodock-WEBSERVER/3B_ServerDocks.py)
- [4B_LSFbatch.py](/Users/jxs794/Documents/autodock-WEBSERVER/4B_LSFbatch.py)

Each of these now loads `hpc_profile.json` from the generated package when present, so later `.lsf` regeneration stays aligned with the selected package profile instead of silently reverting to Joey defaults.

### UI and docs

- [templates/build.html](/Users/jxs794/Documents/autodock-WEBSERVER/templates/build.html)
  - explicit `joey_lsf` mode
  - new `custom_lsf` mode
  - new custom fields for email, notifications, queue, project, workers, memory, walltimes, Conda, Python, Vina, and multiline setup commands

- Updated:
  - [README.md](/Users/jxs794/Documents/autodock-WEBSERVER/README.md)
  - [docs/HEADLESS_API.md](/Users/jxs794/Documents/autodock-WEBSERVER/docs/HEADLESS_API.md)

## Validation

- `python3 -m py_compile app.py packager.py lsf_templates.py hpc_profiles.py 1B_confgen_batch.py 3B_ServerDocks.py 4B_LSFbatch.py`
- Focused packaging/profile unit tests were added and run after implementation.
