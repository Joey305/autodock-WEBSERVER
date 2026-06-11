# Contributing

## Development workflow

1. Fork the repository.
2. Create a feature branch.
3. Create and activate a virtual environment.
4. Install dependencies with `pip install -r requirements.txt`.
5. Make focused, reviewable changes.
6. Run the validation commands before opening a pull request:

```bash
python3 -m unittest discover -s tests -v
python3 -m py_compile app.py packager.py runner_templates.py manage.py
```

## Pull request guidance

- Keep changes scoped to one improvement or bug fix where practical.
- Update `README.md` when setup, behavior, or workflow details change.
- Do not commit `.env`, databases, logs, caches, or generated docking outputs.
- If your change affects scientific workflow outputs, describe the expected effect clearly in the pull request.

## Reporting issues

When opening an issue, include:

- what you were trying to do
- the command or UI step that failed
- the error message
- whether `obabel`, `vina`, and RDKit are installed in your environment
