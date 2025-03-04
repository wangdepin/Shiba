# Shiba - Commands and Style Guidelines

## Build/Test Commands
- Run all tests: `python -m unittest discover test`
- Run a single test file: `python -m unittest test/test_general.py`
- Run a specific test: `python -m unittest test.test_general.TestGeneralModule.test_load_config_success`
- Run Shiba: `python shiba.py -p <threads> config.yaml`
- Run scShiba: `python scshiba.py -p <threads> config.yaml`
- Run SnakeShiba: `snakemake -s snakeshiba.smk --configfile config.yaml --cores <threads> --use-singularity`
- Run SnakeScShiba: `snakemake -s snakescshiba.smk --configfile config.yaml --cores <threads> --use-singularity`

## Code Style
- **Imports**: Standard library first, then third-party, then local modules; sorted alphabetically
- **Docstrings**: Use Google style docstrings with parameter descriptions and return values
- **Functions**: Use snake_case naming convention
- **Classes**: Use PascalCase naming convention
- **Error handling**: Use specific exceptions with descriptive messages
- **Type annotations**: Not currently used, but preferred for new code
- **Formatting**: 4-space indentation (no tabs)
- **Line length**: Aim for 80 characters
- **Comments**: Focus on explaining "why" rather than "what"