# pyure

Research web app for generating PURE cell-free protein expression CRN models with
BioCRNpyler, exporting SBML, and simulating locally with Bioscrape.

Released under the Apache License 2.0.

This pass focuses on local development and testing. The app is packaged so it can
later be deployed as a container, including Azure Container Apps, without changing
the Python entry point.

## Current scope

- Streamlit web UI with a labeled model hierarchy control.
- DNA sequence input and validation.
- Nucleotide-level detailed model generation using the existing
  `pure_most_detailed.py` implementation.
- Lower-detail modular BioCRNpyler builders following the style in
  `pure_mid_detail.py`.
- SBML model download.
- Bioscrape simulation from generated SBML with editable non-zero initial
  conditions.
- Manual species-name plotting after simulation.

The following hierarchy levels are wired into the UI:

- protein expression only
- TX-TL only model
- lumped resources
- nucleotide-level detail

The highest-detail level uses the current script-style detailed PURE CRN model.
The lower-detail levels use BioCRNpyler mixtures in the modular style shown by
`pure_mid_detail.py`: `ExpressionExtract`, `TxTlExtract`, and `BasicPURE`.

## Repository layout

```text
pyure/
  app/
    streamlit_app.py
  src/
    pyure/
      export.py
      exceptions.py
      model_builders.py
      plotting.py
      simulation.py
      validation.py
  tests/
  .github/workflows/ci.yml
  .streamlit/config.toml
  Dockerfile
  PURE_TXTL_initial_values_Final.csv
  pyproject.toml
  requirements.txt
  pure_most_detailed.py
  fMGG_synthesis_parameters_CRN.csv
```

## Local setup with conda

If using the existing `pyure310` environment:

```powershell
conda activate pyure310
cd C:\Users\apand\Box\Research\PURE_research\pyure
python -m pip install -e .
```

Run the app:

```powershell
streamlit run app/streamlit_app.py
```

Open the URL printed by Streamlit, usually `http://localhost:8501`.

Run tests:

```powershell
pytest
```

## Local setup with pip

Bioscrape can require a working C/C++ build toolchain. Conda or Docker is usually
the smoother path on Windows, but a standard Python environment can be used if
the compiler dependencies are available.

```powershell
cd C:\Users\apand\Box\Research\PURE_research\pyure
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
streamlit run app/streamlit_app.py
```

## Docker

Build the image:

```powershell
docker build -t pyure-streamlit .
```

Run the container:

```powershell
docker run --rm -p 8501:8501 pyure-streamlit
```

Then open `http://localhost:8501`.

## Notes for Azure later

- Streamlit is configured to bind to `0.0.0.0`.
- The Dockerfile exposes port `8501`.
- Runtime settings should be added through environment variables rather than
  hard-coded local paths.
- No cloud deployment resources are included yet.
