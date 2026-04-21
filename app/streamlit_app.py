from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import streamlit as st

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pyure.exceptions import (  # noqa: E402
    ModelGenerationError,
    ModelNotImplementedError,
    PyureError,
    SequenceValidationError,
    SimulationError,
)
from pyure.export import download_filename  # noqa: E402
from pyure.initial_conditions import (  # noqa: E402
    dataframe_to_initial_conditions,
    initial_conditions_to_dataframe,
)
from pyure.model_builders import HIERARCHY_OPTIONS, generate_model  # noqa: E402
from pyure.plotting import (  # noqa: E402
    build_species_plot,
    default_species_column,
    parse_species_input,
    validate_requested_species,
)
from pyure.simulation import SimulationSettings, simulate_sbml  # noqa: E402
from pyure.validation import normalize_dna_sequence  # noqa: E402


DEFAULT_SEQUENCE = "GGGATCCCGACTGGCGAGAGCCAGGTAACGAATGGATCCAA"


def reaction_search_terms(raw_query: str) -> list[str]:
    query = raw_query.strip()
    if not query:
        return []
    terms = {query}
    if "_" in query:
        material, name = query.split("_", 1)
        terms.add(f"{material}[{name}]")
        terms.add(name)
    if "[" in query and "]" in query:
        material = query.split("[", 1)[0]
        name = query.split("[", 1)[1].split("]", 1)[0]
        terms.add(f"{material}_{name}")
        terms.add(name)
    return [term for term in terms if term]


def find_reactions_for_species(crn: object, raw_query: str) -> list[str]:
    terms = [term.lower() for term in reaction_search_terms(raw_query)]
    if not terms:
        return []
    reactions = getattr(crn, "reactions", [])
    matches = []
    for reaction in reactions:
        reaction_text = str(reaction)
        reaction_text_lower = reaction_text.lower()
        if any(term in reaction_text_lower for term in terms):
            matches.append(reaction_text)
    return matches


st.set_page_config(
    page_title="pyure",
    page_icon="pyure",
    layout="wide",
)

st.markdown(
    """
    <style>
    .stApp, [data-testid="stAppViewContainer"] {
        background: #ffffff;
        color: #1f2933;
    }
    [data-testid="stHeader"] {background: rgba(255, 255, 255, 0.92);}
    .block-container {padding-top: 1.4rem; max-width: 1280px;}
    h1 {font-size: 2rem; margin-bottom: 0.2rem;}
    .pyure-nav a {margin-right: 1.1rem; color: #246a73; text-decoration: none;}
    .pyure-small {color: #52616b; font-size: 0.92rem; margin-top: 0.35rem;}
    .stButton > button, .stDownloadButton > button {border-radius: 6px;}
    .pyure-search-results {
        font-size: 0.92rem;
        line-height: 1.45;
        color: #25313a;
    }
    .pyure-reaction-results {
        max-height: 280px;
        overflow-y: auto;
        border: 1px solid #d7e0e2;
        border-radius: 6px;
        padding: 0.65rem 0.8rem;
        background: #f8faf9;
        font-family: ui-monospace, SFMono-Regular, Consolas, monospace;
        font-size: 0.86rem;
        line-height: 1.45;
        white-space: pre-wrap;
    }
    .pyure-footer {
        color: #52616b;
        font-size: 0.88rem;
        margin-top: 2rem;
        padding-top: 0.8rem;
        border-top: 1px solid #d7e0e2;
    }
    .pyure-slider-heading {
        margin-top: 1rem;
        margin-bottom: 0.15rem;
    }
    .pyure-slider-title {
        display: block;
        font-size: 1.28rem;
        font-weight: 700;
        color: #1f2933;
        margin-bottom: 0.25rem;
    }
    .pyure-slider-guidance {
        display: flex;
        align-items: center;
        gap: 0.65rem;
        color: #52616b;
        font-size: 0.98rem;
    }
    .pyure-detail-arrow {
        position: relative;
        display: inline-block;
        width: 150px;
        height: 3px;
        background: #246a73;
        border-radius: 999px;
    }
    .pyure-detail-arrow::after {
        content: "";
        position: absolute;
        right: -1px;
        top: 50%;
        width: 11px;
        height: 11px;
        border-top: 3px solid #246a73;
        border-right: 3px solid #246a73;
        transform: translateY(-50%) rotate(45deg);
    }
    div[data-testid="stSlider"] label p {
        font-size: 1.25rem;
        font-weight: 650;
    }
    div[data-testid="stSlider"] div[data-baseweb="slider"] {
        margin-top: 0.3rem;
        padding-top: 1rem;
        padding-bottom: 0.8rem;
    }
    div[data-testid="stSlider"] div[data-baseweb="slider"] > div {
        min-height: 34px;
    }
    div[data-testid="stSlider"] div[data-baseweb="slider"] > div:first-child,
    div[data-testid="stSlider"] div[data-baseweb="slider"] > div:first-child > div,
    div[data-testid="stSlider"] div[data-baseweb="slider"] > div:first-child > div > div {
        height: 18px !important;
        border-radius: 999px !important;
    }
    div[data-testid="stSlider"] div[data-baseweb="slider"] > div:first-child {
        box-shadow: inset 0 0 0 1px rgba(36, 106, 115, 0.18);
    }
    div[data-testid="stSlider"] [role="slider"] {
        width: 32px !important;
        height: 32px !important;
        border: 8px solid #ffffff !important;
        background: #246a73 !important;
        border-radius: 10px !important;
        transform: translateY(1px) rotate(45deg);
        box-shadow: 0 2px 8px rgba(31, 41, 51, 0.25) !important;
    }
    div[data-testid="stSlider"] [data-testid="stTickBar"],
    div[data-testid="stSlider"] [class*="tickBar"] {
        display: none !important;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

st.title("pyure: A python app for modeling and simulation of protein expression in PURE cell-free systems")
st.markdown(
    """
    <div class="pyure-nav">
      <a href="https://github.com/BuildACell/bioCRNpyler" target="_blank">BioCRNpyler</a>
      <a href="https://github.com/biocircuits/bioscrape" target="_blank">Bioscrape</a>
      <a href="https://www.biorxiv.org/content/10.64898/2026.02.22.707325v1.abstract" target="_blank">PURE models</a>
      <a href="https://biocrnpyler.readthedocs.io/en/latest/" target="_blank">BioCRNpyler docs</a>
    </div>
    """,
    unsafe_allow_html=True,
)
st.caption("Generate BioCRNpyler PURE models, export SBML, and run local Bioscrape simulations.")

hierarchy_labels = [label for _, label in HIERARCHY_OPTIONS]
st.markdown(
    """
    <div class="pyure-slider-heading">
      <span class="pyure-slider-title">Choose model hierarchy / level of detail</span>
      <span class="pyure-slider-guidance">
        <span>slide right to increase detail</span>
      </span>
    </div>
    """,
    unsafe_allow_html=True,
)
selected_label = st.select_slider(
    "Model hierarchy / level of detail",
    options=hierarchy_labels,
    value="nucleotide-level detail",
    help="Left to right increases model detail. The right-most level is the existing nucleotide-level detailed PURE model.",
    label_visibility="collapsed",
)
hierarchy_index = hierarchy_labels.index(selected_label)
selected_hierarchy = HIERARCHY_OPTIONS[hierarchy_index][0]

for key, default in {
    "generated_model": None,
    "simulation_df": None,
    "plot_species_text": "",
    "plot_species_submitted": [],
    "initial_conditions_df": None,
}.items():
    st.session_state.setdefault(key, default)

left, right = st.columns([0.42, 0.58], gap="large")

with left:
    st.subheader("Inputs")
    dna_sequence = st.text_area(
        "DNA sequence",
        value=DEFAULT_SEQUENCE,
        height=220,
        help="Whitespace is ignored. The detailed model currently expects A/C/G/T DNA with an ATG start codon.",
    )
    normalized = normalize_dna_sequence(dna_sequence)
    st.markdown(
        (
            f"<div class='pyure-small'>Length: {len(normalized)} nt<br>"
            "Detailed models for long sequences take longer to complete.</div>"
        ),
        unsafe_allow_html=True,
    )

    generate_clicked = st.button("Generate model", type="primary", width="stretch")
    if generate_clicked:
        st.session_state.simulation_df = None
        st.session_state.plot_species_text = ""
        st.session_state.initial_conditions_df = None
        try:
            with st.spinner("Generating BioCRNpyler model..."):
                st.session_state.generated_model = generate_model(
                    selected_hierarchy,
                    dna_sequence,
                )
            st.session_state.initial_conditions_df = initial_conditions_to_dataframe(
                st.session_state.generated_model.initial_conditions
            )
            st.success("Model generated.")
        except ModelNotImplementedError as exc:
            st.session_state.generated_model = None
            st.info(str(exc))
        except (SequenceValidationError, ModelGenerationError, PyureError) as exc:
            st.session_state.generated_model = None
            st.error(str(exc))
        except Exception as exc:
            st.session_state.generated_model = None
            st.error(f"Unexpected model-generation error: {exc}")

    st.divider()
    st.subheader("Simulation")
    sim_col_1, sim_col_2 = st.columns(2)
    with sim_col_1:
        sim_end = st.number_input("End time", min_value=1.0, value=3600.0, step=100.0)
    with sim_col_2:
        sim_points = st.number_input("Time points", min_value=2, value=200, step=25)
    model_for_sim = st.session_state.generated_model
    if model_for_sim is not None:
        if st.session_state.initial_conditions_df is None:
            st.session_state.initial_conditions_df = initial_conditions_to_dataframe(
                model_for_sim.initial_conditions
            )
        st.markdown("**Initial conditions**")
        st.caption("Only non-zero values shown. Edit values and hit Run Simulation.")
        st.session_state.initial_conditions_df = st.data_editor(
            st.session_state.initial_conditions_df,
            hide_index=True,
            width="stretch",
            num_rows="dynamic",
            column_config={
                "Species": st.column_config.TextColumn("Species"),
                "Initial concentration": st.column_config.NumberColumn(
                    "Initial concentration",
                    min_value=0.0,
                    format="%.6g",
                ),
            },
        )

    simulate_clicked = st.button(
        "Run Simulation",
        width="stretch",
        disabled=st.session_state.generated_model is None,
    )
    if simulate_clicked and st.session_state.generated_model is not None:
        try:
            settings = SimulationSettings(end=float(sim_end), points=int(sim_points))
            initial_conditions = dataframe_to_initial_conditions(
                st.session_state.initial_conditions_df
            )
            with st.spinner("Running Bioscrape simulation..."):
                st.session_state.simulation_df = simulate_sbml(
                    st.session_state.generated_model.sbml_xml,
                    settings,
                    initial_conditions=initial_conditions,
                )
            first_species = default_species_column(st.session_state.simulation_df)
            st.session_state.plot_species_text = first_species or ""
            st.session_state.plot_species_submitted = [first_species] if first_species else []
            st.success("Simulation complete.")
        except SimulationError as exc:
            st.error(str(exc))
        except Exception as exc:
            st.error(f"Unexpected simulation error: {exc}")

with right:
    st.subheader("Outputs")
    model = st.session_state.generated_model
    if model is None:
        st.info("Generate an implemented hierarchy level to see summaries, SBML export, and simulation output.")
    else:
        metric_cols = st.columns(4)
        metric_cols[0].metric("Species", model.summary["Species"])
        metric_cols[1].metric("Reactions", model.summary["Reactions"])
        metric_cols[2].metric("DNA nt", model.summary["DNA length"])
        metric_cols[3].metric("Build s", f"{model.build_seconds:.2f}")

        st.markdown("**Model summary**")
        summary_df = pd.DataFrame(
            [{"Field": str(key), "Value": str(value)} for key, value in model.summary.items()]
        )
        st.dataframe(summary_df, hide_index=True, width="stretch")

        st.download_button(
            "Download generated SBML model",
            data=model.sbml_xml,
            file_name=download_filename(model.hierarchy.value),
            mime="application/xml",
            width="stretch",
        )

        species_query = st.text_input(
            "Search species by partial name",
            placeholder="Type part of a species name, for example ATP or mRNA",
        )
        if species_query:
            matches = [
                name
                for name in model.species_names
                if species_query.lower() in name.lower()
            ]
            if matches:
                st.markdown(
                    "<div class='pyure-search-results'>"
                    + ", ".join(matches[:100])
                    + "</div>",
                    unsafe_allow_html=True,
                )
                if len(matches) > 100:
                    st.caption(f"Showing 100 of {len(matches)} matches.")
            else:
                st.info("No species matched that search.")

        reaction_query = st.text_input(
            "Search reactions by species",
            placeholder="Enter a species, for example protein_GFP, dna_target, or ATP",
        )
        if reaction_query:
            reaction_matches = find_reactions_for_species(model.crn, reaction_query)
            if reaction_matches:
                st.caption(f"Found {len(reaction_matches)} associated reaction(s).")
                st.markdown(
                    "<div class='pyure-reaction-results'>"
                    + "\n".join(reaction_matches[:150])
                    + "</div>",
                    unsafe_allow_html=True,
                )
                if len(reaction_matches) > 150:
                    st.caption(f"Showing 150 of {len(reaction_matches)} matches.")
            else:
                st.info("No reactions matched that species.")

    simulation_df = st.session_state.simulation_df
    if simulation_df is not None:
        st.divider()
        st.subheader("Simulation plot")
        with st.form("species_plot_form", clear_on_submit=False):
            species_text = st.text_area(
                "Species to plot (separate multiple species by comma or newline)",
                value=st.session_state.plot_species_text,
                height=90,
                help="Examples: ATP, GTP or one species per line.",
            )
            plot_clicked = st.form_submit_button("Plot", width="stretch")

        if plot_clicked:
            st.session_state.plot_species_text = species_text
            requested_species = parse_species_input(species_text)
            if not requested_species:
                default_species = default_species_column(simulation_df)
                requested_species = [default_species] if default_species else []
            st.session_state.plot_species_submitted = requested_species

        requested_species = st.session_state.plot_species_submitted
        if not requested_species:
            default_species = default_species_column(simulation_df)
            requested_species = [default_species] if default_species else []

        found, missing = validate_requested_species(
            requested_species,
            [str(column) for column in simulation_df.columns],
        )
        if missing:
            st.warning("Species not found: " + ", ".join(missing))
        if found:
            st.plotly_chart(
                build_species_plot(simulation_df, found),
                width="stretch",
            )
        else:
            st.info("Enter at least one species name present in the simulation dataframe.")

        with st.expander("Simulation dataframe"):
            st.dataframe(simulation_df, width="stretch")

st.markdown(
    "<div class='pyure-footer'>pyure is released under the BSD 3-Clause License.</div>",
    unsafe_allow_html=True,
)
