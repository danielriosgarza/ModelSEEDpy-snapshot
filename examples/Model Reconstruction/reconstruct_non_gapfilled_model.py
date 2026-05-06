#!/usr/bin/env python
"""Reference script for reconstructing a non-gapfilled ModelSEEDpy model.

This script intentionally uses ``MSBuilder.build_base_model()`` instead of
``MSBuilder.build()`` or ``MSBuilder.gapfill_model()``. The output is therefore
the annotation-derived draft model before ATP correction or gapfilling.
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
from pathlib import Path
from typing import Optional

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import cobra.io

from modelseedpy.biochem.modelseed_biochem import (
    ModelSEEDBiochem,
    from_github,
    from_local,
)
from modelseedpy.core.msbuilder import MSBuilder
from modelseedpy.core.msgenome import MSGenome
from modelseedpy.core.mstemplate import MSTemplateBuilder
from modelseedpy.core.rast_client import split_annotation
from modelseedpy.helpers import config, get_template


LOGGER = logging.getLogger("reconstruct_non_gapfilled_model")

REFERENCE_FASTA = REPO_ROOT / "genomes" / "GCF_000010425.1_ASM1042v1_protein.faa.gz"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "reconstructed_models"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build an annotation-derived draft metabolic model without running "
            "ATP correction or gapfilling."
        )
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        default=REFERENCE_FASTA,
        help=f"Protein FASTA input. Default: {REFERENCE_FASTA}",
    )
    parser.add_argument(
        "--model-id",
        help="COBRA model ID. Defaults to a sanitized ID derived from --fasta.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help=(
            "Output JSON path. Defaults to reconstructed_models/"
            "<model-id>_non_gapfilled.json."
        ),
    )
    parser.add_argument(
        "--sbml-output",
        type=Path,
        help="Optional SBML output path. Requires python-libsbml.",
    )
    parser.add_argument(
        "--template-id",
        default="template_gram_pos",
        help=(
            "Template to use, such as template_gram_pos, template_gram_neg, "
            "template_core, or auto. Ignored when --template-json is provided. "
            "The reference genome is Bifidobacterium, "
            "so template_gram_pos is used by default."
        ),
    )
    parser.add_argument(
        "--template-json",
        type=Path,
        help=(
            "Load a template JSON directly, for example one generated from a "
            "current ModelSEEDDatabase checkout."
        ),
    )
    biochem_group = parser.add_mutually_exclusive_group()
    biochem_group.add_argument(
        "--biochem-path",
        type=Path,
        help=(
            "Path to a local ModelSEEDDatabase checkout to use for direct "
            "ModelSEED biochemistry lookups."
        ),
    )
    biochem_group.add_argument(
        "--biochem-github-ref",
        help=(
            "ModelSEEDDatabase git ref to load directly from GitHub for "
            "biochemistry lookups, such as master or a commit SHA."
        ),
    )
    parser.add_argument(
        "--annotation-source",
        choices=("rast", "fasta-description", "none"),
        default="rast",
        help=(
            "Use remote RAST annotation, approximate roles from FASTA product "
            "descriptions, or build without adding annotations. Default: rast."
        ),
    )
    parser.add_argument(
        "--index",
        default="0",
        help="Compartment index suffix used by ModelSEEDpy. Default: 0.",
    )
    parser.add_argument(
        "--allow-all-non-gpr-reactions",
        action="store_true",
        help="Include all template universal/spontaneous reactions.",
    )
    parser.add_argument(
        "--classic-biomass",
        action="store_true",
        help="Use classic biomass construction where supported by the template.",
    )
    parser.add_argument(
        "--gc",
        type=float,
        default=0.5,
        help="GC fraction passed to template biomass construction. Default: 0.5.",
    )
    parser.add_argument(
        "--add-atpm",
        action="store_true",
        help="Add ATP maintenance after draft reconstruction. No gapfilling is run.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"),
        help="Logging verbosity. Default: INFO.",
    )
    return parser.parse_args()


def model_id_from_fasta(fasta: Path) -> str:
    name = fasta.name
    for suffix in (".faa.gz", ".fasta.gz", ".fa.gz", ".faa", ".fasta", ".fa"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return re.sub(r"[^A-Za-z0-9_]+", "_", name).strip("_")


def load_template(template_id: str):
    if template_id == "auto":
        return None

    template_data = get_template(template_id)
    return MSTemplateBuilder.from_dict(template_data).build()


def load_template_json(template_path: Path):
    import json

    with open(template_path, "r") as fh:
        template_data = json.load(fh)
    template_data.setdefault("__VERSION__", 1)
    return MSTemplateBuilder.from_dict(template_data).build()


def configure_biochemistry(args: argparse.Namespace) -> None:
    if args.biochem_path:
        biochem_path = args.biochem_path.resolve()
        if not biochem_path.exists():
            raise FileNotFoundError(f"ModelSEEDDatabase path not found: {biochem_path}")
        config.set("biochem", "path", str(biochem_path))
        ModelSEEDBiochem.default_biochemistry = from_local(str(biochem_path))
        LOGGER.info("Using local ModelSEED biochemistry at %s", biochem_path)
    elif args.biochem_github_ref:
        ModelSEEDBiochem.default_biochemistry = from_github(args.biochem_github_ref)
        LOGGER.info(
            "Using ModelSEED biochemistry from GitHub ref %s",
            args.biochem_github_ref,
        )


def product_from_description(description: Optional[str]) -> Optional[str]:
    if not description:
        return None

    product = re.sub(r"\s*\[[^\]]+\]\s*$", "", description).strip()
    product = re.sub(r"^MULTISPECIES:\s*", "", product, flags=re.IGNORECASE)
    product = re.sub(r"\s+", " ", product).strip()
    return product or None


def add_fasta_description_annotations(genome: MSGenome) -> int:
    annotated_features = 0
    for feature in genome.features:
        product = product_from_description(feature.description)
        if not product:
            continue
        for role in split_annotation(product):
            feature.add_ontology_term("RAST", role)
        annotated_features += 1
    return annotated_features


def main() -> None:
    args = parse_args()
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(levelname)s:%(name)s:%(message)s",
    )

    fasta = args.fasta.resolve()
    if not fasta.exists():
        raise FileNotFoundError(f"Protein FASTA not found: {fasta}")

    model_id = args.model_id or model_id_from_fasta(fasta)
    output = args.output or DEFAULT_OUTPUT_DIR / f"{model_id}_non_gapfilled.json"
    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)

    LOGGER.info("Loading genome from %s", fasta)
    genome = MSGenome.from_fasta(str(fasta))

    annotate_with_rast = args.annotation_source == "rast"
    if args.annotation_source == "fasta-description":
        annotated = add_fasta_description_annotations(genome)
        LOGGER.info("Seeded RAST-like annotations from %d FASTA descriptions", annotated)

    configure_biochemistry(args)
    template = (
        load_template_json(args.template_json.resolve())
        if args.template_json
        else load_template(args.template_id)
    )
    builder = MSBuilder(genome, template=template)

    LOGGER.info("Building non-gapfilled draft model")
    model = builder.build_base_model(
        model_id,
        index=args.index,
        allow_all_non_grp_reactions=args.allow_all_non_gpr_reactions,
        annotate_with_rast=annotate_with_rast,
        biomass_classic=args.classic_biomass,
        biomass_gc=args.gc,
        add_reaction_from_rast_annotation=True,
    )

    if args.add_atpm:
        MSBuilder.add_atpm(model)

    cobra.io.save_json_model(model, str(output))

    if args.sbml_output:
        sbml_output = args.sbml_output.resolve()
        sbml_output.parent.mkdir(parents=True, exist_ok=True)
        cobra.io.write_sbml_model(model, str(sbml_output))

    print(f"Model ID: {model.id}")
    print(f"Features: {len(genome.features)}")
    print(f"Template: {builder.template.id} ({builder.template.name})")
    print(f"Reactions: {len(model.reactions)}")
    print(f"Metabolites: {len(model.metabolites)}")
    print(f"Genes: {len(model.genes)}")
    print(f"Objective: {model.objective.expression}")
    print(f"JSON: {output}")
    if args.sbml_output:
        print(f"SBML: {args.sbml_output.resolve()}")
    print("Gapfilling: not run")


if __name__ == "__main__":
    main()
