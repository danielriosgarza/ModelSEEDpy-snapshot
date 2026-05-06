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


LOGGER = logging.getLogger("reconstruct_non_gapfilled_model")

REFERENCE_FASTA = REPO_ROOT / "genomes" / "GCF_000010425.1_ASM1042v1_protein.faa.gz"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "reconstructed_models"
DEFAULT_DATABASE_PATH = REPO_ROOT / "ModelSEEDDatabase"
DEFAULT_TEMPLATE_DIR = DEFAULT_OUTPUT_DIR / "templates"
DEFAULT_FASTA_GLOBS = ("*.faa", "*.faa.gz", "*.fna", "*.fna.gz")
SUPPORTED_FASTA_SUFFIXES = (
    ".faa",
    ".faa.gz",
    ".fa",
    ".fa.gz",
    ".fasta",
    ".fasta.gz",
    ".fna",
    ".fna.gz",
)
UPDATED_TEMPLATE_FILES = {
    "positive": "GramPositive.latest.modeltemplate.json",
    "negative": "GramNegative.latest.modeltemplate.json",
}
GRAM_STATUS_BY_ACCESSION = {
    # GCF_000005845.2_ASM584v2 is Escherichia coli str. K-12 substr. MG1655.
    "GCF_000005845.2": "negative",
    # GCF_000010425.1_ASM1042v1 is Bifidobacterium adolescentis.
    "GCF_000010425.1": "positive",
}
GRAM_STATUS_BY_GENUS = {
    # Curated for the current genomes/ batch. Add genera here or use --gram-map
    # for new batches.
    "bacillus": "positive",
    "bifidobacterium": "positive",
    "enterococcus": "positive",
    "lacticaseibacillus": "positive",
    "lactiplantibacillus": "positive",
    "lactococcus": "positive",
    "latilactobacillus": "positive",
    "lentilactobacillus": "positive",
    "leuconostoc": "positive",
    "levilactobacillus": "positive",
    "liquorilactobacillus": "positive",
    "loigolactobacillus": "positive",
    "pediococcus": "positive",
    "schleiferilactobacillus": "positive",
    "escherichia": "negative",
    "hafnia": "negative",
    "pseudomonas": "negative",
    "rahnella": "negative",
}


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
        "--fasta-dir",
        type=Path,
        help=(
            "Directory of FASTA files to reconstruct as a batch. Supports "
            ".faa, .faa.gz, .fna, and .fna.gz by default. "
            "When set, --model-id, --output, and --sbml-output cannot be used."
        ),
    )
    parser.add_argument(
        "--fasta-glob",
        help=(
            "Comma-separated glob(s) used with --fasta-dir. Default: "
            "*.faa,*.faa.gz,*.fna,*.fna.gz."
        ),
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
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=(
            "Directory for default JSON outputs. Used for single runs when "
            "--output is omitted and for all --fasta-dir batch outputs. "
            f"Default: {DEFAULT_OUTPUT_DIR}."
        ),
    )
    parser.add_argument(
        "--protein-output-dir",
        type=Path,
        help=(
            "Directory for .faa files generated from nucleotide .fna inputs. "
            "Defaults to <output-dir>/protein_fastas."
        ),
    )
    parser.add_argument(
        "--sbml-output",
        type=Path,
        help="Optional SBML output path. Requires python-libsbml.",
    )
    parser.add_argument(
        "--template-id",
        help=(
            "Legacy packaged template to use, such as template_gram_pos, "
            "template_gram_neg, template_core, or auto. Ignored when "
            "--template-json is provided. If omitted, the script uses the "
            "updated local Gram-positive/Gram-negative templates."
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
    parser.add_argument(
        "--updated-template-dir",
        type=Path,
        default=DEFAULT_TEMPLATE_DIR,
        help=(
            "Directory containing updated GramPositive.latest.modeltemplate.json "
            f"and GramNegative.latest.modeltemplate.json. Default: {DEFAULT_TEMPLATE_DIR}"
        ),
    )
    parser.add_argument(
        "--gram-status",
        choices=("auto", "positive", "negative"),
        default="auto",
        help=(
            "Gram status used to select the updated local template when "
            "--template-json and --template-id are omitted. Default: auto."
        ),
    )
    parser.add_argument(
        "--gram-map",
        type=Path,
        help=(
            "Optional TSV/CSV mapping accession, genus, model ID, or FASTA "
            "filename to gram status. Rows are key,status."
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
        "--dry-run-template-selection",
        action="store_true",
        help=(
            "Print model ID, FASTA, inferred gram status, and template path, "
            "then exit without loading genomes, annotating, or reconstructing."
        ),
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
    for suffix in SUPPORTED_FASTA_SUFFIXES:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return re.sub(r"[^A-Za-z0-9_]+", "_", name).strip("_")


def fasta_name_without_sequence_suffix(fasta: Path) -> str:
    name = fasta.name
    for suffix in SUPPORTED_FASTA_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return fasta.stem


def is_supported_fasta(path: Path) -> bool:
    name = path.name.lower()
    return any(name.endswith(suffix) for suffix in SUPPORTED_FASTA_SUFFIXES)


def is_nucleotide_fasta(path: Path) -> bool:
    name = path.name.lower()
    return name.endswith(".fna") or name.endswith(".fna.gz")


def accession_from_fasta(fasta: Path) -> Optional[str]:
    match = re.match(r"^(GC[AF]_\d+\.\d+)", fasta_name_without_sequence_suffix(fasta))
    return match.group(1) if match else None


def genus_from_fasta(fasta: Path) -> Optional[str]:
    name = fasta_name_without_sequence_suffix(fasta)
    parts = name.split("_")
    if len(parts) >= 3 and parts[0] in {"GCA", "GCF"}:
        return parts[2]
    return parts[0] if parts else None


def normalize_gram_status(value: str) -> str:
    normalized = re.sub(r"[\s_-]+", "-", value.strip().lower())
    positive_values = {"p", "pos", "positive", "gram-positive", "grampositive"}
    negative_values = {"n", "neg", "negative", "gram-negative", "gramnegative"}
    if normalized in positive_values:
        return "positive"
    if normalized in negative_values:
        return "negative"
    raise ValueError(f"Unsupported gram status: {value}")


def load_gram_map(path: Optional[Path]) -> dict[str, str]:
    if not path:
        return {}

    mapping = {}
    with open(path, "r") as fh:
        for line_number, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = [field.strip() for field in re.split(r"[\t,]", line, maxsplit=1)]
            if len(fields) != 2:
                raise ValueError(f"Invalid gram map row {line_number}: {line}")
            key, status = fields
            if key.lower() in {"key", "accession", "genus", "model_id", "fasta"}:
                continue
            mapping[key.lower()] = normalize_gram_status(status)
    return mapping


def infer_gram_status(fasta: Path, model_id: str, gram_map: dict[str, str]) -> str:
    accession = accession_from_fasta(fasta)
    genus = genus_from_fasta(fasta)
    lookup_keys = [
        model_id.lower(),
        fasta.name.lower(),
        accession.lower() if accession else None,
        genus.lower() if genus else None,
    ]
    for key in lookup_keys:
        if key and key in gram_map:
            return gram_map[key]
    if accession and accession in GRAM_STATUS_BY_ACCESSION:
        return GRAM_STATUS_BY_ACCESSION[accession]
    if genus and genus.lower() in GRAM_STATUS_BY_GENUS:
        return GRAM_STATUS_BY_GENUS[genus.lower()]

    raise ValueError(
        "Could not infer gram status for "
        f"{fasta.name}. Pass --gram-status positive|negative or add it to --gram-map."
    )


def load_template(template_id: str):
    from modelseedpy.core.mstemplate import MSTemplateBuilder
    from modelseedpy.helpers import get_template

    if template_id == "auto":
        return None

    template_data = get_template(template_id)
    return MSTemplateBuilder.from_dict(template_data).build()


def load_template_json(template_path: Path):
    import json

    from modelseedpy.core.mstemplate import MSTemplateBuilder

    with open(template_path, "r") as fh:
        template_data = json.load(fh)
    template_data.setdefault("__VERSION__", 1)
    return MSTemplateBuilder.from_dict(template_data).build()


def select_template(
    args: argparse.Namespace,
    fasta: Path,
    model_id: str,
    gram_map: dict[str, str],
):
    if args.template_json:
        template_path = args.template_json.resolve()
        return load_template_json(template_path), template_path, None

    if args.template_id:
        return load_template(args.template_id), args.template_id, None

    gram_status = (
        infer_gram_status(fasta, model_id, gram_map)
        if args.gram_status == "auto"
        else args.gram_status
    )
    template_path = (
        args.updated_template_dir.resolve() / UPDATED_TEMPLATE_FILES[gram_status]
    )
    if not template_path.exists():
        raise FileNotFoundError(
            f"Updated {gram_status} template not found: {template_path}. "
            "Build it with examples/Model Reconstruction/"
            "build_template_from_modelseed_database.py."
        )
    return load_template_json(template_path), template_path, gram_status


def configure_biochemistry(args: argparse.Namespace) -> None:
    from modelseedpy.biochem.modelseed_biochem import (
        ModelSEEDBiochem,
        from_github,
        from_local,
    )
    from modelseedpy.helpers import config

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
    elif DEFAULT_DATABASE_PATH.exists():
        config.set("biochem", "path", str(DEFAULT_DATABASE_PATH))
        ModelSEEDBiochem.default_biochemistry = from_local(str(DEFAULT_DATABASE_PATH))
        LOGGER.info(
            "Using default local ModelSEED biochemistry at %s", DEFAULT_DATABASE_PATH
        )


def product_from_description(description: Optional[str]) -> Optional[str]:
    if not description:
        return None

    product = re.sub(r"\s*\[[^\]]+\]\s*$", "", description).strip()
    product = re.sub(r"^MULTISPECIES:\s*", "", product, flags=re.IGNORECASE)
    product = re.sub(r"\s+", " ", product).strip()
    return product or None


def add_fasta_description_annotations(genome: MSGenome) -> int:
    from modelseedpy.core.rast_client import split_annotation

    annotated_features = 0
    for feature in genome.features:
        product = product_from_description(feature.description)
        if not product:
            continue
        for role in split_annotation(product):
            feature.add_ontology_term("RAST", role)
        annotated_features += 1
    return annotated_features


def open_text_sequence_file(path: Path):
    if path.name.lower().endswith(".gz"):
        import gzip

        return gzip.open(path, "rt")
    return open(path, "r")


def iter_fasta_records(path: Path):
    record_id = None
    sequence_lines = []
    with open_text_sequence_file(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if record_id is not None:
                    yield record_id, "".join(sequence_lines)
                record_id = line[1:].split(None, 1)[0]
                sequence_lines = []
            else:
                sequence_lines.append(line)
    if record_id is not None:
        yield record_id, "".join(sequence_lines)


def protein_fasta_dir(args: argparse.Namespace) -> Path:
    if args.protein_output_dir:
        return args.protein_output_dir.resolve()
    return (args.output_dir / "protein_fastas").resolve()


def convert_fna_to_faa(fna_path: Path, faa_path: Path) -> Path:
    try:
        import pyrodigal
    except ImportError as err:
        raise ImportError(
            "Nucleotide FASTA inputs (.fna/.fna.gz) require pyrodigal. "
            "Install it with `uv pip install pyrodigal` or rerun the batch "
            "script with `-InstallMissingPackages`."
        ) from err

    if faa_path.exists() and faa_path.stat().st_mtime >= fna_path.stat().st_mtime:
        LOGGER.info("Using existing translated protein FASTA %s", faa_path)
        return faa_path

    faa_path.parent.mkdir(parents=True, exist_ok=True)
    gene_finder = pyrodigal.GeneFinder(meta=True)
    translated_contigs = 0
    with open(faa_path, "w") as fh:
        for sequence_id, sequence in iter_fasta_records(fna_path):
            if not sequence:
                continue
            genes = gene_finder.find_genes(sequence)
            genes.write_translations(fh, sequence_id=sequence_id)
            translated_contigs += 1

    LOGGER.info(
        "Translated %d contigs from %s to %s",
        translated_contigs,
        fna_path,
        faa_path,
    )
    return faa_path


def prepare_protein_fasta(input_fasta: Path, args: argparse.Namespace) -> Path:
    if not is_nucleotide_fasta(input_fasta):
        return input_fasta

    translated_name = f"{fasta_name_without_sequence_suffix(input_fasta)}.faa"
    return convert_fna_to_faa(input_fasta, protein_fasta_dir(args) / translated_name)


def iter_fasta_inputs(args: argparse.Namespace) -> list[Path]:
    if args.fasta_dir:
        if args.model_id or args.output or args.sbml_output:
            raise ValueError(
                "--model-id, --output, and --sbml-output are only valid for "
                "single-FASTA runs."
            )
        fasta_dir = args.fasta_dir.resolve()
        if not fasta_dir.exists():
            raise FileNotFoundError(f"FASTA directory not found: {fasta_dir}")
        patterns = (
            [pattern.strip() for pattern in args.fasta_glob.split(",")]
            if args.fasta_glob
            else DEFAULT_FASTA_GLOBS
        )
        fasta_files = []
        seen = set()
        for pattern in patterns:
            for fasta in sorted(fasta_dir.glob(pattern)):
                resolved = fasta.resolve()
                if resolved in seen or not is_supported_fasta(fasta):
                    continue
                seen.add(resolved)
                fasta_files.append(fasta)
        if not fasta_files:
            raise FileNotFoundError(
                f"No supported FASTA files matched {','.join(patterns)} in {fasta_dir}"
            )
        return fasta_files

    fasta = args.fasta.resolve()
    if not is_supported_fasta(fasta):
        raise ValueError(f"Unsupported FASTA extension: {fasta.name}")
    return [fasta]


def describe_template_selection(
    args: argparse.Namespace,
    fasta: Path,
    model_id: str,
    gram_map: dict[str, str],
) -> tuple[str, str]:
    if args.template_json:
        return "explicit", str(args.template_json.resolve())
    if args.template_id:
        return "legacy-auto" if args.template_id == "auto" else "legacy", args.template_id

    gram_status = (
        infer_gram_status(fasta, model_id, gram_map)
        if args.gram_status == "auto"
        else args.gram_status
    )
    template_path = (
        args.updated_template_dir.resolve() / UPDATED_TEMPLATE_FILES[gram_status]
    )
    return gram_status, str(template_path)


def reconstruct_one(args: argparse.Namespace, fasta: Path, gram_map: dict[str, str]) -> None:
    import cobra.io

    from modelseedpy.core.msbuilder import MSBuilder
    from modelseedpy.core.msgenome import MSGenome

    fasta = fasta.resolve()
    if not fasta.exists():
        raise FileNotFoundError(f"FASTA input not found: {fasta}")

    model_id = args.model_id or model_id_from_fasta(fasta)
    output = args.output or args.output_dir / f"{model_id}_non_gapfilled.json"
    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)

    protein_fasta = prepare_protein_fasta(fasta, args)
    LOGGER.info("Loading genome from %s", protein_fasta)
    genome = MSGenome.from_fasta(str(protein_fasta))

    annotate_with_rast = args.annotation_source == "rast"
    if args.annotation_source == "fasta-description":
        annotated = add_fasta_description_annotations(genome)
        LOGGER.info("Seeded RAST-like annotations from %d FASTA descriptions", annotated)

    template, template_source, gram_status = select_template(
        args, fasta, model_id, gram_map
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
    print(f"Input FASTA: {fasta}")
    if protein_fasta != fasta:
        print(f"Protein FASTA: {protein_fasta}")
    print(f"Features: {len(genome.features)}")
    if gram_status:
        print(f"Gram status: {gram_status}")
    print(f"Template source: {template_source}")
    print(f"Template: {builder.template.id} ({builder.template.name})")
    print(f"Reactions: {len(model.reactions)}")
    print(f"Metabolites: {len(model.metabolites)}")
    print(f"Genes: {len(model.genes)}")
    print(f"Objective: {model.objective.expression}")
    print(f"JSON: {output}")
    if args.sbml_output:
        print(f"SBML: {args.sbml_output.resolve()}")
    print("Gapfilling: not run")


def main() -> None:
    args = parse_args()
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(levelname)s:%(name)s:%(message)s",
    )

    gram_map = load_gram_map(args.gram_map.resolve() if args.gram_map else None)
    fasta_files = iter_fasta_inputs(args)
    if args.dry_run_template_selection:
        print("model_id\tfasta\tgram_status\ttemplate_source")
        for fasta in fasta_files:
            model_id = args.model_id or model_id_from_fasta(fasta)
            gram_status, template_source = describe_template_selection(
                args, fasta, model_id, gram_map
            )
            print(f"{model_id}\t{fasta.name}\t{gram_status}\t{template_source}")
        return

    configure_biochemistry(args)
    for fasta in fasta_files:
        reconstruct_one(args, fasta, gram_map)


if __name__ == "__main__":
    main()
