#!/usr/bin/env python
"""Build a ModelSEEDpy template JSON from a local ModelSEEDDatabase checkout.

The ModelSEEDDatabase repository stores current templates as TSV files plus the
current master biochemistry JSON. This script converts those source files into
the template JSON shape consumed by ``MSTemplateBuilder.from_dict()``.
"""

from __future__ import annotations

import argparse
import builtins
import contextlib
import io
import json
import logging
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, MutableMapping

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATABASE_PATH = REPO_ROOT / "ModelSEEDDatabase"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "reconstructed_models" / "templates"

LOGGER = logging.getLogger("build_template_from_modelseed_database")

TEMPLATE_NAME_TO_ID = {
    "Core": "CoreModelTemplateV2",
    "GramNegative": "GramNegative.modeltemplate",
    "GramPositive": "GramPositive.modeltemplate",
    "Plant": "Plant.modeltemplate",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate a ModelSEEDpy-compatible template JSON from a local "
            "ModelSEEDDatabase checkout."
        )
    )
    parser.add_argument(
        "--database-path",
        type=Path,
        default=DEFAULT_DATABASE_PATH,
        help=f"Path to ModelSEEDDatabase. Default: {DEFAULT_DATABASE_PATH}",
    )
    parser.add_argument(
        "--template",
        default="GramPositive",
        help=(
            "Template folder under ModelSEEDDatabase/Templates, such as "
            "GramPositive, GramNegative, Core, or Plant. Default: GramPositive."
        ),
    )
    parser.add_argument(
        "--template-id",
        help="Template object ID. Defaults to the historical ID for known templates.",
    )
    parser.add_argument(
        "--name",
        help="Template display name. Defaults to --template-id.",
    )
    parser.add_argument(
        "--domain",
        default="Bacteria",
        help="Template domain value. Default: Bacteria.",
    )
    parser.add_argument(
        "--type",
        default="GenomeScale",
        help="Template type value. Default: GenomeScale.",
    )
    parser.add_argument(
        "--biochem-ref",
        default="local/ModelSEEDDatabase",
        help=(
            "Reference string stored in template['biochemistry_ref']. "
            "Default: local/ModelSEEDDatabase."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        help=(
            "Output template JSON path. Defaults to "
            "reconstructed_models/templates/<template-id>.json."
        ),
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"),
        help="Logging verbosity. Default: INFO.",
    )
    parser.add_argument(
        "--show-template-notices",
        action="store_true",
        help="Show notices printed by the legacy ModelSEEDDatabase template helper.",
    )
    return parser.parse_args()


def normalize_nulls(value: Any) -> Any:
    if value is None:
        return "null"
    if isinstance(value, list):
        return [normalize_nulls(item) for item in value]
    if isinstance(value, dict):
        return {key: normalize_nulls(item) for key, item in value.items()}
    return value


def index_records_by_id(records: Any, label: str) -> Dict[str, MutableMapping[str, Any]]:
    if isinstance(records, dict):
        return {key: normalize_nulls(value) for key, value in records.items()}

    indexed = {}
    for record in records:
        record = normalize_nulls(record)
        record_id = record.get("id")
        if not record_id:
            raise ValueError(f"{label} record is missing an id: {record}")
        indexed[record_id] = record
    return indexed


def load_template_helper(database_path: Path):
    scripts_path = database_path / "Scripts"
    helper_path = scripts_path / "TemplateHelper.py"
    if not helper_path.exists():
        raise FileNotFoundError(f"TemplateHelper.py not found: {helper_path}")

    if not hasattr(builtins, "long"):
        builtins.long = int
    if str(database_path) not in sys.path:
        sys.path.insert(0, str(database_path))

    from Scripts.TemplateHelper import TemplateHelper

    return TemplateHelper


def add_reaction_refs(
    reactions: Iterable[MutableMapping[str, Any]], biochem_ref: str
) -> None:
    for reaction in reactions:
        base_reaction_id = reaction["id"].rsplit("_", 1)[0]
        reaction["reaction_ref"] = f"{biochem_ref}/reactions/id/{base_reaction_id}"


def build_template(args: argparse.Namespace) -> Dict[str, Any]:
    database_path = args.database_path.resolve()
    template_dir = database_path / "Templates" / args.template
    if not template_dir.exists():
        raise FileNotFoundError(f"Template source folder not found: {template_dir}")

    TemplateHelper = load_template_helper(database_path)
    helper = TemplateHelper(
        str(database_path / "Biochemistry" / "compounds.json"),
        str(database_path / "Biochemistry" / "reactions.json"),
    )
    helper.masterCompounds = index_records_by_id(helper.masterCompounds, "compound")
    helper.masterReactions = index_records_by_id(helper.masterReactions, "reaction")

    def read_template_sources() -> None:
        helper.readCompartmentsFile(str(template_dir / "Compartments.tsv"), False)
        helper.readBiomassesFile(
            str(template_dir / "Biomasses.tsv"),
            str(template_dir / "BiomassCompounds.tsv"),
            False,
        )
        helper.readRolesFile(str(database_path / "Annotations" / "Roles.tsv"), False)
        helper.readComplexesFile(
            str(database_path / "Annotations" / "Complexes.tsv"), False
        )
        helper.readReactionsFile(str(template_dir / "Reactions.tsv"), False)

    if args.show_template_notices:
        read_template_sources()
    else:
        notices = io.StringIO()
        with contextlib.redirect_stdout(notices):
            read_template_sources()
        notice_lines = [line for line in notices.getvalue().splitlines() if line]
        if notice_lines:
            LOGGER.info(
                "Suppressed %d template helper notices; rerun with "
                "--show-template-notices to inspect them.",
                len(notice_lines),
            )

    add_reaction_refs(helper.reactions.values(), args.biochem_ref)

    template_id = args.template_id or TEMPLATE_NAME_TO_ID.get(
        args.template, f"{args.template}.modeltemplate"
    )
    name = args.name or template_id

    return {
        "__VERSION__": 1,
        "id": template_id,
        "name": name,
        "type": args.type,
        "domain": args.domain,
        "biochemistry_ref": args.biochem_ref,
        "pathways": [],
        "subsystems": [],
        "compartments": list(helper.compartments.values()),
        "biomasses": list(helper.biomasses.values()),
        "roles": list(helper.roles.values()),
        "complexes": list(helper.complexes.values()),
        "reactions": list(helper.reactions.values()),
        "compounds": list(helper.compounds.values()),
        "compcompounds": list(helper.compCompounds.values()),
    }


def main() -> None:
    args = parse_args()
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(levelname)s:%(name)s:%(message)s",
    )

    template = build_template(args)
    output = args.output
    if output is None:
        output = DEFAULT_OUTPUT_DIR / f"{template['id']}.json"
    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w") as fh:
        json.dump(template, fh, indent=2, sort_keys=True)
        fh.write("\n")

    print(f"Template ID: {template['id']}")
    print(f"Template source: {args.database_path.resolve() / 'Templates' / args.template}")
    print(f"Biochemistry ref: {template['biochemistry_ref']}")
    print(f"Compartments: {len(template['compartments'])}")
    print(f"Biomasses: {len(template['biomasses'])}")
    print(f"Roles: {len(template['roles'])}")
    print(f"Complexes: {len(template['complexes'])}")
    print(f"Reactions: {len(template['reactions'])}")
    print(f"Compounds: {len(template['compounds'])}")
    print(f"Compartment compounds: {len(template['compcompounds'])}")
    print(f"JSON: {output}")


if __name__ == "__main__":
    main()
