"""Inspect the public Gaia TAP service for DR4 schemas and table definitions.

Run after the Gaia DR4 archive tables become available:

    python experiments/gaia/dr4_schema_probe.py > gaia-dr4-schema.json

This deliberately discovers the published schema instead of assuming the DR3
table layout carries over unchanged.
"""

import csv
import json
from io import StringIO

import requests


TAP = "https://gea.esac.esa.int/tap-server/tap/sync"


def query(sql):
    response = requests.get(
        TAP,
        params={"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": sql},
        headers={"User-Agent": "sky-survey/0.1 Gaia-DR4-schema-probe"},
        timeout=(10, 120),
    )
    response.raise_for_status()
    reader = csv.DictReader(StringIO(response.text.lstrip("\ufeff")))
    if not reader.fieldnames:
        raise ValueError("Gaia TAP returned no CSV column names.")
    return list(reader)


def run():
    schemas = query("SELECT schema_name, description FROM TAP_SCHEMA.schemas")
    tables = query(
        "SELECT schema_name, table_name, description "
        "FROM TAP_SCHEMA.tables ORDER BY schema_name, table_name"
    )
    dr4_schemas = [row for row in schemas if "dr4" in row["schema_name"].lower()]
    dr4_tables = [
        row for row in tables
        if "dr4" in row["schema_name"].lower() or "dr4" in row["table_name"].lower()
    ]
    columns = {}
    for table in dr4_tables:
        key = (table["schema_name"], table["table_name"])
        columns[f"{key[0]}.{key[1]}"] = query(
            "SELECT column_name, datatype, arraysize, unit, ucd, description "
            "FROM TAP_SCHEMA.columns "
            f"WHERE schema_name='{key[0]}' AND table_name='{key[1]}' "
            "ORDER BY principal DESC, column_name"
        )
    return {
        "archive": TAP,
        "dr4_schemas": dr4_schemas,
        "dr4_tables": dr4_tables,
        "columns": columns,
        "note": "An empty result means DR4 tables were not discoverable at query time; it is not evidence that DR4 lacks those products.",
    }


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
