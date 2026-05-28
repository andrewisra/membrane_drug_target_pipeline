from __future__ import annotations



import argparse

import csv

import sys

import time

from pathlib import Path

from typing import Iterable



import requests



UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"





def fetch_fasta_by_accession(accession: str, timeout: int = 30) -> str:

    url = UNIPROT_FASTA_URL.format(accession=accession)

    response = requests.get(url, timeout=timeout)

    response.raise_for_status()

    return response.text.strip()





def fetch_accession_table(input_csv: str, output_fasta: str, delay: float = 0.2) -> None:

    """Fetch FASTA sequences for accessions listed in a CSV.

    CSV columns:
    accession,label,organism,function
    """

    output = []

    with open(input_csv, newline="", encoding="utf-8") as f:

        reader = csv.DictReader(f)

        for row in reader:

            accession = row["accession"].strip()

            if not accession:

                continue

            fasta = fetch_fasta_by_accession(accession)

            header, *seq_lines = fasta.splitlines()


            meta = (

                f" label={row.get('label','unknown').strip() or 'unknown'}"

                f" | organism={(row.get('organism','unknown').strip() or 'unknown').replace(' ', '_')}"

                f" | function={(row.get('function','unknown').strip() or 'unknown').replace(' ', '_')}"

            )

            output.append(header + meta)

            output.extend(seq_lines)

            time.sleep(delay)

    Path(output_fasta).parent.mkdir(parents=True, exist_ok=True)

    Path(output_fasta).write_text("\n".join(output) + "\n", encoding="utf-8")





def search_uniprot(query: str, output_tsv: str, size: int = 25) -> None:

    """Search UniProtKB and save a small accession table.

    Example query:
    (organism_id:1280) AND (cc_subcellular_location:membrane) AND reviewed:true
    """

    params = {

        "query": query,

        "format": "tsv",

        "fields": "accession,id,protein_name,organism_name,length,cc_subcellular_location",

        "size": size,

    }

    response = requests.get(UNIPROT_SEARCH_URL, params=params, timeout=60)

    response.raise_for_status()

    Path(output_tsv).parent.mkdir(parents=True, exist_ok=True)

    Path(output_tsv).write_text(response.text, encoding="utf-8")





def main() -> None:

    parser = argparse.ArgumentParser(description="Fetch protein sequences or accession tables from UniProt.")

    sub = parser.add_subparsers(dest="command", required=True)



    p_fetch = sub.add_parser("fetch-accessions", help="Fetch FASTA for accessions listed in CSV.")

    p_fetch.add_argument("--input-csv", default="data/uniprot_accessions_example.csv")

    p_fetch.add_argument("--output-fasta", default="data/uniprot_proteins.fasta")



    p_search = sub.add_parser("search", help="Run a UniProtKB search and save TSV results.")

    p_search.add_argument("--query", required=True)

    p_search.add_argument("--output-tsv", default="data/uniprot_search_results.tsv")

    p_search.add_argument("--size", type=int, default=25)



    args = parser.parse_args()

    if args.command == "fetch-accessions":

        fetch_accession_table(args.input_csv, args.output_fasta)

        print(f"Saved FASTA to {args.output_fasta}")

    elif args.command == "search":

        search_uniprot(args.query, args.output_tsv, args.size)

        print(f"Saved search results to {args.output_tsv}")





if __name__ == "__main__":

    main()

