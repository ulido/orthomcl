from __future__ import annotations
import argparse
from gzip import GzipFile
import pathlib
from typing import NamedTuple
from collections.abc import Mapping

import requests
from tqdm.auto import tqdm
import platformdirs

ORTHOMCL_GROUPS_URL = "https://orthomcl.org/common/downloads/Current_Release/coreGroups_OrthoMCL-CURRENT/GroupsFile.txt.gz"  # noqa: E501


class OrthoEntry(NamedTuple):
    """A single entry in the OrthoMCL database.

    Has three attributes:
     - `group`: The OrthoMCL group name this entry belongs to.
     - `organism`: The entry's / gene's organism.
     - `gene_id`: The entry's / gene's ID name.
    """
    group: str
    organism: str
    gene_id: str


class OrthoEntryCollection(Mapping[str, list[OrthoEntry]]):
    """Collection of OrthoMCL database entries, queryable by organism.
    """
    def __init__(self, entries: list[OrthoEntry]):
        self._entries = entries
        self._organisms = {
            entry[1] for entry in self._entries
        }

    def __getitem__(self, organism: str):
        return [
            OrthoEntry(*entry) for entry in self._entries
            if entry[1] == organism
        ]

    def __contains__(self, organism: str) -> bool:
        return organism in self._organisms

    def __iter__(self):
        return iter(self._organisms)

    def __len__(self):
        return len(self._organisms)


def _download_groups(path: pathlib.Path):
    """Downloads the OrthoMCL database's current release and stores it in
    `path`.
    """
    r = requests.get(ORTHOMCL_GROUPS_URL, stream=True)
    total_size = int(r.headers.get("content-length", 0))
    block_size = 16*1024

    path.parent.mkdir(parents=True, exist_ok=True)

    with (
        tqdm(
            desc="Downloading OrthoMCL groups",
            total=total_size,
            unit="B",
            unit_scale=True,
        ) as progress_bar,
        path.open("wb") as outfile
    ):
        for chunk in r.iter_content(block_size, decode_unicode=False):
            progress_bar.update(len(chunk))
            outfile.write(chunk)


def _load_groups():
    """Loads the cached OrthoMCL database. If it doesn't exist, downloads it"""
    path = (
        pathlib.Path(platformdirs.user_cache_dir()) /
        "orthomcl/GroupsFile.txt.gz"
    )

    if not path.is_file():
        _download_groups(path)

    with GzipFile(path, "rb") as f:
        for line in f:
            yield line


class _OrthoMCL(Mapping[str, OrthoEntry]):
    """Queryable OrthoMCL database.

    Querying the database uses a simple key-value interface. This illustrative
    example queries OrthoMCL for the Trypanosoma brucei orthologs of the
    Leishmania mexicana PF16 gene (ID LmxM.20.1400):
    ```python-repl
    >>> gene_id = "LmxM.20.1400"
    >>> results = OrthoMCL[gene_id]
    >>> results["tbre"]
    [OrthoEntry(group='OG7_0009222', organism='tbre', gene_id='Tb427_010020200')]
    ```

    The database can also be queried by OrthoMCL group name:
    ```python-repl
    >>> results = OrthoMCL.by_group("OG7_0009222")
    >>> results["lmex"]
    [OrthoEntry(group='OG7_0009222', organism='lmex', gene_id='LmxM.20.1400')]
    ```
    """
    def __init__(self):
        self._groups: dict[str, list[OrthoEntry]] = {}
        self._entries: dict[str, OrthoEntry] = {}
        self._initialised: bool = False

    def _initialise(self):
        stream = _load_groups()

        for line in stream:
            self._add_group_from_line(line.decode())

        self._initialised = True

    def _add_group_from_line(self, line: str):
        name, remainder = line.split(":", maxsplit=1)
        group: list[OrthoEntry] = [
            self._entry_from_string(name, entry_string)
            for entry_string in remainder[1:].split(" ")
        ]
        self._add_group(name, group)

    def _add_group(self, name: str, group: list[OrthoEntry]):
        for entry in group:
            self._entries[entry[2]] = entry
        self._groups[name] = group

    @staticmethod
    def _entry_from_string(group: str, string: str):
        organism, gene_id = string.split("|")
        return OrthoEntry(group, organism, gene_id)

    def __getitem__(self, gene_id: str):
        if not self._initialised:
            self._initialise()
        try:
            entry = self._entries[gene_id]
            return OrthoEntryCollection(self._groups[entry[0]])
        except KeyError:
            # If we have no knowledge of the gene, just return an empty list.
            return OrthoEntryCollection([])

    def __iter__(self):
        if not self._initialised:
            self._initialise()
        return iter(self._entries)

    def __len__(self):
        if not self._initialised:
            self._initialise()
        return len(self._entries)

    def by_group(self, group_name: str):
        if not self._initialised:
            self._initialise()
        return OrthoEntryCollection(self._groups[group_name])

    def all_organisms(self):
        if not self._initialised:
            self._initialise()
        return list({
            entry.organism for entry in self._entries.values()
        })
OrthoMCL = _OrthoMCL()  # noqa: E305


def cli():
    parser = argparse.ArgumentParser(
        description=(
            "Find VEuPathDB gene orthologues via the OrthoMCL database."),
    )
    parser.add_argument(
        "gene_id",
        nargs="+",
        type=str,
        help="The gene IDs for which the OrthoMCL database should be queried.",
    )
    parser.add_argument(
        "--organisms",
        "-o",
        type=str,
        help=(
            "Only display orthologues from the given (comma-separated) list "
            "of OrthoMCL organism IDs."
        ),
        default=None,
    )

    args = parser.parse_args()

    if args.organisms is not None:
        organisms: list[str] = [
            organism.strip() for organism in args.organisms.split(",")]
    else:
        organisms = OrthoMCL.all_organisms()

    for geneid in args.gene_id:
        orthologs = OrthoMCL[geneid]
        print(f"{geneid}:")
        for organism in organisms:
            if organism in orthologs:
                geneids = ','.join(
                    gene.gene_id for gene in orthologs[organism]
                    if gene.gene_id != geneid
                )
                print(f"  {organism}: {geneids}")
