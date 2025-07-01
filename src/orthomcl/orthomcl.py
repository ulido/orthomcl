from __future__ import annotations
import argparse
from gzip import GzipFile
import pathlib
from typing import NamedTuple
from collections.abc import Mapping

import requests
from tqdm.auto import tqdm
import platformdirs

ORTHOMCL_ROOT_URL = "https://orthomcl.org/common/downloads/Current_Release/"
ORTHOMCL_GROUPS_URL = (
    ORTHOMCL_ROOT_URL +
    "coreGroups_OrthoMCL-CURRENT/GroupsFile.txt.gz"
)
ORTHOMCL_DEFLINES_URL = (
    ORTHOMCL_ROOT_URL +
    "deflines_OrthoMCL-7.txt.gz"
)


class OrthoEntry(NamedTuple):
    """A single entry in the OrthoMCL database.

    Has three attributes:
     - `group`: The OrthoMCL group name this entry belongs to.
     - `organism`: The entry's / gene's organism.
     - `gene_id`: The entry's / gene's ID name.
     - `description`: The entry's /gene's description.
    """
    group: str
    organism: str
    gene_id: str
    description: str


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


def _download_file(url: str, path: pathlib.Path):
    """Downloads file from `url` and stores it at `path`.
    """
    r = requests.get(url, stream=True)
    total_size = int(r.headers.get("content-length", 0))
    block_size = 16*1024

    path.parent.mkdir(parents=True, exist_ok=True)

    with (
        tqdm(
            desc=f"Downloading {path.stem}",
            total=total_size,
            unit="B",
            unit_scale=True,
        ) as progress_bar,
        path.open("wb") as outfile
    ):
        for chunk in r.iter_content(block_size, decode_unicode=False):
            progress_bar.update(len(chunk))
            outfile.write(chunk)


def _download_groups(path: pathlib.Path):
    """Downloads the OrthoMCL database's current release and stores it in
    `path`.
    """
    _download_file(ORTHOMCL_GROUPS_URL, path)


def _download_deflines(path: pathlib.Path):
    """Downloads the OrthoMCL database's current deflines and stores it in
    `path`.
    """
    _download_file(ORTHOMCL_DEFLINES_URL, path)


def _load_gzipped_file(url: str, rel_path: str):
    """Loads or downloads a cached file and presents it as a line Generator."""
    path = (
        pathlib.Path(platformdirs.user_cache_dir()) /
        rel_path
    )

    if not path.is_file():
        _download_file(url, path)

    with GzipFile(path, "rb") as f:
        for line in f:
            yield line


def _load_groups():
    """Loads the cached OrthoMCL database. If it doesn't exist, downloads it"""
    return _load_gzipped_file(
        ORTHOMCL_GROUPS_URL, "orthomcl/GroupsFile.txt.gz")


def _load_deflines():
    """Loads the cached OrthoMCL deflines file. If it doesn't exist,
    downloads it"""
    return _load_gzipped_file(
        ORTHOMCL_DEFLINES_URL, "orthomcl/deflines.txt.gz")


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
        self._descriptions: dict[str, str] = {}
        self._initialised: bool = False

    def _initialise(self):
        stream = _load_deflines()
        for line in stream:
            org_geneid, _, rest = line.split(maxsplit=2)
            description = rest[::-1].split(maxsplit=1)[1][::-1]
            self._descriptions[org_geneid[1:].decode()] = description.decode()

        stream = _load_groups()
        for line in stream:
            self._add_group_from_line(line.decode())

        self._initialised = True

    def _add_group_from_line(self, line: str):
        name, remainder = line.split(":", maxsplit=1)
        group: list[OrthoEntry] = [
            self._entry_from_string(name, entry_string)
            for entry_string in remainder[1:-1].split(" ")
        ]
        self._add_group(name, group)

    def _add_group(self, name: str, group: list[OrthoEntry]):
        for entry in group:
            self._entries[entry[2]] = entry
        self._groups[name] = group

    def _entry_from_string(self, group: str, string: str):
        organism, gene_id = string.split("|")
        description = self._descriptions[string]
        return OrthoEntry(group, organism, gene_id, description)

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

    def get_entry(self, gene_id: str):
        return self._entries[gene_id]

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
