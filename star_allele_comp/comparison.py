"""
Evaluation module for HLA and KIR data.
"""
from __future__ import annotations
import re
from typing import Any, ClassVar
from itertools import chain
from dataclasses import dataclass, field
from collections import defaultdict
from collections.abc import Iterable, Mapping

import pandas as pd


CohortInput = Mapping[str, Iterable[str]]


class AlleleError(BaseException):
    """
    Exception for invalid allele name.
    """


@dataclass
class Allele:
    """
    Represents an allele as an object with methods
    for querying gene, getting resolution,
    and returning the allele with a specified resolution.

    Note that resolution here is defined as the biological meaning resolution.
    For example, HLA has at most 8 digits but only 4 resolutions (fields, indeed),
    while KIR has 3 resolutions.

    Examples:
        >>> HlaAllele("A").resolution                   # not allow
        >>> HlaAllele("A*").resolution                  # 0
        >>> HlaAllele("A*03:05:07").resolution          # 3
        >>> HlaAllele("A*03:05:07").as_resolution(2)    # "A*03:05"
        >>> KirAllele("2DL1*0010305").as_resolution(0)  # "2DL1"
        >>> KirAllele("2DL1*0010305").as_resolution(1)  # "2DL1*001"
        >>> KirAllele("2DL1*0010305").as_resolution(2)  # "2DL1*00103"
    """

    allele: str = ""
    parts: tuple[str, ...] = field(init=False, compare=False)
    max_resolution: ClassVar[int]
    ignore_suffix: bool = False

    def __format__(self, *args: Any) -> str:
        return self.allele.__format__(*args)

    def __repr__(self) -> str:
        return self.allele

    def __str__(self) -> str:
        return self.allele

    def __post_init__(self) -> None:
        try:
            self.parts = self.split()
        except ValueError:
            raise AlleleError(f"{self.allele} is invalid allele name")

    def split(self) -> tuple[str, ...]:
        """
        The method to split the allele string
        depends on the specific nomenclature used.
        After split, the result will save in a tuple `parts`,
        which is structured as follows: (gene, resolution1_str, resolution2_str, ...)
        """
        raise NotImplementedError

    @classmethod
    def join(cls, parts: Iterable[str]) -> Allele:
        """Joining the splited allele string is also depended on nomeclature"""
        raise NotImplementedError

    @property
    def gene(self) -> str:
        """
        Get the gene of allele.

        Examples:
            >>> HlaAllele("A*03:05").gene
            "A"
        """
        return self.parts[0]

    @property
    def resolution(self) -> int:
        """
        Get the resolution of allele.

        Examples:
            >>> HlaAllele("A*03:05").resolution
            2
        """
        return len(self.parts) - 1

    def trim_resolution(self, resolution: int) -> Allele:
        """
        Trim the resolution of allele.

        Examples:
        >>> HlaAllele("A*03:05:07").trim_resolution(2)    # "A*03:05"
        >>> HlaAllele("A*03:05").trim_resolution(3)       # "A*03:05"
        """
        return self.join(self.parts[: resolution + 1])

    def as_resolution(self, resolution: int) -> Allele | None:
        """
        Trim the resolution of allele.

        Return None if the resolution is larger than the max_resolution.

        Examples:
        >>> HlaAllele("A*03:05:07").as_resolution(2)    # "A*03:05"
        >>> HlaAllele("A*03:05").as_resolution(3)       # None
        """
        if resolution > self.resolution:
            return None
        return self.trim_resolution(resolution)


class HlaAllele(Allele):
    """Define HLA allele"""

    max_resolution: ClassVar[int] = 4

    def split(self) -> tuple[str, ...]:
        """
        The allele should match the regex `(\\w+)\\*([\\w:]*)`.
        And first four fields of the resolution will retained.
        """
        allele_reg = re.match(r"(\w+)\*([\w:]*)", self.allele)
        if not allele_reg:
            raise ValueError
        gene_str, allele_str = allele_reg.groups()
        if not allele_str:
            return (gene_str,)
        part_allele = allele_str.split(":")
        if not all(part_allele):
            raise ValueError
        part_allele = part_allele[:4]
        if (
            self.ignore_suffix
            and part_allele
            and not re.match(r"^\d+$", part_allele[-1])
        ):
            part_allele[-1] = re.findall(r"^\d+", part_allele[-1])[0]
        return (gene_str, *part_allele)

    @classmethod
    def join(cls, parts: Iterable[str]) -> Allele:
        iter_p = iter(parts)
        return cls(next(iter_p) + "*" + ":".join(iter_p))


class KirAllele(Allele):
    """Define KIR allele"""

    max_resolution: ClassVar[int] = 3

    def split(self) -> tuple[str, ...]:
        """
        The allele should match regex `(\\w+)\\*(\\d*)(\\w*)`

        Don't directly append string after the allele
        e.g.
        * `2DL1*00102N` OK -> 001 02n
        * `2DL1*0010203test` not OK -> 001 02 03test
        * `2Dl1*001 this_is_ok` OK -> 001
        """
        allele_reg = re.match(r"(\w+)\*(\d*)(\w*)", self.allele)
        if not allele_reg:
            raise ValueError
        gene_str, allele_str, suffix_str = allele_reg.groups()
        parts = [gene_str]
        for digit in [3, 2, 2]:
            if len(allele_str) >= digit:
                parts.append(allele_str[:digit])
                allele_str = allele_str[digit:]
            elif allele_str:  # 0 < len(allele_str) < digit
                raise ValueError
        if not self.ignore_suffix and suffix_str:
            if len(parts) > 1:
                parts[-1] += suffix_str
            else:
                parts.append(suffix_str)
        return tuple(parts)

    @classmethod
    def join(cls, parts: Iterable[str]) -> Allele:
        iter_p = iter(parts)
        return cls(next(iter_p) + "*" + "".join(iter_p))


@dataclass
class MatchResult:
    """
    A format that saved the comparison of the alleles

    Examples:
       ```
       allele1 = A*03:05:01
       allele2 = A*03:05:02
       match_res = 2
       allele1_res = 3
       gene = "A"
       ```
    """

    allele1: Allele | None = None
    allele2: Allele | None = None
    match_res: int = -1
    match_str: str = ""
    allele1_res: int = field(init=False)
    gene: str = field(init=False)

    def __post_init__(self) -> None:
        if self.allele1:
            self.allele1_res = self.allele1.resolution
            self.gene = self.allele1.gene
        else:  # FP
            self.allele1_res = -1
            assert self.allele2
            self.gene = self.allele2.gene

        # default is -1
        # if self.match_str == "FP" or self.match_str == "FP":
        #     self.match_res  = -1

    def __str__(self) -> str:
        if len(self.match_str) == 1:
            match_style = f"={self.match_str}="
        elif self.match_str == "FP":
            match_style = "FP>"
        elif self.match_str == "FN":
            match_style = "<FN"

        # print string if the type is None, use empty string
        allele1 = self.allele1 or ""
        allele2 = self.allele2 or ""
        return f"{allele1:16s} {match_style} {allele2:16s}"


class CohortResult(dict[str, list[MatchResult]]):
    """
    A class that stores the comparison result of cohort,
    which is indeed `dict[str, list[MatchResult]]`

    The class provides methods to transform the result into dataframe (`to_dataframe`)
    and print it out (`print(result)`)
    """

    def to_dataframe(self) -> pd.DataFrame:
        """
        Transform comparison result of cohort into flatten dataframe.

        Returns:
        ```
                allele1        allele2  match_res match_str  allele1_res   gene   id   method
        0      2DL1*00101     2DL1*00101          2         2            2   2DL1  id1  method1
        1      2DL1*00101     2DL1*00101          2         2            2   2DL1  id1  method1
        2        2DL3*001       2DL3*001          1         1            1   2DL3  id1  method1
        3    2DL4*0030201   2DL4*0030201          3         3            3   2DL4  id1  method1
        ```
        """
        list_df = []
        for sample_id, results in self.items():
            cohort_df = pd.DataFrame(results)
            cohort_df["id"] = sample_id
            cohort_df["allele1"] = [
                i["allele"] if i else None for i in cohort_df["allele1"]
            ]
            cohort_df["allele2"] = [
                i["allele"] if i else None for i in cohort_df["allele2"]
            ]
            list_df.append(cohort_df)
        return pd.concat(list_df)

    def __str__(self) -> str:
        """
        Print result allele by allele

        The input can be from both compare_cohort or compare_method.

        Returns:
        ```
        Sample id1
        A*01:02:03       =3= A*01:02:03
        A*01:01          =2= A*01:01
        ...
        ```
        """
        text = ""
        for sample_id, results in self.items():
            text += f"Sample {sample_id}\n"
            for result in results:
                text += str(result) + "\n"
        return text


class MethodResult(dict[str, CohortResult]):
    """
    A class that stores the comparison result of different methods,
    the structure is `dict[str, CohortResult]`.
    """

    def to_dataframe(self) -> pd.DataFrame:
        """
        Same to CohortResult.to_dataframe
        with an additional `method` column.
        """
        list_df = []
        for method, results in self.items():
            if isinstance(results, dict):
                method_df = results.to_dataframe()
                method_df["method"] = method
                list_df.append(method_df)
        return pd.concat(list_df)

    def __str__(self) -> str:
        """
        Same to CohortResult.__str__(),
        but add `method` indicator.
        """
        text = ""
        for method, results in self.items():
            text += f"Method {method}\n"
            text += str(results)
        return text


def str_to_allele(
    alleles: Iterable[str], allele_type: str, ignore_suffix: bool = False
) -> list[Allele]:
    """Transform list of string to list of Allele object"""
    if allele_type.lower() == "hla":
        return [HlaAllele(allele, ignore_suffix=ignore_suffix) for allele in alleles]
    if allele_type.lower() == "kir":
        return [KirAllele(allele, ignore_suffix=ignore_suffix) for allele in alleles]
    raise ValueError


def compare_gene(
    alleles1: list[Allele], alleles2: list[Allele]
) -> Iterable[MatchResult]:
    """
    Compare two sets of alleles.

    The function attempts to find the best possible pairwise matches
    between the alleles in alleles1 and alleles2.

    Args:
        alleles1: A list of alleles, treated as reference.
        alleles2: Another list of alleles.

    Returns:
        An iterable of MatchResult objects representing the pairwise matches.

    Note:
        All alleles in alleles1 and alleles2 should have the same gene name.
    """
    alleles1 = alleles1.copy()
    alleles2 = alleles2.copy()
    max_resolution = list(chain.from_iterable([alleles1, alleles2]))[0].max_resolution
    # match highest resolution first and then the 2nd highest resolution and ... to gene-level(resolution=0)
    for resolution in range(max_resolution, -1, -1):
        for allele1 in alleles1.copy():
            for allele2 in alleles2.copy():
                allele1_limit = allele1.as_resolution(resolution)
                allele2_limit = allele2.as_resolution(resolution)
                if (
                    allele1_limit and allele1_limit == allele2_limit
                ):  # not None  i.e. occurs in resolution > allele.resolution
                    alleles1.remove(allele1)
                    alleles2.remove(allele2)
                    yield MatchResult(
                        allele1=allele1,
                        allele2=allele2,
                        match_res=resolution,
                        match_str=str(resolution),
                    )
                    break

    # No allele1 can match to allele2 = FN
    for allele1 in alleles1:
        yield MatchResult(
            allele1=allele1,
            match_str="FN",
        )
    # No allele2 can match to allele1 = FP
    for allele2 in alleles2:
        yield MatchResult(
            allele2=allele2,
            match_str="FP",
        )


def group_by_gene(alleles: Iterable[Allele]) -> dict[str, list[Allele]]:
    """
    Group alleles by their gene name

    Returns: A dict with key=gene's name, value=list of alleles
    """
    gene_group = defaultdict(list)
    for allele in alleles:
        gene_group[allele.gene].append(allele)
    return gene_group


def compare_sample(
    sample1: Iterable[Allele], sample2: Iterable[Allele]
) -> Iterable[MatchResult]:
    """
    Compare two set of alleles from two samples

    Args:
        alleles1: A list of alleles, treated as reference.
        alleles2: Another list of alleles.

    """
    sample1_gene_group = group_by_gene(sample1)
    sample2_gene_group = group_by_gene(sample2)
    genes = sample1_gene_group.keys() | sample2_gene_group.keys()
    for gene in sorted(genes):
        # print("Gene", gene)
        result = compare_gene(
            sample1_gene_group.get(gene, []), sample2_gene_group.get(gene, [])
        )
        # for r in result: print(r)
        yield from result


def compare_cohort(
    cohort1: CohortInput,
    cohort2: CohortInput,
    allele_type: str,
    ignore_suffix: bool = False,
) -> CohortResult:
    """
    Compare two cohorts

    Args:
        cohort1: The first cohort, a dictionary where key = sample IDs and values = lists of alleles for the sample.
        cohort2: The second cohort, same format as cohort1.
        allele_type: The type of allele in the cohorts, either "kir" or "hla".
        ignore_suffix: Whether to ignore suffix in the alleles.

    Returns:
        A dict where keys = sample IDs and values = lists of MatchResult objects.

    Note:
        * The sample IDs in cohort1 and cohort2 shall be the same.
        * cohort1 is treated as reference.

    """
    ids = cohort1.keys() | cohort2.keys()
    results_cohort = {}
    for sample_id in sorted(ids):
        results = compare_sample(
            str_to_allele(
                cohort1.get(sample_id, []), allele_type, ignore_suffix=ignore_suffix
            ),
            str_to_allele(
                cohort2.get(sample_id, []), allele_type, ignore_suffix=ignore_suffix
            ),
        )
        results_cohort[sample_id] = list(results)
    return CohortResult(results_cohort)


def compare_method(
    method_cohort_dict: Mapping[str, CohortInput],
    ref_method: str,
    *args: Any,
    **kwargs: Any,
) -> MethodResult:
    """
    Compare cohorts of alleles generated from different methods.

    The `ref_method` cohort is treated as the reference for comparison.

    Args:
        method_cohort_dict: A dictionary where keys = method names and values = Cohort (see compare_cohort).
        ref_method: The name of the reference method.
        *args, **kwargs: Set allele_type, ignore_suffix, See compare_cohort.

    Returns:
        A dictionary where keys = method names and values = CohortResult objects (see compare_cohort).

    Note:
        - The sample IDs in each cohort should be the same.
    """
    if ref_method not in method_cohort_dict:
        raise ValueError(
            f"Reference method {ref_method} not found. Possible methods: {method_cohort_dict.keys()}"
        )

    method_result = {}
    for method in method_cohort_dict.keys():
        method_result[method] = compare_cohort(
            method_cohort_dict[ref_method],
            method_cohort_dict[method],
            *args,
            **kwargs,
        )
    return MethodResult(method_result)
