#!/usr/bin/env python

import os
import sys
import logging
from argparse import ArgumentParser
from typing import Optional, Any

from Bio import SeqIO

args_parser = ArgumentParser(description="combine V-gene, J-gene and CDR3 into one sequence")
args_parser.add_argument("tcrv_gene_name", help="full name of the TCR V-gene (alpha/beta), insert empty string to omit")
args_parser.add_argument("cdr3_sequence", help="sequence of the TCR CDR3 loop (alpha/beta)")
args_parser.add_argument("tcrj_gene_name", help="full name of the TCR J-gene (alpha/beta), insert empty string to omit")
args_parser.add_argument("--debug", "-d", action='store_const', const=True, default=False, help="prints debug information")
args_parser.add_argument("--no-gaps", "-n", action='store_const', const=True, default=False, help="disable IMGT gaps")


_log = logging.getLogger(__name__)


class MisMatch(Exception):
    pass


def get_gene_sequences(gene_name: str, gapped: Optional[bool] = True) -> str:
    """
    Retrieve the sequence of any gene (V/J & alpha/beta) from the fasta files.
    """

    if gene_name.startswith("TRAV"):
        if gapped:
            path = "trav-gaps.fa"
        else:
            path = "trav-gaps.fa"

    elif gene_name.startswith("TRBV"):
        if gapped:
            path = "trbv-gaps.fa"
        else:
            path = "trbv.fa"

    elif gene_name.startswith("TRAJ"):
        path = "traj.fa"

    elif gene_name.startswith("TRBJ"):
        path = "trbj.fa"

    else:
        raise TypeError(f"unknown gene type: {gene_name}")

    path = os.path.join(os.path.dirname(__file__), 'data', path)

    matching = {}
    for record in SeqIO.parse(path, "fasta"):
        record_gene_name = record.id.split('|')[1]
        if record_gene_name == gene_name:
            # exact match
            matching[record_gene_name] = str(record.seq)
            break

        elif record_gene_name.startswith(f"{gene_name}-"):
            matching[record_gene_name] = str(record.seq)

        elif record_gene_name.startswith(f"{gene_name}*"):
            matching[record_gene_name] = str(record.seq)

    return matching


def combine_at_cterm(cdr3_seq: str, j_seq: str) -> str:
    """
    Attempts to combine a J-gene sequence with the C-terminal part of a CDR3 sequence
    and returns the combination.

    The CDR3 sequence is dominant here, meaning that it overrides the J-gene sequence.
    """

    # We must match as much of CDR3 as possible, so try the longest part first.
    # Try multiple offsets in the CDR3 sequence, starting at the most C-terminal part of CDR3
    solutions = {}
    for length in [4,3,2]:
        for offset in [0,1]:

            cterm_seq = cdr3_seq[-(length + offset): len(cdr3_seq) - offset]
            i = j_seq.find(cterm_seq)
            if i == -1:
                continue

            # amount to take off J-sequence
            ndel = i + length + offset

            _log.debug(f"combine at cterm: {cdr3_seq}{'.' * (len(j_seq) - ndel)}")
            _log.debug(f"combine at cterm: {'.' * (len(cdr3_seq) - ndel)}{j_seq}")

            solutions[i] = cdr3_seq + j_seq[ndel:]

    if len(solutions) == 0:
        raise MisMatch(f"cannot match CDR3 {cdr3_seq} with {j_seq}")

    # the most N-terminal solution must be chosen
    return solutions[min(solutions.keys())]


def combine_at_nterm(v_seq: str, cdr3_seq: str) -> str:
    """
    Attempts to combine a V-gene sequence with the N-terminal part of a CDR3 sequence
    and returns the combination.

    The CDR3 sequence is dominant here, meaning that it overrides the V-gene sequence.
    """

    # We must match as much of CDR3 as possible, so try the longest part first.
    # Try multiple offsets in the CDR3 sequence, starting at the most N-terminal part of CDR3
    solutions = {}
    for length in [4,3,2]:
        for offset in [0,1]:

            nterm = cdr3_seq[offset:(length + offset)]
            i = v_seq.rfind(nterm)
            if i == -1:
                continue

            # amount to take off V-sequence
            ndel = len(v_seq) - i + offset

            _log.debug(f"combine at nterm: {'.' * (len(v_seq) - ndel)}{cdr3_seq}")
            _log.debug(f"combine at nterm: {v_seq}{'.' * (len(cdr3_seq) - ndel)}")

            solutions[i] = v_seq[:(i - offset)] + cdr3_seq

    if len(solutions) == 0:
        raise MisMatch(f"cannot match {v_seq} with CDR3 {cdr3_seq}")

    # the most C-terminal match must be chosen
    return solutions[max(solutions.keys())]


def is_valid_name(input_: Any) -> bool:
    """
    Checks whether the name is omitted by the user.
    """

    if type(input_) == str:

        return len(input_.strip()) > 0
    else:
        return False


def get_full_tcr_sequence(
    tcrv_gene_name: Any,
    cdr3_sequence: str,
    tcrj_gene_name: Any,
    gapped: Optional[bool] = True,
) -> str:
    """
    Takes the full genes name as input. Example: TRAV1-1*01
    """

    # start from cdr3
    full_sequence = cdr3_sequence

    if is_valid_name(tcrv_gene_name):

        v_sequences = get_gene_sequences(tcrv_gene_name, gapped)
        if len(v_sequences) == 0:
            raise ValueError(f"no sequences found for {tcrv_gene_name}")

        for v_name, v_sequence in v_sequences.items():
            try:
                full_sequence = combine_at_nterm(v_sequence, full_sequence)
                break
            except MisMatch:
                continue
        else:
            choice_s = "\n".join([f"{n}: {s}" for n, s in v_sequences.items()])
            raise MisMatch(f"cannot match any V-sequence with CDR3 {cdr3_sequence}:\n{choice_s}")

    if is_valid_name(tcrj_gene_name):

        j_sequences = get_gene_sequences(tcrj_gene_name, gapped)
        if len(j_sequences) == 0:
            raise ValueError(f"no sequences found for {tcrj_gene_name}")

        for j_name, j_sequence in j_sequences.items():
            try:
                full_sequence = combine_at_cterm(full_sequence, j_sequence)
                break
            except MisMatch:
                continue
        else:
            choice_s = "\n".join([f"{n}: {s}" for n, s in j_sequences.items()])
            raise MisMatch(f"cannot match any J-sequence with CDR3 {cdr3_sequence}:\n{choice_s}")

    return full_sequence


if __name__ == "__main__":

    args = args_parser.parse_args()
    logging.basicConfig(stream=sys.stdout, level=logging.DEBUG if args.debug else logging.INFO)

    full_sequence = get_full_tcr_sequence(
        args.tcrv_gene_name,
        args.cdr3_sequence,
        args.tcrj_gene_name,
        not args.no_gaps,
    )

    print(full_sequence)
