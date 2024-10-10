#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
import textwrap
from pathlib import Path
from collections import Counter
from typing import Iterator, Dict, List
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "HaidingWANG"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Haiding WANG"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Haiding WANG"
__email__ = "wanghaidingfr@126.com"
__status__ = "Developpement"



def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=Path,
                        default=Path("OTU.fasta"), help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file: Path, minseqlen: int) -> Iterator[str]:
    """Read a compressed fasta and extract all fasta sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :return: A generator object that provides the Fasta sequences (str).
    """
    with gzip.open(amplicon_file, 'rt') as monfich:
        sequence = ""
        for line in monfich:
            line = line.strip()
            if line.startswith('>'):
                if len(sequence) >= minseqlen:
                    yield sequence
                sequence = ""
            else :
                sequence += line
        if sequence and len(sequence) >= minseqlen:  # Vérifier et yield la dernière séquence
            yield sequence



def dereplication_fulllength(amplicon_file: Path, minseqlen: int, mincount: int) -> Iterator[List]:
    """Dereplicate the set of sequence

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :param mincount: (int) Minimum amplicon count
    :return: A generator object that provides a (list)[sequences, count] of sequence with a count >= mincount and a length >= minseqlen.
    """
    sequence_count = {}
    for sequence in read_fasta(amplicon_file, minseqlen):
        sequence_count[sequence] = sequence_count.get(sequence, 0) + 1

    for seq, count in sorted(sequence_count.items(), key=lambda x: x[1], reverse=True):
        if count >= mincount:
            yield [seq, count]


def get_identity(alignment_list: List[str]) -> float:
    """Compute the identity rate between two sequences

    :param alignment_list:  (list) A list of aligned sequences in the format ["SE-QUENCE1", "SE-QUENCE2"]
    :return: (float) The rate of identity between the two sequences.
    """
    matches = sum(1 for seq1, seq2 in zip(alignment_list[0], alignment_list[1]) if seq1 == seq2)
    return (matches / len(alignment_list[0])) * 100

def abundance_greedy_clustering(amplicon_file: Path, minseqlen: int, mincount: int, chunk_size: int, kmer_size: int) -> List:
    """Compute an abundance greedy clustering regarding sequence count and identity.
    Identify OTU sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length.
    :param mincount: (int) Minimum amplicon count.
    :param chunk_size: (int) A fournir mais non utilise cette annee
    :param kmer_size: (int) A fournir mais non utilise cette annee
    :return: (list) A list of all the [OTU (str), count (int)] .
    """
    OTU_bank = []
    
    for sequence, count in dereplication_fulllength(amplicon_file, minseqlen, mincount):
        is_OTU = True  

        # Comparer la séquence actuelle avec toutes les OTUs existantes
        for otu_seq, _ in OTU_bank:
            alignment = nw.global_align(sequence, otu_seq, gap_open=-1, gap_extend=-1, matrix=str(Path(__file__).parent / "MATCH"))
            identity = get_identity(alignment)  # Calcule le pourcentage d'identité
            
            # Si la similarité est supérieure ou égale à 97%, elle n'est pas considérée comme une nouvelle OTU
            if identity >= 97:
                is_OTU = False
                break
        
        if is_OTU:
            OTU_bank.append([sequence, count])  
    
    return OTU_bank




def write_OTU(OTU_list: List, output_file: Path) -> None:
    """Write the OTU sequence in fasta format.

    :param OTU_list: (list) A list of OTU sequences
    :param output_file: (Path) Path to the output file
    """
    with open(output_file, "w") as file:
        for i, (seq, count) in enumerate(OTU_list):
            nb = i + 1
            file.write(f">OTU_{nb} occurrence:{count}\n")
            formatted_OTU_seq = textwrap.fill(seq, width=80)
            file.write(f"{formatted_OTU_seq}\n")
    


#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Votre programme ici
    OTU = abundance_greedy_clustering(args.amplicon_file, args.minseqlen, args.mincount, 0, 0)
    write_OTU(OTU, args.output_file)



if __name__ == '__main__':
    main()
