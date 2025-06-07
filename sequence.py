import re
from Bio import Entrez
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
from pathlib import Path
from collections import deque
from itertools import chain

class Sequence(object):
    _sequence_pairs = []
    gap = -10
    match = 2
    mismatch = -4

    def __init__(self, database, file:str, /, file_format="gb", email="", export=True) -> None:
        """Constructor"""
        Entrez.email = email
        self.database = database
        self.file_format = file_format

        self._add_sequences(file)

        if export:
            self._export()

        match database:
            case "nucleotide":
                self._do_nucleotid_alignment_matrices()
            case "protein":
                self._do_protein_alignment()
            case _:
                return

    def __str__(self):
        """ Imprime los valores de los scores """
        return f"gap: {self.gap}, match: {self.match}, mismatch: {self.mismatch}"

    def _add_sequences(self, file:str) -> None:
        """Agrega una nueva secuencia a la lista de secuencias buscándolas desde la BD

        Args:
            file (str): Nombre del archivo en donde están las secuencias
        """
        with open(file) as f:
            for line in f.readlines():
                sequences = re.split(r'\s+', line)[:2]
                print(f"Getting sequences {sequences}")

                with Entrez.efetch(db=self.database, id=sequences, rettype=self.file_format) as handle:
                    data = [handle for handle in SeqIO.parse(
                        handle, self.file_format)]
                    self._sequence_pairs.append(deque(data))

    def _export(self) -> None:
        """Exporta todas las secuencias"""
        folder = Path("sequences")
        folder.mkdir(exist_ok=True)
        
        for seq in list(chain(*self._sequence_pairs)):
            filepath = folder.joinpath(f"{seq.id}.{self.file_format}")
            with open(filepath, "w") as output_handle:
                SeqIO.write(seq, output_handle, self.file_format)

    def _do_protein_alignment(self) -> None:
        """ Realiza el alineamiento de las secuencias para proteinas """
        self._process_pairs(self._aligner, ["global", "local"])

    def _do_nucleotid_alignment_matrices(self) -> None:
        """ Realiza el alineamiento de las secuencias para nucleotidos utilizando BLOSUM62 y PAM250 """
        self._process_pairs(self._alignment_substitution,["BLOSUM62", "PAM250"])
        self._process_pairs(self._aligner, ["global", "local"])

    def _process_pairs(self, method, modes):
        """ Procesa los pares """
        for pair in self._sequence_pairs:
            print(f"{pair[0].description=}\n{pair[1].description=}")
            for mode in modes:
                method(pair[0], pair[1], mode)

    def print_sequences(self) -> None:
        """Imprime todas las secuencias agregadas"""
        for pair in self._sequence_pairs:
            for sequence in pair:
                print(sequence)

    def _alignment_substitution(self, seq1, seq2, type):
        """ Implementación PairwiseAligner para nucleotidos"""
        aligner = PairwiseAligner()
        aligner.substitution_matrix = substitution_matrices.load(type)
        alignment = aligner.align(seq1.seq, seq2.seq)
        self._get_resume(alignment, seq1.id, seq2.id, type)

    def _aligner(self, seq1, seq2, type:str):
        """ Implementación"""
        aligner = PairwiseAligner()
        aligner.mode = type
        aligner.match_score = self.match
        aligner.mismatch_score = self.mismatch
        aligner.open_gap_score = self.gap
        aligner.extend_gap_score = self.gap
        alignment = aligner.align(seq1.seq, seq2.seq)
        self._get_resume(alignment, seq1.id, seq2.id, type)

    def _get_resume(self, alignment, seq1, seq2, type) -> tuple:
        st = str(alignment[0]).count("|")
        smax = max([*alignment[0].coordinates[0],*alignment[0].coordinates[1]])
        identity = (st/smax)*100
        print(f"Sequences: {seq1, seq2} {type=}, score: {alignment.score} {identity=}%")
