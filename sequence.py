from Bio import Entrez
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices


class Sequence(object):
    _sequence_pairs = []
    gap = -10
    match = 2
    mismatch = -4

    def __init__(self, database, file_format="gb", email="") -> None:
        """Constructor"""
        Entrez.email = email
        self.database = database
        self.file_format = file_format
        print(f"gap: {self.gap}, match: {
              self.match}, mismatch: {self.mismatch}\n")

    def add_sequence(self, seqs: list) -> None:
        """
            Agrega una nueva secuencia a la lista de secuencias

            @return: None
        """
        print(f"Getting sequences {seqs}")
        data = self._get_data(accession_number=seqs)
        self._sequence_pairs.append(data)

    def export(self) -> None:
        """Exporta todas las secuencias"""
        for pair in self._sequence_pairs:
            for seq in pair:
                with open(f"sequences/{seq.id}.{self.file_format}", "w") as output_handle:
                    SeqIO.write(seq, output_handle, self.file_format)

    def do_alignment(self):
        for pair in self._sequence_pairs:
            print(f"{pair[0].description=}\n{pair[1].description=}")
            self._aligner(pair[0], pair[1], "global")
            self._aligner(pair[0], pair[1], "local")
            print()

    def do_alignment_matrices(self):
        for pair in self._sequence_pairs:
            print(f"{pair[0].description=}\n{pair[1].description=}")
            self._alignment_sustitution(pair[0], pair[1], "BLOSUM62")
            self._alignment_sustitution(pair[0], pair[1], "PAM250")
            self._aligner(pair[0], pair[1], "global")
            self._aligner(pair[0], pair[1], "local")
            print()

    def print_sequences(self) -> None:
        """Imprime todas las secuencias agregadas"""
        for pair in self._sequence_pairs:
            for sequence in pair:
                print(sequence)

    def _get_data(self, accession_number):
        """
            Obtiene la secuencia desde GenBank

            @param database: base de datos de la secuencia
            @param accession_number: numero de acceso de la secuencia

            @return: list con la secuencia o None si no se encuentra
        """
        with Entrez.efetch(db=self.database, id=accession_number,
                           rettype=self.file_format) as handle:
            return [handle for handle in SeqIO.parse(handle, self.file_format)]

    def _alignment_sustitution(self, seq1, seq2, type):
        aligner = PairwiseAligner()
        aligner.substitution_matrix = substitution_matrices.load(type)
        alignment = aligner.align(seq1.seq, seq2.seq)
        self._get_resume(alignment, seq1.id, seq2.id, type)

    def _aligner(self, seq1, seq2, type):
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
        smax = max([*alignment[0].coordinates[0],
                   *alignment[0].coordinates[1]])
        identity = (st/smax)*100
        print(
            f"Sequences: {seq1, seq2} {type=} score: {alignment.score} {identity=}%")
