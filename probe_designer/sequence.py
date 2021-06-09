from .thermodynamics import calc_gc, calc_tm, calc_deltaG, ThermoParams
from .utils import complement, run_blast, FastaReader

REF = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta"


class Sequence:
    thermo = ThermoParams()

    def __init__(self, sequence, thermo=thermo, reference=REF):
        """Class to characterize a sequence of DNA

        Parameters
        ----------
        sequence : str
            Sequence of ACTG bases
        thermo : ThermoParams object, optional
            Thermodynamic parameters of hybridization (see ThermoParams class)
        reference : path, optional
            Reference file used for running BLAST
        """
        self.validate(sequence)
        self.sequence = sequence
        self.thermo = thermo
        self.reference = reference

    def __len__(self):
        return len(self.sequence)

    def __str__(self):
        return self.sequence

    def __repr__(self):
        args = ", ".join(f"{k}={v}" for k, v in self.__dict__.items())
        return f"{self.__class__.__name__}({args})"

    def validate(self, sequence):
        if len(sequence) == 0:
            raise ValueError("Sequence length must be greater than 0")
        else:
            accepted_bases = {"A", "T", "G", "C"}
            for base in sequence:
                if base not in accepted_bases:
                    raise ValueError("Sequence must contain only ACTG")

    @property
    def gc(self):
        return calc_gc(self.sequence)

    @property
    def tm(self):
        t = self.thermo
        return calc_tm(
            self.sequence, dnac1=t.dnac1, dnac2=t.dnac2, Na=t.Na, Mg=t.Mg
        )

    @property
    def dg(self):
        t = self.thermo
        return calc_deltaG(self.sequence, temp=t.temp, paraG=t.paraG)

    def blast(self):
        blast_df = run_blast(
            self.sequence,
            self.reference,
            qcov_hsp_perc=40,
            perc_identity=90,
            threads=1,
        )
        return blast_df.to_dict("records")


class Probe(Sequence):
    thermo = ThermoParams()

    def __init__(
        self,
        chrom,
        alt_position,
        length,
        ref_base,
        alt_base,
        reference,
        thermo=thermo,
        design_ref=False,
        same_strand=False,
    ):
        """Class for creating a probe. Inherits thermodynamic attributes from
        `Sequence` class.

        Parameters
        ----------
        chrom : str
            Chromosome
        alt_position : int
            Positon of ALT base (center of probe)
        length : int
            Length of probe in bp
        ref_base : str
            REF base at alt_position
        alt_base : str
            ALT base at alt_position
        reference : path or FastaReader object
            Reference genome used for fetching sequence/BLAST
        thermo : ThermoParams object, optional
            Thermodynamic parameters of hybridization
        design_ref : bool, optional
            If true, the probe will contain the REF base instead of ALT base
        same_strand : bool, optional
            If true, design probe on + strand. If false, when REF base = C,
            probe will be designed for reverse strand. We have found better
            mutant discrimination when doing so.
        """
        self.chrom = chrom
        self.alt_position = alt_position - 1
        self.length = length
        self.ref_base = ref_base
        self.alt_base = alt_base
        self.fasta, self.reference = self.create_reference(reference)
        self.thermo = thermo
        self.design_ref = design_ref
        self.same_strand = same_strand
        self.strand = "+"
        self.offset = 0
        self.sequence = self.create_sequence()
        self.blast_count = float("nan")

    @staticmethod
    def create_reference(reference):
        """ Check if reference is a file or FastaReader object
        and return FastaReader object and reference file path """
        fasta = reference
        if isinstance(fasta, FastaReader):
            reference = fasta.fasta_path
        elif isinstance(fasta, str):
            reference = fasta
            fasta = FastaReader(fasta)
        else:
            raise ValueError("Fasta must be file path or FastaReader object")
        return fasta, reference

    def create_sequence(self):
        """ Create sequence from coordinates using given Fasta. """
        chrom = self.chrom
        start = self.start
        end = self.end
        if abs(self.offset) >= 0.5 * self.length:
            err = (
                f"Length of {self.length} and offset of {self.offset} "
                "are incompatible - sequence cannot be created"
            )
            raise ValueError(err)
        probe_sequence = self.fasta.get(chrom, start, end).upper()
        assert probe_sequence[self.rel_alt_pos] == self.ref_base
        if not self.design_ref:
            probe_sequence = list(probe_sequence)
            probe_sequence[self.rel_alt_pos] = self.alt_base
            probe_sequence = "".join(probe_sequence)
        if not self.same_strand:
            if self.ref_base == "C":
                probe_sequence = complement(probe_sequence, reverse=True)
                self.strand = "-"
        return probe_sequence

    def get_filtered_blast(self):
        """ Many BLAST hits will not result in pull down. We filter
        those hits here using melting temperature """
        t = self.thermo
        blast_hits = self.blast()

        filtered_hits = []
        for blast_hit in blast_hits:
            probe = blast_hit["qseq"]
            hit = blast_hit["sseq"]
            tm_hit = calc_tm(
                probe,
                dnac1=t.dnac1,
                dnac2=t.dnac2,
                Na=t.Na,
                Mg=t.Mg,
                c_seq=complement(hit),
            )
            if tm_hit >= t.temp:
                filtered_hits.append(blast_hit)
        self.unfiltered_blast = blast_hits
        self.blast_hits = filtered_hits
        self.blast_count = len(filtered_hits)
        return filtered_hits

    def set_offset(self, offset):
        self.offset = offset
        self.sequence = self.create_sequence()

    def set_length(self, length):
        self.length = length
        self.sequence = self.create_sequence()

    @property
    def rel_alt_pos(self):
        """ Property that returns the alt_position within probe sequence """
        return self.length // 2 + self.offset

    @property
    def start(self):
        return self.alt_position - self.rel_alt_pos

    @property
    def end(self):
        return self.start + self.length
