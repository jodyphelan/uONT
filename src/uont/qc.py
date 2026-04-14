import pysam
from .types import FullPath 


def get_sequence_name_from_filename(fasta_file):

    if fasta_file.endswith(".fasta"):
        return fasta_file.split("/")[-1][:-6]
    elif fasta_file.endswith(".fa"):
        return fasta_file.split("/")[-1][:-3]
    elif fasta_file.endswith(".fna"):
        return fasta_file.split("/")[-1][:-4]
    else:
        return None

class Fasta:
    def __init__(self, fasta_file: FullPath, sample_id: str = None):
        self.fasta_file = fasta_file
        self.fasta = pysam.FastaFile(fasta_file)
        self.all_sequence = "".join([self.fasta.fetch(contig) for contig in self.fasta.references])
        self.sample_id = sample_id if sample_id else get_sequence_name_from_filename(fasta_file)

    def fetch(self, contig, start=None, end=None) -> str:
        return self.fasta.fetch(contig, start, end)
    
    def n50(self) -> int:
        lengths = [len(self.fasta.fetch(contig)) for contig in self.fasta.references]
        lengths.sort(reverse=True)
        total_length = sum(lengths)
        cumsum = 0
        for length in lengths:
            cumsum += length
            if cumsum >= total_length / 2:
                return length
        
        return 0
    
    def l90(self) -> int:
        lengths = [len(self.fasta.fetch(contig)) for contig in self.fasta.references]
        lengths.sort(reverse=True)
        total_length = sum(lengths)
        cumsum = 0
        for i, length in enumerate(lengths):
            cumsum += length
            if cumsum >= total_length * 0.9:
                return i + 1
        
        return len(lengths)
    
    def num_contigs(self, min_contig_length=0) -> int:
        return sum(
            1 for contig in self.fasta.references 
            if len(self.fasta.fetch(contig)) >= min_contig_length
        )
    
    def gc_content(self) -> float:
        gc_count = self.all_sequence.count("G") + self.all_sequence.count("C")
        return gc_count / len(self.all_sequence) if len(self.all_sequence) > 0 else 0
    
    def Ns_per_kb(self) -> float:
        n_count = self.all_sequence.count("N")
        return n_count / (len(self.all_sequence) / 1000) if len(self.all_sequence) > 0 else 0
    
    def qc_metrics(self, min_contig_length: int = 0) -> dict:

        metrics = {
            "length": len(self.all_sequence),
            "gc_content": self.gc_content(),
            "Ns_per_kb": self.Ns_per_kb(),
            "num_contigs": self.num_contigs(min_contig_length),
            "n50": self.n50(),
            "l90": self.l90()
        }
        
        return metrics
    
    def write_qc_report(self, output_file, min_contig_length: int =0):
        metrics = self.qc_metrics(min_contig_length)
        with open(output_file, "w") as f:
            f.write("SampleID\tLength\tGC_Content\tNs_per_kb\tNum_Contigs\tN50\tL90\n")
            f.write(f"{self.sample_id}\t{metrics['length']}\t{metrics['gc_content']:.4f}\t{metrics['Ns_per_kb']:.4f}\t{metrics['num_contigs']}\t{metrics['n50']}\t{metrics['l90']}\n")