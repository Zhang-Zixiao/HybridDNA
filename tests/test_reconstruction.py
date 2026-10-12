from pathlib import Path
import random
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import encode
import decode
from reconstruction import (ReconstructionConfig, cluster_reads, minimizers,
                            poa_consensus, reconstruct_records, UnionFind)


def noisy_copies(sequence, count=7):
    result = []
    for i in range(count):
        bases = list(sequence)
        # Different reads have independent error positions; none is a clean copy.
        position = 50 + 3*i
        bases[position] = "A" if bases[position] != "A" else "C"
        bases.insert(85+3*i, "G")
        del bases[120+3*i]
        result.append((f"unlabelled{i}","".join(bases)))
    return result


class ReconstructionTests(unittest.TestCase):
    def test_minimizers_short_and_shifted(self):
        self.assertEqual(minimizers("ACG",6,10),set())
        sequence = "ACGTCGATCGAGCTACGATCGTACGATC"
        self.assertTrue(minimizers(sequence) & minimizers("T"+sequence))

    def test_poa_indels(self):
        rng = random.Random(9)
        sequence = "".join(rng.choice("ACGT") for _ in range(200))
        self.assertEqual(poa_consensus([s for _,s in noisy_copies(sequence)]),sequence)

    def test_union_cannot_bridge_conflicting_addresses(self):
        union = UnionFind([("Polar",(1,1,"Polar")),(None,None),
                           ("Polar",(2,1,"Polar"))])
        self.assertTrue(union.union(0,1))
        self.assertFalse(union.union(1,2))

    def test_invalid_length_is_not_padded(self):
        blocks,length = encode.make_outer_blocks("0"*128,.96,.08,.12)
        records,_ = encode.build_records(blocks,length,True,.4,.6,3)
        shortened = records[0][1][:90]+records[0][1][91:]
        restored,stats = reconstruct_records([("anonymous",shortened)]*5)
        self.assertEqual(restored,[])
        self.assertEqual(stats["rejected_clusters"],1)

    def test_headers_do_not_control_clusters(self):
        rng = random.Random(8)
        data = bytes(rng.randrange(256) for _ in range(80))
        blocks,length = encode.make_outer_blocks(encode.bytes_to_bits(data),.96,.08,.12)
        records,_ = encode.build_records(blocks,length,True,.4,.6,3)
        reads = [read for _,sequence in records for read in noisy_copies(sequence)]
        groups = cluster_reads(reads)
        relabelled = [(f"img_0_frag_0_read_{i}",s) for i,(_,s) in enumerate(reads)]
        self.assertEqual(groups,cluster_reads(relabelled))
        restored,stats = reconstruct_records(reads)
        self.assertEqual(stats["validated_consensus"],len(records))
        self.assertEqual({s for _,s in restored},{s for _,s in records})

    def test_reconstruction_cli_pipeline(self):
        rng = random.Random(11)
        data = bytes(80) + bytes(rng.randrange(256) for _ in range(80))
        blocks,length = encode.make_outer_blocks(encode.bytes_to_bits(data),.96,.08,.12)
        records,_ = encode.build_records(blocks,length,True,.4,.6,3)
        reads = [read for _,sequence in records for read in noisy_copies(sequence)]
        with tempfile.TemporaryDirectory() as folder:
            source,output = Path(folder)/"reads.fasta",Path(folder)/"recovered.bin"
            encode.write_fasta(reads,source)
            arguments = ["decode.py","-i",str(source),"-o",str(output),"--output-bytes",
                         "--bit-length",str(length),"--reconstruct"]
            with patch.object(sys,"argv",arguments):
                decode.main()
            self.assertEqual(output.read_bytes(),data)


if __name__ == "__main__":
    unittest.main()
