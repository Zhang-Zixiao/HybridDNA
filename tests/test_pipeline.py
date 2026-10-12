import random
import tempfile
import unittest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import encode
import decode


def flip(word, positions):
    bits = list(word)
    for index in positions:
        bits[index] = "1" if bits[index] == "0" else "0"
    return "".join(bits)


class HybridDNATests(unittest.TestCase):
    def test_payload_rs_correction(self):
        blocks, length = encode.make_outer_blocks("0"*128, .96, .08, .12)
        records, _ = encode.build_records(blocks, length, False, .4, .6, 3)
        sequence = records[0][1]
        position = decode.INDEX_RS_NT + decode.MARKER_NT + 5
        replacement = "A" if sequence[position] != "A" else "C"
        mutated = sequence[:position] + replacement + sequence[position+1:]
        self.assertEqual(decode.decode_sequence(mutated), decode.decode_sequence(sequence))

    def test_polar_scl_correction(self):
        rng = random.Random(42)
        for _ in range(3):
            message = "".join(rng.choice("01") for _ in range(128))
            word = encode.polar_encode(message)
            self.assertEqual(decode.polar_decode(word, .02), message)
            self.assertEqual(decode.polar_decode(flip(word, [17]), .02), message)

    def test_ldpc_bp_correction(self):
        rng = random.Random(42)
        _, generator = encode.ldpc_matrices()
        message = "".join(rng.choice("01") for _ in range(640))
        word = encode.ldpc_encode(message, generator)
        self.assertEqual(decode.ldpc_decode(word), message)
        self.assertEqual(decode.ldpc_decode(flip(word, [17, 81, 239])), message)

    def test_mixed_end_to_end(self):
        rng = random.Random(7)
        data = bytes(80) + bytes(rng.randrange(256) for _ in range(80)) + b"DNA"
        bits = encode.bytes_to_bits(data)
        blocks, length = encode.make_outer_blocks(bits, .96, .08, .12)
        self.assertEqual({block.method for block in blocks}, {"Polar", "LDPC"})
        records, _ = encode.build_records(blocks, length, True, .4, .6, 3)
        fragments = [decode.decode_sequence(sequence) for _, sequence in records]
        self.assertTrue(all(fragment is not None for fragment in fragments))
        restored = []
        for outer_index in range(1, len(blocks)+1):
            group = sorted([f for f in fragments if f.outer_index == outer_index],
                           key=lambda f: f.inner_index)
            restored.append(decode.decode_outer("".join(f.payload for f in group),
                                                group[0].method, .02, 4, 300))
        self.assertEqual(decode.bits_to_bytes("".join(restored)[:length]), data)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"reads.fasta"
            encode.write_fasta(records, path)
            loaded, header_length = decode.read_sequences(path)
            self.assertEqual(header_length, length)
            self.assertEqual(len(loaded), len(records))
            fastq = Path(directory)/"reads.fastq"
            fastq.write_text("".join(f"@{name}\n{sequence}\n+\n{'I'*len(sequence)}\n"
                                     for name,sequence in records), encoding="ascii")
            self.assertEqual(decode.read_sequences(fastq)[1], length)


if __name__ == "__main__":
    unittest.main()
