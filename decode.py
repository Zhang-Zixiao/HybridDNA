"""
HybridDNA decoder.

Pipeline:
1. Read FASTA/FASTQ sequences.
   Optionally cluster reads and reconstruct POA consensus sequences.
2. Locate method markers and split index/payload regions.
3. RS-decode index and payload, then reverse pseudo-random masks.
4. Reassemble outer codewords by hierarchical index.
5. Decode each codeword with the marker-selected Polar or LDPC decoder.
6. Concatenate recovered blocks and trim to the original bit length.
"""

from __future__ import annotations

import argparse
import hashlib
import re
import time
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
from Bio import SeqIO
from pyldpc import make_ldpc
from reedsolo import RSCodec, ReedSolomonError
from outer_codes import information_mask, scl_decode, ldpc_decode as bp_outer_decode


POLAR_DATA_BITS = 128
POLAR_CODE_BITS = 256
LDPC_DATA_BITS = 640
LDPC_CODE_BITS = 1024
LDPC_DV = 3
LDPC_DC = 8
LDPC_SEED = 42
OUTER_INDEX_BITS = 28
INNER_INDEX_BITS = 4
INNER_PAYLOAD_BITS = 256
INDEX_RS_BYTES = 1
DATA_RS_BYTES = 2
PRIMER_FWD = "CCTGCAGAGTAGCATGTCAG"
PRIMER_REV = "CTGACACTGATGCATCCGAT"
INDEX_RS_BITS = OUTER_INDEX_BITS + INNER_INDEX_BITS + INDEX_RS_BYTES * 8
DATA_RS_BITS = INNER_PAYLOAD_BITS + DATA_RS_BYTES * 8
INDEX_RS_NT = INDEX_RS_BITS // 2
DATA_RS_NT = DATA_RS_BITS // 2
MARKER_NT = 4
EXPECTED_INNER_FRAGMENTS = {
    "Polar": POLAR_CODE_BITS // INNER_PAYLOAD_BITS,
    "LDPC": LDPC_CODE_BITS // INNER_PAYLOAD_BITS,
}

MARKERS = {
    "LDPC": {
        0: "ACCT",
        1: "ACGT",
        2: "AGCT",
        3: "AGGT",
    },
    "Polar": {
        0: "TCCA",
        1: "TCGA",
        2: "TGCA",
        3: "TGGA",
    },
}

MARKERS_REVERSE = {
    marker: (method, tag)
    for method, method_markers in MARKERS.items()
    for tag, marker in method_markers.items()
}

BASE_TO_BITS = {"A": "00", "C": "01", "G": "10", "T": "11"}


@dataclass(frozen=True)
class DecodedInnerBlock:
    outer_index: int
    inner_index: int
    method: str
    payload: str


def pseudo_random_mask(length: int, namespace: str, tag: int) -> str:
    out = []
    counter = 0
    while len(out) < length:
        seed = f"HybridDNA:{namespace}:{tag}:{counter}".encode("ascii")
        digest = hashlib.sha256(seed).digest()
        out.extend(f"{byte:08b}" for byte in digest)
        counter += 1
    return "".join(out)[:length]


def xor_bits(left: str, right: str) -> str:
    return "".join("1" if a != b else "0" for a, b in zip(left, right))


def dna_to_bits(sequence: str) -> str:
    try:
        return "".join(BASE_TO_BITS[base] for base in sequence)
    except KeyError as exc:
        raise ValueError(f"Unsupported DNA base: {exc.args[0]}") from exc


def rs_decode(bits: str, rs_bytes: int) -> str:
    if len(bits) % 8 != 0:
        raise ValueError("RS input length must be byte-aligned.")
    payload = bytes(int(bits[i : i + 8], 2) for i in range(0, len(bits), 8))
    decoded = RSCodec(rs_bytes).decode(payload)
    if isinstance(decoded, tuple):
        decoded_payload = decoded[0]
    else:
        decoded_payload = decoded
    return "".join(f"{byte:08b}" for byte in decoded_payload)


def polar_decode(codeword: str, bit_error_probability: float, list_size: int = 4) -> str:
    return scl_decode(codeword, information_mask(POLAR_CODE_BITS, POLAR_DATA_BITS),
                      bit_error_probability, list_size)


@lru_cache(maxsize=1)
def ldpc_matrices():
    return make_ldpc(LDPC_CODE_BITS, LDPC_DV, LDPC_DC, systematic=True, seed=LDPC_SEED)


def ldpc_decode(codeword: str, bit_error_probability: float = .02, maxiter: int = 300) -> str:
    parity_check, generator = ldpc_matrices()
    return bp_outer_decode(codeword, parity_check, generator, LDPC_DATA_BITS,
                           bit_error_probability, maxiter)


def decode_outer(codeword, method, probability, list_size, maxiter):
    if method == "Polar":
        return polar_decode(codeword, probability, list_size)
    if method == "LDPC":
        return ldpc_decode(codeword, probability, maxiter)
    raise ValueError(f"Unsupported outer code: {method}")



def detect_format(input_path: Path) -> str:
    suffix = input_path.suffix.lower()
    if suffix in {".fastq", ".fq"}:
        return "fastq"
    return "fasta"


def strip_primers(sequence: str) -> str:
    sequence = sequence.upper().strip()
    if sequence.startswith(PRIMER_FWD) and sequence.endswith(PRIMER_REV):
        return sequence[len(PRIMER_FWD) : -len(PRIMER_REV)]
    return sequence


def find_marker(sequence: str) -> tuple[str, int, int] | None:
    expected_pos = INDEX_RS_NT
    exact = sequence[expected_pos : expected_pos + MARKER_NT]
    if exact in MARKERS_REVERSE:
        method, tag = MARKERS_REVERSE[exact]
        return method, tag, expected_pos

    search_start = max(0, expected_pos - 2)
    search_end = min(len(sequence), expected_pos + MARKER_NT + 3)
    candidates = []
    for marker, (method, tag) in MARKERS_REVERSE.items():
        offset = sequence[search_start:search_end].find(marker)
        if offset >= 0:
            pos = search_start + offset
            candidates.append((abs(pos - expected_pos), method, tag, pos))
    if candidates:
        _, method, tag, pos = sorted(candidates, key=lambda item: item[0])[0]
        return method, tag, pos
    return None


def decode_sequence(sequence: str) -> DecodedInnerBlock | None:
    sequence = strip_primers(sequence)
    marker_hit = find_marker(sequence)
    if marker_hit is None:
        return None

    method, tag, marker_pos = marker_hit
    index_region = sequence[:marker_pos]
    payload_region = sequence[marker_pos + MARKER_NT : marker_pos + MARKER_NT + DATA_RS_NT]
    if len(index_region) != INDEX_RS_NT or len(payload_region) != DATA_RS_NT:
        return None

    try:
        index_bits_rs = dna_to_bits(index_region)
        payload_bits_rs = dna_to_bits(payload_region)
        index_bits_masked = rs_decode(index_bits_rs, INDEX_RS_BYTES)
        payload_bits_unmasked = payload_bits_rs if tag == 0 else xor_bits(
            payload_bits_rs,
            pseudo_random_mask(len(payload_bits_rs), "payload", tag),
        )
        payload_bits = rs_decode(payload_bits_unmasked, DATA_RS_BYTES)
    except (ValueError, ReedSolomonError):
        return None

    index_bits = xor_bits(index_bits_masked, pseudo_random_mask(len(index_bits_masked), "index", 0))

    outer_index = int(index_bits[:OUTER_INDEX_BITS], 2)
    inner_index = int(index_bits[OUTER_INDEX_BITS : OUTER_INDEX_BITS + INNER_INDEX_BITS], 2)
    return DecodedInnerBlock(outer_index, inner_index, method, payload_bits)


def read_sequences(input_path: Path) -> tuple[list[tuple[str, str]], int | None]:
    records = []
    bit_length = None
    for record in SeqIO.parse(str(input_path), detect_format(input_path)):
        description = record.description
        if bit_length is None:
            match = re.search(r"(?:^|\|)bits=(\d+)(?:\||$)", description)
            if match:
                bit_length = int(match.group(1))
        records.append((description, str(record.seq)))
    return records, bit_length


def bits_to_bytes(bits: str) -> bytes:
    if len(bits) % 8:
        bits += "0" * (8 - len(bits) % 8)
    return bytes(int(bits[i : i + 8], 2) for i in range(0, len(bits), 8))


def write_output(bits: str, output_path: Path, as_bytes: bool) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if as_bytes:
        output_path.write_bytes(bits_to_bytes(bits))
    else:
        output_path.write_text(bits, encoding="ascii")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Decode a HybridDNA FASTA/FASTQ file.")
    parser.add_argument("-i", "--input", required=True, type=Path, help="Input FASTA/FASTQ file.")
    parser.add_argument("-o", "--output", required=True, type=Path, help="Output recovered bit-string or bytes.")
    parser.add_argument("--bit-length", type=int, help="Original bit length. Overrides FASTA header metadata.")
    parser.add_argument("--output-bytes", action="store_true", help="Write recovered bytes instead of a 0/1 text bit-string.")
    parser.add_argument("--bit-error-probability", type=float, default=0.02,
                        help="Residual binary error probability after reconstruction and RS.")
    parser.add_argument("--polar-list-size", type=int, default=4)
    parser.add_argument("--ldpc-maxiter", type=int, default=300)
    parser.add_argument("--reconstruct", action="store_true",
                        help="Cluster noisy reads and generate POA consensus before decoding.")
    parser.add_argument("--read-length", type=int, choices=[160, 200], default=200)
    parser.add_argument("--reconstruction-max-distance", type=int, default=20)
    parser.add_argument("--poa-executable", type=Path,
                        help="Optional compiled tools/poa_consensus executable.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    start = time.perf_counter()
    raw_records, header_bit_length = read_sequences(args.input)
    original_bit_length = args.bit_length if args.bit_length is not None else header_bit_length
    if args.reconstruct:
        from reconstruction import ReconstructionConfig, reconstruct_records
        raw_records, reconstruction_stats = reconstruct_records(
            raw_records, ReconstructionConfig(expected_length=args.read_length,
                max_distance=args.reconstruction_max_distance), args.poa_executable)
        print(f"Reconstruction: {reconstruction_stats}")

    decoded_blocks = []
    for _, sequence in raw_records:
        block = decode_sequence(sequence)
        if block is not None:
            decoded_blocks.append(block)

    grouped: dict[int, list[DecodedInnerBlock]] = {}
    for block in decoded_blocks:
        grouped.setdefault(block.outer_index, []).append(block)

    recovered_blocks = []
    skipped_outer_blocks = 0
    if original_bit_length is None:
        raise ValueError("Original bit length is required in headers or --bit-length")
    if grouped and sorted(grouped) != list(range(1, max(grouped) + 1)):
        raise ValueError("Missing outer indices; refusing to concatenate shifted data")
    for outer_index in sorted(grouped):
        by_inner: dict[int, DecodedInnerBlock] = {}
        for item in grouped[outer_index]:
            by_inner.setdefault(item.inner_index, item)
        inner_blocks = [by_inner[index] for index in sorted(by_inner)]
        method = inner_blocks[0].method
        if any(item.method != method for item in grouped[outer_index]):
            raise ValueError("Conflicting methods for an outer index")
        for index in by_inner:
            if any(item.payload != by_inner[index].payload for item in grouped[outer_index]
                   if item.inner_index == index):
                raise ValueError("Conflicting duplicate fragments; reconstruct reads first")
        expected_fragments = EXPECTED_INNER_FRAGMENTS[method]
        if any(index not in by_inner for index in range(1, expected_fragments + 1)):
            raise ValueError("Incomplete outer codeword; refusing partial file output")
        outer_codeword = "".join(by_inner[index].payload for index in range(1, expected_fragments + 1))
        recovered_blocks.append(
            decode_outer(
                outer_codeword,
                method,
                args.bit_error_probability,
                args.polar_list_size,
                args.ldpc_maxiter,
            )
        )

    recovered_bits = "".join(recovered_blocks)
    if len(recovered_bits) < original_bit_length:
        raise ValueError("Insufficient recovered bits; refusing partial file output")
    if original_bit_length is not None:
        recovered_bits = recovered_bits[:original_bit_length]

    write_output(recovered_bits, args.output, args.output_bytes)
    elapsed = time.perf_counter() - start
    print(f"Input sequences: {len(raw_records)}")
    print(f"Decoded inner blocks: {len(decoded_blocks)}")
    print(f"Recovered outer blocks: {len(recovered_blocks)}")
    print(f"Skipped incomplete outer blocks: {skipped_outer_blocks}")
    print(f"Recovered bits: {len(recovered_bits)}")
    print(f"Decoding time: {elapsed:.4f} s")


if __name__ == "__main__":
    main()
