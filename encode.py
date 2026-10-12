"""
HybridDNA encoder.

Pipeline:
1. Convert an input file to a binary stream and split it into fixed-size blocks.
2. Classify each block by bit-distribution features.
3. Encode distribution-biased blocks with Polar codes and random-like blocks with LDPC.
4. Protect hierarchical index and payload with RS inner codes.
5. Apply pseudo-random DNA transcoding and method markers, then write FASTA.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import time
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
from pyldpc import make_ldpc, utils
from reedsolo import RSCodec
from outer_codes import information_mask


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

DNA_MAP = {"00": "A", "01": "C", "10": "G", "11": "T"}


@dataclass(frozen=True)
class DistributionFeatures:
    entropy: float
    balance: float
    adjacent_equal_ratio: float
    method: str


@dataclass(frozen=True)
class OuterBlock:
    bits: str
    method: str
    features: DistributionFeatures


def bytes_to_bits(data: bytes) -> str:
    return "".join(f"{byte:08b}" for byte in data)


def bits_to_int_list(bits: str) -> list[int]:
    return [1 if bit == "1" else 0 for bit in bits]


def int_array_to_bits(values: np.ndarray | list[int]) -> str:
    return "".join(str(int(value)) for value in values)


def pad_bits(bits: str, block_size: int) -> str:
    if len(bits) >= block_size:
        return bits
    return bits + ("0" * (block_size - len(bits)))


def make_outer_blocks(
    bits: str,
    entropy_threshold: float,
    balance_threshold: float,
    correlation_threshold: float,
) -> tuple[list[OuterBlock], int]:
    original_bit_length = len(bits)
    blocks: list[OuterBlock] = []

    for start in range(0, original_bit_length, LDPC_DATA_BITS):
        analysis_bits = bits[start : start + LDPC_DATA_BITS]
        features = classify_distribution(
            analysis_bits,
            entropy_threshold,
            balance_threshold,
            correlation_threshold,
        )

        if features.method == "LDPC":
            blocks.append(OuterBlock(pad_bits(analysis_bits, LDPC_DATA_BITS), "LDPC", features))
            continue

        for sub_start in range(0, len(analysis_bits), POLAR_DATA_BITS):
            polar_bits = analysis_bits[sub_start : sub_start + POLAR_DATA_BITS]
            polar_features = classify_distribution(
                polar_bits,
                entropy_threshold,
                balance_threshold,
                correlation_threshold,
            )
            blocks.append(OuterBlock(pad_bits(polar_bits, POLAR_DATA_BITS), "Polar", polar_features))

    return blocks, original_bit_length


def binary_entropy(one_ratio: float) -> float:
    if one_ratio <= 0.0 or one_ratio >= 1.0:
        return 0.0
    return -(one_ratio * math.log2(one_ratio) + (1.0 - one_ratio) * math.log2(1.0 - one_ratio))


def classify_distribution(
    block: str,
    entropy_threshold: float,
    balance_threshold: float,
    correlation_threshold: float,
) -> DistributionFeatures:
    ones = block.count("1")
    one_ratio = ones / len(block)
    entropy = binary_entropy(one_ratio)
    balance = abs(one_ratio - 0.5)
    if len(block) > 1:
        adjacent_equal_ratio = sum(1 for i in range(len(block) - 1) if block[i] == block[i + 1]) / (len(block) - 1)
    else:
        adjacent_equal_ratio = 1.0

    is_biased = (
        entropy < entropy_threshold
        or balance > balance_threshold
        or abs(adjacent_equal_ratio - 0.5) > correlation_threshold
    )
    method = "Polar" if is_biased else "LDPC"
    return DistributionFeatures(entropy, balance, adjacent_equal_ratio, method)


def polar_info_index(code_bits: int, data_bits: int) -> np.ndarray:
    return information_mask(code_bits, data_bits)



def polar_transform(bits: np.ndarray) -> np.ndarray:
    encoded = bits.copy()
    step = 1
    n = len(encoded)
    while step < n:
        for start in range(0, n, step * 2):
            left = slice(start, start + step)
            right = slice(start + step, start + 2 * step)
            encoded[left] ^= encoded[right]
        step *= 2
    return encoded


def polar_encode(block: str) -> str:
    if len(block) != POLAR_DATA_BITS:
        raise ValueError(f"Polar block must be {POLAR_DATA_BITS} bits.")
    info_index = polar_info_index(POLAR_CODE_BITS, POLAR_DATA_BITS)
    u = np.zeros(POLAR_CODE_BITS, dtype=np.int8)
    u[info_index] = np.array(bits_to_int_list(block), dtype=np.int8)
    return int_array_to_bits(polar_transform(u))


@lru_cache(maxsize=1)
def ldpc_matrices() -> tuple[np.ndarray, np.ndarray]:
    return make_ldpc(LDPC_CODE_BITS, LDPC_DV, LDPC_DC, systematic=True, seed=LDPC_SEED)


def ldpc_encode(block: str, generator: np.ndarray) -> str:
    if len(block) != LDPC_DATA_BITS:
        raise ValueError(f"LDPC block must be {LDPC_DATA_BITS} bits.")
    message_len = generator.shape[1]
    message = np.array(bits_to_int_list(block + ("0" * (message_len - len(block)))), dtype=np.int8)
    codeword = utils.binaryproduct(generator, message) % 2
    return int_array_to_bits(codeword)


def encode_outer(block: str, method: str, generator: np.ndarray) -> str:
    if method == "Polar":
        return polar_encode(block)
    if method == "LDPC":
        return ldpc_encode(block, generator)
    raise ValueError(f"Unsupported outer code: {method}")


def rs_encode(bits: str, rs_bytes: int) -> str:
    if len(bits) % 8 != 0:
        raise ValueError("RS input length must be byte-aligned.")
    payload = bytes(int(bits[i : i + 8], 2) for i in range(0, len(bits), 8))
    encoded = RSCodec(rs_bytes).encode(payload)
    return "".join(f"{byte:08b}" for byte in encoded)


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


def binary_to_dna(bits: str) -> str:
    if len(bits) % 2:
        bits += "0"
    return "".join(DNA_MAP[bits[i : i + 2]] for i in range(0, len(bits), 2))


def gc_ratio(sequence: str) -> float:
    return (sequence.count("G") + sequence.count("C")) / len(sequence)


def max_homopolymer(sequence: str) -> int:
    if not sequence:
        return 0
    longest = 1
    current = 1
    for previous, base in zip(sequence, sequence[1:]):
        if base == previous:
            current += 1
            longest = max(longest, current)
        else:
            current = 1
    return longest


def passes_dna_constraints(sequence: str, min_gc: float, max_gc: float, max_run: int) -> bool:
    gc = gc_ratio(sequence)
    return min_gc <= gc <= max_gc and max_homopolymer(sequence) <= max_run


def randomize_payload(bits: str, method: str, min_gc: float, max_gc: float, max_run: int) -> tuple[str, int]:
    for tag in range(4):
        candidate = bits if tag == 0 else xor_bits(bits, pseudo_random_mask(len(bits), "payload", tag))
        marker = MARKERS[method][tag]
        dna_payload = binary_to_dna(candidate)
        if passes_dna_constraints(marker + dna_payload, min_gc, max_gc, max_run):
            return candidate, tag
    return xor_bits(bits, pseudo_random_mask(len(bits), "payload", 3)), 3


def build_records(
    blocks: list[OuterBlock],
    original_bit_length: int,
    include_primers: bool,
    min_gc: float,
    max_gc: float,
    max_run: int,
) -> tuple[list[tuple[str, str]], dict[str, int]]:
    _, generator = ldpc_matrices()
    records = []
    stats = {"Polar": 0, "LDPC": 0}

    for outer_index, block in enumerate(blocks, start=1):
        stats[block.method] += 1
        outer_codeword = encode_outer(block.bits, block.method, generator)

        for inner_index, start in enumerate(range(0, len(outer_codeword), INNER_PAYLOAD_BITS), start=1):
            payload = outer_codeword[start : start + INNER_PAYLOAD_BITS]
            if len(payload) < INNER_PAYLOAD_BITS:
                payload += "0" * (INNER_PAYLOAD_BITS - len(payload))

            index_bits = f"{outer_index:0{OUTER_INDEX_BITS}b}{inner_index:0{INNER_INDEX_BITS}b}"
            index_rs = rs_encode(xor_bits(index_bits, pseudo_random_mask(len(index_bits), "index", 0)), INDEX_RS_BYTES)
            data_rs = rs_encode(payload, DATA_RS_BYTES)
            randomized_payload, tag = randomize_payload(data_rs, block.method, min_gc, max_gc, max_run)

            dna = binary_to_dna(index_rs) + MARKERS[block.method][tag] + binary_to_dna(randomized_payload)
            if include_primers:
                dna = PRIMER_FWD + dna + PRIMER_REV

            record_id = (
                f"seq{len(records) + 1}|outer={outer_index}|inner={inner_index}|"
                f"method={block.method}|bits={original_bit_length}|"
                f"outer_bits={len(block.bits)}|"
                f"entropy={block.features.entropy:.4f}|balance={block.features.balance:.4f}|"
                f"corr={block.features.adjacent_equal_ratio:.4f}"
            )
            records.append((record_id, dna))

    return records, stats


def write_fasta(records: list[tuple[str, str]], output_path: Path) -> None:
    with output_path.open("w", encoding="ascii") as handle:
        for record_id, sequence in records:
            handle.write(f">{record_id}\n")
            handle.write(f"{sequence}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Encode a file with the HybridDNA adaptive Polar/LDPC strategy.")
    parser.add_argument("-i", "--input", required=True, type=Path, help="Input binary file.")
    parser.add_argument("-o", "--output", required=True, type=Path, help="Output FASTA file.")
    parser.add_argument("--no-primers", action="store_true", help="Do not add universal primers.")
    parser.add_argument("--entropy-threshold", type=float, default=0.96)
    parser.add_argument("--balance-threshold", type=float, default=0.08)
    parser.add_argument("--correlation-threshold", type=float, default=0.12)
    parser.add_argument("--min-gc", type=float, default=0.40)
    parser.add_argument("--max-gc", type=float, default=0.60)
    parser.add_argument("--max-homopolymer", type=int, default=3)
    parser.add_argument("--strict-constraints", action="store_true",
                        help="Refuse output if any complete strand fails GC/run constraints.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    start = time.perf_counter()
    input_bytes = args.input.read_bytes()
    input_bits = bytes_to_bits(input_bytes)
    blocks, original_bit_length = make_outer_blocks(
        input_bits,
        args.entropy_threshold,
        args.balance_threshold,
        args.correlation_threshold,
    )
    records, stats = build_records(
        blocks,
        original_bit_length,
        include_primers=not args.no_primers,
        min_gc=args.min_gc,
        max_gc=args.max_gc,
        max_run=args.max_homopolymer,
    )
    violations = sum(not passes_dna_constraints(sequence, args.min_gc, args.max_gc,
                                               args.max_homopolymer)
                     for _, sequence in records)
    if violations and args.strict_constraints:
        raise ValueError(f"{violations} strands fail constraints; no FASTA written")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_fasta(records, args.output)
    elapsed = time.perf_counter() - start
    print(f"Input bits: {original_bit_length}")
    print(f"Outer blocks: {len(blocks)}")
    print(f"DNA sequences: {len(records)}")
    print(f"Polar blocks: {stats['Polar']}")
    print(f"LDPC blocks: {stats['LDPC']}")
    print(f"Complete-strand biochemical constraint violations: {violations}")
    print(f"Encoding time: {elapsed:.4f} s")


if __name__ == "__main__":
    main()
