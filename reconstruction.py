"""Sequence-only minimizer clustering and medoid-seeded SPOA reconstruction."""
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
import argparse
import shutil
import subprocess
import tempfile

from Levenshtein import distance
from reedsolo import ReedSolomonError


@dataclass(frozen=True)
class ReconstructionConfig:
    expected_length: int = 200
    length_tolerance: int = 25
    max_distance: int = 20
    k: int = 6
    window: int = 10
    candidates: int = 20


def minimizers(sequence, k=6, window=10):
    if k < 1 or window < k:
        raise ValueError("Minimizer window must be at least k")
    kmers = [sequence[i:i+k] for i in range(len(sequence)-k+1)]
    if not kmers:
        return set()
    width = window-k+1
    if len(kmers) <= width:
        return {min(kmers)}
    return {min(kmers[i:i+width]) for i in range(len(kmers)-width+1)}


def structural_hint(sequence, expected_length):
    from decode import (MARKERS_REVERSE, INDEX_RS_NT, PRIMER_FWD,
                        rs_decode, dna_to_bits, xor_bits, pseudo_random_mask,
                        OUTER_INDEX_BITS, EXPECTED_INNER_FRAGMENTS)
    offset = len(PRIMER_FWD) if expected_length == 200 else 0
    expected = offset + INDEX_RS_NT
    hits = []
    for position in range(max(offset, expected-5), min(len(sequence)-3, expected+6)):
        for marker, (method, tag) in MARKERS_REVERSE.items():
            mismatches = sum(a != b for a,b in zip(marker, sequence[position:position+4]))
            if mismatches <= 1:
                hits.append((mismatches, abs(position-expected), position, method, tag))
    if not hits:
        return None, None
    hits.sort()
    best = [hit for hit in hits if hit[:2] == hits[0][:2]]
    methods = {hit[3] for hit in best}
    method = next(iter(methods)) if len(methods) == 1 else None
    address = None
    if hits[0][0] == 0 and len(best) == 1 and hits[0][2]-offset == INDEX_RS_NT:
        try:
            masked = rs_decode(dna_to_bits(sequence[offset:expected]), 1)
            bits = xor_bits(masked, pseudo_random_mask(len(masked), "index", 0))
            outer, inner = int(bits[:OUTER_INDEX_BITS],2), int(bits[OUTER_INDEX_BITS:],2)
            if outer > 0 and 1 <= inner <= EXPECTED_INNER_FRAGMENTS[method]:
                address = (outer, inner, method)
        except (ValueError, ReedSolomonError):
            pass
    return method, address


def compatible(left, right):
    lm, la = left
    rm, ra = right
    return not ((lm is not None and rm is not None and lm != rm) or
                (la is not None and ra is not None and la != ra))


class UnionFind:
    def __init__(self, hints):
        self.parent = list(range(len(hints)))
        self.hints = list(hints)

    def find(self, node):
        while self.parent[node] != node:
            self.parent[node] = self.parent[self.parent[node]]
            node = self.parent[node]
        return node

    def union(self, left, right):
        left, right = self.find(left), self.find(right)
        if left == right or not compatible(self.hints[left], self.hints[right]):
            return False
        lm, la = self.hints[left]
        rm, ra = self.hints[right]
        self.parent[right] = left
        self.hints[left] = (lm or rm, la or ra)
        return True


def cluster_reads(records, config=ReconstructionConfig()):
    if config.expected_length not in (160,200) or config.max_distance < 0:
        raise ValueError("Expected length must be 160 or 200 and distance nonnegative")
    counts = Counter(sequence.upper() for _,sequence in records
                     if set(sequence.upper()) <= set("ACGT") and
                     abs(len(sequence)-config.expected_length) <= config.length_tolerance)
    centers, hints, groups = [], [], []
    seed_index = defaultdict(set)
    def seeds(sequence):
        # Shared universal primers must not dominate candidate ranking.
        region = sequence[24:-20] if config.expected_length == 200 else sequence[4:]
        return minimizers(region, config.k, config.window)
    for sequence in sorted(counts, key=lambda s:(-counts[s],s)):
        hint = structural_hint(sequence, config.expected_length)
        hits = Counter(cid for seed in seeds(sequence) for cid in seed_index[seed])
        matches = []
        for cid,_ in hits.most_common(config.candidates):
            if compatible(hint,hints[cid]):
                value = distance(sequence, centers[cid])
                if value <= config.max_distance:
                    matches.append((value,cid))
        if matches:
            cid = min(matches)[1]
            method,address = hints[cid]
            hints[cid] = (method or hint[0], address or hint[1])
            groups[cid].extend([sequence]*counts[sequence])
        else:
            cid = len(centers)
            centers.append(sequence)
            hints.append(hint)
            groups.append([sequence]*counts[sequence])
            for seed in seeds(sequence):
                seed_index[seed].add(cid)
    uf = UnionFind(hints)
    for cid,center in enumerate(centers):
        hits = Counter(other for seed in seeds(center) for other in seed_index[seed]
                       if other < cid)
        for other,_ in hits.most_common(config.candidates):
            if compatible(uf.hints[uf.find(cid)],uf.hints[uf.find(other)]) and \
                    distance(center,centers[other]) <= config.max_distance:
                uf.union(other,cid)
    result = defaultdict(list)
    for cid,reads in enumerate(groups):
        result[uf.find(cid)].extend(reads)
    return dict(result)


def medoid(reads):
    # Retain frequency weights when choosing among a bounded set of candidates.
    counts = Counter(reads)
    sample = sorted(counts, key=lambda s:(-counts[s],s))[:20]
    return min(sample, key=lambda s:(sum(distance(s,r)*count for r,count in counts.items()),s))


def poa_consensus(reads, executable=None):
    if not reads:
        raise ValueError("POA needs at least one read")
    first = medoid(reads)
    ordered = list(reads)
    ordered.remove(first)
    ordered.insert(0,first)
    local = Path(__file__).parent/"tools"/("poa_consensus.exe" if __import__("os").name == "nt" else "poa_consensus")
    if executable or local.exists():
        command = str(executable or local)
        result = subprocess.run([command],input="\n".join(ordered)+"\n",text=True,
                                capture_output=True,check=True,timeout=120)
        consensus = result.stdout.strip()
    else:
        try:
            import spoa
        except ImportError:
            program = shutil.which("spoa")
            if not program:
                raise RuntimeError("Install pyspoa or build tools/poa_consensus; see --help") from None
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory)/"reads.fasta"
                path.write_text("".join(f">read{i}\n{s}\n" for i,s in enumerate(ordered)),encoding="ascii")
                result = subprocess.run([program,"-l","1","-r","0",str(path)],
                                        capture_output=True,text=True,check=True,timeout=120)
                consensus = "".join(line.strip() for line in result.stdout.splitlines()
                                    if line and not line.startswith(">"))
        else:
            consensus,_ = spoa.poa(ordered,algorithm=1,m=5,n=-4,g=-8,e=-6)
    if not consensus or set(consensus) - set("ACGT"):
        raise ValueError("POA returned an invalid DNA sequence")
    return consensus


def reconstruct_records(records, config=ReconstructionConfig(), executable=None):
    from decode import decode_sequence
    clusters = cluster_reads(records,config)
    reconstructed = []
    for cid,reads in clusters.items():
        consensus = poa_consensus(reads,executable)
        if len(consensus) != config.expected_length or decode_sequence(consensus) is None:
            continue
        reconstructed.append((f"cluster{cid}|coverage={len(reads)}",consensus))
    return reconstructed,{"input_reads":len(records),"clusters":len(clusters),
                           "validated_consensus":len(reconstructed),
                           "rejected_clusters":len(clusters)-len(reconstructed)}


def main():
    from decode import read_sequences
    from encode import write_fasta
    parser = argparse.ArgumentParser(description="Minimizer clustering and SPOA consensus reconstruction")
    parser.add_argument("-i","--input",type=Path,required=True)
    parser.add_argument("-o","--output",type=Path,required=True)
    parser.add_argument("--bit-length",type=int)
    parser.add_argument("--read-length",type=int,choices=[160,200],default=200)
    parser.add_argument("--max-distance",type=int,default=20)
    parser.add_argument("--poa-executable",type=Path,help="Optional compiled tools/poa_consensus executable")
    args = parser.parse_args()
    records,length = read_sequences(args.input)
    if args.bit_length is not None:
        length = args.bit_length
    restored,stats = reconstruct_records(records,ReconstructionConfig(
        expected_length=args.read_length,max_distance=args.max_distance),args.poa_executable)
    if length is not None:
        restored = [(name+f"|bits={length}",sequence) for name,sequence in restored]
    args.output.parent.mkdir(parents=True,exist_ok=True)
    write_fasta(restored,args.output)
    print(stats)


if __name__ == "__main__":
    main()
