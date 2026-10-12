"""Polar SCL and sparse LDPC BP decoding for the HybridDNA outer layer."""
import math

import numpy as np
from pyldpc import decode as bp_decode, utils


def transform(bits):
    out = np.asarray(bits, dtype=np.int8).copy()
    step = 1
    while step < len(out):
        view = out.reshape(-1, 2 * step)
        view[:, :step] ^= view[:, step:]
        step *= 2
    return out


def information_mask(n, k):
    if n <= 0 or n & (n - 1) or not 0 < k <= n:
        raise ValueError("Polar requires power-of-two N and 0 < K <= N")
    weights = [sum(((i >> j) & 1) * 2 ** (j / 4)
                   for j in range(int(math.log2(n)))) for i in range(n)]
    mask = np.zeros(n, dtype=bool)
    mask[np.argsort(weights, kind="stable")[-k:]] = True
    return mask


def channel_llrs(codeword, probability):
    if not 0 < probability < .5:
        raise ValueError("Residual bit error probability must be between 0 and 0.5")
    if any(bit not in "01" for bit in codeword):
        raise ValueError("Codeword must contain binary digits")
    bits = np.array([int(bit) for bit in codeword], dtype=np.int8)
    return (1 - 2 * bits.astype(float)) * math.log((1 - probability) / probability)


def bit_llr(channel, prefix, index):
    if len(channel) == 1:
        return float(channel[0])
    half = len(channel) // 2
    a, b = channel[:half], channel[half:]
    if index < half:
        check = np.logaddexp(0, a + b) - np.logaddexp(a, b)
        return bit_llr(check, prefix[:half], index)
    partial = transform(prefix[:half])
    return bit_llr(b + (1 - 2 * partial) * a, prefix[half:], index - half)


def scl_decode(codeword, mask, probability=.02, list_size=4):
    """Exact log-domain path metric; no CRC or clean inverse shortcut."""
    if list_size < 1 or len(codeword) != len(mask):
        raise ValueError("Invalid list size or Polar codeword length")
    channel = channel_llrs(codeword, probability)
    paths = [(0., np.zeros(len(mask), dtype=np.int8))]
    for index in range(len(mask)):
        candidates = []
        for metric, prefix in paths:
            llr = bit_llr(channel, prefix, index)
            for bit in ([0, 1] if mask[index] else [0]):
                updated = prefix.copy()
                updated[index] = bit
                metric_next = metric + float(np.logaddexp(0, -(1 - 2 * bit) * llr))
                candidates.append((metric_next, updated))
        paths = sorted(candidates, key=lambda item: item[0])[:list_size]
    return "".join(str(int(bit)) for bit in paths[0][1][mask])


def ldpc_decode(codeword, parity_check, generator, k, probability=.02, maxiter=300):
    if len(codeword) != parity_check.shape[1] or maxiter < 1:
        raise ValueError("Invalid LDPC codeword length or iteration limit")
    # pyldpc uses LLR=2*y/variance. At snr=0 variance is one.
    observation = channel_llrs(codeword, probability) / 2
    decoded = bp_decode(parity_check, observation, snr=0, maxiter=maxiter)
    if not utils.incode(parity_check, decoded):
        raise ValueError("LDPC BP did not converge to a valid codeword")
    dimension = generator.shape[1]
    if not np.array_equal(generator[:dimension], np.eye(dimension, dtype=int)):
        raise ValueError("A systematic LDPC generator is required")
    message = decoded[:dimension]
    if np.any(message[k:]):
        raise ValueError("Invalid shortened LDPC information positions")
    return "".join(str(int(bit)) for bit in message[:k])
