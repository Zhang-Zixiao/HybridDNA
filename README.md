# HybridDNA

HybridDNA is an adaptive hybrid coding framework for DNA data storage. It
classifies binary blocks by their statistical characteristics and assigns
Polar or LDPC outer codes, followed by RS inner protection, hierarchical
indexing and pseudo-random DNA transcoding. Data retrieval combines minimizer
clustering, partial order alignment (POA) consensus reconstruction and
concatenated decoding using Polar SCL and sparse LDPC BP.

## Code Structure

```text
HybridDNA/
|-- encode.py                 # Adaptive encoding and DNA sequence generation
|-- decode.py                 # Reconstruction-aided decoding and file recovery
|-- outer_codes.py            # Polar SCL and LDPC BP
|-- reconstruction.py         # Minimizers, clustering and POA consensus
|-- requirements.txt          # Python dependencies
|-- tools/
|   |-- CMakeLists.txt         # SPOA backend build
|   `-- poa_consensus.cpp      # SPOA consensus interface
|-- tests/
|   |-- test_pipeline.py      # Coding and end-to-end tests
|   `-- test_reconstruction.py # Clustering and noisy-read reconstruction tests
`-- datasets/Dataset/         # Binary, text/document, image, audio and video inputs
```

## Installation

Python 3.11 is recommended.

```bash
python -m pip install -r requirements.txt
```

POA uses [SPOA](https://github.com/rvaser/spoa), through pyspoa on Linux/macOS
or the included C++ interface. On Windows, build the interface with CMake,
Git and a C++17 compiler (for example, MinGW):

```bash
cmake -S tools -B tools/build -G "MinGW Makefiles"
cmake --build tools/build --parallel
```

## Workflow

Encode a file and recover it from DNA sequences:

```bash
python encode.py -i "datasets/Dataset/binaries/run-kolibri-distr.sh" -o results/encoded.fasta
python decode.py -i results/encoded.fasta -o results/recovered.sh --output-bytes
```

Reconstruct noisy FASTA/FASTQ reads and decode in one workflow:

```bash
python decode.py -i results/noisy_reads.fastq -o results/recovered.bin --reconstruct --output-bytes
```

Reconstruction can also run separately:

```bash
python reconstruction.py -i results/noisy_reads.fastq -o results/consensus.fasta
python decode.py -i results/consensus.fasta -o results/recovered.bin --output-bytes
```

Command-line options are available through `--help`.

## Tests

```bash
python -m unittest discover -s tests -v
```
