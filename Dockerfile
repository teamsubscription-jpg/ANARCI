# ANARCI image for RunPod Serverless
#
# Build:  docker build -t anarci-serverless .
# Local test:
#         docker run --rm anarci-serverless python -u handler.py --test_input '{"input": {"sequence": "EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSA"}}'

# Python 3.11 still ships distutils, which setup.py imports
FROM python:3.11-slim-bookworm

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

# HMMER3 (hmmscan/hmmbuild/hmmpress) is required to build and search the HMMs
RUN apt-get update \
    && apt-get install -y --no-install-recommends hmmer \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace/ANARCI
COPY . .

# setup.py's install step downloads germlines from IMGT and builds the HMMs (needs muscle on PATH).
# It only prints pipeline errors, so finish with a smoke test that fails the build if numbering is broken.
RUN pip install biopython runpod "setuptools<70" \
    && chmod +x bin/muscle bin/ANARCI \
    && cp bin/muscle /usr/local/bin/muscle \
    && python setup.py install \
    && rm -rf build build_pipeline/IMGT_sequence_files build_pipeline/muscle_alignments build_pipeline/curated_alignments build_pipeline/HMMs \
    && ANARCI -i EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSA | grep -q "^H "

# Start the RunPod serverless worker
CMD ["python", "-u", "handler.py"]
