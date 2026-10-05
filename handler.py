"""
RunPod serverless handler for ANARCI (antibody / TCR numbering).

Example job input:
    {"input": {"sequence": "EVQLQQSGAEVVRSGASVKLSCTAS..."}}
    {"input": {"sequences": [{"id": "H", "sequence": "EVQL..."}, {"id": "L", "sequence": "DIVM..."}],
               "scheme": "kabat", "assign_germline": true, "allowed_species": ["human", "mouse"]}}
    {"input": {"fasta": ">seq1\nEVQL...\n>seq2\nDIVM..."}}

`sequences` may also be a plain list of strings. Optional: scheme (imgt, chothia, kabat, martin, aho, wolfguy),
restrict (list of chain types from H, K, L, A, B, G, D), bit_score_threshold (default 80).
Gap positions ('-') are left out of the returned numbering.
"""
import os

import runpod

from anarci import run_anarci
from anarci.anarci import validate_sequence

SCHEMES = ('imgt', 'chothia', 'kabat', 'martin', 'aho', 'wolfguy', 'i', 'c', 'k', 'm', 'a', 'w')
CHAIN_TYPES = ('H', 'K', 'L', 'A', 'B', 'G', 'D')
MAX_SEQUENCES = int(os.environ.get('ANARCI_MAX_SEQUENCES', '10000'))


def parse_fasta(text):
    sequences, name, chunks = [], None, []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith('>'):
            if name is not None:
                sequences.append((name, ''.join(chunks)))
            name, chunks = line[1:].strip() or f'seq{len(sequences) + 1}', []
        elif line:
            if name is None:
                raise ValueError("FASTA text must start with a '>' header line")
            chunks.append(line)
    if name is not None:
        sequences.append((name, ''.join(chunks)))
    return sequences


def get_sequences(job_input):
    given = [k for k in ('sequence', 'sequences', 'fasta') if job_input.get(k)]
    if len(given) != 1:
        raise ValueError("Provide exactly one of 'sequence', 'sequences' or 'fasta'")
    key = given[0]
    value = job_input[key]

    if key == 'sequence':
        if not isinstance(value, str):
            raise ValueError("'sequence' must be a string")
        sequences = [('seq1', value)]
    elif key == 'fasta':
        if not isinstance(value, str):
            raise ValueError("'fasta' must be a string")
        sequences = parse_fasta(value)
    else:
        if not isinstance(value, list):
            raise ValueError("'sequences' must be a list")
        sequences = []
        for i, item in enumerate(value, 1):
            if isinstance(item, str):
                sequences.append((f'seq{i}', item))
            elif isinstance(item, dict) and isinstance(item.get('sequence'), str):
                sequences.append((str(item.get('id') or f'seq{i}'), item['sequence']))
            else:
                raise ValueError("Each entry of 'sequences' must be a string or {'id': ..., 'sequence': ...}")

    sequences = [(name, ''.join(seq.split()).upper()) for name, seq in sequences]
    if not sequences:
        raise ValueError('No sequences given')
    if len(sequences) > MAX_SEQUENCES:
        raise ValueError(f'At most {MAX_SEQUENCES} sequences per job')
    return sequences


def format_domain(numbering, start, end, details):
    domain = {
        'chain_type': details.get('chain_type'),
        'species': details.get('species'),
        'evalue': details.get('evalue'),
        'bitscore': details.get('bitscore'),
        'scheme': details.get('scheme'),
        'start': start,
        'end': end,
        'numbering': [
            {'position': pos, 'insertion': ins.strip(), 'residue': aa}
            for (pos, ins), aa in numbering if aa != '-'
        ],
    }
    germlines = details.get('germlines')
    if germlines:
        domain['germlines'] = {
            gene: {'species': hit[0][0], 'gene': hit[0][1], 'identity': hit[1]}
            for gene, hit in germlines.items() if hit and hit[0]
        }
    return domain


def handler(job):
    job_input = job.get('input') or {}
    try:
        sequences = get_sequences(job_input)

        scheme = str(job_input.get('scheme', 'imgt')).lower()
        if scheme not in SCHEMES:
            raise ValueError(f"Invalid scheme: {scheme}. Choose from imgt, chothia, kabat, martin, aho, wolfguy.")

        allow = set(job_input.get('restrict') or CHAIN_TYPES)
        if not allow <= set(CHAIN_TYPES):
            raise ValueError(f"'restrict' must be a subset of {', '.join(CHAIN_TYPES)}")
        # Chothia, Kabat, Martin and Wolfguy are only defined for antibodies (same as the ANARCI CLI)
        if scheme in ('chothia', 'kabat', 'martin', 'wolfguy', 'c', 'k', 'm', 'w'):
            allow &= {'H', 'K', 'L'}

        kwargs = {
            'scheme': scheme,
            'allow': allow,
            'assign_germline': bool(job_input.get('assign_germline', False)),
            'bit_score_threshold': float(job_input.get('bit_score_threshold', 80)),
        }
        if job_input.get('allowed_species'):
            kwargs['allowed_species'] = list(job_input['allowed_species'])

        for name, seq in sequences:
            validate_sequence(seq)

        ncpu = max(1, min(os.cpu_count() or 1, len(sequences)))
        _, numbered, details, _ = run_anarci(sequences, ncpu=ncpu, **kwargs)
    except (ValueError, TypeError, AssertionError) as exc:
        return {'error': str(exc)}

    results = []
    for (name, seq), domains, domain_details in zip(sequences, numbered, details):
        results.append({
            'id': name,
            'sequence': seq,
            'domains': [
                format_domain(numbering, start, end, d)
                for (numbering, start, end), d in zip(domains or [], domain_details or [])
            ],
        })
    return {'results': results}


if __name__ == '__main__':
    runpod.serverless.start({'handler': handler})
