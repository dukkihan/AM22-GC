#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
calculate_repair.py — DNA repair pathway completeness from DRAM annotations
===========================================================================

WHAT THIS DOES
--------------
DRAM (`DRAM.py distill`) scores KEGG *modules* and electron transport chain
complexes, and writes them to `product.tsv`. The four canonical DNA repair
systems are KEGG *pathways* (ko03410 BER, ko03420 NER, ko03430 MMR,
ko03440 HR), not modules, so DRAM does not evaluate them at all.

This script fills that gap. It reads the same `annotations.tsv` that DRAM
produced, applies curated step definitions for the four repair pathways, and
writes per-genome completeness so the result can be displayed alongside
DRAM's own output.

USAGE
-----
    python calculate_repair.py                          # uses the default paths below
    python calculate_repair.py -i path/to/annotations.tsv -o results_dir

Requires: pandas. Nothing else.

INPUT
-----
`annotations.tsv` from `DRAM.py annotate`. Only two of its columns are used:
    fasta   genome (bin) name
    ko_id   KEGG orthologue assignment(s) for that gene, comma-separated

OUTPUT (three files)
--------------------
1. repair_completeness.tsv   one row per genome, two columns per pathway
     <pathway>          fraction: steps recovered / steps defined
     <pathway> (Pass)   boolean:  diagnostic steps recovered

2. repair_steps.tsv         one row per genome x pathway x step, with the
                            KOs that satisfied it. Use this for a supplementary
                            table and to audit any single cell of a figure.

3. repair_diagnostics.txt   which KOs and which steps were never recovered from
                            ANY genome in the dataset. Read this before
                            interpreting the numbers — see WHY below.

HOW A STEP IS SCORED
--------------------
Each pathway is a list of steps. A step carries a logic type:

    OR      satisfied by any one of the listed KOs
            (alternative orthologues doing the same job)
    AND     satisfied only when every listed KO is present
            (obligate subunits of one complex, e.g. UvrA+UvrB+UvrC)
    GROUPS  satisfied when any one group is complete, AND within each group
            (alternative complete systems, e.g. RecBCD *or* RecFOR *or* MRN)

The two output columns then differ in how they combine the steps:

    fraction  = (number of steps satisfied) / (number of steps defined)
    Pass      = the pathway's diagnostic steps are satisfied

COMPARABILITY WITH DRAM — READ THIS BEFORE PUTTING BOTH IN ONE FIGURE
---------------------------------------------------------------------
DRAM's *module completeness* is `num_steps_present / num_steps`. The **fraction**
columns here use the same definition, so the two are comparable in kind.

DRAM's *functional presence* call (the TRUE/FALSE cells of `product.tsv`) is
`np.all(...)` over every gene set it lists for that function — each set
satisfied by any one of its alternative orthologues. It therefore requires
100% of the listed sets. The **Pass** columns here require only the diagnostic
steps, which is a weaker criterion. The two are NOT equivalent.

If a figure shows both, its legend must state which criterion each group of
columns uses. Preferring the fraction columns avoids the problem entirely.

WHY THE KO LISTS LOOK THE WAY THEY DO
-------------------------------------
Repair enzymes are not conserved between Bacteria and Archaea. If a step is
defined with bacterial orthologues only, every archaeal genome fails it for
reasons that have nothing to do with its repair capacity, and the resulting
completeness reads domain rather than biology. Each step below therefore
carries the counterparts of both domains:

    recombinase      RecA (K03553)            | RadA (K04483), RAD51 (K04482)
    DNA ligase       LigA, NAD+ (K01972)      | LIG1, ATP-dependent (K10747)
    repair synthesis PolA (K02335),            | archaeal PolB/PolD
                     PolIII alpha (K02337)     |   (K02319, K02322, K02323)
    NER incision     UvrABC (K03701-03)       | XPB+XPD (K10843, K10844)
    HJ resolution    RuvC, RecG               | Hjc/Hje (K03552)
    DSB processing   RecBCD, RecFOR           | Mre11+Rad50 (K10865, K10866)
    mismatch sensing MutS+MutL (K03555/72)    | also NucS (K07503), the
                                              |   MutS/MutL-independent
                                              |   mismatch endonuclease

Two pitfalls this avoids, both worth checking in any dataset:

  * A step nobody can satisfy silently caps every score. MutH (K03573) is
    restricted to enterobacteria; including it as an MMR step holds every other
    bacterium at 5/6 regardless of its actual repair capacity. It is not used
    here. `repair_diagnostics.txt` reports any step with this problem.

  * Under OR logic, one broadly conserved accessory protein hands out a free
    pass. The DNA polymerase III holoenzyme is a clear case: its beta clamp
    (K02338) and clamp loader (K02339-45) are essential in every bacterium, so
    listing them beside the catalytic alpha subunit makes the polymerase step
    unconditionally true. Only catalytic subunits are listed here.

  * K04485 is annotated "radA, sms" but is a *bacterial* protein and not a
    recombinase. The archaeal RecA/Rad51 counterpart is K04483. Do not
    substitute one for the other.

Every KO can be verified at https://www.kegg.jp/entry/<KO>.
"""

import argparse
import os
import sys

try:
    import pandas as pd
except ImportError:
    sys.exit("pandas is required:  pip install pandas")


# =============================================================================
# Pathway definitions
# =============================================================================
REPAIR_RULES = {
    'Mismatch repair (MMR)': {
        'steps': {
            # Mismatch sensing. Two mutually independent systems:
            # canonical MutS/MutL, or NucS in archaea and actinobacteria.
            'Step1_Recognition':   {'type': 'GROUPS',
                                    'groups': [['K03555', 'K03572'],   # mutS + mutL
                                               ['K07503']]},           # nucS
            'Step2_Helicase':      {'type': 'OR',
                                    'kos': ['K03657']},                # uvrD/pcrA
            'Step3_Exonuclease':   {'type': 'OR',
                                    'kos': ['K01141',                  # sbcB/exoI
                                            'K03601', 'K03602',        # xseA, xseB
                                            'K07462',                  # recJ
                                            'K10857']},                # exoX
            'Step4_Polymerase':    {'type': 'OR',
                                    'kos': ['K02337',                  # dnaE  (PolIII alpha)
                                            'K03763',                  # polC  (PolIII alpha, Gram+)
                                            'K02319',                  # archaeal pol
                                            'K02322', 'K02323']},      # polD large, small
            'Step5_Ligation':      {'type': 'OR',
                                    'kos': ['K01972',                  # ligA  (NAD+)
                                            'K10747']},                # LIG1  (ATP)
        },
        # Diagnostic: mismatch has to be recognised, and the nick has to be sealed.
        'binary_pass': lambda passed: ('Step1_Recognition' in passed)
                                      and ('Step5_Ligation' in passed),
    },

    'Base excision repair (BER)': {
        'steps': {
            # Glycosylases are split by lesion class. They are not
            # interchangeable: a genome able to excise oxidised bases is not
            # thereby able to excise uracil. Deamination (cytosine -> uracil) is
            # the dominant lesion in ancient and long-buried DNA, so keeping
            # that class separate is what makes this pathway informative.
            'Step1a_Uracil_Glycosylase':     {'type': 'OR',
                                              'kos': ['K03648',        # UNG/UDG
                                                      'K03649',        # mug
                                                      'K21929']},      # udg
            'Step1b_Oxidative_Glycosylase':  {'type': 'OR',
                                              'kos': ['K03575',        # mutY
                                                      'K05522',        # nei
                                                      'K10563',        # mutM/fpg
                                                      'K10773']},      # nth
            'Step1c_Alkylation_Glycosylase': {'type': 'OR',
                                              'kos': ['K01246',        # tag
                                                      'K01247',        # alkA
                                                      'K13529']},      # ada-alkA
            'Step2_AP_Endonuclease':         {'type': 'OR',
                                              'kos': ['K01142',        # xthA  (ExoIII class)
                                                      'K01151',        # nfo   (EndoIV class)
                                                      'K10771']},      # APEX1
            'Step3_Synthesis':               {'type': 'OR',
                                              'kos': ['K02335',        # polA
                                                      'K02319',
                                                      'K02322', 'K02323']},
            'Step4_Ligation':                {'type': 'OR',
                                              'kos': ['K01972', 'K10747']},
        },
        # Diagnostic: excise a damaged base, then cut the abasic site.
        'binary_pass': lambda passed: any(s in passed for s in
                                          ['Step1a_Uracil_Glycosylase',
                                           'Step1b_Oxidative_Glycosylase',
                                           'Step1c_Alkylation_Glycosylase'])
                                      and ('Step2_AP_Endonuclease' in passed),
    },

    'Nucleotide excision repair (NER)': {
        'steps': {
            # Bacterial excinuclease, or the archaeal/eukaryotic XPB-XPD system.
            'Step1_Damage_Incision':        {'type': 'GROUPS',
                                             'groups': [['K03701', 'K03702', 'K03703'],
                                                        ['K10843', 'K10844']]},
            'Step2_Helicase':               {'type': 'OR',
                                             'kos': ['K03657']},
            'Step3_Synthesis':              {'type': 'OR',
                                             'kos': ['K02335', 'K02319',
                                                     'K02322', 'K02323']},
            'Step4_Ligation':               {'type': 'OR',
                                             'kos': ['K01972', 'K10747']},
            # Transcription-coupled repair. Mfd has no archaeal counterpart, so
            # this step is bacteria-only by construction; it is kept because it
            # is a real part of the pathway, and excluded from the diagnostic.
            'Step5_Transcription_Coupling': {'type': 'OR',
                                             'kos': ['K03723']},       # mfd
        },
        # Diagnostic: incise the lesion, then complete the patch.
        'binary_pass': lambda passed: ('Step1_Damage_Incision' in passed)
                                      and (('Step3_Synthesis' in passed)
                                           or ('Step4_Ligation' in passed)),
    },

    'Homologous recombination (HR)': {
        'steps': {
            'Step1_Synapsis':         {'type': 'OR',
                                       'kos': ['K03553',               # recA
                                               'K04483',               # radA  (NOT K04485)
                                               'K04482']},             # RAD51
            # Three alternative end-processing systems. Each group is required
            # in full: a single RecFOR component on its own is not a pathway.
            'Step2_End_Processing':   {'type': 'GROUPS',
                                       'groups': [['K03582', 'K03583'],           # recB + recC
                                                  ['K03629', 'K03584', 'K06187'], # recF + recO + recR
                                                  ['K10865', 'K10866']]},         # mre11 + rad50
            # RuvAB has no archaeal counterpart with a KEGG orthologue, so this
            # step is bacteria-only and is excluded from the diagnostic.
            'Step3_Branch_Migration': {'type': 'AND',
                                       'kos': ['K03550', 'K03551']},   # ruvA + ruvB
            'Step4_Resolution':       {'type': 'OR',
                                       'kos': ['K01159',               # ruvC
                                               'K03655',               # recG
                                               'K03552']},             # hjc/hje
        },
        # Diagnostic: strand invasion, and resolution of the joint molecule.
        'binary_pass': lambda passed: ('Step1_Synapsis' in passed)
                                      and ('Step4_Resolution' in passed),
    },
}


# =============================================================================
# Scoring
# =============================================================================
def step_kos(rule):
    """Every KO named by a step, whatever its logic type."""
    if rule['type'] == 'GROUPS':
        return [k for group in rule['groups'] for k in group]
    return list(rule['kos'])


def step_satisfied(rule, genome_kos):
    t = rule['type']
    if t == 'AND':
        return all(k in genome_kos for k in rule['kos'])
    if t == 'GROUPS':
        return any(all(k in genome_kos for k in group) for group in rule['groups'])
    if t == 'OR':
        return any(k in genome_kos for k in rule['kos'])
    raise ValueError("unknown step type: %r" % t)


def genome_sort_key(name):
    """Sort MAG2 before MAG10 when names end in a number."""
    digits = ''.join(c for c in str(name) if c.isdigit())
    return (0, int(digits)) if digits else (1, str(name))


def load_annotations(path):
    if not os.path.exists(path):
        sys.exit("annotations file not found: %s\n"
                 "Point at it with -i, e.g.  -i DRAM_annotation/annotations.tsv" % path)
    df = pd.read_csv(path, sep='\t', low_memory=False)
    for col in ('fasta', 'ko_id'):
        if col not in df.columns:
            sys.exit("column '%s' missing from %s — is this a DRAM annotations.tsv?"
                     % (col, path))
    genome_kos = {}
    for genome, ko_field in zip(df['fasta'], df['ko_id']):
        g = str(genome).replace('.fasta', '').replace('.fa', '').strip()
        if not g or g == 'nan':
            continue
        target = genome_kos.setdefault(g, set())
        if pd.isna(ko_field):
            continue
        for k in str(ko_field).split(','):
            k = k.strip()
            if k:
                target.add(k)
    if not genome_kos:
        sys.exit("no genomes parsed from %s" % path)
    return genome_kos


def main():
    ap = argparse.ArgumentParser(
        description="DNA repair pathway completeness from DRAM annotations.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument('-i', '--annotations', default='DRAM_annotation/annotations.tsv',
                    help="annotations.tsv written by DRAM.py annotate")
    ap.add_argument('-o', '--outdir', default='.', help="directory for the three output files")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    print("Reading %s ..." % args.annotations)
    genome_kos = load_annotations(args.annotations)
    genomes = sorted(genome_kos, key=genome_sort_key)
    print("  %d genomes, %d distinct KOs" % (len(genomes), len(set().union(*genome_kos.values()))))

    summary_rows, step_rows = [], []
    for g in genomes:
        kos = genome_kos[g]
        row = {'genome': g}
        for pathway, spec in REPAIR_RULES.items():
            steps = spec['steps']
            passed = [s for s, r in steps.items() if step_satisfied(r, kos)]
            row[pathway] = round(len(passed) / len(steps), 4)
            row['%s (Pass)' % pathway] = spec['binary_pass'](passed)
            for s, r in steps.items():
                step_rows.append({
                    'genome': g, 'pathway': pathway, 'step': s,
                    'logic': r['type'], 'present': s in passed,
                    'KOs_defined': ';'.join(step_kos(r)),
                    'KOs_found': ';'.join(k for k in step_kos(r) if k in kos),
                })
        summary_rows.append(row)

    summary = pd.DataFrame(summary_rows)
    steps_df = pd.DataFrame(step_rows)
    f_sum = os.path.join(args.outdir, 'repair_completeness.tsv')
    f_stp = os.path.join(args.outdir, 'repair_steps.tsv')
    f_dia = os.path.join(args.outdir, 'repair_diagnostics.txt')
    summary.to_csv(f_sum, sep='\t', index=False)
    steps_df.to_csv(f_stp, sep='\t', index=False)

    # ---- diagnostics: what this dataset's annotation could not see ----------
    all_kos = set().union(*genome_kos.values())
    lines = ["DIAGNOSTICS — %s" % args.annotations,
             "%d genomes scored.\n" % len(genomes),
             "KOs named by the rules but never assigned to any genome.",
             "These contribute nothing to the scores. If one is the only",
             "counterpart a taxonomic group has for its step, that group cannot",
             "reach full completeness and the step must not be read as loss.\n"]
    missing = False
    for pathway, spec in REPAIR_RULES.items():
        for s, r in spec['steps'].items():
            absent = [k for k in step_kos(r) if k not in all_kos]
            if absent:
                missing = True
                lines.append("  %-34s %-30s %s" % (pathway, s, ' '.join(absent)))
    if not missing:
        lines.append("  (none — every KO in the rules was assigned at least once)")

    lines += ["", "Steps never satisfied by ANY genome.",
              "Such a step caps every genome's fraction without carrying",
              "information. Check whether its KOs exist in your lineages at all",
              "before keeping it.\n"]
    dead = False
    for pathway, spec in REPAIR_RULES.items():
        for s, r in spec['steps'].items():
            if not any(step_satisfied(r, genome_kos[g]) for g in genomes):
                dead = True
                lines.append("  %-34s %s" % (pathway, s))
    if not dead:
        lines.append("  (none)")

    lines += ["", "Steps satisfied by every genome.",
              "These cannot discriminate between genomes. Verify the step is not",
              "being carried by a broadly essential accessory protein.\n"]
    flat = False
    for pathway, spec in REPAIR_RULES.items():
        for s, r in spec['steps'].items():
            if all(step_satisfied(r, genome_kos[g]) for g in genomes):
                flat = True
                lines.append("  %-34s %s" % (pathway, s))
    if not flat:
        lines.append("  (none)")

    with open(f_dia, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')

    # ---- console summary ---------------------------------------------------
    print("\nCompleteness (fraction of steps recovered; * = diagnostic Pass)")
    header = "%-14s" % 'genome'
    order = list(REPAIR_RULES)
    for p in order:
        header += "%14s" % p.split('(')[-1].rstrip(')')
    print(header)
    for row in summary_rows:
        line = "%-14s" % row['genome']
        for p in order:
            line += "%13s%s" % (row[p], '*' if row['%s (Pass)' % p] else ' ')
        print(line)

    print("\nWrote:\n  %s\n  %s\n  %s" % (f_sum, f_stp, f_dia))
    print("\nRead repair_diagnostics.txt before interpreting these numbers.")


if __name__ == '__main__':
    main()
