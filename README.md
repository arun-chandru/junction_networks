# Pólya inequalities for compact planar junction networks

**R. Arun Chandru**  
2 October 2026

## Abstract

We prove the Dirichlet Pólya inequality for every finite network formed from unit square chambers centered on the square lattice of spacing three and centered passages of length two and width one half along any selected nearest-neighbour edges. All six local chamber types are included, in particular four-port junctions. There is no periodicity, connectivity, density, or size condition on the selected finite graph. The bound holds for every eigenvalue, counted with multiplicity, with the actual physical area of the domain. Fully occupied grids give connected non-tiling examples with arbitrarily many holes and area-normalized diameter below the square root of six. The proof combines mixed-boundary comparisons, quantitative boundary-moment estimates with complete infinite-complement bounds, a phase-uniform fractional-part estimate, and exact finite arithmetic certificates. A complementary result allows independent port displacements along rectangular chamber chains.

## An accessible account

Imagine a drum membrane made from square rooms connected by short rectangular passages. The membrane is fixed along every exposed boundary, including the boundaries of its holes. Its vibration frequencies depend on the geometry. The eigenvalues studied in the paper are the squares of those frequencies, after the usual physical constants have been normalized.

At high frequencies, the area of a drum determines the leading number of vibration modes. Pólya's conjecture asks whether that leading area expression is already an upper bound at every frequency. This is more demanding than a statement that becomes correct only at sufficiently high frequencies, or a bound on an average of many eigenvalues.

The paper establishes that every-frequency conclusion for a precisely specified family. Place unit square rooms at any finite selection of sites three units apart on a square grid, and choose any of the possible connections between adjacent rooms. A chosen connection is a passage two units long and one-half unit wide. Missing connections, branches, four-way junctions, cycles, and disconnected pieces are all allowed. If there are n rooms and m passages, the physical area is exactly n+m. The kth Dirichlet eigenvalue is at least 4 pi k/(n+m). The number of rooms is not bounded by the verification: the same small catalogue of local comparisons proves the result for every finite network in this family.

The difficult part is that cutting a network into rooms and passages changes its boundary conditions and can create too many low-frequency modes in the comparison. The proof retains the fixed boundary segments around each opening to control these extra modes. Boundary integrals give quantitative constraints on whole groups of eigenfunctions, including repeated eigenvalues. Every omitted high-frequency component is bounded analytically; a finite matrix is not mistaken for the complete operator. A separate argument controls the entire high-energy tail, and exact arithmetic checks the remaining bounded range, including the endpoints at which a chosen counting bound changes.

The computations use integer and rational arithmetic with outward-rounded intervals. They do not treat approximate eigenvalues from a numerical solver as rigorous lower bounds. The accompanying programs verify the finite arithmetic parts of the argument; the variational, trace, and counting arguments are proved in the manuscript. This is not a formal proof in a proof assistant or a claim of completed external peer review.

The manuscript also gives explicit nonuniform extensions. In coordinates where a room has side two, a width-one passage can have length four or any length at least seven. More generally, a passage of width w and length L is allowed whenever 0 < w <= 1, L >= 9, and Lw >= 4. Its window may be displaced from the face midpoint by up to (1-w)/2. Different incident edges can have different widths, lengths and displacements, subject to the actual geometric requirements of straight connections and no unintended intersections. A concrete nonuniform rectangular-grid construction realizes these choices. The physical area is exactly 4 times the room count plus the sum of the passage areas Lw; no extra terminal rooms are needed.

For example, with side-two rooms and passages of length nine, widths may vary independently between 4/9 and 1. A width-one-half window may shift by one-quarter unit. In the original side-one scale, that is a displacement of one-eighth of the room side. Another admitted range keeps each passage area equal to four while letting its width tend to zero and its length grow. The bounded normalized-diameter statement belongs to the original compact grid family, not to all these length-varying extensions.

The geometric restrictions matter. The result does not cover all planar drums or arbitrary width perturbations of the original shortest passages. A further theorem treats a different chain geometry where independent displacements are allowed. None of these results is presented as a proof of the unrestricted Pólya conjecture.

## Files

- `r_arun_chandru_Pólya inequalities for compact planar junction networks.pdf`: complete paper, figures, references, and exact-certificate appendices.
- `manuscript.tex`: the entire self-contained LaTeX source. All figures are drawn within this file using TikZ; no external images or bibliography database are needed.
- `verify_plj.py`: entry point for the complete arithmetic verification, including the lower-degree cells and displaced chains.
- `verify_degree_four_tail.py`: centered-block integral bounds and the complete four-port separated-event calculation.
- `verify_degree_four_boundary.py`: outward-rounded small-matrix boundary certificates and coverage of the corrected finite energy range.
- `verify_degree_four_odd33.py`: the additional exact odd-parity boundary certificate needed at one four-port threshold.
- `verify_nonuniform.py`: exact scalar constants for the analytic variable-width, variable-length, and displaced-window extension.
- `verification_results.txt`: output of a successful complete verification run on this release.

Keep the five Python files together. No research notes, external datasets, downloaded spectra, or unpublished manuscripts are inputs to the proof or the programs.

## Compile the manuscript

A current TeX Live or MiKTeX installation with the standard AMS, Latin Modern, geometry, TikZ/PGF, microtype, xurl, and hyperref packages is sufficient. From this directory, run:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error manuscript.tex
```

Alternatively, run `pdflatex manuscript.tex` repeatedly until references and the table of contents are stable (normally two or three passes). Tectonic also supports the source:

```sh
tectonic manuscript.tex
```

There is no BibTeX or Biber step. LaTeX shell escape is not required. For an arXiv source submission, upload `manuscript.tex`; the Python files and this README may be retained with the repository/verification supplement according to the submission service's rules.

## Run all exact checks

Use Python 3.8 or later; only its standard library is required:

```sh
python -B -O verify_plj.py
```

The `-B` flag prevents bytecode-cache files. The `-O` flag deliberately demonstrates that acceptance does not depend on removable Python assertions: every required comparison raises an exception on failure. A successful complete run finishes with an `all_certificates` record reporting `passed: true`. Any exception or nonzero exit status is a failed or incomplete verification, not an accepted result.

The largest stage uses two unsigned histograms, about 128 MiB together, plus the Python interpreter and modest overhead. Stages run sequentially and use one process; no GPU, network connection, paid service, or external numerical library is needed. Runtime depends on the machine. The finite limits are part of the proof, so shortening a scan does not reproduce the full certificate.

The scripts verify all listed spectral-event multiplicities, artificial reset events, rational constants, integral lower sums, interval matrix pivots, and the complete finite/infinite handoff. Their analytic interpretation is specified in the paper's appendices. Floating exploratory calculations are not included in the release or used by the acceptance checks.
