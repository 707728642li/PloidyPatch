# Public repository boundary

PloidyPatch keeps the software distribution separate from unpublished research
materials. The public software repository may contain source code, tests,
tutorials, small synthetic examples, packaging recipes, and software-release
documentation. It must not contain manuscript drafts, journal-submission
records, unpublished figures, biological result trees, reviewer packets, or
author-administration files.

The following workspace paths are intentionally excluded from Git and from
release archives:

- `manuscript/`
- `results/`
- manuscript-, submission-, reviewer-, and paper-figure build artifacts
- private author and journal-submission metadata

Before changing repository visibility to public, verify all of the following:

1. `git status --ignored` shows prepublication material only as ignored local
   content.
2. `git ls-tree -r --name-only HEAD` contains no `manuscript/` or `results/`
   entries and no private submission metadata.
3. The wheel and source archive contain only the intended software payload.
4. Repository history and release attachments have been audited separately;
   deleting a path from the current tree does not remove it from older commits.
5. PyPI and Bioconda publication is performed only after the repository and
   release URLs are intentionally public and stable.

The local scientific workspace is authoritative for unpublished research
materials. Removing those paths from Git tracking must never delete the local
files.
