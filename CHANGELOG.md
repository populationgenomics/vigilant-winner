# Changelog

All notable changes to this project will be documented in this file.

Suggested headings per release (as appropriate) are:

* `Added` for new features.
* `Changed` for changes in existing functionality.
* `Deprecated` for soon-to-be removed features.
* `Removed` for now removed features.
* `Fixed` for any bug fixes.
* `Security` in case of vulnerabilities.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

[6.0.0] - 2024-10-21

### Added

* MakePhenopackets stage, which generates a Phenopacket file from Metamist

### Changed

* Renamed `GenerateSeqrFile` to `MinimiseOutputForSeqr`
* Renamed `vcf_to_mt` to `VcfToMt`
* The 'extended PED' file format is no longer used - instead we use a normal 6-col PED file, and a separate file with
  phenotypic data and external IDs for each participant in a GA4GH-compliant Cohort/Phenopacket format
  * An adapter remains in [talos/CPG](src/talos/ConvertPedToPhenopackets.py) to convert the old format to the new format

### Deprecated

* GeneratePED has been deleted

[5.0.0] - 2024-07-08

### Changed

* Throwing all the source code in a `/src` folder
* Principal scripts are renamed, and accessible as command-line entrypoints

1. GeneratePED
2. GeneratePanelData
3. QueryPanelapp
4. RunHailFiltering
5. RunHailFilteringSV
6. ValidateMOI
7. CreateTalosHTML

* extensive removing of the pre-commit exclusions, and consequent reformatting to pass new tests
* documentation updated

[4.0.0] - 2024-06-09

### Changed

- Name AIP -> Talos
- Pre-commit tooling changed again (linting with Ruff instead of Black)
- Removes use of QoB, and relies on 'local' Hail Query operations
- Drops Peddy in favour of Peds for a more minimal install
- Splits metamist querying into a separate script, with all other downstream operations using an extended PED format (6
  cols + ext IDs + HPO IDs)

MOI testing

- Uses the altered pedigree representation for MOI
- Each variant is interpreted relative to the immediate Trio (i.e. a quad is interpreted as 2 separate trios, not as a
  single variant/variants solving all presentations in a family)

### Added

- An adapter to take a VEP-annotated VCF and restructure in Hail to match the VEP-from-JSON annotated version we
  natively expect from the seqr loader pipeline

### Removed

- Everything to do with VEP... This does not contain a vep installation, and no longer determines how/if VEP is to be
  run. Except for ClinVar reanalysis - still in the process of dog-fooding that from ClinvArbitration instead of
  reimplementing here
- Removes even an example of analysis outcome registration from this repository, not even we use that functionality...

[3.2.2] - 2024-04-12

### Changed

Huge re-working of the pre-commit tooling, resulting in reduced overall line count

[3.0.0] - 2023-12-08

This was a MASSIVE refactor affecting almost all parts of the codebase. Instead of characterising
data in various places using type hinting, I've now introducted Pydantic models to formally represent
the data used in each part of the application (see reanalysis/models.py). This has the benefit of
making the code more readable, and also allows for the use of Pydantic's validation and parsing.

A key note is that all serialised/deserialised data is now done through a Pydantic validation layer,
forcing the data to conform to the model. This may require some objects (e.g. state/history files)
to be reformatted, as during this work I identified some inconsistencies in the data.

### Added

* Data models for every aspect of the analysis
* Pydantic validation layer for all serialised/deserialised data
* Specifically, all variants inherit from VariantCommon, and Structural and Small Variants each have
  a fundamentally different object representation, with a common interface but some internal methods
  specific to each variant type. This is extensible (e.g. further inheritance for STRs)

### Changed

* All minimial-case representations of Variants or Reportable events used in unit tests are now
  represented as full Pydantic models for the corresponding data type
* So many little things, this was a +3000/-2000 line change

[2.1.0] - 2023-11-20

This change involved the introduction of SV data, and the co-processing of SV and Small Variant
data in a unified analysis process. The SV data format expected is based on the GATK-SV pipeline
and annotation.

This also includes the addition of the `sv1` category, which represents LoF structural variants

### Added

* SV data is now processed alongside small variants
* `sv1` category specific to SV data

### Changed

* Some of the AbstractVariant loading is modified to allow for SV data (e.g. defaults for depths)

[2.0.0 & 2.0.1] - 2023-10-26

This bump is less substantial than the version number suggests, but in hindsight the previous
version should probably have been a major bump due to wholesale changes in how input is provided.

### Added

* CategoryBoolean6 (AlphaMissense) is now active, and should be implemented in a backwards
  compatible way (no failure if annotation is absent)
* VCF export includes `am_class` and `am_pathogenicity` (defaulting to empty String)
* If these fields are not missing, they should render in the report's variant drop-down drawer

### Changed

* The Hail Labelling stage also takes a dataset argument (used to generate a dataset-specific
  path for the temporary files/checkpoints)
* Removed a second config file from the test framework (consolidate all fields into 1 file)
* MyPy and Ruff are implemented as replacement linting tools. Stricter, and a billion times faster

### Deprecated

* Flake8 and PyLint are removed as linting tools

[1.2.0] - 2023-10-16

### Changed

* a lot. Like... a lot
* All the Cohort-specific content in the configuration files has been removed from the main config file. For CPG usage
  we have some sensitivities around the open publication of the cohorts and details of the analysis carried out. As a
  precaution (although this config does not contain sample-level data) we have removed the cohort-specific content from
  this repository and moved it to a private repository.
* Structurally the cohort-specific region of the config file has changed. This is to better facilitate running of Talos
  through our main analysis pipeline, requiring less manual intervention to run exomes and genomes separately. The
  utils methods `get_cohort_config` and `get_cohort_seq_type_conf` have been updated to reflect this.
    1. `config.cohorts.COHORT` contains details general to a whole cohort - blacklisted genes, additional panels, solved
       cases
    2. `config.cohorts.COHORT.GENOME` contains details specific to the genome analysis of a cohort - the Seqr project
       name, mapping between individual and seqr family IDs, and any labelled variants for this group
    3. `config.cohorts.COHORT.EXOME` is the counterpart to `.GENOME`
* Hail Labelling/Filtering now takes an output path, to be more easily included in a Hail Batch workflow
* The run preparation and HPO-panel script now use GQL queries, which drops a few methods
* The external_lookup file is no longer used - the same Internal-External ID mapping is noq stored alongside phenotype
  data in the HPO panel file.
* `--plink` argument to interpretation_runner is renamed to `--pedigree`
* Methods to reduce a cohort by a fixed percentage have been deleted - was useful as a potential test/train set, but has
  not been used recently

[1.1.5] - 2023-08-28

### Changed

* reanalysis/interpretation_runner.py creates two reports; the first report is all standard content from the run, the
  second report is variants new within this latest run
* Relevant interfaces with the output file registration and report hunter scripts are altered
* The report hunter script is altered to create two separate index pages; one for the standard reports, one to only look
  for the latest-variant report

[1.1.4] - 2023-06-28

### Changed

* Updated to use the latest models in Metamist/sample metadata DB
* Queries updated to use GraphQL

[1.1.3] - 2023-06-22

### Added

* Adds a helper script which exports a file formatted specifically for ingestion into Seqr

[1.1.2] - 2023-05-31

### Changed

* Adds a number of changes to how ClinVar pathogenic variants are processed
* See [this PR](https://github.com/populationgenomics/talos/pull/305)

[1.1.1] - 2023-05-08

### Added

* See [this PR](https://github.com/populationgenomics/talos/pull/293)
* Adds the PM5 category, and the VEP content to create the relevant annotation data

[1.0.1] - 2023-04-01

### Added

* Added this changelog
* Allow for addition of non-Metamist metadata registrars, selectable in config
* Added Variant Tagging in report, adopted from MW@MSFT's fork - #271
* Add family ID into top level variant table

### Changed

* Separated out metadata registration into a separate file & job
* File registration is set to `always_run`, if indicated in config
* De Novo calling (hail method) can have `max_parent_ab` ratio overridden in config
* Establishes a method-contract for altering metadata registrars and HTML scripts
