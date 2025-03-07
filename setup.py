"""
setup.py for the talos package
"""

from setuptools import find_packages, setup


with open('README.md', encoding='utf-8') as handle:
    readme = handle.read()


def read_reqs(filename: str) -> list[str]:
    """
    Read requirements from a file, return as a list
    TODO eventually split out the requirements vs. the setup content

    Args:
        filename (str): the requirements file to parse

    Returns:
        list[str]: the requirements
    """
    with open(filename, encoding='utf-8') as filehandler:
        return [line.strip() for line in filehandler if line.strip() and not line.startswith('#')]


setup(
    name='talos',
    description='Centre for Population Genomics Variant Prioritisation',
    long_description=readme,
    version='6.6.3',
    author='Matthew Welland, CPG',
    author_email='matthew.welland@populationgenomics.org.au, cas.simons@populationgenomics.org.au',
    package_data={'talos': ['templates/*.jinja', 'example_config.toml']},
    url='https://github.com/populationgenomics/talos',
    license='MIT',
    classifiers=[
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    packages=find_packages('src'),
    package_dir={'': 'src'},
    include_package_data=True,
    install_requires=read_reqs('requirements.txt'),
    extras_require={
        'test': read_reqs('requirements-dev.txt'),
        'cpg': read_reqs('requirements-cpg.txt'),
    },
    entry_points={
        'console_scripts': [
            # ----------------------------------------------------------------------------------------
            # these scripts are used during the annotation & data formatting pre-workflow
            # ----------------------------------------------------------------------------------------
            # Parse an Ensembl GFF3 file, output a BED file of gene regions
            'CreateRoiFromGff3 = talos.annotation_scripts.CreateRoiFromGff3:cli_main',
            # Parse AlphaMissense into Hail Table
            'ParseAlphaMissenseIntoHt = talos.annotation_scripts.ParseAlphaMissenseIntoHt:cli_main',
            # Parse MANE summary file into JSON
            'ParseManeIntoJson = talos.annotation_scripts.ParseManeIntoJson:cli_main',
            # Reformats an annotated VCF to a MatrixTable, adding additional annotations
            'ReformatVcfToMt = talos.annotation_scripts.ReformatVcfToMt:cli_main',
            # ----------------------------------------------------------------------------------------
            # these scripts are specific to internal CPG functions
            # ----------------------------------------------------------------------------------------
            # Scans Metamist for published reports, collects into an index page
            'BuildReportIndexPage = talos.cpg_internal_scripts.BuildReportIndexPage:main',
            # Create cohort phenopackets from Metamist
            'MakePhenopackets = talos.cpg_internal_scripts.MakePhenopackets:cli_main',
            # Parse the summary JSON, generating a file for ingestion by Seqr
            'MinimiseOutputForSeqr = talos.cpg_internal_scripts.MinimiseOutputForSeqr:cli_main',
            # ----------------------------------------------------------------------------------------
            # these scripts are part of the core Talos workflow
            # ----------------------------------------------------------------------------------------
            # Combine multiple separate Exomiser result files into a single Hail Table/JSON
            'AggregateExomiserVariantTsvs = talos.AggregateExomiserVariantTsvs:cli_main',
            # turns the SVDB TSV into a Hail Table
            'ConvertSpliceVarDb = talos.ConvertSpliceVarDb:cli_main',
            # Builds cohort phenopackets from a pedigree with HPO terms
            'ConvertPedToPhenopackets = talos.ConvertPedToPhenopackets:cli_main',
            # use the HPO terms to select panels for this analysis
            'GeneratePanelData = talos.GeneratePanelData:cli_main',
            # query PanelApp for those selected panels
            'QueryPanelapp = talos.QueryPanelapp:cli_main',
            # use API queries to find the gene symbol for each gene ID
            'FindGeneSymbolMap = talos.FindGeneSymbolMap:cli_main',
            # Filter and label a small-variant MatrixTable
            'RunHailFiltering = talos.RunHailFiltering:cli_main',
            # Filter and label a SV MatrixTable
            'RunHailFilteringSV = talos.RunHailFilteringSV:cli_main',
            # Run each of the category-labelled variants through MOI filters
            'ValidateMOI = talos.ValidateMOI:cli_main',
            # catch variants which have strong phenotypic matches
            'HPOFlagging = talos.HPOFlagging:cli_main',
            # Generate a HTML file from the report summary
            'CreateTalosHTML = talos.CreateTalosHTML:cli_main',
            # Generate a brief summary of an existing report (variant count, family size, etc.)
            'SummariseReport = talos.SummariseReport:cli_main',
        ],
    },
)
