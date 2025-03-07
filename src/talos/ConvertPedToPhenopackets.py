"""
script for converting the 'extended PED format' to a regular PED file and Phenopackets file
"""

from argparse import ArgumentParser
from os import makedirs
from os.path import exists
from pathlib import Path

import phenopackets.schema.v2 as pps2
from google.protobuf.json_format import MessageToJson
from peds import open_ped


# map the integer reported sex values to the enum
reported_sex_map = {1: pps2.Sex.MALE, 2: pps2.Sex.FEMALE}


def cli_main():
    parser = ArgumentParser()
    parser.add_argument('--input', help='path to the extended PED file')
    parser.add_argument('--output', help='stem path to write phenopacket and new PED file')
    args = parser.parse_args()
    main(ped_file=args.input, output=args.output)


def main(ped_file: str, output: str) -> None:
    """
    Create a Phenopacket from a pedigree file

    Args:
        ped_file (str): path to the extended PED file
        output (str): stem path to write phenopacket and new PED file
    """
    extended_ped = open_ped(ped_file)
    minimal_ped: list[str] = []

    # create the cohort shell
    cohort = pps2.Cohort(
        id='unknown',
        description=f'Phenotypic data from {ped_file}',
        members=[],
        meta_data=pps2.MetaData(
            created_by='Talos',
            resources=[
                pps2.Resource(
                    id='hp',
                    name='Human Phenotype Ontology',
                    url='http://www.human-phenotype-ontology.org',
                    version='2024-08-13',
                    namespace_prefix='HP',
                    iri_prefix='http://purl.obolibrary.org/obo/HP_',
                ),
            ],
        ),
    )

    for family in extended_ped:
        for participant in family:
            # create a Phenopacket for this individual, and add it to the cohort
            cohort.members.append(
                pps2.Phenopacket(
                    # primary ID is the internal CPG ID
                    id=participant.id,
                    subject=pps2.Individual(
                        # ID is the Family ID; Cohort model doesn't explicitly accommodate Family ID, but we need it
                        id=participant.family,
                        # alternate_ids captures the external ID, and is optional
                        alternate_ids=[participant.data[0]] if participant.data else None,
                        sex=reported_sex_map.get(participant.sex, pps2.Sex.UNKNOWN_SEX),
                    ),
                    # if HPO terms were provided, add them
                    phenotypic_features=[
                        pps2.PhenotypicFeature(type=pps2.OntologyClass(id=hpo)) for hpo in participant.data[1:]
                    ],
                ),
            )
            # save the minimal pedigree entry
            participant.data = []
            minimal_ped.append(str(participant))

    output_dir = Path(output).parent
    # make a parent directory for the output, if required
    if not exists(output_dir):
        makedirs(output_dir)

    # save the phenopacket JSON
    with open(f'{output}_phenopackets.json', 'w', encoding='utf-8') as handle:
        handle.write(MessageToJson(cohort))

    # save the pedigree
    with open(f'{output}_pedigree.ped', 'w', encoding='utf-8') as handle:
        handle.write(''.join(minimal_ped))


if __name__ == '__main__':
    cli_main()
