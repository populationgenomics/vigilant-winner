"""
Participant HPO ~ panel matching stage
Designed to operate on an extended pedigree format (see docs)
- read participants and assc. HPOs
- query for panels and HPOs
- write a new file containing participant-panel matches
"""

import re
from argparse import ArgumentParser
from collections import defaultdict

import phenopackets.schema.v2 as pps2
from google.protobuf.json_format import ParseDict
from networkx import dfs_successors
from obonet import read_obo

from talos.config import config_retrieve
from talos.models import ParticipantHPOPanels, PhenoPacketHpo, PhenotypeMatchedPanels
from talos.static_values import get_logger
from talos.utils import get_json_response, read_json_from_path

HPO_RE = re.compile(r'HP:\d+')

PANELAPP_HARD_CODED_DEFAULT = 'https://panelapp-aus.org/api/v1/panels'
PANELAPP_HARD_CODED_BASE_PANEL = 137

try:
    PANELS_ENDPOINT = config_retrieve(['GeneratePanelData', 'panelapp'], PANELAPP_HARD_CODED_DEFAULT)
    DEFAULT_PANEL = config_retrieve(['GeneratePanelData', 'default_panel'], 137)
except KeyError:
    get_logger(__file__).warning('Config environment variable TALOS_CONFIG not set, falling back to Aussie PanelApp')
    PANELS_ENDPOINT = PANELAPP_HARD_CODED_DEFAULT
    DEFAULT_PANEL = PANELAPP_HARD_CODED_BASE_PANEL


def get_panels(endpoint: str = PANELS_ENDPOINT) -> dict[str, set[int]]:
    """
    query panelapp, and collect panels by HPO term

    Args:
        endpoint (str): URL for panels

    Returns:
        dict: {HPO_Term: [panel_ID, panel_ID],}
    """

    panels_by_hpo = defaultdict(set)

    while True:
        endpoint_data = get_json_response(endpoint)
        for panel in endpoint_data['results']:
            # can be split over multiple strings
            relevant_disorders = ' '.join(panel['relevant_disorders'] or [])
            for match in re.findall(HPO_RE, relevant_disorders):
                panels_by_hpo[match].add(int(panel['id']))

        # cycle through additional pages
        if endpoint := endpoint_data['next']:
            continue
        break

    return dict(panels_by_hpo)


def set_up_cohort_pmp(cohort: pps2.Cohort) -> tuple[PhenotypeMatchedPanels, set[str]]:
    """
    read the extended pedigree file, pull out family details and HPO terms

    Args:
        cohort (str): GA4GH Cohort/PhenoPacket file

    Returns:
        PhenotypeMatchedPanels object with Cohort details, set of all HPO terms
    """

    hpo_dict = PhenotypeMatchedPanels()
    all_hpos: set[str] = set()

    for member in cohort.members:
        hpo_dict.samples[member.id] = ParticipantHPOPanels(
            external_id=member.subject.alternate_ids[0] if member.subject.alternate_ids else member.id,
            family_id=member.subject.id,
            hpo_terms=[PhenoPacketHpo(id=hpo.type.id, label=hpo.type.label) for hpo in member.phenotypic_features],
            panels={DEFAULT_PANEL},
        )
        all_hpos.update({hp.type.id for hp in member.phenotypic_features})

    return hpo_dict, all_hpos


def match_hpos_to_panels(hpo_panel_map: dict[str, set[int]], hpo_file: str, all_hpos: set[str]) -> dict[str, set[int]]:
    """
    take the HPO terms from the participant metadata, and match to panels

    Args:
        hpo_panel_map (dict): panel IDs to all related panels
        hpo_file (str): path to an obo file containing HPO tree
        all_hpos (set[str]): collection of all unique HPO terms

    Returns:
        a dictionary linking all HPO terms to a corresponding set of Panel IDs
    """

    hpo_graph = read_obo(hpo_file, ignore_obsolete=False)

    hpo_to_panels = defaultdict(set)
    for hpo in all_hpos:
        # identify all HPO terms back to the ontology root
        successor_hpo_terms = set(dfs_successors(hpo_graph, hpo))

        for hpo_term in successor_hpo_terms:
            if hpo_term in hpo_panel_map:
                hpo_to_panels[hpo].update(hpo_panel_map[hpo_term])

    return hpo_to_panels


def match_participants_to_panels(participant_hpos: PhenotypeMatchedPanels, hpo_panels: dict):
    """
    take the two maps of Participants: HPOs, and HPO: Panels
    blend the two to find panels per participant

    For each participant, find any HPO terms which were matched to panels
    for each matched term, add the panel(s) to the participant's private set

    Args:
        participant_hpos (PhenotypeMatchedPanels): CPG ID to phenotype details
        hpo_panels (dict): lookup of panels per HPO term
    """

    for party_data in participant_hpos.samples.values():
        for hpo_term in party_data.hpo_terms:
            if panel_list := hpo_panels.get(hpo_term.id):
                # add relevant panels for this participant
                party_data.panels.update(panel_list)
                # and add to the collection of all panels
                participant_hpos.all_panels.update(panel_list)


def cli_main():
    get_logger(__file__).info('Starting HPO~Panel matching')
    parser = ArgumentParser()
    parser.add_argument('--input', help='GA4GH Cohort/PhenoPacket File')
    parser.add_argument('--output', help='Path to write PhenotypeMatchedPanels to (JSON)')
    parser.add_argument('--hpo', help='Local copy of HPO obo file', required=True)
    args = parser.parse_args()
    main(ga4gh_cohort_file=args.input, panel_out=args.output, hpo_file=args.hpo)


def main(ga4gh_cohort_file: str, panel_out: str, hpo_file: str):
    """
    query PanelApp - get all panels and their assc. HPOs
    read Cohort/PhenoPacket file
    associate each participant with panels
    write a PhenotypeMatchedPanels instance to a local file

    Args:
        ga4gh_cohort_file (str): path to GA4GH Cohort/PhenoPacket file
        panel_out (str): where to write PhenotypeMatchedPanels file
        hpo_file (str): path to an obo file containing HPO tree
    """

    # read the Cohort/PhenoPacket file as a JSON
    ga4gh_cohort = ParseDict(read_json_from_path(ga4gh_cohort_file), pps2.Cohort())

    # query PanelApp to get all PanelIDs, and their relevant HPO terms
    panels_by_hpo = get_panels()

    # build the PhenotypeMatchedPanels object from the Cohort data
    # also collect each unique HPO term in the Cohort
    pmp_dict, all_hpos = set_up_cohort_pmp(cohort=ga4gh_cohort)

    # match HPO terms to panel IDs
    # returns a lookup of each HPO term in the cohort, and panels it is associated with
    hpo_to_panels = match_hpos_to_panels(hpo_panel_map=panels_by_hpo, all_hpos=all_hpos, hpo_file=hpo_file)

    match_participants_to_panels(pmp_dict, hpo_to_panels)

    # validate the object
    valid_pheno_dict = PhenotypeMatchedPanels.model_validate(pmp_dict)

    # validate and write using pydantic
    if panel_out:
        with open(panel_out, 'w', encoding='utf-8') as handle:
            handle.write(valid_pheno_dict.model_dump_json(indent=4))


if __name__ == '__main__':
    cli_main()
