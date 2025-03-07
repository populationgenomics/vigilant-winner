"""
tests for the PanelApp parser
"""

from copy import deepcopy

from talos.models import PanelApp, PanelDetail
from talos.QueryPanelapp import get_best_moi, get_panel, parse_panel_activity

empty_gene_dict = PanelApp(genes={})


def test_activity_parser(panel_activities):
    """
    check that we correctly parse the activities JSON
    """

    activity_dict = parse_panel_activity(panel_activities)

    # this should be absent
    assert 'NOT_GENE' not in activity_dict

    assert 'GENE1' in activity_dict
    assert activity_dict['GENE1'].strftime('%Y-%m-%d') == '2022-02-01'

    assert 'GENE2' in activity_dict
    assert activity_dict['GENE2'].strftime('%Y-%m-%d') == '2024-04-25'

    assert 'GENE3' in activity_dict
    assert activity_dict['GENE3'].strftime('%Y-%m-%d') == '2023-08-15'


def test_panel_query(httpx_mock, latest_mendeliome):
    """check that the default parsing delivers correct data"""

    httpx_mock.add_response(url='https://panelapp-aus.org/api/v1/panels/137/', json=latest_mendeliome)
    httpx_mock.add_response(url='https://panelapp-aus.org/api/v1/panels/137/activities/', json=[])

    gd = deepcopy(empty_gene_dict)
    get_panel(gd, forbidden_genes=set())
    assert gd.genes['ENSG00ABCD'].all_moi == {'biallelic'}
    assert gd.genes['ENSG00ABCD'].panels == {137}
    assert gd.genes['ENSG00EFGH'].all_moi == {'monoallelic'}


def test_panel_query_removal(httpx_mock, latest_mendeliome):
    """check that the default parsing delivers correct data"""

    httpx_mock.add_response(url='https://panelapp-aus.org/api/v1/panels/137/', json=latest_mendeliome)
    httpx_mock.add_response(url='https://panelapp-aus.org/api/v1/panels/137/activities/', json=[])

    gd = deepcopy(empty_gene_dict)
    get_panel(gd, blacklist=['ENSG00EFGH'])
    assert gd.genes['ENSG00ABCD'].all_moi == {'biallelic'}
    assert gd.genes['ENSG00ABCD'].panels == {137}
    assert 'ENSG00EFGH' not in gd.genes


def test_panel_query_forbidden(httpx_mock, latest_mendeliome):
    """check that the default parsing delivers correct data"""

    httpx_mock.add_response(url='https://panelapp-aus.org/api/v1/panels/137/', json=latest_mendeliome)
    httpx_mock.add_response(url='https://panelapp-aus.org/api/v1/panels/137/activities/', json=[])

    gd = deepcopy(empty_gene_dict)
    get_panel(gd, forbidden_genes={'ENSG00EFGH'})
    assert gd.genes['ENSG00ABCD'].all_moi == {'biallelic'}
    assert gd.genes['ENSG00ABCD'].all_moi == {'biallelic'}
    assert gd.genes['ENSG00ABCD'].panels == {137}
    assert 'ENSG00EFGH' not in gd.genes


def test_panel_query_removal_2(httpx_mock, latest_mendeliome):
    """check skipping by symbol works as well. This test has a weirdly long runtime"""

    httpx_mock.add_response(url='https://panelapp-aus.org/api/v1/panels/137/', json=latest_mendeliome)
    httpx_mock.add_response(url='https://panelapp-aus.org/api/v1/panels/137/activities/', json=[])

    gd = deepcopy(empty_gene_dict)
    get_panel(gd, blacklist=['EFGH'])
    assert gd.genes['ENSG00ABCD'].all_moi == {'biallelic'}
    assert gd.genes['ENSG00ABCD'].panels == {137}
    assert 'ENSG00EFGH' not in gd.genes


def test_panel_query_forbidden_2(latest_mendeliome, httpx_mock):
    """check skipping by symbol works as well"""

    httpx_mock.add_response(url='https://panelapp-aus.org/api/v1/panels/137/', json=latest_mendeliome)
    httpx_mock.add_response(url='https://panelapp-aus.org/api/v1/panels/137/activities/', json=[])

    gd = deepcopy(empty_gene_dict)
    get_panel(gd, forbidden_genes={'EFGH'})
    assert gd.genes['ENSG00ABCD'].all_moi == {'biallelic'}
    assert gd.genes['ENSG00ABCD'].panels == {137}
    assert 'ENSG00EFGH' not in gd.genes


def test_panel_query_addition(latest_incidentalome, httpx_mock):
    """
    check that the default parsing delivers correct data
    oof, this was a tricky one
    """

    # assumed data we already gathered
    gd = PanelApp(
        metadata=[{'version': '0.11088', 'name': 'Mendeliome', 'id': 137}],
        genes={
            'ENSG00ABCD': {
                'symbol': 'ABCD',
                'all_moi': {'monoallelic'},
                'new': [],
                'panels': [137],
            },
            'ENSG00IJKL': {
                'symbol': 'IJKL',
                'all_moi': {'both'},
                'new': [137],
                'panels': [123, 137],
            },
        },
    )
    httpx_mock.add_response(url='https://panelapp-aus.org/api/v1/panels/126/', json=latest_incidentalome)
    httpx_mock.add_response(url='https://panelapp-aus.org/api/v1/panels/126/activities/', json=[])

    # should query for and integrate the incidentalome content
    get_panel(gd, panel_id=126)
    assert gd.genes['ENSG00ABCD'].all_moi == {'monoallelic', 'biallelic'}
    assert gd.genes['ENSG00ABCD'].panels == {137, 126}
    assert gd.genes['ENSG00IJKL'].all_moi == {'both'}
    assert gd.genes['ENSG00IJKL'].panels == {123, 137}
    assert 'ENSG00EFGH' not in gd.genes


def test_get_best_moi_empty():
    """
    check that the MOI summary works
    """

    d = {'ensg1': PanelDetail(all_moi=set(), chrom='1', symbol='ensg1')}
    get_best_moi(d)
    assert d['ensg1'].moi == 'Biallelic'

    d = {'ensg1': PanelDetail(all_moi=set(), chrom='X', symbol='ensgX')}
    get_best_moi(d)
    assert d['ensg1'].moi == 'Hemi_Bi_In_Female'


def test_get_best_moi_mono():
    """
    check that the MOI summary works
    """

    d = {'ensg1': PanelDetail(all_moi={'monoallelic'}, chrom='1', symbol='ensg1')}
    get_best_moi(d)
    assert d['ensg1'].moi == 'Monoallelic'


def test_get_best_moi_mono_and_biallelic():
    """
    check that the MOI summary works
    """

    d = {'ensg1': PanelDetail(all_moi={'monoallelic', 'biallelic'}, chrom='1', symbol='ensg1')}
    get_best_moi(d)
    assert d['ensg1'].moi == 'Mono_And_Biallelic'


def test_get_best_moi_1():
    """
    check that the MOI summary works
    """

    d = {
        'ensg1': PanelDetail(
            all_moi={'Monoallelic', 'Biallelic', 'both'},
            chrom='1',
            symbol='ensg1',
        ),
    }
    get_best_moi(d)
    assert d['ensg1'].moi == 'Mono_And_Biallelic'


def test_get_best_moi_x():
    """
    check that the MOI summary works
    """

    d = {
        'ensg1': PanelDetail(
            all_moi={'x-linked biallelic', 'x-linked'},
            chrom='X',
            symbol='ensgX',
        ),
    }
    get_best_moi(d)
    assert d['ensg1'].moi == 'Hemi_Mono_In_Female'
