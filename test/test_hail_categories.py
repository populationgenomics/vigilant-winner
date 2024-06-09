"""
unit testing collection for the hail MT methods
"""

import pandas as pd
import pytest

import hail as hl

from talos.hail_filter_and_label import (
    annotate_category_1,
    annotate_category_2,
    annotate_category_3,
    annotate_category_5,
    annotate_category_6,
    annotate_talos_clinvar,
    filter_to_categorised,
    filter_to_population_rare,
    green_and_new_from_panelapp,
    split_rows_by_gene_and_filter_to_green,
)
from talos.models import PanelApp
from test.test_utils import ONE_EXPECTED, TWO_EXPECTED, ZERO_EXPECTED

category_1_keys = ['locus', 'clinvar_talos_strong']
category_2_keys = ['locus', 'clinvar_talos', 'cadd', 'revel', 'geneIds', 'consequence_terms']
category_3_keys = ['locus', 'clinvar_talos', 'lof', 'consequence_terms']
hl_locus = hl.Locus(contig='chr1', position=1, reference_genome='GRCh38')


@pytest.mark.parametrize('value,classified', [(0, 0), (1, 1), (2, 0)])
def test_class_1_assignment(value, classified, make_a_mt):
    """
    use some fake annotations, apply to the single fake variant
    check that the classification process works as expected based
    on the provided annotations
    """
    anno_matrix = make_a_mt.annotate_rows(
        info=make_a_mt.info.annotate(
            clinvar_talos_strong=value,
        ),
    )

    anno_matrix = annotate_category_1(anno_matrix)
    assert anno_matrix.info.categoryboolean1.collect() == [classified]


@pytest.mark.parametrize(
    'clinvar_talos,c6,gene_id,consequence_terms,classified',
    [
        (0, 0, 'GREEN', 'missense', ZERO_EXPECTED),
        (1, 1, 'RED', 'frameshift_variant', ZERO_EXPECTED),
        (1, 0, 'GREEN', 'missense', ONE_EXPECTED),
        (1, 1, 'GREEN', 'frameshift_variant', ONE_EXPECTED),
        (0, 1, 'GREEN', 'synonymous', ONE_EXPECTED),
        (0, 1, 'GREEN', 'synonymous', ONE_EXPECTED),
    ],
)
def test_cat_2_assignment(clinvar_talos, c6, gene_id, consequence_terms, classified, make_a_mt):
    """
    use some fake annotations, apply to the single fake variant
    Args:
        clinvar_talos ():
        c6 ():
        gene_id ():
        consequence_terms ():
        classified ():
        make_a_mt ():
    """

    anno_matrix = make_a_mt.annotate_rows(
        geneIds=gene_id,
        info=make_a_mt.info.annotate(clinvar_talos=clinvar_talos, categoryboolean6=c6),
        vep=hl.Struct(
            transcript_consequences=hl.array([hl.Struct(consequence_terms=hl.set([consequence_terms]))]),
        ),
    )

    anno_matrix = annotate_category_2(anno_matrix, new_genes=hl.set(['GREEN']))
    assert anno_matrix.info.categoryboolean2.collect() == [classified]


@pytest.mark.parametrize(
    'clinvar_talos,loftee,consequence_terms,classified',
    [
        (0, 'hc', 'frameshift_variant', ZERO_EXPECTED),
        (0, 'HC', 'frameshift_variant', ONE_EXPECTED),
        (1, 'lc', 'frameshift_variant', ONE_EXPECTED),
        (1, hl.missing(hl.tstr), 'frameshift_variant', ONE_EXPECTED),
    ],
)
def test_class_3_assignment(clinvar_talos, loftee, consequence_terms, classified, make_a_mt):
    """

    Args:
        clinvar_talos ():
        loftee ():
        consequence_terms ():
        classified ():
        make_a_mt ():
    """

    anno_matrix = make_a_mt.annotate_rows(
        info=make_a_mt.info.annotate(clinvar_talos=clinvar_talos),
        vep=hl.Struct(
            transcript_consequences=hl.array(
                [
                    hl.Struct(consequence_terms=hl.set([consequence_terms]), lof=loftee),
                ],
            ),
        ),
    )

    anno_matrix = annotate_category_3(anno_matrix)
    assert anno_matrix.info.categoryboolean3.collect() == [classified]


@pytest.mark.parametrize(
    'spliceai_score,flag',
    [(0.1, 0), (0.11, 0), (0.3, 0), (0.49, 0), (0.5, 1), (0.69, 1), (0.9, 1)],
)
def test_category_5_assignment(spliceai_score: float, flag: int, make_a_mt):
    """

    Args:
        spliceai_score ():
        flag ():
        make_a_mt ():
    """

    matrix = make_a_mt.annotate_rows(info=make_a_mt.info.annotate(splice_ai_delta=spliceai_score))
    matrix = annotate_category_5(matrix)
    assert matrix.info.categoryboolean5.collect() == [flag]


@pytest.mark.parametrize(
    'am_class,classified',
    [('likely_pathogenic', 1), ('not_pathogenic', 0), ('', 0), (hl.missing('tstr'), 0)],
)
def test_class_6_assignment(am_class, classified, make_a_mt):
    """

    Args:
        am_class ():
        classified ():
        make_a_mt ():
    """

    anno_matrix = make_a_mt.annotate_rows(
        vep=hl.Struct(transcript_consequences=hl.array([hl.Struct(am_class=am_class)])),
    )

    anno_matrix = annotate_category_6(anno_matrix)
    anno_matrix.rows().show()
    assert anno_matrix.info.categoryboolean6.collect() == [classified]


def annotate_c6_missing(make_a_mt, caplog):
    """
    test what happens if the am_class attribute is missing

    Args:
        make_a_mt ():
    """
    anno_matrix = make_a_mt.annotate_rows(
        vep=hl.Struct(transcript_consequences=hl.array([hl.Struct(not_am='a value')])),
    )

    anno_matrix = annotate_category_6(anno_matrix)
    assert anno_matrix.info.categoryboolean6.collect() == [0]
    assert 'AlphaMissense class not found, skipping annotation' in caplog.text


def test_green_and_new_from_panelapp():
    """
    TODO make a proper object
    check that the set expressions from panelapp data are correct
    this is collection of ENSG names from panelapp
    2 set expressions, one for all genes, one for new genes only
    """
    mendeliome_data = {
        'ENSG00ABCD': {'new': [1], 'symbol': 'ABCD'},
        'ENSG00EFGH': {'new': [], 'symbol': 'EFHG'},
        'ENSG00IJKL': {'new': [2], 'symbol': 'IJKL'},
    }
    mendeliome = PanelApp.model_validate({'genes': mendeliome_data})
    green_expression, new_expression = green_and_new_from_panelapp(mendeliome)

    # check types
    assert isinstance(green_expression, hl.SetExpression)
    assert isinstance(new_expression, hl.SetExpression)

    # check content by collecting
    assert sorted(green_expression.collect()[0]) == ['ENSG00ABCD', 'ENSG00EFGH', 'ENSG00IJKL']
    assert new_expression.collect()[0] == {'ENSG00ABCD', 'ENSG00IJKL'}


@pytest.mark.parametrize(
    'exomes,genomes,clinvar,length',
    [(0, 0, 0, 1), (1.0, 0, 0, 0), (1.0, 0, 1, 1), (0.0001, 0.0001, 0, 1), (0.0001, 0.0001, 1, 1)],
)
def test_filter_rows_for_rare(exomes, genomes, clinvar, length, make_a_mt):
    """

    Args:
        exomes ():
        genomes ():
        clinvar ():
        length ():
        make_a_mt ():
    """
    anno_matrix = make_a_mt.annotate_rows(
        info=make_a_mt.info.annotate(gnomad_ex_af=exomes, gnomad_af=genomes, clinvar_talos=clinvar),
    )
    matrix = filter_to_population_rare(anno_matrix)
    assert matrix.count_rows() == length


@pytest.mark.parametrize(
    'gene_ids,length',
    [
        ({'not_green'}, 0),
        ({'green'}, 1),
        ({'gene'}, 1),
        ({'gene', 'not_green'}, 1),
        ({'green', 'gene'}, 2),
        ({hl.missing(t=hl.tstr)}, 0),
    ],
)
def test_filter_to_green_genes_and_split(gene_ids, length, make_a_mt):
    """

    Args:
        gene_ids ():
        length ():
        make_a_mt ():
    """
    green_genes = hl.literal({'green', 'gene'})
    anno_matrix = make_a_mt.annotate_rows(
        geneIds=hl.literal(gene_ids),
        vep=hl.Struct(
            transcript_consequences=hl.array([hl.Struct(gene_id='gene', biotype='protein_coding', mane_select='')]),
        ),
    )
    matrix = split_rows_by_gene_and_filter_to_green(anno_matrix, green_genes)
    assert matrix.count_rows() == length


def test_filter_to_green_genes_and_split__consequence(make_a_mt):
    """

    Args:
        make_a_mt ():
    """

    green_genes = hl.literal({'green'})
    anno_matrix = make_a_mt.annotate_rows(
        geneIds=green_genes,
        vep=hl.Struct(
            transcript_consequences=hl.array(
                [
                    hl.Struct(gene_id='green', biotype='protein_coding', mane_select=''),
                    hl.Struct(gene_id='green', biotype='batman', mane_select='NM_Bane'),
                    hl.Struct(gene_id='green', biotype='non_coding', mane_select=''),
                    hl.Struct(gene_id='NOT_GREEN', biotype='protein_coding', mane_select=''),
                ],
            ),
        ),
    )
    matrix = split_rows_by_gene_and_filter_to_green(anno_matrix, green_genes)
    assert matrix.count_rows() == 1
    matrix = matrix.filter_rows(hl.len(matrix.vep.transcript_consequences) == TWO_EXPECTED)
    assert matrix.count_rows() == 1


@pytest.mark.parametrize(
    'one,two,three,four,five,six,pm5,length',
    [
        (0, 0, 0, 'missing', 0, 0, 'missing', 0),
        (1, 0, 0, 'missing', 0, 0, 'missing', 1),
        (0, 1, 0, 'missing', 0, 0, 'missing', 1),
        (0, 0, 1, 'missing', 0, 0, 'missing', 1),
        (0, 0, 0, 'present', 0, 0, 'missing', 1),
        (0, 0, 0, 'missing', 1, 0, 'missing', 1),
        (0, 0, 0, 'missing', 0, 1, 'missing', 1),
        (0, 0, 0, 'missing', 0, 0, 'present', 1),
        (0, 1, 1, 'missing', 0, 0, 'missing', 1),
    ],
)
def test_filter_to_classified(one, two, three, four, five, six, pm5, length, make_a_mt):
    """

    Args:
        one argument per category
        make_a_mt (): a template matrix table
    """
    anno_matrix = make_a_mt.annotate_rows(
        info=make_a_mt.info.annotate(
            categoryboolean1=one,
            categoryboolean2=two,
            categoryboolean3=three,
            categorysample4=four,
            categoryboolean5=five,
            categoryboolean6=six,
            categorydetailsPM5=pm5,
        ),
    )
    matrix = filter_to_categorised(anno_matrix)
    assert matrix.count_rows() == length


def test_talos_clinvar_default(make_a_mt):
    """
    no private annotations applied
    Args:
        make_a_mt (hl.MatrixTable):
    """

    mt = annotate_talos_clinvar(make_a_mt)
    assert mt.count_rows() == ONE_EXPECTED
    assert not [x for x in mt.info.clinvar_talos.collect() if x == 1]
    assert not [x for x in mt.info.clinvar_talos_strong.collect() if x == 1]


@pytest.mark.parametrize(
    'rating,stars,rows,regular,strong',
    [
        ('benign', 0, 1, 0, 0),
        ('benign', 1, 0, 0, 0),
        ('other', 7, 1, 0, 0),
        ('pathogenic', 0, 1, 1, 0),
        ('pathogenic', 1, 1, 1, 1),
    ],
)
def test_annotate_talos_clinvar(rating, stars, rows, regular, strong, tmp_path, make_a_mt):
    """
    Test intention
    - take a VCF of two variants w/default clinvar annotations
    - create a single variant annotation table with each run
    - apply the parametrized annotations to the table
    """

    # make into a data frame
    table = hl.Table.from_pandas(
        pd.DataFrame(
            [
                {
                    'locus': hl.Locus(contig='chr1', position=12345),
                    'alleles': ['A', 'G'],
                    'clinical_significance': rating,
                    'gold_stars': stars,
                    'allele_id': 1,
                },
            ],
        ),
        key=['locus', 'alleles'],
    )

    table_path = str(tmp_path / 'anno.ht')
    table.write(table_path)

    returned_table = annotate_talos_clinvar(make_a_mt, clinvar=table_path)
    assert returned_table.count_rows() == rows
    assert len([x for x in returned_table.info.clinvar_talos.collect() if x == 1]) == regular
    assert len([x for x in returned_table.info.clinvar_talos_strong.collect() if x == 1]) == strong
