"""
code for lifting over models from 1.1.0 to 1.2.0
"""

TRIO_SIZE = 3


def resultdata(data_dict: dict) -> dict:
    """
    Lift over ResultData from 1.1.0 to 1.2.0
    Requires the adjustment of ResultData.results.ParticipantResults.metadata.ParticipantMeta
    adding accurate values for `family_size` and `part_of_trio`
    """
    # check we're upgrading the right version
    # this could be from any prior version
    if not data_dict['version'] < '1.2.0':
        raise AssertionError(f'This method cannot upgrade from {data_dict["version"]}')

    for _sample, content in data_dict['results'].items():
        sample_meta = content['metadata']
        family_size = len(sample_meta['members'])
        # absolute guess here, we don't currently record enough info to know
        part_of_trio = family_size == TRIO_SIZE

        sample_meta.update({'family_size': family_size, 'part_of_trio': part_of_trio})

    data_dict['version'] = '1.2.0'
    return data_dict
