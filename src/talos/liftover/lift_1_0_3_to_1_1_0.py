"""
code for lifting over models from 1.0.3 to 1.1.0
"""


def resultdata(data_dict: dict) -> dict:
    """
    Lift over ResultData from 1.0.3 to 1.1.0
    Requires the adjustment of ResultData.results.ParticipantResults.metadata.ParticipantMeta.panel_*
    """
    # check we're upgrading the right version
    # this could be from any prior version
    if not data_dict['version'] < '1.1.0':
        raise AssertionError(f'This method cannot upgrade from {data_dict["version"]}')

    for _sample, content in data_dict['results'].items():
        sample_meta = content['metadata']
        assert all(key in sample_meta for key in ['panel_ids', 'panel_names'])
        sample_meta['panel_details'] = dict(
            zip(
                sample_meta.pop('panel_ids'),
                sample_meta.pop('panel_names'),
                strict=True,
            ),
        )

        for variant in content['variants']:
            panels = variant['panels']
            panels['forced'] = {pid: 'UNKNOWN' for pid in panels.get('forced', [])}
            panels['matched'] = {pid: 'UNKNOWN' for pid in panels['matched']}

    data_dict['version'] = '1.1.0'
    return data_dict
