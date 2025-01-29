"""
Provides access to config variables. Complete theft of the CPG config module, minimised
https://github.com/populationgenomics/cpg-utils/blob/main/cpg_utils/config.py
"""

from os import environ
from typing import Any

import toml

# We use these globals for lazy initialization, but pylint doesn't like that.
# pylint: disable=global-statement, invalid-name
CONFIG_TYPE = dict[str, Any]
_config: CONFIG_TYPE | None = None  # Cached config, initialized lazily.
ENV_VAR: str = 'TALOS_CONFIG'


class ConfigError(Exception):
    """
    Error retrieving keys from config.
    """


class Unsupplied:
    pass


def config_retrieve(key: list[str] | str, default: Any | None = Unsupplied, config_path: str | None = ENV_VAR) -> Any:
    """
    Retrieve key from config, assuming nested key specified as a list of strings.

    >> config_retrieve(['workflow', 'access_level'], config={'workflow': {'access_level': 'test'}})
    'test'

    >> config_retrieve(['workflow', 'access_level'], config={}, default='default')
    'default'

    Allow None as default value
    >> config_retrieve(['key1', 'key2', 'key3'], config={}, default=None) is None
    True
    """

    global _config
    if _config is None:  # Lazily initialize the config.
        # when using default, extract from environment variable
        if config_path == ENV_VAR:
            config_path = environ[ENV_VAR]
            if config_path is None:
                raise ValueError(f'Set the environment variable {ENV_VAR}, or provide a specfic config path')
        elif config_path is None:
            raise ValueError(f'Supply a non-default file path, or export to the environment variable {ENV_VAR}')

        with open(config_path) as f:
            _config = toml.loads(f.read())

    if isinstance(key, str):
        key = [key]

    if not key:
        raise ValueError('Key cannot be empty')

    d = _config
    for idx, k in enumerate(key):
        if k not in d:
            if default is Unsupplied:
                message = f'Key "{k}" not found in {d}'
                if idx > 0:
                    key_bits = ' -> '.join(key[: idx + 1])
                    message += f' (path: {key_bits})'

                raise ConfigError(message)
            return default

        d = d[k]

    return d


def config_check(key: list[str], expected_type: type | tuple[type]) -> list[str]:
    """
    take a path to a config entry, and one or more expected types
    return a list of Strings:
        - if the value is present in the config dict, but the wrong type, explain
        - if the keys are not present in the config, explain where the key was absent
        - if the key(s) lead to a value, and the type is correct, return an empty list
    Args:
        key (list[str]): the keys for each layer in the config dict, leading to a value to test
        expected_type (type | tuple[type]): Type(s) we accept for this config value
    Returns:
        a list of the faults in the config search & type check, can be empty
    """

    try:
        value = config_retrieve(key)
        if isinstance(value, expected_type):
            return []
        config_keys = ' -> '.join(key)
        actual_type = type(value)
        return [f'config path {config_keys} was {actual_type}, expected {expected_type}']

    except ConfigError as ce:
        return [str(ce)]
