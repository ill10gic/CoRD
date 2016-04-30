from configparser import ConfigParser

import click

from ..modelrun import modelrun_series


@click.command()
@click.option('--config', type=click.File('r'), default=None)
@click.option('--data_dir')
@click.option('--initial_vegetation_file')
@click.option('--vegzone_map')
@click.option('--ripcas_required_data')
@click.option('--peak_flows_file')
@click.option('--geometry_file')
@click.option('--streambed_roughness')
@click.option('--streambed_slope')
@click.option('--dflow_run_fun')
@click.option('--log_f')
@click.option('--hs_sync', is_flag=True)
@click.option('--hs_username', default=None)
@click.option('--hs_password', default=None)
def cli(config, data_dir, initial_vegetation_file, vegzone_map,
        ripcas_required_data, peak_flows_file, geometry_file,
        streambed_roughness, streambed_slope, dflow_run_fun, log_f,
        hs_sync, hs_username, hs_password):

    if config is not None:
        # load these values from the config file
        cfg = load_args_from_config(config)
        modelrun_series(
            cfg['data_dir'],
            cfg['initial_vegetation_file'],
            cfg['vegzone_map'],
            cfg['ripcas_required_data'],
            cfg['peak_flows_file'],
            cfg['geometry_file'],
            cfg['streambed_roughness'],
            cfg['streambed_slope'],
            cfg['dflow_run_fun'],
            cfg['log_f'],
            cfg['hs_sync'],
            cfg['hs_username'],
            cfg['hs_password']
        )

    else:
        modelrun_series(data_dir, initial_vegetation_file, vegzone_map,
                        ripcas_required_data, peak_flows_file, geometry_file,
                        streambed_roughness, streambed_slope, dflow_run_fun,
                        log_f, hs_sync, hs_username, hs_password)


def load_args_from_config(config_file):
    """
    Load a config file using ConfigParser and add defaults for missing values

    Arguments:
        config_file (str): path to configuration file

    Returns
        (dict) dictionary of kwargs ready for modelrun_series
    """

    cfg = ConfigParser(inline_comment_prefixes='#')
    cfg.read(config_file)

    gen = cfg['General']
    if gen['log_f'] == u'':
        gen['log_f'] = None

    if gen['dflow_run_fun'] == u'':
        gen['dflow_run_fun'] = None

    hs = cfg['HydroShare']

    if hs['sync_hydroshare'] == u'False':
        hs['sync_hydroshare'] = False
        hs['hs_username'] = None
        hs['hs_password'] = None

    d = dict(gen)
    d.update(dict(hs))

    return d
