from configparser import ConfigParser

import click

from ..modelrun import modelrun_series


@click.group()
@click.option('--debug', is_flag=True, default=False)
@click.pass_context
def cli(ctx, debug):

    ctx.obj = dict(DEBUG=debug)


def CPE():
    return click.Path(exists=True)


@cli.command()
@click.option('--data-dir', prompt=True)
@click.option('--initial-veg-map', prompt=True, type=CPE())
@click.option('--vegzone-map', prompt=True, type=CPE())
@click.option('--ripcas-required-data', prompt=True, type=CPE())
@click.option('--peak-flows-file', prompt=True, type=CPE())
@click.option('--geometry-file', prompt=True, type=CPE())
@click.option('--streambed-roughness', prompt=True, type=click.FLOAT)
@click.option('--streambed-slope', prompt=True, type=click.FLOAT)
@click.option('--dflow-run-fun', default=None)
@click.option('--logfile', default=None)
@click.pass_context
def interactive(ctx, data_dir, initial_veg_map, vegzone_map,
                ripcas_required_data, peak_flows_file, geometry_file,
                streambed_roughness, streambed_slope, dflow_run_fun, logfile):

    modelrun_series(data_dir, initial_veg_map, vegzone_map,
                    ripcas_required_data, peak_flows_file, geometry_file,
                    streambed_roughness, streambed_slope, dflow_run_fun,
                    logfile, ctx.obj['DEBUG'])


@cli.command()
@click.argument('config_file', type=CPE())
@click.pass_context
def from_config(ctx, config_file):

    cfg = load_args_from_config(config_file)

    modelrun_series(
        cfg['data_dir'],
        cfg['initial_vegetation_map'],
        cfg['vegzone_map'],
        cfg['ripcas_required_data'],
        cfg['peak_flows_file'],
        cfg['geometry_file'],
        cfg['streambed_roughness'],
        cfg['streambed_slope'],
        cfg['dflow_run_fun'],
        cfg['log_f'],
        ctx.obj['DEBUG']
    )


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

    gen = dict(cfg['General'])
    if gen['log_f'] == u'':
        gen['log_f'] = None

    if gen['dflow_run_fun'] == u'':
        gen['dflow_run_fun'] = None

    gen['streambed_roughness'] = float(gen['streambed_roughness'])
    gen['streambed_slope'] = float(gen['streambed_slope'])

    hs = dict(cfg['HydroShare'])
    if hs['sync_hydroshare'] == u'False':
        hs['sync_hydroshare'] = False
        hs['hs_username'] = None
        hs['hs_password'] = None

    ret = gen
    ret.update(hs)

    return ret
