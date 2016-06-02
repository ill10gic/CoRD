from configparser import ConfigParser
from hs_restclient import HydroShare, HydroShareAuthBasic

import click
import glob
import os
import shutil

from ..modelrun import modelrun_series


@click.group()
@click.option('--debug', is_flag=True, default=False)
@click.option('--logfile', default=None)
@click.pass_context
def cli(ctx, debug, logfile):

    ctx.obj = dict(DEBUG=debug, LOGFILE=logfile)


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
@click.pass_context
def interactive(ctx, data_dir, initial_veg_map, vegzone_map,
                ripcas_required_data, peak_flows_file, geometry_file,
                streambed_roughness, streambed_slope, dflow_run_fun):
    """Run CoRD interactively or with options"""

    modelrun_series(data_dir, initial_veg_map, vegzone_map,
                    ripcas_required_data, peak_flows_file, geometry_file,
                    streambed_roughness, streambed_slope, dflow_run_fun,
                    ctx.obj['LOGFILE'], ctx.obj['DEBUG'])


@cli.command()
@click.argument('config_file', type=CPE())
@click.pass_context
def from_config(ctx, config_file):
    """Run CoRD with params from <config_file>"""

    cfg = load_args_from_config(config_file)

    ctxlog = ctx.obj['LOGFILE']
    logfile = ctxlog if ctxlog is not None else cfg['log_f']

    modelrun_series(
        cfg['data_dir'],
        cfg['initial_vegetation_map'],
        cfg['vegzone_map'],
        cfg['veg_roughness_shearres_lookup'],
        cfg['peak_flows_file'],
        cfg['geometry_file'],
        cfg['streambed_roughness'],
        cfg['streambed_slope'],
        cfg['dflow_run_fun'],
        logfile,
        ctx.obj['DEBUG']
    )


@cli.command()
@click.option('--username', prompt=True)
@click.option('--password', prompt=True)
@click.option('--modelrun-dir', prompt=True)
@click.option('--include-shear-nc', type=click.BOOL, default=False)
@click.option('--resource-title', prompt=True)
@click.option('-k', '--keyword', multiple=True, default=None)
@click.pass_context
def post_hs(ctx, username, password, modelrun_dir, include_shear_nc,
            resource_title, keyword):
    """Post the model run data to HydroShare"""
    # iterate over files and folders of interest, adding files to resource
    # RipCAS files are the vegetation .asc; TODO XXX include base RipCAS XXX TODO

    export_dir = os.path.join(modelrun_dir, 'export')
    veg_export_dir = os.path.join(export_dir, 'vegetation')

    if os.path.isdir(export_dir):
        shutil.rmtree(export_dir)

    os.mkdir(export_dir)
    os.mkdir(veg_export_dir)

    veg_pattern = os.path.join(modelrun_dir, 'ripcas-*', 'vegetation.asc')
    for tstep, veg_map in enumerate(glob.glob(veg_pattern)):

        veg_map_path = os.path.join(
            export_dir, 'vegetation', 'vegetation-%s.asc' % tstep
        )
        shutil.copy(veg_map, veg_map_path)

    shutil.make_archive(veg_export_dir, 'zip', veg_export_dir)

    # connect
    hs = HydroShare(
        auth=HydroShareAuthBasic(username=username, password=password)
    )

    # create new resource
    rtype = 'GenericResource'

    r_id = hs.createResource(
       rtype, resource_title, keywords=keyword  # , abstract=abstract
    )

    print "adding vegmap archive file %s to resource %s" % (veg_map_path, r_id)
    hs.addResourceFile(r_id, os.path.join(export_dir, 'vegetation.zip'))

    inputs_dir = os.path.join(modelrun_dir, 'inputs')
    inputs_export_basename = os.path.join(export_dir, 'inputs')
    shutil.make_archive(inputs_export_basename, 'zip', inputs_dir)


    print "adding inputs archive file %s to resource %s" % (veg_map_path, r_id)
    hs.addResourceFile(
        r_id, os.path.join(export_dir, 'inputs.zip')
    )

    shear_export_dir = os.path.join(export_dir, 'shear')
    os.mkdir(shear_export_dir)

    shear_pattern = os.path.join(modelrun_dir, 'dflow-*', 'shear_out.asc')
    for tstep, shear_map in enumerate(glob.glob(shear_pattern)):

        shear_map_path = os.path.join(
            export_dir, 'shear', 'shear-%s.asc' % tstep
        )
        shutil.copy(shear_map, shear_map_path)

    shutil.make_archive(shear_export_dir, 'zip', shear_export_dir)

    print "adding shear archive file %s to resource %s" % (veg_map_path, r_id)
    hs.addResourceFile(
        r_id, os.path.join(export_dir, 'shear.zip')
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

    # general config options
    gen = dict(cfg['General'])
    if gen['log_f'] == u'':
        gen['log_f'] = None

    if gen['dflow_run_fun'] == u'':
        gen['dflow_run_fun'] = None

    curdir = os.path.dirname(__file__)
    if gen['initial_vegetation_map'] == u'':
        gen['initial_vegetation_map'] = os.path.join(
            curdir, '..', 'data', 'ripcas_inputs', 'vegclass_2z.asc'
        )

    if gen['vegzone_map'] == u'':
        gen['vegzone_map'] = os.path.join(
            curdir, '..', 'data', 'ripcas_inputs', 'zonemap_2z.asc'
        )

    if gen['veg_roughness_shearres_lookup'] == u'':
        gen['veg_roughness_shearres_lookup'] = os.path.join(
            curdir, '..', 'data', 'ripcas_inputs',
            'veg_roughness_shearres.xlsx'
        )

    if gen['geometry_file'] == u'':
        gen['geometry_file'] = os.path.join(
            curdir, '..', 'data', 'dflow_inputs', 'DBC_geometry.xyz'
        )

    if gen['peak_flows_file'] == u'':
        raise RuntimeError('PEAK_FLOWS_FILE must be defined in ' + config_file)

    if gen['streambed_roughness'] == u'':
        raise RuntimeError(
            'STREAMBED_ROUGHNESS must be defined in ' + config_file
        )

    if gen['streambed_slope'] == u'':
        raise RuntimeError('STREAMBED_SLOPE must be defined in ' + config_file)

    gen['streambed_roughness'] = float(gen['streambed_roughness'])
    gen['streambed_slope'] = float(gen['streambed_slope'])

    # hydroshare config options
    hs = dict(cfg['HydroShare'])
    if hs['sync_hydroshare'] == u'False':
        hs['sync_hydroshare'] = False
        hs['hs_username'] = None
        hs['hs_password'] = None

    ret = gen
    ret.update(hs)

    return ret
