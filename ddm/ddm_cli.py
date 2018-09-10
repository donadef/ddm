# -*- coding: utf-8 -*-

"""Console script for ddm."""
import click
from ddm.ddm_run import DDM


@click.command()
@click.option('--config',
              help='Your config file. See the docs for an example: https://ddm.readthedocs.io ',
              # prompt=True,
              type=click.Path(exists=True, allow_dash=True),
              required=True)
@click.option('--complex',
              help='The pbd file of the complex you want to study.',
              type=click.Path(exists=True, allow_dash=True),
              required=True)
@click.version_option(version='0.1.0')
def main(config, complex):
    """Console script for ddm."""
    click.echo(f"dG of binding calculated for the complex {complex} using the specified {config} file.")
    ddm = DDM(config, complex)
    ddm.perform_ddm()
    return 0


if __name__ == "__main__":
    main()
