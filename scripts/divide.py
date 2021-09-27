import click
import os

@click.command()
@click.option("-i", "origin", required=True, type=str)
@click.option("-o", "target", required=True, type=str)
@click.option("-n", "number", required=True, type=int)
def divide(origin, target, number):
    file_names = 
