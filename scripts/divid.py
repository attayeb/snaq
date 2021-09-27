import click
import os
from os.path import join


@click.command()
@click.argument('folder')
@click.argument('target')
@click.argument('number')
def divide(folder, target, number):
    files = os.listdir(folder)
    first = [x for x in files if "_R1_" in x][:int(number)]
    second = [x.replace("_R1_", "_R2_") for x in first]
    os.mkdir(target)
    [os.rename(join(folder, x), join(target, x)) for x in first]
    [os.rename(join(folder, x), join(target, x)) for x in second]
    print("finished")

if __name__ == "__main__":
    divide()

