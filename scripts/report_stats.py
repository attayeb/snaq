import click
from reportlab.platypus import SimpleDocTemplate, Image
from reportlab.lib.pagesizes import letter

@click.command()
@click.option("--inp")
@click.option("--outp")
def make_report(inp, outp):
    story = []
    doc = SimpleDocTemplate(outp,pagesize=letter,
                        rightMargin=72,leftMargin=72,
                        topMargin=72,bottomMargin=18)
    for plot_file in inp.split(","):
            story.append(Image(plot_file))

    doc.build()

if __name__ == "__main__":
    make_report()
