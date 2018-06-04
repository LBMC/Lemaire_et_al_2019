#!/usr/bin/python3.5

# coding: utf8

"""
Description:
    This script aims, for each project to show the number of regulated exon contained in the same gene.

"""

# imports
import plotly.graph_objs as go
import plotly
import os
import figure_producer


def get_gene_id(cnx, id_project, regulation):
    """
    Get every gene_id for every exon ``regulation`` in ``id_project``.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param id_project: (int) a project id
    :param regulation: (string) up or down
    :return: (list of int) the number of regulated exon contained in the same gene.
    """
    gene_id_dic = {}
    if regulation == "up":
        regulation = ">= 0.1"
    else:
        regulation = "<= -0.1"
    cursor = cnx.cursor()
    query = """SELECT gene_id
               FROM ase_event
               WHERE id_project = %s
               AND delta_psi %s
               AND pvalue_glm_cor <= 0.05""" % (id_project, regulation)
    cursor.execute(query)
    res = cursor.fetchall()
    if len(res) == 0:
            query = """SELECT gene_id
               FROM ase_event
               WHERE id_project = %s
               AND delta_psi %s
               AND pvalue <= 0.05""" % (id_project, regulation)
            cursor.execute(query)
            res = cursor.fetchall()
    # filling the dictionary
    for gene_id in res:
        if gene_id not in gene_id_dic.keys():
            gene_id_dic[gene_id] = 1
        else:
            gene_id_dic[gene_id] += 1
    # getting the list of the number of regulated exon by gene
    exon_reg_nb = sorted(list(gene_id_dic.values()))
    return exon_reg_nb


def graphic_maker(exon_reg_nb, project_name, regulation, output):
    """
    Create a graphic that shows the number of regulated exon contained in the same gene for the project ``project_name``

    :param exon_reg_nb: (list of int) the number of regulated exon contained in the same gene.
    :param project_name: (string) the name of the project
    :param regulation: (string) up or down
    :param output: (string) path where the output will be created
    """
    abscissa = list(range(1, 6))
    ordinate = []
    count = 0
    for x in abscissa[:-1]:
        count += exon_reg_nb.count(x) * x
        ordinate.append(float(exon_reg_nb.count(x) * x) / sum(exon_reg_nb) * 100)
    ordinate.append(float(sum(exon_reg_nb)-count) / sum(exon_reg_nb) * 100)
    trace = go.Bar(x=[" 1 ", " 2 ", " 3 " , " 4 ", "5+"], y=ordinate)
    layout = go.Layout(
        title='%s-regulated exons for the project %s\n Number of %s-exons : %s'
              % (regulation, project_name, regulation, sum(exon_reg_nb)),
        yaxis=dict(
            title="Regulated exons(%)",
            showgrid=True,
            zeroline=True,
            autotick=True,
            gridcolor='rgb(255, 255, 255)',
            gridwidth=1,
            zerolinecolor='rgb(255, 255, 255)',
            zerolinewidth=2,
        ),
        margin=dict(
            l=40,
            r=30,
            b=150,
            t=100,
        ),
        xaxis=dict(
            title="Regulated exons sharing the same gene",
            type="category",
            dtick=1
        ),
        paper_bgcolor='rgb(243, 243, 243)',
        plot_bgcolor='rgb(243, 243, 243)',
        showlegend=True
    )

    fig = go.Figure(data=[trace], layout=layout)
    plotly.offline.plot(fig, filename="%s%s_%s_exons_figure.html" % (output, project_name, regulation),
                        auto_open=False)


def main():
    """
    Produce every barplot graphics.
    """
    output = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/result/barplot_of_regulated_gene/"
    # If the output directory does not exist, then we create it !
    if not os.path.isdir(output):
        os.mkdir(output)
    seddb = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/data/sed.db"
    regs = ["up", "down"]
    cnx = figure_producer.connexion(seddb)
    id_projects, name_projects = figure_producer.get_interest_project(cnx)
    for i in range(len(id_projects)):
        for regulation in regs:
            exon_reg_nb = get_gene_id(cnx, id_projects[i], regulation)
            graphic_maker(exon_reg_nb, name_projects[i], regulation, output)

if __name__ == "__main__":
    main()