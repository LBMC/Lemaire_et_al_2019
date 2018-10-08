#!/usr/bin/env python3


# -*- coding utf-8 -*-

import union_dataset_function
import sqlite3
import os
import figure_producer
import numpy as np
import plotly.graph_objs as go
import plotly
import rpy2.robjects as robj
import rpy2.robjects.vectors as v
import group_factor


def group_factors():
    """
    Define some group of factor
    :return: (2 list of strings and one dictionary):
        1. list of factors that down-regulate at_rich exons
        2. list of factors that down-regulate gc_rich exons
        3. dictionary that contains (1 - the list of exons regulated by U1 factors, 2- the list of exons regulated \
    """
    dic_sf = {"U1": group_factor.u1_factors, "U2": group_factor.u2_factors, "chromatin_factors": group_factor.chromatin_factors}
    return group_factor.at_rich_down, group_factor.gc_rich_down, dic_sf


def get_exons_values(cnx, sf_list, target_column, regulation, common_exon_2_remove=None):
    """
    Return the values of target_column in every`\
    `regulation`` exons regulated by a splicing factor in (one or multiple) cell lines.

    :param cnx: (sqlite3 connexion object) allow connexion to sed database
    :param sf_list:  (list of string) the list of splicing factor studied
    :param target_column: (string) the value for which we want to get the median value for the ``regulation`` \
    exon.
    :param regulation: (list of string) up or down or up + down
    :param common_exon_2_remove: (list of int  and int) list of exons
    :return: (list of list of float) each sublist corresponds to the value of `` target_column`` for \
        every exons regulated by a splicing factor
    """
    exon_list = []
    for sf_name in sf_list:
        exon_list += union_dataset_function.get_events_4_a_sl(cnx, sf_name, regulation)
    exon_list = union_dataset_function.washing_events_all(exon_list)
    print(len(exon_list))
    if common_exon_2_remove is not None:
        new_exon_list = [exon for exon in exon_list if exon not in common_exon_2_remove]
        exon_list = new_exon_list
        print("After removing common exons %s" % len(exon_list))
    values = figure_producer.get_redundant_list_of_value(cnx, exon_list, target_column)
    return values


def get_common_exon(cnx, sf_list1, sf_list2, regulation):
    """
    Return the values of target_column in every`\
    `regulation`` exons regulated by a splicing factor in (one or multiple) cell lines.

    :param cnx: (sqlite3 connexion object) allow connexion to sed database
    :param sf_list1:  (list of string) the 1st list of splicing factor studied
    :param sf_list2:  (list of string) the 2nd list of splicing factor studied
    :param regulation: (string) the regulation wanted
    :return: (list of list of int)  list opf common exon regulated by the list sf_list1 and sf_list2
    """
    exon_list1 = []
    for sf_name in sf_list1:
        exon_list1 += union_dataset_function.get_events_4_a_sl(cnx, sf_name, regulation)
    exon_list1 = union_dataset_function.washing_events_all(exon_list1)
    exon_list2 = []
    for sf_name in sf_list2:
        exon_list2 += union_dataset_function.get_events_4_a_sl(cnx, sf_name, regulation)
    exon_list2 = union_dataset_function.washing_events_all(exon_list2)
    common_exon_list = [exon for exon in exon_list1 if exon in exon_list2]
    return common_exon_list


def get_control_exon_information(cnx, exon_type, target_column):
    """
    :param cnx: (sqlite3 object) allow connection to sed database
    :param exon_type: (string) the type of control exon we want to use
    :param target_column:  (string) the column of interest
    :return:
        * result: (list of tuple) every information about control exons
        * names: (list of string) the name of every column in sed table
    """
    cursor = cnx.cursor()
    if exon_type != "ALL":
        query = """SELECT {}
                   FROM sed
                   WHERE exon_type LIKE '%{}%'""".format(target_column, exon_type)
    else:
        query = """SELECT {}
                   FROM sed
                """.format(target_column)
    cursor.execute(query)
    result = cursor.fetchall()
    # turn tuple into list
    values = []
    for val in result:
        values.append(val[0])
    return values


def mann_withney_test_r(list_values1, list_values2):
    wicox = robj.r("""

    function(x, y){
        test = wilcox.test(x,y, alternative='two.sided', correct=F)
        return(test$p.value)
    }

                   """)
    pval = float(wicox(v.FloatVector(list_values1), v.FloatVector(list_values2))[0])
    return pval


def create_figure(list_values, list_name, target_name, output, regulation, exon_lvl=True):
    """
    Create a figure showing the values of ``target_name`` for every exon list regulated by ``list_name`` factor

    :param list_values: (list of list of float) each sublist corresponds to the value of  ``target_name`` \
    for every factor in ``list_name``
    :param list_name: (list of string) list of name of different factor studied
    :param target_name: (string) the name of the target column
    :param output: (string) path where the output_file will be created
    :param regulation: (string) up or down
    """
    data = []
    color_list = ['#1f77b4', '#2ca02c', '#1f77b4', '#2ca02c', '#ff0000',
                  '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    pval_at_cg = mann_withney_test_r(list_values[0], list_values[1])
    pval_u1_u2 = mann_withney_test_r(list_values[2], list_values[3])
    pval_cg_ctrl = mann_withney_test_r(list_values[0], list_values[4])
    pval_at_ctrl = mann_withney_test_r(list_values[1], list_values[4])
    for i in range(len(list_values)):
        data.append({"y": list_values[i], "type": "violin",
                     "name": list_name[i], "visible": True, "line": {"color": color_list[i]},
                     "box": {"visible": True}, "meanline": {"visible": True}})

    layout = go.Layout(
        title='%s of %s exons regulated by different factors'

              '<br> mann whitney test GC exons vs CTRL exons (two sided) : p = %.2E'
              '<br> mann whitney test AT exons vs CTRL exons (two sided): p = %.2E'
              % (target_name, regulation, pval_cg_ctrl, pval_at_ctrl),
        yaxis=dict(
            autorange=True,
            showgrid=True,
            zeroline=True,
            autotick=True,
            title=target_name,
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
        paper_bgcolor='rgb(243, 243, 243)',
        plot_bgcolor='rgb(243, 243, 243)',
        showlegend=True
    )

    fig = {"data": data, "layout": layout}
    if exon_lvl:
        plotly.offline.plot(fig, filename="%s%s_%s_exons_lvl.html" % (output, target_name, regulation),
                        auto_open=False, validate=False)
    else:
        plotly.offline.plot(fig, filename="%s%s_%s_sf_lvl.html" % (output, target_name, regulation),
                        auto_open=False, validate=False)


def get_median_values_sf(cnx, list_sf, target_column, regulation, common_exons_2_remove=None):
    """
    Get the list of median values of target_column for each exons regulated by the splicing \
    factor listed in ``list_sf``
    :param cnx: (sqlite3 connect object) sed database
    :param list_sf:  (string) the list of splicing factors
    :param target_column: (string) the feature for which we want to extract the median values for the sf \
    in ``sf_list``
    :param regulation: (string) down-regulated exons
    :param common_exons_2_remove: (
    :return: (list of float) list of median values of each sf
    """
    list_median_values = []
    for sf in list_sf:
        exon_list = union_dataset_function.get_events_4_a_sl(cnx, sf, regulation)
        if common_exons_2_remove is not None:
            exon_list = [exon for exon in exon_list if exon not in common_exons_2_remove]
        values = figure_producer.get_redundant_list_of_value(cnx, exon_list, target_column)
        list_median_values.append(np.median(values))
    return list_median_values



def main():
    regulation = "down"
    exon_type = "CCE"
    seddb = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/data/sed.db"
    output = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/result/barplot_force_gt_at_exon/"
    if not os.path.isdir(output):
        os.mkdir(output)
    cnx = sqlite3.connect(seddb)
    at_rich_down, gc_rich_down, dic_sf = group_factors()
    at_rich_down = list(at_rich_down)
    gc_rich_down = list(gc_rich_down)
    print(len(at_rich_down))
    print(len(gc_rich_down))

    common_exon_list = get_common_exon(cnx, gc_rich_down, at_rich_down, regulation)
    common_spliceosome = get_common_exon(cnx, dic_sf["U1"], dic_sf["U2"], regulation)
    print("Common exons AT- GC : %s" % len(common_exon_list))
    print("Common spliceosome : %s" % len(common_spliceosome))
    targets = ["force_donor", "force_acceptor"]
    name_values_list = ["GC_sf_exon", "AT_sf_exons", "U1_exons", "U2_exons", exon_type]
    name_values = {"GC_sf_exon": gc_rich_down, "AT_sf_exons":at_rich_down, "U1_exons": dic_sf["U1"], "U2_exons":dic_sf["U2"], exon_type: None}
    for target_column in targets:
        value_list_exon = []
        value_list_sf = []
        for factor in name_values_list:
            if factor == "CCE":
                print("   CCE exons - %s recovery" % target_column)
                value_list_exon.append(get_control_exon_information(cnx, exon_type, target_column))
                value_list_sf.append(get_control_exon_information(cnx, exon_type, target_column))
            else:
                print("%s graphics..." % target_column)
                print("   %s factors - %s recovery" % (factor, target_column))
                if "GC" in factor or "AT" in factor:
                    value_list_exon.append(get_exons_values(cnx, name_values[factor], target_column, regulation, common_exon_list))
                    value_list_sf.append(get_median_values_sf(cnx, name_values[factor], target_column, regulation, common_exon_list))
                elif "U1" in factor or "U2" in factor:
                    value_list_exon.append(get_exons_values(cnx, name_values[factor], target_column, regulation, common_spliceosome))
                    value_list_sf.append(get_median_values_sf(cnx, name_values[factor], target_column, regulation, common_spliceosome))
                    if  "U1" in factor or "U2" in factor:
                        print("----")
                        print(name_values[factor])
                        print(value_list_sf[-1])
                        print("---")
                else:
                    value_list_exon.append(get_exons_values(cnx, name_values[factor], target_column, regulation))
                    value_list_sf.append(get_median_values_sf(cnx, name_values[factor], target_column, regulation))
                print("   " + str(value_list_sf[-1]))
        print(name_values_list)
        print(value_list_sf[:-1])
        create_figure(value_list_exon, name_values_list, target_column, output, regulation, exon_lvl=True)
        create_figure(value_list_sf, name_values_list, target_column, output, regulation, exon_lvl=False)


if __name__ == "__main__":
    main()

