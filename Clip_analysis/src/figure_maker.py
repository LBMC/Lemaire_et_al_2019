#!/usr/bin/env python3

# -*- coding: utf-8 -*-

import rpy2.robjects as robj
import rpy2.robjects.vectors as v
import warnings
from rpy2.rinterface import RRuntimeWarning
import plotly.figure_factory as ff
import plotly
import numpy as np


def web_logo_creator(sequence_list, sequence_name, output):
    """
    :param sequence_list: (tuple of strings) - list of sequences
    :param sequence_name: (string) name of the sequence
    :param output: (string) the folder where the results will be created
    """
    warnings.filterwarnings("ignore", category=RRuntimeWarning)
    weblogo_maker = robj.r("""
    library("ggplot2")
    library("ggseqlogo")

    function(mys_seq, name_file, mytitle, size){
        s1 = 15
        cs1 = make_col_scheme(chars=c('A','T','G','C', 'R', 'Y', 'W', 'S', 'K', 'M'), groups=c('g1','g2','g3','g4','g5', 'g6', 'g7', 'g8', 'g9', 'g10'),cols=c('limegreen','brown1','gold','dodgerblue3','darkorange', "brown1", "limegreen", "dodgerblue3", "darkorchid3", "dodgerblue3"), name='custom1')

        p1 = ggseqlogo(mys_seq,  method = "bit", col_scheme=cs1, namespace = c('A','T','G','C', 'R', 'Y', 'W', 'S', 'K', 'M')) + theme_logo() + scale_x_discrete(limits = as.character(seq(1,size, by=1)), labels = as.character(seq(1,size, by=2)), breaks = as.character(seq(1, size, by=2))) + theme(axis.title.y=element_text(size=s1+25), legend.position="none")
        p1 = p1 + ggtitle(mytitle) +  theme(plot.title = element_text(hjust = 0.5))


        p1 = p1 + theme(axis.text=element_text(size=s1 + 25), plot.title = element_text(size=s1 + 30))
        p1 = p1 + scale_y_discrete(limits = c(0, 0.5, 1), labels = as.character(seq(0,1, length=3)), breaks = as.character(seq(0,1, length=3)), expand = c(0,0.05))
        #p1 = p1 + ylim(0,1)
        png(file=paste(name_file,"_weblogo.png", sep=""),height=149 * 2,width=52 * size * 2 )
        print(p1)
        dev.off()
    }
    """)
    weblogo_maker(v.StrVector(sequence_list), v.StrVector([output + sequence_name]), v.StrVector([sequence_name]), v.IntVector([len(sequence_list[0])]))


def make_density_hist_ctrl(val_median_nt_frequency, ctrl_nt_freq, nt, exon_type, name_project, output, nb_iteration):
    """
    Create a density plot showing the deltapsi in ``name_ctrl`` of exon exon _list ``nt`` \
    and ``name_list2``
    :param val_median_nt_frequency: (float) the mean value of dpsi for exon regulatedby ``nt`` (dpsi value taken \
    from ``name_ctrl``
    :param ctrl_nt_freq:(list of float) list of delta_psi in ``name_ctrl`` of ctrl exon ``name_list2``
    in ``name_ctrl``
    :param nt: (string) the name of the list of exon that have the list of delta_psi ``dpsi_list_1`` \
    in ``name_ctrl``
    :param exon_type: (string) the name of the list of exon that have the list of delta_psi ``dpsi_list_2`` \
    in ``name_ctrl``
    :param name_project: (string) the project for which we have the dpsi ``dpsi_list_1`` and ``dpsi_list_2``
    :param output: (string) path where figure will be created
    :param nb_iteration: (int) the number of iteration
    """
    print(val_median_nt_frequency)
    print(ctrl_nt_freq[0:5])
    print("mean : %s" % np.mean(ctrl_nt_freq))
    pval = min(sum(np.array(ctrl_nt_freq) >= val_median_nt_frequency) / len(ctrl_nt_freq),
               sum(np.array(ctrl_nt_freq) <= val_median_nt_frequency) / len(ctrl_nt_freq))
    hist_data = [ctrl_nt_freq]
    group_labels = [exon_type]
    fig = ff.create_distplot(hist_data, group_labels,
                             curve_type='kde', show_rug =False, show_hist=True,
                             colors=['#1f77b4', '#2ca02c'])
    main_title = 'Distribution of %s frequency in %s between %s and %s exon_list<br> iteration %s - pvalue = %s' % (nt, name_project, name_project, exon_type, nb_iteration, pval)

    x_title = '%s frequency' % nt
    y_title = 'data frequency'
    figname = '%s%s_density_%s.html' % (output, nt, name_project)
    fig['layout'].update(title=main_title)
    fig['layout'].update(xaxis=dict(title=x_title))
    fig['layout'].update(yaxis=dict(title=y_title))
    fig['layout'].update(shapes=[dict(type="line", x0=val_median_nt_frequency, y0=0, x1=val_median_nt_frequency, y1=1, line=dict(color="red", width=2))])
    plotly.offline.plot(fig, filename=figname, auto_open=False)


def make_control_figures(median_test_dic, median_ctrl_dic, output, exon_type, name_project, nb_iteration):
    """
    Create control figure for every wanted nucleotides.

    :param median_test_dic: (dictionary of float) links each nucleotide to their frequencies in \
    the test peaks.
    :param median_ctrl_dic: (dictionary of list of float) links each nucleotide to their frequencies in \
    the ``nb_iteration`` control set created
    :param output: (string) path where the results will be created
    :param exon_type: (string) control exon used
    :param name_project: (string) the name of the project studied
    :param nb_iteration: (int) the number of iteration used to create the control sets of exons
    """
    for nt in median_test_dic.keys():
        print(nt)
        make_density_hist_ctrl(median_test_dic[nt], median_ctrl_dic[nt], nt, exon_type, name_project, output,
                               nb_iteration)