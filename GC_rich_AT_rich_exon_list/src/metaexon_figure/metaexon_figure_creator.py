from matplotlib import pyplot as plt
import sys
import os


def sum_dic(vector_dic, nt):
    """

    :param vector_dic: (dic of list of 150 floats), a list (key of the dic) where each float gives the \
    frequencies of a nucleotide. within a windows of x nucleotides in a particular position (center of the window)
    in a sequence of interest. The position of the float value in the list corresponds to
    the position of the sequence of interest where the frequencies of nucleotides is calculated
    :param nt: (string) a nucleotide in iupac code
    :return: (list of 150 float) the frequencies of the interest nucleotide nt within a window centered at each \
    position of the interest region.
    """
    iupac = {'Y': ['C', 'T'], 'R': ['A', 'G'], 'W': ['A', 'T'], 'S': ['G', 'C'], 'K': ['T', 'G'], 'M': ['C', 'A'],
             'D': ['A', 'G', 'T'], 'V': ['A', 'C', 'G'], 'H': ['A', 'C', 'T'], 'B': ['C', 'G', 'T']}
    res_list = []
    if nt in ['A', 'C', 'G', 'T']:
        return vector_dic[nt]
    else:
        if len(iupac[nt]) == 2:
            for i in range(len(vector_dic['A'])):
                res_list.append(vector_dic[iupac[nt][0]][i] + vector_dic[iupac[nt][1]][i])
        else:
            for i in range(len(vector_dic['A'])):
                res_list.append(vector_dic[iupac[nt][0]][i] + vector_dic[iupac[nt][1]][i] + vector_dic[iupac[nt][2]][i])
        return res_list


def select_modules(exon_type):
    """
    Select the control exons dictionary to load
    :param exon_type: (string) the control exon that we want to use.
    :return: (2 lists of floats) the dictionaries of lists of 150 floats
    """
    dir_path = os.path.dirname(os.path.realpath(__file__))
    dir_path += "/control_dictionaries/"
    sys.path.insert(0, dir_path)
    cce_mod = __import__("%s_metagene_windowsed" % exon_type)
    cce_5p = cce_mod.final_res_5p
    cce_3p = cce_mod.final_res_3p
    return cce_5p, cce_3p


def get_ordinate_value(list_tuple):
    """
    :param list_tuple: (a tuple containing lists or None values)
    :return: (2 int) the maximun and the minimum values of the values in ``list_tuple``
    """
    super_list = []
    for cur_list in list_tuple:
        if isinstance(cur_list, list):
            super_list += cur_list
    max_ord = max(super_list)
    max_ord += max_ord * 0.05
    min_ord = min(super_list)
    return max_ord, min_ord


def make_metagene_graphics(output, p5_vector_list, p3_vector_list, nt,
                           set_name, list_name, color_list, exon_type, legend=False):
    """

    :param output: (string) the path where the graphic will be created
    :param p5_vector_list:  (list dic of list of 150 floats), a dict where each list gives the frequencies of \
     nucleotides in the 5p_sequence of the up-regulated exon set. The position of a dictionary in the list corresponds \
     to the position of 5p_sequence where the frequencies of nucleotides is calculated. Each list corresponds to a \
     different list of exons
    :param p3_vector_list: (list dic of list of 150 floats), a dict where each list gives the frequencies of nucleotides
    in the 3p_sequence of the up-regulated exon. The position of a dictionary in the list corresponds to the position
    of 3p_sequence where the frequencies of nucleotides is calculated. Each list corresponds to a different list of \
    exons.
    :param nt: (string) the nt  for which me want to create the metaexon
    :param set_name: (string) the end of the title within the graphic
    :param list_name: (list of strings string) the names of the exons list given to produce the vectors \
    ``p5_vector_list`` and ``p3_vector_list``
    :param color_list: (list of strings) list of colors that will be use for each exon list we want to display on \
    the metaexon graphics
    :param exon_type: (string) the type of  control exons used.
    :param legend: (boolean) True to draw the exon + display the legned, False else
    """
    # charging the dictionaries we need
    # path we want to use...
    if len(p5_vector_list) != len(p3_vector_list):
        print("ERROR : p5_vector_list and p3_vector_list must have the same length !")
        exit(1)
    if None in p5_vector_list or None in p3_vector_list:
        print("Error : None value found in p5_vector_list or p3_vector_list")
        exit(1)
    cce_5p, cce_3p = \
        select_modules(exon_type)

    # --------- sequence 5p lists --------------------
    nt_list_5p = []
    for my_dic_vector in p5_vector_list:
        nt_list_5p.append(sum_dic(my_dic_vector, nt))
    nt_cce_5p = sum_dic(cce_5p, nt)
    # --------------- sequence 3p list ------------------
    nt_list_3p = []
    for my_dic_vector in p3_vector_list:
        nt_list_3p.append(sum_dic(my_dic_vector, nt))
    nt_cce_3p = sum_dic(cce_3p, nt)

    # -------------------- getting ordinated limit values ---------------------
    super_tuple = nt_list_5p + nt_list_3p + [nt_cce_5p] + [nt_cce_3p]
    max_ord, min_ord = get_ordinate_value(super_tuple)

    fig = plt.figure(figsize=(11, 4))
    ax = plt.subplot2grid((1, 2), (0, 0))
    ax.plot(range(1, 151), nt_cce_5p, color="#FF0000", marker=None, linewidth=2.5, label="CCE")
    for i in range(len(nt_list_5p)):
        ax.plot(range(1, 151), nt_list_5p[i], color=color_list[i], marker=None, linewidth=2, label=list_name[i])

    # drawing exons
    if legend:
        ax.plot([0, 100], [max_ord - max_ord * 0.05, max_ord - max_ord * 0.05], 'k-', lw=2)  # intron line
        ax.plot([100, 100], [max_ord - max_ord * 0.1, max_ord], 'k-', lw=2)  # exon vertical line
        # exon horizontal lines
        ax.plot([100, 150], [max_ord - max_ord * 0.1, max_ord - max_ord * 0.1], 'k-', lw=2)
        ax.plot([100, 150], [max_ord, max_ord], 'k-', lw=2)
        ax.plot([100, 100], [min_ord, max_ord], 'k--', lw=0.5)
        ax.text(50, max_ord - max_ord * 0.03, 'Intron')
        ax.text(125, max_ord - max_ord * 0.05, 'Exon')

        ax.legend(loc="upper left", shadow=True, numpoints=1, fontsize='xx-small')

    ax.set_title(u"5' region of " + set_name + u" exons\nnucleotide : " + nt, fontsize=11)

    # ax.set_xlabel(u"Nombre de " + str(unit_name) + u" dans les exons" )
    ax.set_ylabel(u"Nucleotide proportion " + nt)
    ax.set_ylim([min_ord, max_ord])
    fig.add_subplot(ax)

    # --------------- sequence 3p -----------------------
    ax2 = plt.subplot2grid((1, 2), (0, 1))
    ax2.plot(range(1, 151), nt_cce_3p, color="#FF0000", marker=None, linewidth=2.5, label="CCE")
    for i in range(len(nt_list_3p)):
        ax2.plot(range(1, 151), nt_list_3p[i], color=color_list[i], marker=None, linewidth=2, label=list_name[i])

    # drawing exons
    if legend:
        ax2.plot([50, 150], [max_ord - max_ord * 0.05, max_ord - max_ord * 0.05], 'k-', lw=2)  # intron line
        ax2.plot([50, 50], [max_ord - max_ord * 0.1, max_ord], 'k-', lw=2)  # exon vertical line
        # exon horizontal lines
        ax2.plot([0, 50], [max_ord - max_ord * 0.1, max_ord - max_ord * 0.1], 'k-', lw=2)
        ax2.plot([0, 50], [max_ord, max_ord], 'k-', lw=2)
        ax2.plot([50, 50], [min_ord, max_ord], 'k--', lw=0.5)
        ax2.text(75, max_ord - max_ord * 0.03, 'Intron')
        ax2.text(25, max_ord - max_ord * 0.05, 'Exon')

        ax2.legend(loc="upper left", shadow=True, numpoints=1, fontsize='xx-small')
    ax2.set_title(u"3' region of " + set_name + u" exons\nnucleotide : " + nt, fontsize=11)
    ax2.set_ylim([min_ord, max_ord])
    fig.add_subplot(ax2)
    plt.tight_layout()
    #  plt.savefig(output + set_name + "_metaexon_figure_" + str(nt) + ".png", bbox_inches='tight')
    plt.savefig(output + set_name + "_metaexon_figure_" + str(nt) + ".pdf", bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()
