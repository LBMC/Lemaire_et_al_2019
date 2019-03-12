import rpy2.robjects as robj
import rpy2.robjects.vectors as v
import warnings
from rpy2.rinterface import RRuntimeWarning
import exon_class
import sqlite3
import os


def get_control_exon_information(cnx, exon_type):
    """
    Get the gene symbol, the gene id and the position of every ``exon_type`` exons in fasterDB.

    :param cnx: (sqlite3 object) allow connection to fasterdb database
    :param exon_type: (string) the type of control exon we want to use
    :return:
        * result: (list of tuple) every information about control exons
        * names: (list of string) the name of every column in sed table
    """
    cursor = cnx.cursor()
    if exon_type != "ALL":
        query = """SELECT t2.official_symbol, t1.id_gene, t1.pos_on_gene
                   FROM exons t1, genes t2
                   WHERE t1.exon_type LIKE '%{}%'
                   AND t1.id_gene = t2.id""".format(exon_type)
    else:
        query = """SELECT t2.official_symbol, t1.id_gene, t1.pos_on_gene
                   FROM exons t1, genes t2
                   AND t1.id_gene = t2.id
                """
    cursor.execute(query)
    result = cursor.fetchall()
    # turn tuple into list
    nresult = []
    for exon in result:
        nresult.append(list(exon))
    return nresult


def seq_tranformer(sequence, transfo_type):
    """

    :param sequence: (string) an interest nucleotids sequence
    :param transfo_type: (string) YR or WS or KM
    :return: the sequence with only YR letter, WS letters, KM letters
    """
    if transfo_type == "YR":
        sequence = sequence.replace("C", "Y").replace("T", "Y").replace("A", "R").replace("G", "R").replace("N","X")
    elif transfo_type == "WS":
        sequence = sequence.replace("C", "S").replace("T", "W").replace("A", "W").replace("G", "S").replace("N","X")
    # elif transfo_type == "KM":
    #     sequence = sequence.replace("C", "M").replace("T", "K").replace("A", "M").replace("G", "K").replace("N","X")
    return sequence


def create_3p_sequence_instance(cur_exon_list):
    """
    :param cur_exon_list: (ListExon instance) an exon list
    :return:
    1 - instance_list : (list of list of Seq instance) it is composed of
        instance_5p : list of sequence_5p (Seq instance) - (sequence 5p = 100 nt before the exon of interest and 50
        first nt of the exon) - nt can be A, T, G, C
        instance_3p : list of sequence_3p (Seq instance) - (sequence 3p =  50 last nt of the exon and 100 nt after this
        exon) nt can be A, T, G, C
        instance_5p_yr : list of sequence_5p  - (sequence 5p = 100 nt before the exon of interest and 50
        first nt of the exon) - (same sequence that instance_5p but only with letter Y and R)
        instance_3p_yr : list of sequence_3p - (sequence 3p =  50 last nt of the exon and 100 nt after
        this exon) (same sequence that instance_5p but only with letter Y and R)
        instance_5p_ws : list of sequence_5p - (sequence 5p = 100 nt before the exon of interest and 50
        first nt of the exon) - (same sequence that instance_5p but only with letter W and S)
        instance_3p_ws : list of sequence_3p  - (sequence 3p =  50 last nt of the exon and 100 nt after
        this exon) (same sequence that instance_5p but only with letter W and S)
        instance_5p_km : list of sequence_5p  - (sequence 5p = 100 nt before the exon of interest and 50
        first nt of the exon) - (same sequence that instance_5p but only with letter K and M)
        instance_3p_km : list of sequence_3p  - (sequence 3p =  50 last nt of the exon and 100 nt after
        this exon) (same sequence that instance_5p but only with letter K and M)
    2 - instance_name : (list of string) the name of each list of sequences contained in instance_list

    """
    seq_3p = []
    seq_3p_yr = []
    seq_3p_ws = []
    for exon in cur_exon_list:
        if not isinstance(exon, list):
            if exon.upstream_intron is not None and exon.upstream_intron.sequence_proxi is not None:
                seq_3p.append((exon.upstream_intron.sequence_proxi))
                seq_3p_yr.append(seq_tranformer(seq_3p[-1], "YR"))
                seq_3p_ws.append(seq_tranformer(seq_3p[-1], "WS"))

    seq_3p = v.StrVector(seq_3p)
    seq_3p_yr = v.StrVector(seq_3p_yr)
    seq_3p_ws = v.StrVector(seq_3p_ws)

    sequence_list = (seq_3p, seq_3p_yr, seq_3p_ws)
    sequence_name = ( "3p", "3p_yr", "3p_ws")
    return sequence_list, sequence_name


def web_logo_creator(sequence_list, sequence_name, output, regulation, name_exon):
    """
    :param sequence_list: (tuple of  3 list of strings) - each list in the tuple corresponds to a list of sequence
    :param sequence_name: (list of string) each string identifies on list of sequence in sequence_list
    :param output: (string) the folder where the results will be created
    :param regulation: (string) the name of the regulation
    """
    warnings.filterwarnings("ignore", category=RRuntimeWarning)
    weblogo_maker = robj.r("""
    library("ggplot2")
    library("ggseqlogo")

    function(mys_seq, name_file, mytitle){
        s1 = 15
        cs1 = make_col_scheme(chars=c('A','T','G','C', 'R', 'Y', 'W', 'S', 'K', 'M'), groups=c('g1','g2','g3','g4','g5', 'g6', 'g7', 'g8', 'g9', 'g10'),cols=c('limegreen','brown1','gold','dodgerblue3','darkorange', "brown1", "limegreen", "dodgerblue3", "darkorchid3", "dodgerblue3"), name='custom1')

        p1 = ggseqlogo(mys_seq,  method = "probability", col_scheme=cs1, namespace = c('A','T','G','C', 'R', 'Y', 'W', 'S', 'K', 'M')) +  annotate('text', x=-0.2, y=0.5, label="5\'", size = s1 + 5) + annotate('text', x=100+1.5, y=0.5, label="3\'", size = s1 + 5) + theme_logo() + scale_x_discrete(limits = as.character(seq(-100,-1, by=1)), labels = as.character(seq(-100,-1, by=5)), breaks = as.character(seq(-100, -1, by=5))) + theme(axis.title.y=element_text(size=s1+25), legend.position="none")
        p1 = p1 + ggtitle(mytitle) +  theme(plot.title = element_text(hjust = 0.5))


        p1 = p1 + theme(axis.text=element_text(size=s1 + 25), plot.title = element_text(size=s1 + 30))
        p1 = p1 + scale_y_discrete(limits = c(0, 0.5, 1), labels = as.character(seq(0,1, length=3)), breaks = as.character(seq(0,1, length=3)), expand = c(0,0.05))
        #p1 = p1 + ylim(0,1)
        png(file=name_file,height=159 * 2,width=3305 * 2 * 0.75 )
        print(p1)
        dev.off()
    }
    """)
    for i in range(len(sequence_list)):

        if len(sequence_name[i].split("_")) > 1:
            res = "_" + sequence_name[i].split("_")[1]
        else:
            res = ""

        full_name = output + name_exon + "_exon_" + regulation + "_sequence_" + res + ".png"
        name_file =  regulation + " " + name_exon + " exon - sequence "
        weblogo_maker(sequence_list[i], v.StrVector([full_name]), v.StrVector([name_file]))


def extract_exons_from_file(filename):
    """
    Get the exons in ``filename``
    :param filename: (string) a file containg exons
    :return: (list of list of 2 int) exon_list
    """
    exon_list = []
    with open(filename, "r") as in_file:
        line = in_file.readline()
        while line:
            line = line.split("\t")
            exon_list.append([line[0], int(line[0]), int(line[1])])
            line = in_file.readline()
    return exon_list


def main():

    exon_class.set_debug(0)
    exon_type = "CCE"
    path = os.path.realpath(os.path.dirname(__file__)).replace("src/weblogo_and_pyrimidine_content", "result/")
    output = os.path.realpath(os.path.dirname(__file__)).replace("src/weblogo_and_pyrimidine_content", "result/boxplot_CT_content/")
    if not os.path.isdir(output):
        os.mkdir(output)
    seddb = os.path.realpath(os.path.dirname(__file__)).replace("src/weblogo_and_pyrimidine_content", "data/sed.db")
    fasterdb = os.path.realpath(os.path.dirname(__file__)).replace("src/weblogo_and_pyrimidine_content", "data/fasterDB_lite.db")
    cnx = sqlite3.connect(seddb)
    cnx_fasterdb = sqlite3.connect(fasterdb)
    at_pure_file = "%sAT_rich_exons" % path
    gc_pure_file = "%sGC_rich_exons" % path
    at_all_file = "%sAT_rich_with_intersection_exons" % path
    gc_all_file = "%sGC_rich_with_intersection_exons" % path
    name_file = ["GC_all_exons", "AT_all_exons", "GC_pure_exons", "AT_pure_exons", exon_type]
    list_file = [gc_all_file, at_all_file, gc_pure_file, at_pure_file, None]
    for i in range(len(list_file)):
        print("Working on %s" % name_file[i])
        if name_file[i] == exon_type:
            exon_list = get_control_exon_information(cnx_fasterdb, exon_type)
        else:
            exon_list = extract_exons_from_file(list_file[i])
        print("   ---> Getting upstream sequences")
        object_exon_list = [exon_class.ExonClass(cnx_fasterdb, exon[0], exon[1], exon[2]) for exon in exon_list]
        print("   ---> Creating sequence")
        sequence_list, sequence_name = create_3p_sequence_instance(object_exon_list)
        print("   ---> Creating weblogos")
        web_logo_creator(sequence_list, sequence_name, output, "down", name_file[i])
    cnx.close()
    cnx_fasterdb.close()

if __name__ == "__main__":
    main()
