from dna_features_viewer import BiopythonTranslator
from matplotlib import pyplot as plt


def draw_region_by_coordinates(gff_file, start_, end_, name, prokka=False, gff_without_fasta=None):
    """
    Function produces schematic plot of DNA region specified with coordinates;
    visualizes only gene feature types, so one have to change them depending on color
    adds lables by their Name
    This code is based on https://github.com/rybinaanya/O-antigens/blob/main/operon_visualization/draw_operon_from_gff.py
    :param gff_file: path to GFF annotation file
    :param start_: start coordinate of the desired region
    :param end_: end coordinate of the desired region
    :param prokka: True or False; specifies type of GFF annotation
    :param gff_without_fasta: default None, if Prokka it's used to create new gff
    :return: None
    """
    
    picture_name = name
    
    color_map = {
        "WZZ": "#0E9594",
        "WZX": "#F7B05B",
        "WZM": "#7C77B9",
        "TDP": "#E3BAC6",
        "UDP": "#646F4B",
        "GDP": "#775253",
        "RFAL": "#982649",
        "UDPN": "#406E8E",
        "TRANS": "#9EE493"
        }
    
    if prokka:
        with open(gff_file, 'r') as gff_with_fasta:
            with open(gff_without_fasta, 'w') as new_gff:
                for line in gff_with_fasta:
                    if line.startswith("##FASTA"):
                        break
                    new_gff.write(line)

        # draw operon by new gff file
        translator = BiopythonTranslator(
        features_filters=(lambda f: f.type not in ["gene", "source"],),
        features_properties=lambda f: {"color": color_map.get(f.type, "white")},
        )
        translator.label_fields = ["Name"]
        graphic_record = translator.translate_record(gff_without_fasta)
        operon = graphic_record.crop((start_, end_))
        #operon.plot(figure_width=10, elevate_outline_annotations=False)
        
        plt.figure(figsize=((10,8)))
        operon.plot(figure_width=10, elevate_outline_annotations=False)
        plt.savefig("%s.pdf" % picture_name)
        plt.show()
    else:
        #biopython_translator = BiopythonTranslator()
        biopython_translator = BiopythonTranslator(features_filters=(lambda f: f.type not in ["gene", "source"],),
        features_properties=lambda f: {"color": color_map.get(f.type, "white")},
        )
        #biopython_translator.ignored_features_types = ['CDS']
        #biopython_translator.label_fields = ["gene", "product"]
        biopython_translator.label_fields = ["Name"]
        graphic_record = biopython_translator.translate_record(gff_file)
    
        graphic_record.crop((start_, end_)).plot(figure_width=10, elevate_outline_annotations=False)
        plt.show()

draw_region_by_coordinates(gff_file, 6117090, 6145406, name, prokka=True, gff_without_fasta="no_fasta.gff")