# create a class specifying feature colours
# make feature_list and colour_dict so that each feature name and colour is only added once to legend if more than one share the same name
from itertools import cycle

from dna_features_viewer import BiopythonTranslator

feature_list = []
colour_dict = {}
# set probel_labels on or off here
probe_box = False
probe_labels = False
# turn on or off open chromatin peak annotations
open_chromatin = False
# set font size
font_size = 9


class MyCustomTranslator(BiopythonTranslator):
    """Custom translator implementing the following theme:
    -Colour promoter in pale green
    -colour exons in dark grey
    -colour introns in light grey
    -colour TFs from colour palette"""

    # import colour blind palette
    # colour palette from mkweb.bcgsc.ca/colorblind/palettes.mhtml
    # df = pd.read_csv("colourblind_palette.csv", header=0)
    # #convert to floats
    # floats=df.divide(255)
    # #make sets of each row to get the red, green blue colours
    # CB_colour_palette = floats.apply(tuple, axis = 1)
    # #make into df
    # df = CB_colour_palette.to_frame()
    # #save file
    # df.to_csv('../../data/TRAMP/colour_list')

    # #turn into a list of colour sets
    # list_colours = list(CB_colour_palette)
    # colour_list = list_colours
    # colour_list = ['#88CCEE', '#44AA99', '#117733', '#332288', '#DDCC77', '#999933','#CC6677', '#882255', '#AA4499', '#DDDDDD']
    colour_list = [
        "#2271B2",  # ARF9/18
        "#d55e00",  # ANAC032
        "#228833",  # NLP6/7
        "#3DB7E9",  # ANR1
        "#e69f00",  # DREB26
        "#F748A5",  # TGA1
        # '#f0e442',
        "#900000",
    ]
    # make colour iterator
    colour_iterator = cycle(colour_list)
    # change colour cycle
    # mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=list_colours)
    # count = -1
    # set colour count

    def compute_feature_color(self, feature):
        """set colour of each feature"""

        if feature.type == "promoter":
            return "#F5F5F5"
        elif feature.type == "gene_upstream":
            return "#DDCC77"  # 2f4f4f"#dark slate grey
        elif feature.type == "mRNA_upstream":
            return "#DDCC77"  # dark slate grey
        elif feature.type == "exon_upstream":
            return "#DDCC77"  # dark slate grey
        elif feature.type == "exon":
            return "#635147"  # umber
        elif feature.type == "gene":
            return "#F5F5F5"
            # return (169,169,169)
        elif feature.type == "intron":
            return "lightgrey"
            # return (211,211,211)
        elif feature.type == "5'UTR":
            return "c4aead"  # silver pink
        elif feature.type == "start_codon":
            return "black"
        elif feature.type == "TRAMP_probe":
            return "#FFFFFF00"  # transparent
        #     if feature.qualifiers.get("label")[0] == "NLP7#7":
        #         col = "#e69f00"
        #     elif feature.qualifiers.get("label")[0] == "NLP7#9":
        #         col = "#d55e00"
        #     elif feature.qualifiers.get("label")[0] == "NLP7#10":
        #         col = "#2271B2"
        #     elif feature.qualifiers.get("label")[0] == "ANAC032#4":
        #         col = "#2271B2"
        #     elif feature.qualifiers.get("label")[0] == "ANAC032#8":
        #         col = "#F748A5"
        #     elif feature.qualifiers.get("label")[0] == "ANAC032#3":
        #         col = "#000000"

        #     else:
        #         pass

        #     return col

        elif feature.type == "root_open_chromatin":
            if open_chromatin is True:
                # transparency 25%
                return "#a52a2a40"

            else:
                pass

        elif feature.type == "shoot_open_chromatin":
            if open_chromatin is True:
                # transparency 25%
                return "#00808040"
            else:
                pass

        elif feature.type == "TFBS":
            if feature.qualifiers.get("label")[0] in colour_dict.keys():
                col = colour_dict[feature.qualifiers.get("label")[0]]
            else:
                col = next(self.colour_iterator)
                colour_dict[feature.qualifiers.get("label")[0]] = col

            return col
        else:
            return "black"

    def compute_feature_box_linewidth(self, feature):
        """change feature box linewidth"""
        if feature.type == "start_codon":
            return 1
        if feature.type == "TRAMP_probe_tested":
            if probe_labels is True:
                print("probe labels = True")
                # if prob labels = True, add probe labels
                return 1
        else:
            return 0

    def compute_feature_thickness(self, feature):
        """change feature vertical thickness"""
        if open_chromatin is True:
            if feature.type == "root_open_chromatin":
                return 700
            if feature.type == "shoot_open_chromatin":
                return 0.1
        else:
            return 14

    def compute_feature_linewidth(self, feature):
        """change linewidth of feature's arrow/rectangle"""
        # remove border from certain features
        if feature.type == "gene_upstream":
            return 0
        if feature.type == "start_codon":
            return 1
        elif feature.type == "mRNA_upstream":
            return 0
        elif feature.type == "exon_upstream":
            return 0
        elif feature.type == "misc_RNA_upstream":
            return 0
        elif feature.type == "exon":
            return 0
        elif feature.type == "intron":
            return 0
        elif feature.type == "5'UTR":
            return 0
        elif feature.type == "TFBS":
            return 0
        elif feature.type == "TRAMP_probe":
            if probe_box is True:
                # if probe box = True, add probe box
                return 1
            elif probe_box is False:
                return 0
        # open chromatin features
        if open_chromatin is True:
            if feature.type == "root_open_chromatin":
                return 0
            if feature.type == "shoot_open_chromatin":
                return 0
        return 0

    # this doesn't seem to work - so I forced it to work when calling this class instead (see make_graphic_record() function within make plot function)
    # def compute_feature_linecolor(self, feature):
    #     if feature.type=="root_open_chromatin":
    #         return "brown"
    #     if feature.type=="shoot_open_chromatin":
    #         return "teal"

    #     if feature.type == "TRAMP_probe_tested":
    #         if probe_labels is True:
    #             if feature.qualifiers.get("label")[0] == "NLP7#7":
    #                 col = "#e69f00"
    #             elif feature.qualifiers.get("label")[0] == "NLP7#9":
    #                 col = "#d55e00"
    #             elif feature.qualifiers.get("label")[0] == "NLP7#10":
    #                 col = "#2271B2"
    #             elif feature.qualifiers.get("label")[0] == "ANAC032#4":
    #                 col = "#2271B2"
    #             elif feature.qualifiers.get("label")[0] == "ANAC032#8":
    #                 col = "#F748A5"
    #             elif feature.qualifiers.get("label")[0] == "ANAC032#3":
    #                 col = "#000000"
    #         else:
    #             pass

    # #open chromatin features
    # if open_chromatin is True:

    #     else:
    #     #         pass
    # def compute_linecolor(self, feature):
    #     if feature.type=="root_open_chromatin":
    #         return "brown"
    #     if feature.type=="shoot_open_chromatin":
    #         return "teal"

    def compute_feature_box_color(self, feature):
        """change colour of feature box border"""
        if feature.type == "TRAMP_probe_tested":
            if feature.qualifiers.get("label")[0] == "NLP7#7":
                # col = colour_dict["NLP6/7"]
                col = "#228833"
            elif feature.qualifiers.get("label")[0] == "NLP7#9":
                # col = colour_dict["DREB26"]
                col = "#e69f00"
            elif feature.qualifiers.get("label")[0] == "NLP7#10":
                # col = colour_dict["ANAC032*"]
                col = "#d55e00"
            elif feature.qualifiers.get("label")[0] == "ANAC032#4":
                # col = colour_dict["ANAC032*"]
                col = "#d55e00"
            elif feature.qualifiers.get("label")[0] == "ANAC032#8":
                # col = colour_dict["ARF9/18*"]
                col = "#2271B2"
            elif feature.qualifiers.get("label")[0] == "ANAC032#3":
                # col = colour_dict["TGA1"]
                col = "#F748A5"
            else:
                return "black"

            return col
        if feature.type == "root_open_chromatin":
            return "brown"
        if feature.type == "shoot_open_chromatin":
            return "teal"
            # 'ARF9/18*': '#2271B2', 'ANR1*': '#3DB7E9', 'TGA1': '#F748A5', 'ANAC032*': '#d55e00', 'DREB26': '#e69f00', 'NLP6/7': '#228833'

    def compute_feature_open_left(self):
        """set to true if feature does not end on the left"""
        return False

    def compute_feature_open_right(self):
        """set to true if feature does not end on the right"""
        return False

    def compute_feature_label(self, feature):
        """Remove most feature labels"""
        if feature.type == "start_codon":
            return "ATG"
        if feature.type == "TRAMP_probe_tested":
            if probe_labels is True:
                return feature.qualifiers.get("label")[0]
        else:
            pass

        # return super().compute_feature_label(feature)

    # def compute_feature_min_y_height_of_text_line(self, feature):
    #     return 0.1

    def compute_feature_fontdict(self, feature):
        """change label font to arial, fontsize to 9"""
        if feature.type == "TRAMP_probe_tested":
            # if certain label, align to the right
            if feature.qualifiers.get("label")[0] == "ANAC032#3":
                return dict(
                    family="sans-serif", fontsize=font_size, ha="right"
                )
            # elif feature.qualifiers.get("label")[0] == "ANAC032#8":
            #     return dict(family='sans-serif',fontsize=10, ha='left')
            if feature.qualifiers.get("label")[0] == "ANAC032#8":
                return dict(family="sans-serif", fontsize=font_size)
            else:
                return dict(family="sans-serif", fontsize=font_size)

        else:
            return dict(family="sans-serif", fontsize=font_size)
        # return dict([('family','sans-serif'),('sans-serif','Arial'),('size',10)])

    # make feature_list so that each feature name is only added once if more than one share the same name
    # feature_list = []
    def compute_feature_legend_text(self, feature):
        """add legend if feature label has not been added to legend already"""
        if feature.type == "promoter":
            pass
        # elif feature.type=='exon':
        #     pass
        # elif feature.type=='intron':
        #     pass
        # elif feature.type=="5'UTR":
        #     pass
        # elif feature.qualifiers.get("label")[0] in self.feature_list:
        #     pass
        elif feature.qualifiers.get("label")[0] in feature_list:
            pass
        else:

            feature_list.append(feature.qualifiers.get("label")[0])

            # feature_list.append(feature.qualifiers.get("label")[0])
            return feature.qualifiers.get("label")[0]

    def compute_filtered_features(self, features):
        """Displays feature if feature is not one of the following"""
        if probe_labels is True:
            return [
                feature
                for feature in features
                if (feature.type != "TRAMP_probe")
                and (feature.type != "none")
                and (feature.type != "DHS")
                and (feature.type != "misc_feature")
                and (feature.type != "primer")
                # and (feature.type != "gene")
                and (feature.type != "mRNA")
                and (feature.type != "CDS")
                and (feature.type != "source")
                and (feature.type != "misc_RNA")
                and (feature.type != "EMSA_probe_long")
                and (
                    feature.qualifiers.get("label")[0] != "ARID5_ARID6"
                    and (feature.qualifiers.get("label")[0] != "ARDI5_ARID6")
                    and (
                        "Translation" not in feature.qualifiers.get("label")[0]
                    )
                    # return feature only if it is one of the following:
                    and (
                        (feature.type == "TFBS")
                        # or (feature.type == "root_open_chromatin")
                        # or (feature.type == "shoot_open_chromatin")
                        or (feature.type == "promoter")
                        or (feature.type == "exon")
                        or (feature.type == "intron")
                        or (feature.type == "start_codon")
                        # or (feature.type == "gene")
                        or (feature.type == "exon_upstream")
                        or (feature.type == "mRNA_upstream")
                        or (feature.type == "gene_upstream")
                        or (feature.type == "misc_RNA_upstream")
                        or (feature.type == "5'UTR")
                        or (feature.type == "TRAMP_probe_tested")
                        or (feature.type == "TRAMP_probe")
                    )
                )
            ]

        else:
            return [
                feature
                for feature in features
                if (feature.type != "TRAMP_probe_tested")
                # and (feature.type != "TRAMP_probe")
                and (feature.type != "none")
                and (feature.type != "DHS")
                and (feature.type != "misc_feature")
                and (feature.type != "primer")
                # and (feature.type != "gene")
                and (feature.type != "mRNA")
                and (feature.type != "CDS")
                and (feature.type != "source")
                and (feature.type != "misc_RNA")
                and (feature.type != "EMSA_probe_long")
                and (
                    feature.qualifiers.get("label")[0] != "ARID5_ARID6"
                    and (feature.qualifiers.get("label")[0] != "ARDI5_ARID6")
                    and (
                        "Translation" not in feature.qualifiers.get("label")[0]
                    )
                    # return feature only if it is one of the following:
                    and (
                        (feature.type == "TFBS")
                        # or (feature.type == "root_open_chromatin")
                        # or (feature.type == "shoot_open_chromatin")
                        or (feature.type == "promoter")
                        or (feature.type == "exon")
                        or (feature.type == "intron")
                        or (feature.type == "start_codon")
                        # or (feature.type == "gene")
                        or (feature.type == "exon_upstream")
                        or (feature.type == "mRNA_upstream")
                        or (feature.type == "gene_upstream")
                        or (feature.type == "misc_RNA_upstream")
                        or (feature.type == "5'UTR")
                        or (feature.type == "TRAMP_probe")
                    )
                )
            ]
